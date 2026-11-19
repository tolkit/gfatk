//! Resolve a repeat-rich assembly graph into a circular genome using PacBio
//! HiFi read-path evidence (GAF), rather than the local-coverage greedy walk
//! used by [`crate::linear`].
//!
//! The approach automates what a human does when curating a tangled
//! mitogenome graph in Bandage: look at real read depth to work out how many
//! times each segment truly appears (its "copy number"), and combine that
//! with real read-observed adjacencies to decide which neighbours belong
//! together. Concretely:
//!
//!   - each segment's **budget** (how many times it may appear) comes from
//!     its GAF-derived read depth relative to the graph's baseline (typical
//!     single-copy) depth
//!   - **candidate edges** come from three progressively weaker evidence
//!     tiers: an exact `(prev, segment, next)` triple seen in a single read,
//!     a plain pairwise adjacency seen in a read's path, or (as a last
//!     resort) the bare GFA link topology with no read support at all
//!   - the actual selection is posed as a small integer program: choose the
//!     set of candidate edges that maximises coverage (weighted towards
//!     annotated genes, when a GFF is supplied), subject to each segment's
//!     budget and flow conservation at every state -- which guarantees the
//!     selected edges decompose cleanly into one or more closed circuits,
//!     with no dangling ends
//!
//! Segments left out entirely, but which have real read evidence connecting
//! them to both sides of the resolved circuit, are reported separately as
//! recombination-bubble arms: a real, alternate local conformation that lost
//! out to the dominant path, not a data gap.

use anyhow::{Context, Result};
use gfa::gfa::Orientation;
use good_lp::{
    variable, Expression, ProblemVariables, Solution, SolutionStatus, SolverModel, Variable,
    WithTimeLimit,
};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::gaf::{read_gaf_evidence, GafRecord};
use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;

type Seg = Vec<u8>;
type State = (Seg, Orientation);
/// Read support for each `(prev, segment) -> next` triple; see
/// [`build_triple_lookup`].
type TripleLookup = HashMap<(Seg, Orientation, Seg, Orientation), HashMap<State, u32>>;

/// A candidate transition between two (segment, orientation) states, with
/// the strength of evidence backing it.
#[derive(Debug, Clone)]
struct CandidateEdge {
    from: State,
    to: State,
    /// 0 = exact triple context, 1 = plain pairwise adjacency, 2 = bare GFA
    /// topology with no read support at all. Lower is stronger.
    tier: u8,
    reads: u32,
}

const MIN_EDGE_READS: u32 = 5;
const GENE_WEIGHT: f64 = 20.0;
const BASE_WEIGHT: f64 = 1.0;

fn flip(o: Orientation) -> Orientation {
    match o {
        Orientation::Forward => Orientation::Backward,
        Orientation::Backward => Orientation::Forward,
    }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => *other,
        })
        .collect()
}

fn orient_char(o: Orientation) -> char {
    match o {
        Orientation::Forward => '+',
        Orientation::Backward => '-',
    }
}

/// Per-segment read depth (aligned bases within that segment / its length),
/// attributing each read's aligned span to the specific segments of its
/// path it actually overlaps, not just segments the path merely touches.
///
/// A GAF's path coordinates count the bases two linked segments share once,
/// so each segment starts its link's overlap before the previous one ends.
/// (Laying segments end to end instead puts every later segment too far
/// along, by all the overlaps before it, and credits bases to the wrong
/// segment.) Reads lying within a single segment count too: leaving them
/// out, as versions up to 0.6.3 did, measures only reads that cross a
/// junction, which a long segment has few of for its length - segments over
/// 30 kb read at 0.4x the depth of the rest. Only reads aligned over most
/// of their length count (`counts_for_depth`, `single_segment_bases`; see
/// [`crate::gaf::read_gaf_evidence`]).
fn compute_depth(
    records: &[GafRecord],
    single_segment_bases: &HashMap<Seg, u64>,
    seg_len: &HashMap<Seg, usize>,
    overlaps: &HashMap<(State, State), usize>,
) -> HashMap<Seg, f64> {
    let mut aligned_bases: HashMap<Seg, u64> = single_segment_bases.clone();

    for rec in records.iter().filter(|r| r.counts_for_depth) {
        let mut offset: usize = 0;
        for (i, (seg, _o)) in rec.path.iter().enumerate() {
            let Some(&len) = seg_len.get(seg) else {
                continue;
            };
            if i > 0 {
                let shared = overlaps
                    .get(&(rec.path[i - 1].clone(), rec.path[i].clone()))
                    .copied()
                    .unwrap_or(0);
                offset = offset.saturating_sub(shared);
            }
            let seg_start = offset;
            let seg_end = offset + len;
            let ov_start = seg_start.max(rec.path_start);
            let ov_end = seg_end.min(rec.path_end);
            if ov_end > ov_start {
                *aligned_bases.entry(seg.clone()).or_insert(0) += (ov_end - ov_start) as u64;
            }
            offset += len;
        }
    }

    seg_len
        .iter()
        .map(|(seg, &len)| {
            let bases = aligned_bases.get(seg).copied().unwrap_or(0) as f64;
            (seg.clone(), bases / len as f64)
        })
        .collect()
}

/// The depth of single-copy sequence, and each segment's estimated copy
/// number relative to it, floored at 1 for segments with too little
/// evidence to say otherwise.
///
/// The baseline is the depth at the middle base of the graph: the median
/// over segments with reads, each weighted by its length. Most of an
/// organelle genome is single-copy, so most *bases* are at single-copy
/// depth - which is not true of most *segments* in a graph with many short
/// repeat or low-coverage pieces, and those would set an unweighted median.
fn estimate_budget(depth: &HashMap<Seg, f64>, seg_len: &HashMap<Seg, usize>) -> (f64, HashMap<Seg, u32>) {
    let mut values: Vec<(f64, usize, &Seg)> = depth
        .iter()
        .filter(|(_, d)| **d > 0.0)
        .map(|(seg, &d)| (d, seg_len.get(seg).copied().unwrap_or(1).max(1), seg))
        .collect();
    values.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap().then_with(|| a.2.cmp(b.2)));
    let total: usize = values.iter().map(|v| v.1).sum();
    let mut baseline = 1.0;
    let mut seen = 0usize;
    for (d, len, _) in &values {
        seen += len;
        if 2 * seen >= total {
            baseline = *d;
            break;
        }
    }

    let budget = depth
        .iter()
        .map(|(seg, &d)| {
            let ratio = d / baseline;
            let b = ratio.round().max(1.0) as u32;
            (seg.clone(), b)
        })
        .collect();

    (baseline, budget)
}

/// The most specific evidence tier: an exact `(prev, segment, next)` triple
/// actually seen in a read, canonicalised so a read sequenced from either
/// strand of the same molecule counts as the same evidence.
fn build_triple_lookup(records: &[GafRecord]) -> TripleLookup {
    let mut raw: HashMap<(Seg, Orientation, Seg, Orientation, Seg, Orientation), u32> =
        HashMap::new();
    for rec in records {
        let path = &rec.path;
        if path.len() < 3 {
            continue;
        }
        for i in 1..path.len() - 1 {
            let (prev_seg, prev_o) = path[i - 1].clone();
            let (seg, o) = path[i].clone();
            let (next_seg, next_o) = path[i + 1].clone();
            *raw.entry((prev_seg, prev_o, seg, o, next_seg, next_o))
                .or_insert(0) += 1;
        }
    }

    let mut canon: HashMap<(Seg, Orientation, Seg, Orientation, Seg, Orientation), u32> =
        HashMap::new();
    let mut seen: HashSet<(Seg, Orientation, Seg, Orientation, Seg, Orientation)> = HashSet::new();
    for (key, &cnt) in raw.iter() {
        if seen.contains(key) {
            continue;
        }
        let (a, ao, seg, so, b, bo) = key.clone();
        let rev = (b, flip(bo), seg, flip(so), a, flip(ao));
        seen.insert(key.clone());
        seen.insert(rev.clone());
        let rev_cnt = raw.get(&rev).copied().unwrap_or(0);
        let canon_key = if *key <= rev { key.clone() } else { rev };
        canon.insert(canon_key, cnt + rev_cnt);
    }

    let mut lookup: TripleLookup = HashMap::new();
    for ((a, ao, seg, so, b, bo), cnt) in canon {
        if cnt < MIN_EDGE_READS {
            continue;
        }
        lookup
            .entry((a.clone(), ao, seg.clone(), so))
            .or_default()
            .insert((b.clone(), bo), cnt);
        // the same molecule read from the other strand walks
        // b' -> seg' -> a', so *every* state is flipped -- including `a`.
        // Leaving `a` unflipped registers a transition `seg' -> a` that
        // no GFA link licenses.
        lookup
            .entry((b, flip(bo), seg, flip(so)))
            .or_default()
            .insert((a, flip(ao)), cnt);
    }
    lookup
}

/// A broader, weaker tier: every consecutive segment pair actually seen in
/// any read's path, regardless of what precedes it. This is what recovers
/// real adjacencies whose exact-triple registration is missing purely
/// because reads happen to start right at that segment.
fn build_edge_support(records: &[GafRecord]) -> HashMap<State, HashMap<State, u32>> {
    let mut edge_support: HashMap<State, HashMap<State, u32>> = HashMap::new();
    for rec in records {
        for w in rec.path.windows(2) {
            *edge_support
                .entry(w[0].clone())
                .or_default()
                .entry(w[1].clone())
                .or_insert(0) += 1;
        }
    }
    // A read from the other strand of the same molecule walks the same
    // junction as v' -> u': pool the two before thresholding, so support
    // split between strands (3 reads one way, 4 the other) still counts,
    // and both traversals of a junction carry the same support.
    let mut pooled: HashMap<State, HashMap<State, u32>> = HashMap::new();
    for (u, m) in &edge_support {
        for (v, &cnt) in m {
            let (ru, rv) = ((v.0.clone(), flip(v.1)), (u.0.clone(), flip(u.1)));
            let mirrored = (&ru, &rv) != (u, v);
            let other = if mirrored {
                edge_support.get(&ru).and_then(|m| m.get(&rv)).copied().unwrap_or(0)
            } else {
                0 // its own reverse complement: one junction, counted once
            };
            pooled.entry(u.clone()).or_default().insert(v.clone(), cnt + other);
            if mirrored {
                pooled.entry(ru).or_default().insert(rv, cnt + other);
            }
        }
    }
    pooled.retain(|_, m| {
        m.retain(|_, &mut cnt| cnt >= MIN_EDGE_READS);
        !m.is_empty()
    });
    pooled
}

/// The bare GFA topology. A link `L A oA B oB` licenses two traversals:
/// `A(oA) -> B(oB)` as written, and its reverse-complement equivalent
/// `B(flip oB) -> A(flip oA)` -- both are registered.
fn build_links_index(
    links: &[(Seg, Orientation, Seg, Orientation)],
) -> HashMap<State, HashMap<State, u32>> {
    let mut idx: HashMap<State, HashMap<State, u32>> = HashMap::new();
    for (fs, fo, ts, to) in links {
        let fwd = ((fs.clone(), *fo), (ts.clone(), *to));
        let rev = ((ts.clone(), flip(*to)), (fs.clone(), flip(*fo)));
        // a link that is its own reverse complement (e.g. `L x + x -`)
        // is still only one link
        let traversals = if fwd == rev { vec![fwd] } else { vec![fwd, rev] };
        for (u, v) in traversals {
            *idx.entry(u).or_default().entry(v).or_insert(0) += 1;
        }
    }
    idx
}

fn consider(best: &mut HashMap<(State, State), (u8, u32)>, u: State, v: State, cnt: u32, tier: u8) {
    let key = (u, v);
    let better = match best.get(&key) {
        None => true,
        Some(&(t, c)) => tier < t || (tier == t && cnt > c),
    };
    if better {
        best.insert(key, (tier, cnt));
    }
}

/// Union every evidence tier into one candidate-edge pool, keeping the best
/// (strongest tier, then highest read count) evidence for each distinct
/// (from, to) state transition.
///
/// Read-derived transitions (tiers 0 and 1) are only accepted when the GFA
/// itself licenses them: a read can't walk a junction that doesn't exist in
/// the graph it was aligned to, so one that claims to (a GAF aligned to a
/// different graph, or a parsing slip) must not reach the solver. Returns
/// the edges alongside the number of distinct read-derived transitions
/// rejected this way.
fn build_candidate_edges(
    triple_lookup: &TripleLookup,
    edge_support: &HashMap<State, HashMap<State, u32>>,
    links_index: &HashMap<State, HashMap<State, u32>>,
) -> (Vec<CandidateEdge>, usize) {
    let mut best: HashMap<(State, State), (u8, u32)> = HashMap::new();
    let mut unlicensed: HashSet<(State, State)> = HashSet::new();
    let licensed = |u: &State, v: &State| links_index.get(u).is_some_and(|m| m.contains_key(v));

    let mut consider_read_edge = |best: &mut HashMap<(State, State), (u8, u32)>,
                                  u: State,
                                  v: State,
                                  cnt: u32,
                                  tier: u8| {
        if licensed(&u, &v) {
            consider(best, u, v, cnt, tier);
        } else {
            unlicensed.insert((u, v));
        }
    };

    for ((a, ao, seg, so), ctr) in triple_lookup {
        for (&(ref b, bo), &cnt) in ctr {
            consider_read_edge(&mut best, (a.clone(), *ao), (seg.clone(), *so), cnt, 0);
            consider_read_edge(&mut best, (seg.clone(), *so), (b.clone(), bo), cnt, 0);
        }
    }
    for (state, ctr) in edge_support {
        for (next, &cnt) in ctr {
            consider_read_edge(&mut best, state.clone(), next.clone(), cnt, 1);
        }
    }
    for (state, ctr) in links_index {
        for (next, &cnt) in ctr {
            consider(&mut best, state.clone(), next.clone(), cnt, 2);
        }
    }

    let mut edges: Vec<CandidateEdge> = best
        .into_iter()
        .map(|((from, to), (tier, reads))| CandidateEdge {
            from,
            to,
            tier,
            reads,
        })
        .collect();
    // a fixed order: everything downstream (the solver's variables and
    // constraints, tie-breaks between equally supported edges) follows it,
    // so the same input always resolves the same way
    edges.sort_by(|a, b| (&a.from, &a.to).cmp(&(&b.from, &b.to)));
    (edges, unlicensed.len())
}

/// Solve the integer program: pick the subset of candidate edges maximising
/// segment coverage (weighted towards genes), subject to each segment's
/// depth budget and flow conservation at every state. Given a `root` state,
/// an additional single-commodity-flow constraint forces the result into
/// one circuit through it rather than however many disjoint ones flow
/// conservation alone happens to allow -- see the comment further down for
/// why that's needed. This is the harder problem to solve, and on a few inputs proved
/// too slow for microlp within a reasonable time budget, so the caller
/// retries without it on failure rather than leaving the tool hanging.
///
/// Returns the selected edges and whether the solver proved them optimal;
/// if it ran out of time instead, they are only the best found so far.
fn solve_ilp(
    edges: &[CandidateEdge],
    budget: &HashMap<Seg, u32>,
    segments: &HashMap<Seg, Vec<u8>>,
    gene_segs: &HashSet<Seg>,
    root: Option<&State>,
    time_limit: f64,
) -> Result<(Vec<(State, State)>, bool)> {
    // states and segments in sorted order, so the model is built the same
    // way on every run
    let mut states: BTreeSet<State> = BTreeSet::new();
    for e in edges {
        states.insert(e.from.clone());
        states.insert(e.to.clone());
    }

    let mut vars = ProblemVariables::new();

    let mut x: HashMap<(State, State), Variable> = HashMap::new();
    for e in edges {
        x.insert(
            (e.from.clone(), e.to.clone()),
            vars.add(variable().binary()),
        );
    }
    let mut y: BTreeMap<Seg, Variable> = BTreeMap::new();
    for seg in segments.keys().collect::<BTreeSet<_>>() {
        y.insert(seg.clone(), vars.add(variable().binary()));
    }

    let mut incoming: HashMap<State, Vec<Variable>> = HashMap::new();
    let mut outgoing: HashMap<State, Vec<Variable>> = HashMap::new();
    for e in edges {
        let var = x[&(e.from.clone(), e.to.clone())];
        outgoing.entry(e.from.clone()).or_default().push(var);
        incoming.entry(e.to.clone()).or_default().push(var);
    }

    // Flow conservation alone (added below) only guarantees the selected
    // edges decompose into *some* set of closed circuits -- nothing stops
    // the solver picking several disjoint ones instead of merging them,
    // even when a real, well-supported bridging edge exists (confirmed
    // empirically: every extra circuit turned out to connect directly to
    // the main one with strong read evidence). A connected graph where
    // every vertex has equal in/out degree is guaranteed to have a single
    // Eulerian circuit, so adding a connectivity requirement on top of the
    // existing flow conservation is enough to force one circuit: require
    // every visited state be reachable from a root *state* via the
    // selected edges. This is a
    // standard single-commodity-flow subtour-elimination technique -- a
    // second set of continuous flow variables piggybacks on the same edges
    // (gated by whether x_e is selected) purely to prove reachability.
    //
    // This must be bookkept per *state* (segment, orientation), not per
    // segment: an earlier version aggregated flow across both orientations
    // of a segment, which let a repeat segment used with budget > 1 "pay"
    // for two otherwise-disconnected loops out of the same aggregate
    // consumption -- e.g. segment X used as X+ in one loop and X- in a
    // completely separate loop could each be individually "reachable"
    // because the solver was free to route all the flow through whichever
    // orientation was convenient, without either loop actually connecting
    // to the other. Requiring each *state's* own selected in-degree
    // (its actual usage, which can exceed 1 if budget allows the same
    // orientation to be revisited) to be drawn as flow closes that gap:
    // every individual visit must trace back to the root, not just "some
    // orientation of this segment".
    //
    // The root is forced into the circuit, so it has to be a state the
    // best circuit really passes through: the caller takes it from the
    // solution without this constraint (see `choose_root`). Leaving the
    // choice to the solver (a root variable per state) is exact but, with
    // microlp, too slow to be usable.
    let root_state: Option<&State> = root.filter(|r| states.contains(*r));
    let mut flow: HashMap<(State, State), Variable> = HashMap::new();
    if root_state.is_some() {
        for e in edges {
            flow.insert(
                (e.from.clone(), e.to.clone()),
                vars.add(variable().min(0.0)),
            );
        }
    }

    let total_reads: f64 = edges.iter().map(|e| e.reads as f64).sum::<f64>().max(1.0);

    let mut objective = Expression::from(0.0);
    for (seg, &yvar) in &y {
        let w = if gene_segs.contains(seg) {
            GENE_WEIGHT
        } else {
            BASE_WEIGHT
        };
        objective += w * yvar;
    }
    for e in edges {
        let var = x[&(e.from.clone(), e.to.clone())];
        objective += (e.reads as f64 / total_reads) * 0.01 * var;
    }

    let mut model = vars
        .maximise(objective)
        .using(good_lp::microlp)
        .with_time_limit(time_limit);

    for s in &states {
        let inc: Expression = incoming
            .get(s)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .sum();
        let out: Expression = outgoing
            .get(s)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .sum();
        model.add_constraint((inc - out).eq(0.0));
    }

    for (seg, &yvar) in &y {
        let mut entering: Vec<Variable> = Vec::new();
        if let Some(v) = incoming.get(&(seg.clone(), Orientation::Forward)) {
            entering.extend(v.iter().copied());
        }
        if let Some(v) = incoming.get(&(seg.clone(), Orientation::Backward)) {
            entering.extend(v.iter().copied());
        }
        if entering.is_empty() {
            model.add_constraint((1.0 * yvar).eq(0.0));
        } else {
            let entering_sum: Expression = entering.into_iter().sum();
            model.add_constraint((1.0 * yvar).leq(entering_sum.clone()));
            let b = *budget.get(seg).unwrap_or(&1) as f64;
            model.add_constraint(entering_sum.leq(b));
        }
    }

    if let Some(root_state) = root_state {
        // total connectivity supply can never exceed the sum of every
        // segment's budget (each unit of usage anywhere costs at most one
        // unit of flow along the path that reaches it), so this is a
        // valid, if generous, capacity bound
        let capacity = budget.values().map(|b| *b as f64).sum::<f64>().max(1.0);

        let mut flow_in: HashMap<State, Vec<Variable>> = HashMap::new();
        let mut flow_out: HashMap<State, Vec<Variable>> = HashMap::new();
        for e in edges {
            let key = (e.from.clone(), e.to.clone());
            let fvar = flow[&key];
            model.add_constraint((1.0 * fvar).leq(capacity * x[&key]));
            flow_out.entry(e.from.clone()).or_default().push(fvar);
            flow_in.entry(e.to.clone()).or_default().push(fvar);
        }

        // force the root state to actually be visited, rather than being a
        // phantom source that never appears in the solution
        let root_usage: Expression = incoming
            .get(root_state)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .sum();
        model.add_constraint(root_usage.geq(1.0));

        // every non-root state must draw exactly its own usage (its actual
        // selected in-degree, which can exceed 1 for a revisited state) in
        // connectivity flow, so a state can only be "reachable" by really
        // being connected to the root through selected edges
        let mut others_usage = Expression::from(0.0);
        for st in &states {
            if st == root_state {
                continue;
            }
            if let Some(inc) = incoming.get(st) {
                others_usage += inc.iter().cloned().sum::<Expression>();
            }
        }

        for st in &states {
            let inflow: Expression = flow_in
                .get(st)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .sum();
            let outflow: Expression = flow_out
                .get(st)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .sum();
            if st == root_state {
                // the root supplies exactly enough connectivity flow for
                // every other visited state to draw its own usage from it
                model.add_constraint((inflow - outflow + others_usage.clone()).eq(0.0));
            } else {
                let usage: Expression = incoming
                    .get(st)
                    .cloned()
                    .unwrap_or_default()
                    .into_iter()
                    .sum();
                model.add_constraint((inflow - outflow - usage).eq(0.0));
            }
        }
    }

    let solution = model.solve().context("ILP solve failed")?;
    let optimal = matches!(solution.status(), SolutionStatus::Optimal);

    let selected = edges
        .iter()
        .filter(|e| solution.value(x[&(e.from.clone(), e.to.clone())]) > 0.5)
        .map(|e| (e.from.clone(), e.to.clone()))
        .collect();
    Ok((selected, optimal))
}

/// The state to force a single circuit through: the best-connected state
/// (most candidate edges in and out; the first in sorted order on a tie)
/// of the largest circuit found without the connectivity constraint.
/// Versions up to 0.6.2 took the best-connected state of the whole graph,
/// whether or not a good circuit passes through it - and where the best
/// one doesn't, forcing it in cost segments.
fn choose_root(circuits: &[Vec<State>], edges: &[CandidateEdge]) -> Option<State> {
    let largest = circuits
        .iter()
        .max_by(|a, b| a.len().cmp(&b.len()).then_with(|| b.cmp(a)))?;
    let mut degree: HashMap<&State, usize> = HashMap::new();
    for e in edges {
        *degree.entry(&e.from).or_insert(0) += 1;
        *degree.entry(&e.to).or_insert(0) += 1;
    }
    let mut states: Vec<&State> = largest.iter().collect();
    states.sort();
    states.dedup();
    let mut root: Option<&State> = None;
    for s in states {
        let d = |s: &State| degree.get(s).copied().unwrap_or(0);
        if root.is_none_or(|r| d(s) > d(r)) {
            root = Some(s);
        }
    }
    root.cloned()
}

/// How many reads walk `a -> s -> b` (0 if fewer than [`MIN_EDGE_READS`]).
fn triple_reads(triples: &TripleLookup, a: &State, s: &State, b: &State) -> u32 {
    triples
        .get(&(a.0.clone(), a.1, s.0.clone(), s.1))
        .and_then(|m| m.get(b))
        .copied()
        .unwrap_or(0)
}

/// Pair each of a state's `k` incoming edges with one of its `k` outgoing
/// edges so the pairs' total `support[i][j]` is as large as possible:
/// `result[i]` is the outgoing edge for incoming edge `i`. Exhaustive up to
/// 8 visits (the first best pairing in lexicographic order wins a tie),
/// best-pair-first beyond that.
fn pair_through_state(support: &[Vec<u32>]) -> Vec<usize> {
    let k = support.len();
    if k <= 8 {
        fn search(
            i: usize,
            support: &[Vec<u32>],
            used: &mut Vec<bool>,
            current: &mut Vec<usize>,
            total: u64,
            best: &mut Option<(u64, Vec<usize>)>,
        ) {
            if i == support.len() {
                if best.as_ref().is_none_or(|(b, _)| total > *b) {
                    *best = Some((total, current.clone()));
                }
                return;
            }
            for j in 0..support.len() {
                if !used[j] {
                    used[j] = true;
                    current.push(j);
                    search(i + 1, support, used, current, total + support[i][j] as u64, best);
                    current.pop();
                    used[j] = false;
                }
            }
        }
        let mut best = None;
        search(0, support, &mut vec![false; k], &mut Vec::with_capacity(k), 0, &mut best);
        return best.map(|(_, p)| p).unwrap_or_default();
    }
    let mut pairs: Vec<(u32, usize, usize)> = (0..k)
        .flat_map(|i| (0..k).map(move |j| (i, j)))
        .map(|(i, j)| (support[i][j], i, j))
        .collect();
    pairs.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let mut result = vec![usize::MAX; k];
    let mut used = vec![false; k];
    for (_, i, j) in pairs {
        if result[i] == usize::MAX && !used[j] {
            result[i] = j;
            used[j] = true;
        }
    }
    result
}

/// Split a flow-conserving set of edges into closed circuits, one per
/// connected group of edges. Each returned circuit lists the state every
/// one of its edges leaves from, in order; the closing return-to-start edge
/// is implied. A state the edges pass through more than once (a repeat)
/// appears that many times.
///
/// Which way a circuit leaves a repeat after each entry is decided by the
/// reads: at every state the incoming and outgoing edges are paired to
/// maximise the support of the `(prev, state, next)` triples walked. That
/// can leave several closed walks through the same state; they are joined
/// (by swapping the exits of two entries, one from each walk) where it
/// costs the least triple support, until no two circuits share a state.
/// The result depends only on the edges and the triples, not on the order
/// they are given in.
fn decompose_circuits(
    selected: &[(State, State)],
    triples: &TripleLookup,
) -> Result<Vec<Vec<State>>> {
    let mut edges: Vec<&(State, State)> = selected.iter().collect();
    edges.sort();
    let n = edges.len();

    let mut ins: BTreeMap<&State, Vec<usize>> = BTreeMap::new();
    let mut outs: BTreeMap<&State, Vec<usize>> = BTreeMap::new();
    for (e, (u, v)) in edges.iter().map(|e| (&e.0, &e.1)).enumerate() {
        outs.entry(u).or_default().push(e);
        ins.entry(v).or_default().push(e);
    }
    let support = |e_in: usize, e_out: usize| -> u32 {
        triple_reads(triples, &edges[e_in].0, &edges[e_in].1, &edges[e_out].1)
    };

    // succ[e]: the edge the circuit takes after edge e
    let mut succ = vec![usize::MAX; n];
    let states: BTreeSet<&State> = ins.keys().chain(outs.keys()).copied().collect();
    for state in states {
        let empty = Vec::new();
        let (e_in, e_out) = (ins.get(state).unwrap_or(&empty), outs.get(state).unwrap_or(&empty));
        if e_in.len() != e_out.len() {
            anyhow::bail!(
                "selected edges do not balance at {} ({} in, {} out)",
                state_str(state),
                e_in.len(),
                e_out.len()
            );
        }
        let matrix: Vec<Vec<u32>> = e_in
            .iter()
            .map(|&i| e_out.iter().map(|&o| support(i, o)).collect())
            .collect();
        for (i, j) in pair_through_state(&matrix).into_iter().enumerate() {
            succ[e_in[i]] = e_out[j];
        }
    }

    let label = |succ: &[usize]| -> Vec<usize> {
        let mut cycle = vec![usize::MAX; n];
        let mut next_id = 0;
        for start in 0..n {
            if cycle[start] != usize::MAX {
                continue;
            }
            let mut e = start;
            while cycle[e] == usize::MAX {
                cycle[e] = next_id;
                e = succ[e];
            }
            next_id += 1;
        }
        cycle
    };

    // join closed walks that share a state, cheapest join first
    let mut cycle = label(&succ);
    loop {
        let mut best: Option<(i64, usize, usize)> = None;
        for e_in in ins.values() {
            for (x, &a) in e_in.iter().enumerate() {
                for &b in &e_in[x + 1..] {
                    if cycle[a] == cycle[b] {
                        continue;
                    }
                    let kept = support(a, succ[a]) as i64 + support(b, succ[b]) as i64;
                    let swapped = support(a, succ[b]) as i64 + support(b, succ[a]) as i64;
                    let loss = kept - swapped;
                    if best.is_none_or(|(l, _, _)| loss < l) {
                        best = Some((loss, a, b));
                    }
                }
            }
        }
        let Some((_, a, b)) = best else { break };
        succ.swap(a, b);
        cycle = label(&succ);
    }

    let mut circuits: Vec<Vec<State>> = Vec::new();
    let mut done = vec![false; n];
    for start in 0..n {
        if done[start] {
            continue;
        }
        let mut circuit = Vec::new();
        let mut e = start;
        while !done[e] {
            done[e] = true;
            circuit.push(edges[e].0.clone());
            e = succ[e];
        }
        circuits.push(circuit);
    }
    Ok(circuits)
}

/// One fixed way of writing a circuit, whichever state it was found from
/// and whichever strand it was walked on: of every rotation of the circuit
/// and of its reverse complement, the one that sorts first. So it starts at
/// its first segment by name, forward if that segment is ever walked
/// forward.
fn canonical_circuit(circuit: &[State]) -> Vec<State> {
    let n = circuit.len();
    let rc: Vec<State> = circuit.iter().rev().map(|(s, o)| (s.clone(), flip(*o))).collect();
    let mut best: Option<(&[State], usize)> = None;
    for strand in [circuit, &rc[..]] {
        for k in 0..n {
            let candidate = strand[k..].iter().chain(&strand[..k]);
            let better = match best {
                None => true,
                Some((b, bk)) => candidate.clone().lt(b[bk..].iter().chain(&b[..bk])),
            };
            if better {
                best = Some((strand, k));
            }
        }
    }
    match best {
        Some((strand, k)) => strand[k..].iter().chain(&strand[..k]).cloned().collect(),
        None => Vec::new(),
    }
}

/// For every pair of separately-resolved circuits, the strongest real
/// (non-bare-topology) candidate edge directly connecting them, if any.
///
/// Flow conservation alone only guarantees the ILP's selected edges
/// decompose into *some* set of closed circuits -- nothing in that
/// constraint forces them to be a single circuit, even when the underlying
/// GFA graph is one connected component. This tells you which case you're
/// in for a given pair: a real, read-supported bridge that budget
/// contention forced apart (the split is a modelling limitation), versus no
/// direct evidence linking them at all (the split may reflect genuinely
/// separate structure).
/// Try to splice `satellite` (a whole closed circuit) into `main` (another
/// closed circuit, state-disjoint from it) using real read evidence alone --
/// no solver involved. A satellite is a self-contained loop, so grafting it
/// into `main` only needs *two* new real edges that share a single junction
/// state `r` in the satellite: one from some `p` in `main` into `r` (entry),
/// and one from `r`'s immediate predecessor in the satellite's own cycle
/// back out to `p`'s original successor `q` in `main` (exit). Walking
/// `main` up to `p`, across into the satellite at `r`, all the way around
/// the satellite's existing (already-valid) internal edges back to that
/// predecessor, then out to `q`, reuses every edge either circuit already
/// had except the one `main` edge (p, q) being replaced -- so this can
/// never violate a segment's budget (no state is visited any more or fewer
/// times than before) and always yields a genuinely valid single circuit,
/// unlike trusting a single best-evidence edge in isolation.
fn try_splice(main: &[State], satellite: &[State], edges: &[CandidateEdge]) -> Option<Vec<State>> {
    let mut real_edge: HashSet<(State, State)> = HashSet::new();
    for e in edges {
        if e.tier >= 2 {
            continue; // bare topology only, not real evidence
        }
        real_edge.insert((e.from.clone(), e.to.clone()));
    }

    let n = main.len();
    let m = satellite.len();
    if n == 0 || m == 0 {
        return None;
    }

    for i in 0..n {
        let p = &main[i];
        let q = &main[(i + 1) % n];
        for k in 0..m {
            let r = &satellite[k];
            let r_prev = &satellite[(k + m - 1) % m];
            if real_edge.contains(&(p.clone(), r.clone()))
                && real_edge.contains(&(r_prev.clone(), q.clone()))
            {
                let rotated: Vec<State> = satellite[k..]
                    .iter()
                    .chain(satellite[..k].iter())
                    .cloned()
                    .collect();
                let mut merged = Vec::with_capacity(n + m);
                merged.extend_from_slice(&main[..=i]);
                merged.extend(rotated);
                merged.extend_from_slice(&main[(i + 1)..]);
                return Some(merged);
            }
        }
    }
    None
}

/// Repeatedly try to splice every circuit but the largest into it, using
/// [`try_splice`]. Order matters only for which circuit plays the role of
/// `main` at each step -- always the current largest -- so a satellite that
/// only attaches to another satellite (not yet-merged into main) still gets
/// picked up once that satellite has itself been absorbed. Returns the
/// possibly-reduced set of circuits, largest first, plus how many merges
/// were made.
fn splice_circuits(
    mut circuits: Vec<Vec<State>>,
    edges: &[CandidateEdge],
) -> (Vec<Vec<State>>, usize) {
    let mut merges = 0;
    loop {
        circuits.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));
        if circuits.len() < 2 {
            break;
        }
        let mut spliced_any = false;
        let mut i = 0;
        while i < circuits.len() {
            if i == 0 {
                i += 1;
                continue;
            }
            if let Some(merged) = try_splice(&circuits[0], &circuits[i], edges) {
                circuits[0] = merged;
                circuits.remove(i);
                merges += 1;
                spliced_any = true;
                // restart the sweep against the newly-grown main circuit
                break;
            }
            i += 1;
        }
        if !spliced_any {
            break;
        }
    }
    (circuits, merges)
}

/// Explain precisely why [`try_splice`] could not graft `satellite` into
/// `main`, for stderr reporting. A valid graft needs one real edge into some
/// junction state `r` in the satellite (entry) and one real edge out of
/// `r`'s own immediate predecessor in the satellite's cycle (exit) -- this
/// distinguishes that specific, narrow requirement from the coarser "is
/// there any evidence at all" check, since it's the actual reason a
/// single-direction bridge (however well-supported) isn't enough on its
/// own.
fn explain_unspliceable(main: &[State], satellite: &[State], edges: &[CandidateEdge]) -> String {
    let mut real_edge: HashMap<(State, State), u32> = HashMap::new();
    for e in edges {
        if e.tier >= 2 {
            continue; // bare topology only, not real evidence
        }
        let key = (e.from.clone(), e.to.clone());
        let best = real_edge.entry(key).or_insert(0);
        if e.reads > *best {
            *best = e.reads;
        }
    }

    let m = satellite.len();
    let mut best_entry: Option<(State, State, u32)> = None; // (p, r, reads)
    let mut best_exit: Option<(State, State, u32)> = None; // (r_prev, q, reads)

    for k in 0..m {
        let r = &satellite[k];
        let r_prev = &satellite[(k + m - 1) % m];

        let entry_here = main
            .iter()
            .filter_map(|p| {
                real_edge
                    .get(&(p.clone(), r.clone()))
                    .map(|&reads| (p.clone(), reads))
            })
            .max_by_key(|(_, reads)| *reads);
        let exit_here = main
            .iter()
            .filter_map(|q| {
                real_edge
                    .get(&(r_prev.clone(), q.clone()))
                    .map(|&reads| (q.clone(), reads))
            })
            .max_by_key(|(_, reads)| *reads);

        if let Some((p, reads)) = &entry_here {
            if best_entry.as_ref().is_none_or(|(_, _, r0)| reads > r0) {
                best_entry = Some((p.clone(), r.clone(), *reads));
            }
        }
        if let Some((q, reads)) = &exit_here {
            if best_exit.as_ref().is_none_or(|(_, _, r0)| reads > r0) {
                best_exit = Some((r_prev.clone(), q.clone(), *reads));
            }
        }
    }

    match (&best_entry, &best_exit) {
        (None, None) => "no real (non-bare-topology) read evidence connects it to the main \
             circuit in either direction"
            .to_string(),
        (Some((p, r, reads)), None) => format!(
            "real entry evidence only ({} -> {}, ~{reads} reads): the satellite's own state \
             immediately before {} has no real edge back out to the main circuit, so grafting \
             would mean either skipping part of the satellite's own loop or revisiting a \
             segment -- neither is safe",
            state_str(p),
            state_str(r),
            state_str(r)
        ),
        (None, Some((r_prev, q, reads))) => format!(
            "real exit evidence only ({} -> {}, ~{reads} reads): nothing in the main circuit \
             has a real edge into the satellite state immediately before {}, so there's no \
             safe entry point to pair with it",
            state_str(r_prev),
            state_str(q),
            state_str(r_prev)
        ),
        (Some((p, r, er)), Some((r_prev, q, xr))) => format!(
            "entry and exit evidence exist but through different junctions ({} -> {}, ~{er} \
             reads entering vs {} -> {}, ~{xr} reads leaving): grafting through mismatched \
             junctions would skip or double-visit part of the satellite's own loop, so this \
             can't be reconciled into one valid circuit without more/different read evidence",
            state_str(p),
            state_str(r),
            state_str(r_prev),
            state_str(q)
        ),
    }
}

fn find_inter_circuit_links(
    circuits: &[Vec<State>],
    edges: &[CandidateEdge],
) -> Vec<(usize, usize, Option<(State, State, u8, u32)>)> {
    let circuit_state_sets: Vec<HashSet<State>> = circuits
        .iter()
        .map(|c| c.iter().cloned().collect())
        .collect();

    let mut results = Vec::new();
    for i in 0..circuits.len() {
        for j in (i + 1)..circuits.len() {
            let mut best: Option<(State, State, u8, u32)> = None;
            for e in edges {
                if e.tier >= 2 {
                    continue; // bare topology only -- every segment is graph-connected by definition, not informative here
                }
                let forward = circuit_state_sets[i].contains(&e.from)
                    && circuit_state_sets[j].contains(&e.to);
                let backward = circuit_state_sets[j].contains(&e.from)
                    && circuit_state_sets[i].contains(&e.to);
                if forward || backward {
                    let better = match &best {
                        None => true,
                        Some((_, _, t, c)) => {
                            (e.tier, std::cmp::Reverse(e.reads)) < (*t, std::cmp::Reverse(*c))
                        }
                    };
                    if better {
                        best = Some((e.from.clone(), e.to.clone(), e.tier, e.reads));
                    }
                }
            }
            results.push((i, j, best));
        }
    }
    results
}

/// A segment excluded from every resolved circuit, but with real read
/// evidence connecting it to both sides of the circuit -- a recombination
/// bubble's minor arm, not a data gap.
struct BubbleArm {
    segment: Seg,
    /// The strand the arm is walked on between `entry` and `exit`.
    orientation: Orientation,
    entry: State,
    exit: State,
    support_reads: u32,
}

fn detect_bubbles(
    edges: &[CandidateEdge],
    excluded: &HashSet<Seg>,
    circuit_states: &HashSet<State>,
) -> Vec<BubbleArm> {
    // Keyed by the arm's *state*: an entry into x+ and an exit from x- are
    // not one walk through x, so the two must agree on the strand.
    let mut into: HashMap<State, (State, State, u32)> = HashMap::new();
    let mut out_of: HashMap<State, (State, State, u32)> = HashMap::new();

    for e in edges {
        if e.tier >= 2 {
            continue; // bare topology only, not real read evidence
        }
        if excluded.contains(&e.to.0) && circuit_states.contains(&e.from) {
            let entry = into
                .entry(e.to.clone())
                .or_insert((e.from.clone(), e.to.clone(), 0));
            if e.reads > entry.2 {
                *entry = (e.from.clone(), e.to.clone(), e.reads);
            }
        }
        if excluded.contains(&e.from.0) && circuit_states.contains(&e.to) {
            let entry = out_of
                .entry(e.from.clone())
                .or_insert((e.from.clone(), e.to.clone(), 0));
            if e.reads > entry.2 {
                *entry = (e.from.clone(), e.to.clone(), e.reads);
            }
        }
    }

    let mut bubbles = Vec::new();
    let mut excluded: Vec<&Seg> = excluded.iter().collect();
    excluded.sort();
    for seg in excluded {
        // the better-supported strand; forward on a tie (the other strand
        // is the same arm seen from the circuit's reverse complement)
        let mut best: Option<BubbleArm> = None;
        for orientation in [Orientation::Forward, Orientation::Backward] {
            let state = (seg.clone(), orientation);
            if let (Some(in_edge), Some(out_edge)) = (into.get(&state), out_of.get(&state)) {
                let support_reads = in_edge.2.min(out_edge.2);
                if best.as_ref().is_none_or(|b| support_reads > b.support_reads) {
                    best = Some(BubbleArm {
                        segment: seg.clone(),
                        orientation,
                        entry: in_edge.0.clone(),
                        exit: out_edge.1.clone(),
                        support_reads,
                    });
                }
            }
        }
        bubbles.extend(best);
    }
    bubbles
}

fn parse_gff_gene_segments(path: &PathBuf) -> Result<HashMap<Seg, Vec<String>>> {
    let file = File::open(path).with_context(|| format!("Failed to open GFF file: {path:?}"))?;
    let reader = BufReader::new(file);
    let mut out: HashMap<Seg, Vec<String>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[2] != "gene" {
            continue;
        }
        let mut target_name = None;
        let mut coverage = None;
        for kv in fields[8].split(';') {
            if let Some((k, v)) = kv.split_once('=') {
                match k {
                    "target_name" => target_name = Some(v.to_string()),
                    "coverage" => coverage = Some(v.to_string()),
                    _ => {}
                }
            }
        }
        if coverage.as_deref() == Some("full") {
            if let Some(name) = target_name {
                out.entry(fields[0].as_bytes().to_vec())
                    .or_default()
                    .push(name);
            }
        }
    }
    Ok(out)
}

/// Parse an oatk-style BED annotation directly -- `seq_name align_from
/// align_to gene_name score_capped_at_1000 strand`, tab-separated, `seq_name`
/// already referencing the GFA's own segment IDs -- so no GFF conversion
/// step is needed between oatk's own gene calls and `gfatk resolve`. Every
/// row counts as a gene call regardless of score: unlike GFF's
/// `coverage=full` flag, oatk's score isn't a clean complete/partial
/// signal, so filtering on it risks silently dropping real single-copy
/// genes that just don't happen to reach the cap.
fn parse_bed_gene_segments(path: &PathBuf) -> Result<HashMap<Seg, Vec<String>>> {
    let file = File::open(path).with_context(|| format!("Failed to open BED file: {path:?}"))?;
    let reader = BufReader::new(file);
    let mut out: HashMap<Seg, Vec<String>> = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        let seg = fields[0].as_bytes().to_vec();
        let gene_name = fields[3];
        out.entry(seg).or_default().push(gene_name.to_string());
    }
    Ok(out)
}

fn state_str(s: &State) -> String {
    format!("{}{}", String::from_utf8_lossy(&s.0), orient_char(s.1))
}

fn write_fasta_record(header: &str, seq: &[u8]) {
    println!(">{header}");
    for chunk in seq.chunks(70) {
        println!("{}", String::from_utf8_lossy(chunk));
    }
}

/// The overlap, in bases, a GFA link's CIGAR declares between its two
/// segments. `*` (overlap not given) counts as no overlap.
fn link_overlap(cigar: &[u8]) -> Result<usize> {
    if cigar == b"*" {
        return Ok(0);
    }
    utils::parse_cigar(cigar)
}

/// Overlap length for every traversal the GFA's links license: each link
/// as written, and its reverse-complement equivalent (the same pairs
/// [`build_links_index`] registers).
fn build_overlaps(
    links: &[(Seg, Orientation, Seg, Orientation)],
    cigars: &[Vec<u8>],
) -> Result<HashMap<(State, State), usize>> {
    let mut overlaps = HashMap::new();
    for ((fs, fo, ts, to), cigar) in links.iter().zip(cigars) {
        let overlap = link_overlap(cigar).with_context(|| {
            format!(
                "link {}{} -> {}{}",
                String::from_utf8_lossy(fs),
                orient_char(*fo),
                String::from_utf8_lossy(ts),
                orient_char(*to)
            )
        })?;
        overlaps.insert(((fs.clone(), *fo), (ts.clone(), *to)), overlap);
        overlaps.insert(((ts.clone(), flip(*to)), (fs.clone(), flip(*fo))), overlap);
    }
    Ok(overlaps)
}

/// The sequence a walk through the graph spells out. Two linked segments
/// share the bases of their link's overlap, so each segment after the first
/// is written without its first `overlap` bases -- as `gfatk path` and
/// `gfatk linear` do.
///
/// A circuit (`circular`) also has a closing junction, from its last
/// segment back to its first. That overlap is trimmed off the start of the
/// first segment, so the record begins `overlap` bases into it. This is
/// where oatk puts the origin of a circular contig, so the same `path=`
/// maps to the same coordinates in either tool's output.
fn path_sequence(
    path: &[State],
    segments: &HashMap<Seg, Vec<u8>>,
    overlaps: &HashMap<(State, State), usize>,
    circular: bool,
) -> Result<Vec<u8>> {
    let mut seq = Vec::new();
    for (i, state) in path.iter().enumerate() {
        let prev = match i {
            0 if circular => path.last(),
            0 => None,
            _ => Some(&path[i - 1]),
        };
        let trim = match prev {
            Some(prev) => *overlaps
                .get(&(prev.clone(), state.clone()))
                .with_context(|| {
                    format!(
                        "no GFA link for the transition {} -> {}",
                        state_str(prev),
                        state_str(state)
                    )
                })?,
            None => 0,
        };
        let (seg, o) = state;
        let s = &segments[seg];
        anyhow::ensure!(
            trim <= s.len(),
            "the link {} -> {} overlaps by {trim} bases, more than the {} in {}",
            state_str(prev.unwrap()),
            state_str(state),
            s.len(),
            String::from_utf8_lossy(seg)
        );
        match o {
            Orientation::Forward => seq.extend_from_slice(&s[trim..]),
            Orientation::Backward => seq.extend_from_slice(&revcomp(s)[trim..]),
        }
    }
    Ok(seq)
}

/// A segment's sequence as read in the given orientation.
fn oriented(state: &State, segments: &HashMap<Seg, Vec<u8>>) -> Vec<u8> {
    match state.1 {
        Orientation::Forward => segments[&state.0].clone(),
        Orientation::Backward => revcomp(&segments[&state.0]),
    }
}

/// Junctions of a circuit where the two segments do not read the same over
/// their link's overlap, as `(from, to, overlap, differing bases)`. A GFA
/// overlap of `41M` promises 41 aligned bases, not 41 identical ones, and
/// [`path_sequence`] keeps the earlier segment's copy of them.
fn inexact_junctions(
    circuit: &[State],
    segments: &HashMap<Seg, Vec<u8>>,
    overlaps: &HashMap<(State, State), usize>,
) -> Vec<(State, State, usize, usize)> {
    let mut out = Vec::new();
    for (i, to) in circuit.iter().enumerate() {
        let from = &circuit[(i + circuit.len() - 1) % circuit.len()];
        let Some(&overlap) = overlaps.get(&(from.clone(), to.clone())) else {
            continue;
        };
        let (a, b) = (oriented(from, segments), oriented(to, segments));
        if overlap == 0 || overlap > a.len() || overlap > b.len() {
            continue;
        }
        let differing = a[a.len() - overlap..]
            .iter()
            .zip(&b[..overlap])
            .filter(|(x, y)| !x.eq_ignore_ascii_case(y))
            .count();
        if differing > 0 {
            out.push((from.clone(), to.clone(), overlap, differing));
        }
    }
    out
}

/// Levenshtein edit distance, space-optimised to two rolling rows. Guards
/// against a pathologically large comparison (bubble segments are normally
/// small, but a main-path detour spanning several segments could in
/// principle be large) by refusing to run past `max_cells` total DP cells,
/// returning `None` rather than eating memory/time on an unbounded input.
fn levenshtein(a: &[u8], b: &[u8], max_cells: usize) -> Option<usize> {
    if a.len().saturating_mul(b.len()) > max_cells {
        return None;
    }
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0usize; b.len() + 1];
    for (i, &ac) in a.iter().enumerate() {
        cur[0] = i + 1;
        for (j, &bc) in b.iter().enumerate() {
            let cost = if ac == bc { 0 } else { 1 };
            cur[j + 1] = (prev[j + 1] + 1).min(cur[j] + 1).min(prev[j] + cost);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    Some(prev[b.len()])
}

// No real hop cap: a repeat segment's nearest recurrence of the matching
// exit state is often tens to hundreds of positions away in a tangled,
// repeat-heavy tour purely because of how many times other repeats get
// revisited in between -- that's not the same as the bubble being
// structurally distant. `find_main_detour` already bounds this to the
// circuit's own length, so usize::MAX just means "search the whole
// circuit". The real cost control is `LEVENSHTEIN_MAX_CELLS` below, since
// that's what actually scales with how much sequence a long detour pulls
// in.
const MAX_DETOUR_HOPS: usize = usize::MAX;
// A single DP cell is trivial work, so this comfortably covers even a
// small bubble compared against a detour spanning most of a several-
// hundred-kb circuit while still refusing a truly unbounded comparison.
const LEVENSHTEIN_MAX_CELLS: usize = 2_000_000_000;
// A genuine local (two-allele-style) bubble's detour is roughly the same
// size as the bubble arm itself. When the nearest matching detour is many
// times larger, it's evidence of something else going on structurally
// (e.g. the same evidence that leaves circuits un-spliceable), and an edit
// distance dominated by that length gap reads as "very different" no
// matter what the actual sequence looks like -- not a useful signal.
const MAX_DETOUR_SIZE_RATIO: usize = 5;

/// How a bubble arm's sequence compares to whatever the main path actually
/// does between the same entry/exit points -- purely informational, never
/// affects what gets included or written.
enum BubbleDiff {
    /// A genuinely local, size-comparable detour was found and diffed.
    Compared {
        edit_distance: usize,
        pct_identity: f64,
        compared_bp: usize,
    },
    /// A detour exists, but it's many times larger than the bubble arm
    /// itself -- not the same local two-allele structure this metric is
    /// meant for, so diffing it would just report the length gap dressed
    /// up as sequence divergence.
    DetourTooLarge { detour_bp: usize, ratio: usize },
    /// A size-comparable detour was found, but the comparison itself was
    /// too large to run within `LEVENSHTEIN_MAX_CELLS`.
    TooLargeToDiff,
    /// No occurrence of the bubble's entry state reaches its exit state
    /// within the resolved circuit(s) at all (e.g. the entry/exit states
    /// sit in different, un-merged circuits).
    NoDetourFound,
}

impl BubbleDiff {
    /// (stderr summary suffix, FASTA header suffix)
    fn describe(&self) -> (String, String) {
        match self {
            BubbleDiff::Compared {
                edit_distance,
                pct_identity,
                compared_bp,
            } => (
                format!(
                    ", vs main path: {edit_distance}bp edit distance, {pct_identity:.1}% \
                     identity over {compared_bp}bp"
                ),
                format!(":edit_distance={edit_distance}:pct_identity={pct_identity:.1}"),
            ),
            BubbleDiff::DetourTooLarge { detour_bp, ratio } => (
                format!(
                    ", vs main path: nearest matching detour is {detour_bp}bp ({ratio}x this \
                     bubble's own length) -- too large relative to the bubble to be a \
                     meaningful local comparison, not diffed"
                ),
                format!(":main_path_comparison=detour_too_large:detour_bp={detour_bp}"),
            ),
            BubbleDiff::TooLargeToDiff => (
                ", vs main path: too large to diff, skipped".to_string(),
                ":main_path_comparison=too_large".to_string(),
            ),
            BubbleDiff::NoDetourFound => (
                ", vs main path: no comparable detour found within local structure".to_string(),
                ":main_path_comparison=unavailable".to_string(),
            ),
        }
    }
}

fn diff_bubble_arm(bubble_seq: &[u8], detour_seq: &[u8]) -> BubbleDiff {
    let bubble_len = bubble_seq.len().max(1);
    if detour_seq.len() > bubble_len * MAX_DETOUR_SIZE_RATIO {
        return BubbleDiff::DetourTooLarge {
            detour_bp: detour_seq.len(),
            ratio: detour_seq.len() / bubble_len,
        };
    }
    match levenshtein(bubble_seq, detour_seq, LEVENSHTEIN_MAX_CELLS) {
        Some(dist) => {
            let longer = bubble_seq.len().max(detour_seq.len()).max(1);
            let pct_identity = 100.0 * (1.0 - dist as f64 / longer as f64);
            BubbleDiff::Compared {
                edit_distance: dist,
                pct_identity,
                compared_bp: longer,
            }
        }
        None => BubbleDiff::TooLargeToDiff,
    }
}

/// Find the shortest run of states a resolved circuit actually uses between
/// `entry` and `exit` -- i.e. what the main path does instead of a bubble
/// arm at that same graph position -- searching every occurrence of `entry`
/// across all circuits (a repeat segment can appear more than once) up to
/// `max_hops` forward steps. Returns the states strictly between `entry`
/// and `exit` (empty if `exit` immediately follows `entry`: the main path
/// has no sequence there at all, i.e. the bubble is a pure insertion
/// relative to it). `None` means no occurrence of `entry` reaches `exit`
/// within the hop bound -- the local structure is too tangled (usually
/// repeat-heavy) to isolate a single comparable detour.
fn find_main_detour(
    circuits: &[Vec<State>],
    entry: &State,
    exit: &State,
    max_hops: usize,
) -> Option<Vec<State>> {
    let mut best: Option<Vec<State>> = None;
    for circuit in circuits {
        let n = circuit.len();
        if n == 0 {
            continue;
        }
        for (i, s) in circuit.iter().enumerate() {
            if s != entry {
                continue;
            }
            for hop in 1..=max_hops.min(n) {
                let j = (i + hop) % n;
                if &circuit[j] == exit {
                    let detour: Vec<State> =
                        (1..hop).map(|h| circuit[(i + h) % n].clone()).collect();
                    if best.as_ref().is_none_or(|b| detour.len() < b.len()) {
                        best = Some(detour);
                    }
                    break;
                }
            }
        }
    }
    best
}

/// The main entry point for `gfatk resolve`.
pub fn resolve(matches: &clap::ArgMatches) -> Result<()> {
    let gfa_file = matches.get_one::<PathBuf>("GFA");
    let gaf_file = matches
        .get_one::<PathBuf>("gaf")
        .context("`--gaf` is required for `gfatk resolve`")?;
    let gff_file = matches.get_one::<PathBuf>("gff");
    let bed_file = matches.get_one::<PathBuf>("bed");
    let min_mapq = *matches
        .get_one::<u32>("min-mapq")
        .expect("defaulted by clap");
    let min_identity = *matches
        .get_one::<f64>("min-identity")
        .expect("defaulted by clap");
    let min_read_aligned = *matches
        .get_one::<f64>("min-read-aligned")
        .expect("defaulted by clap");
    let time_limit = *matches
        .get_one::<f64>("time-limit")
        .expect("defaulted by clap");
    let no_bubble_fasta = matches.get_flag("no-bubble-fasta");

    let gfa: GFAtk = match gfa_file {
        Some(f) => GFAtk(load_gfa(f)?),
        None => match utils::is_stdin() {
            true => GFAtk(load_gfa_stdin(std::io::stdin().lock())?),
            false => anyhow::bail!("No input from STDIN. Run `gfatk resolve -h` for help."),
        },
    };

    let mut segments: HashMap<Seg, Vec<u8>> = HashMap::new();
    for seg in &gfa.0.segments {
        segments.insert(seg.name.clone(), seg.sequence.clone());
    }
    let seg_len: HashMap<Seg, usize> = segments.iter().map(|(k, v)| (k.clone(), v.len())).collect();

    let links: Vec<(Seg, Orientation, Seg, Orientation)> = gfa
        .0
        .links
        .iter()
        .map(|l| {
            (
                l.from_segment.clone(),
                l.from_orient,
                l.to_segment.clone(),
                l.to_orient,
            )
        })
        .collect();

    eprintln!("[+]\tReading GAF alignments from {gaf_file:?}");
    let evidence = read_gaf_evidence(gaf_file, min_mapq, min_identity, min_read_aligned)?;
    let records = evidence.records;
    eprintln!(
        "[+]\t{} usable read-path records after quality filtering",
        records.len()
    );
    eprintln!(
        "[+]\t{} of {} reads are aligned over at least {:.0}% of their length and count for depth",
        evidence.n_reads_for_depth,
        evidence.n_reads,
        100.0 * min_read_aligned
    );


    let triple_lookup = build_triple_lookup(&records);
    let edge_support = build_edge_support(&records);
    let links_index = build_links_index(&links);
    let cigars: Vec<Vec<u8>> = gfa.0.links.iter().map(|l| l.overlap.clone()).collect();
    let overlaps = build_overlaps(&links, &cigars)?;
    let depth = compute_depth(&records, &evidence.single_segment_bases, &seg_len, &overlaps);
    let (baseline, budget) = estimate_budget(&depth, &seg_len);
    eprintln!("[+]\tBaseline (single-copy) read depth: {baseline:.1}x");
    let (edges, n_unlicensed) = build_candidate_edges(&triple_lookup, &edge_support, &links_index);
    eprintln!(
        "[+]\t{} candidate edges from read/graph evidence",
        edges.len()
    );
    if n_unlicensed > 0 {
        eprintln!(
            "[-]\t{n_unlicensed} read-derived transition(s) not licensed by any GFA link were discarded -- was the GAF aligned to this exact GFA?"
        );
    }

    let mut gene_segs: HashSet<Seg> = HashSet::new();
    if let Some(path) = gff_file {
        let genes = parse_gff_gene_segments(path)?;
        eprintln!(
            "[+]\t{} segments carry a high-confidence annotated gene (from {path:?})",
            genes.len()
        );
        gene_segs.extend(genes.into_keys());
    }
    if let Some(path) = bed_file {
        let genes = parse_bed_gene_segments(path)?;
        eprintln!(
            "[+]\t{} segments carry a high-confidence annotated gene (from {path:?})",
            genes.len()
        );
        gene_segs.extend(genes.into_keys());
    }

    // First without the single-circuit constraint: quick, and it shows
    // which part of the graph the best solution uses, so the constrained
    // solve can be rooted there.
    let (free_selected, free_optimal) = solve_ilp(&edges, &budget, &segments, &gene_segs, None, time_limit)?;
    let free_circuits = decompose_circuits(&free_selected, &triple_lookup)?;
    let root = choose_root(&free_circuits, &edges);
    let (selected, optimal) = if free_circuits.len() <= 1 {
        // already one circuit: the constraint has nothing to add
        (free_selected, free_optimal)
    } else {
        match solve_ilp(&edges, &budget, &segments, &gene_segs, root.as_ref(), time_limit) {
            Ok(solved) => solved,
            Err(e) => {
                eprintln!(
                    "[-]\tCould not solve for a single connected circuit within the time budget ({e}); \
                     falling back to the unconstrained solve, which may report multiple circuits."
                );
                (free_selected, free_optimal)
            }
        }
    };
    if !optimal {
        eprintln!(
            "[-]\tThe solver stopped at its time limit ({time_limit} s) before proving its answer \
             the best one: the circuit(s) below are the best found so far, and may differ \
             between runs. Raise --time-limit to let it finish."
        );
    }
    let circuits = {
        let c = decompose_circuits(&selected, &triple_lookup)?;
        let n_before = c.len();
        let (spliced, merges) = splice_circuits(c, &edges);
        if merges > 0 {
            eprintln!(
                "[+]\t{merges} satellite circuit(s) grafted into a larger one using real \
                 entry/exit read evidence ({n_before} circuit(s) -> {}), with no solver \
                 involved -- this can never be wrong when it succeeds, only unable to find \
                 a valid splice point.",
                spliced.len()
            );
        }
        let mut canonical: Vec<Vec<State>> = spliced.iter().map(|c| canonical_circuit(c)).collect();
        canonical.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a.cmp(b)));
        canonical
    };

    let included: HashSet<Seg> = circuits.iter().flatten().map(|(s, _)| s.clone()).collect();
    let circuit_states: HashSet<State> = circuits.iter().flatten().cloned().collect();
    let excluded: HashSet<Seg> = segments
        .keys()
        .filter(|s| !included.contains(*s))
        .cloned()
        .collect();

    eprintln!(
        "[+]\tResolved {} circuit(s), {} of {} segments included",
        circuits.len(),
        included.len(),
        segments.len()
    );

    for (idx, circuit) in circuits.iter().enumerate() {
        let seq = path_sequence(circuit, &segments, &overlaps, true)?;
        let path_str = circuit.iter().map(state_str).collect::<Vec<_>>().join(",");
        eprintln!(
            "[+]\tCircuit {}: {} segments, {} bp, path={}",
            idx + 1,
            circuit.len(),
            seq.len(),
            path_str
        );
        let inexact = inexact_junctions(circuit, &segments, &overlaps);
        if !inexact.is_empty() {
            let shown: Vec<String> = inexact
                .iter()
                .take(5)
                .map(|(from, to, overlap, differing)| {
                    format!(
                        "{} -> {} ({differing} of {overlap} bp differ)",
                        state_str(from),
                        state_str(to)
                    )
                })
                .collect();
            eprintln!(
                "[-]\tCircuit {}: at {} junction(s) the two segments differ within their link's \
                 overlap; the earlier segment's bases are kept there: {}{}",
                idx + 1,
                inexact.len(),
                shown.join(", "),
                if inexact.len() > shown.len() { ", ..." } else { "" }
            );
        }
        let species_tag = gfa_file
            .and_then(|p| p.file_stem())
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "gfatk_resolve".to_string());
        let header = format!(
            "{species_tag}_circuit{}:n_segments={}:bp={}:circular=true:path={}",
            idx + 1,
            circuit.len(),
            seq.len(),
            path_str
        );
        write_fasta_record(&header, &seq);
    }

    if circuits.len() > 1 {
        eprintln!(
            "[-]\t{} separate circuits remain after automatic grafting -- flow conservation alone \
             guarantees a valid decomposition into closed loops, not that they merge into one, and \
             the read evidence available did not support a safe graft for the rest. Checking for \
             direct evidence between each pair:",
            circuits.len()
        );
        for (i, j, link) in find_inter_circuit_links(&circuits, &edges) {
            match link {
                Some((from, to, tier, reads)) => {
                    let tier_name = if tier == 0 {
                        "exact read-path"
                    } else {
                        "pairwise read"
                    };
                    eprintln!(
                        "[-]\t  Circuit {} <-> Circuit {}: real {} evidence between {} and {} (~{} reads).",
                        i + 1,
                        j + 1,
                        tier_name,
                        state_str(&from),
                        state_str(&to),
                        reads
                    );
                }
                None => {
                    eprintln!(
                        "[-]\t  Circuit {} <-> Circuit {}: no direct read evidence found connecting them \
                         (only via bare graph topology, if at all) -- consistent with genuinely separate structure.",
                        i + 1,
                        j + 1
                    );
                }
            }
        }
        eprintln!("[-]\tWhy the automatic graft (2 real edges through a matching junction) couldn't absorb these into Circuit 1:");
        for (idx, satellite) in circuits.iter().enumerate().skip(1) {
            eprintln!(
                "[-]\t  Circuit {}: {}",
                idx + 1,
                explain_unspliceable(&circuits[0], satellite, &edges)
            );
        }
    }

    let bubbles = detect_bubbles(&edges, &excluded, &circuit_states);
    if !bubbles.is_empty() {
        let fasta_note = if no_bubble_fasta {
            "not written to the FASTA output (--no-bubble-fasta), but not treated as missing either"
        } else {
            "reported as separate records below, not treated as missing"
        };
        eprintln!(
            "[-]\t{} segment(s) excluded from the resolved circuit(s) are recombination-bubble arms \
             (real read evidence connects them to the main circuit on both sides, but the dominant \
             conformation won out) -- {fasta_note}.",
            bubbles.len()
        );
        for b in &bubbles {
            let genes = gene_segs.contains(&b.segment);
            // the arm as it is walked from entry to exit
            let seq = &oriented(&(b.segment.clone(), b.orientation), &segments);

            // How different is this bubble arm from whatever the main
            // path actually does between the same entry/exit points?
            // Purely informational -- doesn't affect what gets written.
            let diff = match find_main_detour(&circuits, &b.entry, &b.exit, MAX_DETOUR_HOPS) {
                Some(detour) => diff_bubble_arm(seq, &path_sequence(&detour, &segments, &overlaps, false)?),
                None => BubbleDiff::NoDetourFound,
            };
            let (diff_summary, diff_header) = diff.describe();

            eprintln!(
                "[-]\t  {} (between {} and {}, ~{} supporting reads{}{})",
                state_str(&(b.segment.clone(), b.orientation)),
                state_str(&b.entry),
                state_str(&b.exit),
                b.support_reads,
                if genes {
                    ", carries an annotated gene"
                } else {
                    ""
                },
                diff_summary
            );
            if !no_bubble_fasta {
                // "_to_" rather than a bare "-" separator: state_str()
                // already ends in '+' or '-' for the orientation, so
                // "u34--u22-" reads ambiguously -- "u34-_to_u22-" doesn't.
                let header = format!(
                    "bubble_arm:{}:between={}_to_{}:support_reads={}:strand={}{}",
                    String::from_utf8_lossy(&b.segment),
                    state_str(&b.entry),
                    state_str(&b.exit),
                    b.support_reads,
                    orient_char(b.orientation),
                    diff_header
                );
                write_fasta_record(&header, seq);
            }
        }
    }

    let mut truly_missing: Vec<&Seg> = excluded
        .iter()
        .filter(|s| !bubbles.iter().any(|b| &b.segment == *s))
        .collect();
    if !truly_missing.is_empty() {
        truly_missing.sort();
        eprintln!(
            "[-]\t{} segment(s) excluded with no bubble evidence either (genuinely unresolved, \
             not just an alternate conformation) -- written out below as separate, unplaced \
             records rather than dropped, the way an assembler would report an unplaced contig.",
            truly_missing.len()
        );
        for seg in &truly_missing {
            let seq = &segments[*seg];
            eprintln!("[-]\t  {} ({}bp)", String::from_utf8_lossy(seg), seq.len());
            let header = format!(
                "unresolved_segment:{}:bp={}",
                String::from_utf8_lossy(seg),
                seq.len()
            );
            write_fasta_record(&header, seq);
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_bed_gene_segments_counts_every_row_regardless_of_score() {
        let bed = "#seq_name\talign_from\talign_to\tgene_name\tscore_capped_at_1000\tstrand\n\
                    u2\t6877\t7492\tccmFc\t755\t-\n\
                    u5\t1\t939\tatp1\t1000\t+\n\
                    u5\t36404\t36477\trps10\t95\t-\n";
        let path = std::env::temp_dir().join("gfatk_test_parse_bed_gene_segments_no_filter.bed");
        std::fs::write(&path, bed).unwrap();
        let genes = parse_bed_gene_segments(&path).unwrap();
        std::fs::remove_file(&path).ok();

        assert_eq!(genes[&b"u2".to_vec()], vec!["ccmFc".to_string()]);
        assert_eq!(
            genes[&b"u5".to_vec()],
            vec!["atp1".to_string(), "rps10".to_string()]
        );
    }

    #[test]
    fn test_flip() {
        assert_eq!(flip(Orientation::Forward), Orientation::Backward);
        assert_eq!(flip(Orientation::Backward), Orientation::Forward);
    }

    #[test]
    fn test_revcomp() {
        assert_eq!(revcomp(b"ACGT"), b"ACGT");
        assert_eq!(revcomp(b"AACCGGTT"), b"AACCGGTT");
        assert_eq!(revcomp(b"GATTACA"), b"TGTAATC");
    }

    #[test]
    fn test_estimate_budget_single_copy_baseline() {
        let mut depth = HashMap::new();
        depth.insert(b"a".to_vec(), 100.0);
        depth.insert(b"b".to_vec(), 100.0);
        depth.insert(b"c".to_vec(), 100.0);
        depth.insert(b"repeat".to_vec(), 300.0);
        let len: HashMap<Seg, usize> = depth.keys().map(|k| (k.clone(), 1000)).collect();
        let (baseline, budget) = estimate_budget(&depth, &len);
        assert!((baseline - 100.0).abs() < 1e-9);
        assert_eq!(budget[&b"a".to_vec()], 1);
        assert_eq!(budget[&b"repeat".to_vec()], 3);
    }

    #[test]
    fn test_estimate_budget_floors_at_one() {
        // a segment with very little read support should still get a
        // budget of at least 1, never 0 -- it's real, just underrepresented
        let mut depth = HashMap::new();
        depth.insert(b"a".to_vec(), 100.0);
        depth.insert(b"b".to_vec(), 100.0);
        depth.insert(b"weak".to_vec(), 5.0);
        let len: HashMap<Seg, usize> = depth.keys().map(|k| (k.clone(), 1000)).collect();
        let (_, budget) = estimate_budget(&depth, &len);
        assert_eq!(budget[&b"weak".to_vec()], 1);
    }

    #[test]
    fn test_estimate_budget_baseline_is_weighted_by_length() {
        // one long single-copy segment and four short two-copy repeats:
        // most segments are at 200x, most bases at 100x
        let mut depth = HashMap::new();
        let mut len = HashMap::new();
        depth.insert(b"unique".to_vec(), 100.0);
        len.insert(b"unique".to_vec(), 100_000);
        for name in ["r1", "r2", "r3", "r4"] {
            depth.insert(name.as_bytes().to_vec(), 200.0);
            len.insert(name.as_bytes().to_vec(), 2_000);
        }
        let (baseline, budget) = estimate_budget(&depth, &len);
        assert!((baseline - 100.0).abs() < 1e-9);
        assert_eq!(budget[&b"unique".to_vec()], 1);
        assert_eq!(budget[&b"r1".to_vec()], 2);
    }

    #[test]
    fn test_compute_depth_uses_trimmed_path_coordinates_and_single_segment_reads() {
        use Orientation::Forward as F;
        // a (10 bp) and b (10 bp) overlap by 4: the path a,b is 16 long and
        // b covers path positions 6..16
        let seg_len: HashMap<Seg, usize> = [(b"a".to_vec(), 10), (b"b".to_vec(), 10)].into_iter().collect();
        let mut overlaps = HashMap::new();
        overlaps.insert((s("a", F), s("b", F)), 4usize);
        let spanning = GafRecord { path: vec![s("a", F), s("b", F)], path_start: 8, path_end: 16, mapq: 60, identity: None, counts_for_depth: true };
        let depth = compute_depth(std::slice::from_ref(&spanning), &HashMap::new(), &seg_len, &overlaps);
        // a: 8..10 = 2 bases; b: 8..16 = 8 bases (laid end to end, b would
        // start at 10 and get only 6)
        assert_eq!(depth[&b"a".to_vec()], 0.2);
        assert_eq!(depth[&b"b".to_vec()], 0.8);
        // a read inside b alone (5 aligned bases) adds to b's depth
        let within: HashMap<Seg, u64> = [(b"b".to_vec(), 5)].into_iter().collect();
        let depth = compute_depth(std::slice::from_ref(&spanning), &within, &seg_len, &overlaps);
        assert_eq!(depth[&b"b".to_vec()], 1.3);
        // a junction-crossing read that is mostly unaligned is not depth
        let partial = GafRecord { counts_for_depth: false, ..spanning };
        let depth = compute_depth(&[partial], &HashMap::new(), &seg_len, &overlaps);
        assert_eq!(depth[&b"b".to_vec()], 0.0);
    }

    fn s(name: &str, o: Orientation) -> State {
        (name.as_bytes().to_vec(), o)
    }

    #[test]
    fn test_decompose_circuits_single_triangle() {
        use Orientation::Forward as F;
        let selected = vec![
            (s("a", F), s("b", F)),
            (s("b", F), s("c", F)),
            (s("c", F), s("a", F)),
        ];
        let circuits = decompose_circuits(&selected, &HashMap::new()).unwrap();
        assert_eq!(circuits.len(), 1);
        // the closing edge back to the start must not appear as a duplicate
        // final element
        assert_eq!(circuits[0].len(), 3);
        // Hierholzer's may start from any of the three nodes (HashMap
        // iteration order is unspecified) -- what matters is the cycle is
        // a valid rotation of a -> b -> c -> a.
        let expected = [s("a", F), s("b", F), s("c", F)];
        let start_idx = circuits[0]
            .iter()
            .position(|st| st == &expected[0])
            .expect("a must appear in the circuit");
        let rotated: Vec<_> = (0..3)
            .map(|i| circuits[0][(start_idx + i) % 3].clone())
            .collect();
        assert_eq!(rotated, expected);
    }

    #[test]
    fn test_decompose_circuits_two_disjoint_loops() {
        use Orientation::Forward as F;
        let selected = vec![
            (s("a", F), s("b", F)),
            (s("b", F), s("a", F)),
            (s("x", F), s("y", F)),
            (s("y", F), s("x", F)),
        ];
        let mut circuits = decompose_circuits(&selected, &HashMap::new()).unwrap();
        circuits.sort_by_key(|c| c.len());
        assert_eq!(circuits.len(), 2);
        assert_eq!(circuits[0].len(), 2);
        assert_eq!(circuits[1].len(), 2);
    }

    #[test]
    fn test_decompose_circuits_shared_repeat_node() {
        // a budget-2 repeat node visited via two distinct real junctions
        // should decompose into a single circuit that passes through it
        // twice, not two separate broken pieces
        use Orientation::Forward as F;
        let selected = vec![
            (s("a", F), s("repeat", F)),
            (s("repeat", F), s("b", F)),
            (s("b", F), s("repeat", F)),
            (s("repeat", F), s("a", F)),
        ];
        let circuits = decompose_circuits(&selected, &HashMap::new()).unwrap();
        assert_eq!(circuits.len(), 1);
        assert_eq!(circuits[0].len(), 4);
        let repeat_visits = circuits[0].iter().filter(|st| st.0 == b"repeat").count();
        assert_eq!(repeat_visits, 2);
    }

    #[test]
    fn test_decompose_circuits_keeps_a_single_segment_circuit() {
        // a circular unitig linked to itself is a circuit of one state
        use Orientation::Forward as F;
        let selected = vec![(s("a", F), s("a", F))];
        let circuits = decompose_circuits(&selected, &HashMap::new()).unwrap();
        assert_eq!(circuits, vec![vec![s("a", F)]]);
    }

    #[test]
    fn test_decompose_circuits_rejects_unbalanced_edges() {
        use Orientation::Forward as F;
        let selected = vec![(s("a", F), s("b", F))];
        assert!(decompose_circuits(&selected, &HashMap::new()).is_err());
    }

    /// A repeat `r` walked three times, between unique segments: every
    /// pairing of entries with exits gives closed walks, and only the reads
    /// say which arrangement is the real one.
    fn three_visit_repeat() -> Vec<(State, State)> {
        use Orientation::Forward as F;
        let mut selected = Vec::new();
        for (a, b) in [("a", "b"), ("c", "d"), ("e", "f")] {
            selected.push((s(a, F), s("r", F)));
            selected.push((s("r", F), s(b, F)));
        }
        // b -> c, d -> e, f -> a close the circle a r b c r d e r f
        for (a, b) in [("b", "c"), ("d", "e"), ("f", "a")] {
            selected.push((s(a, F), s(b, F)));
        }
        selected
    }

    fn follows(circuit: &[State], a: &str, via: &str, b: &str) -> bool {
        let n = circuit.len();
        (0..n).any(|i| {
            circuit[i].0 == a.as_bytes()
                && circuit[(i + 1) % n].0 == via.as_bytes()
                && circuit[(i + 2) % n].0 == b.as_bytes()
        })
    }

    #[test]
    fn test_decompose_circuits_pairs_a_repeat_by_read_triples() {
        use Orientation::Forward as F;
        let selected = three_visit_repeat();
        // reads say a-r-d, c-r-f, e-r-b: still one circuit (a r d e r b c r f)
        let mut walked = Vec::new();
        for (a, b) in [("a", "d"), ("c", "f"), ("e", "b")] {
            walked.extend(reads(&[s(a, F), s("r", F), s(b, F)], 10));
        }
        let triples = build_triple_lookup(&walked);
        let circuits = decompose_circuits(&selected, &triples).unwrap();
        assert_eq!(circuits.len(), 1);
        assert_eq!(circuits[0].len(), 9);
        for (a, b) in [("a", "d"), ("c", "f"), ("e", "b")] {
            assert!(follows(&circuits[0], a, "r", b), "{a} r {b} in {:?}", circuits[0]);
        }
    }

    #[test]
    fn test_decompose_circuits_joins_walks_where_reads_object_least() {
        use Orientation::Forward as F;
        let selected = three_visit_repeat();
        // reads say a-r-b (strongly), c-r-d and e-r-f (weakly): followed to
        // the letter that is three separate walks only if b, d and f led
        // back to a, c and e - here they lead on round the circle, so it is
        // already one circuit and nothing has to give
        let mut walked = reads(&[s("a", F), s("r", F), s("b", F)], 50);
        walked.extend(reads(&[s("c", F), s("r", F), s("d", F)], 6));
        walked.extend(reads(&[s("e", F), s("r", F), s("f", F)], 6));
        let triples = build_triple_lookup(&walked);
        let circuits = decompose_circuits(&selected, &triples).unwrap();
        assert_eq!(circuits.len(), 1);
        for (a, b) in [("a", "b"), ("c", "d"), ("e", "f")] {
            assert!(follows(&circuits[0], a, "r", b));
        }

        // now reads that would close three separate loops (b -> r -> a etc.
        // can't be: use a figure where they do): a r a', i.e. exits pair
        // back to their own entry's loop
        let loops = vec![
            (s("x", F), s("r", F)),
            (s("r", F), s("x", F)),
            (s("y", F), s("r", F)),
            (s("r", F), s("y", F)),
            (s("z", F), s("r", F)),
            (s("r", F), s("z", F)),
        ];
        let mut walked = reads(&[s("x", F), s("r", F), s("x", F)], 50);
        walked.extend(reads(&[s("y", F), s("r", F), s("y", F)], 20));
        walked.extend(reads(&[s("z", F), s("r", F), s("z", F)], 6));
        let triples = build_triple_lookup(&walked);
        let circuits = decompose_circuits(&loops, &triples).unwrap();
        // one circuit all the same: the three loops share r
        assert_eq!(circuits.len(), 1);
        assert_eq!(circuits[0].len(), 6);
    }

    #[test]
    fn test_decompose_circuits_ignores_the_order_edges_come_in() {
        use Orientation::Forward as F;
        let selected = three_visit_repeat();
        let walked = reads(&[s("a", F), s("r", F), s("d", F)], 10);
        let triples = build_triple_lookup(&walked);
        let expected = decompose_circuits(&selected, &triples).unwrap();
        let mut shuffled = selected.clone();
        shuffled.reverse();
        shuffled.swap(1, 5);
        assert_eq!(decompose_circuits(&shuffled, &triples).unwrap(), expected);
    }

    #[test]
    fn test_canonical_circuit_is_the_same_from_any_start_and_strand() {
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let circuit = vec![s("c", F), s("a", B), s("b", F), s("a", F)];
        let expected = vec![s("a", F), s("c", F), s("a", B), s("b", F)];
        assert_eq!(canonical_circuit(&circuit), expected);
        let rotated = vec![s("b", F), s("a", F), s("c", F), s("a", B)];
        assert_eq!(canonical_circuit(&rotated), expected);
        let other_strand: Vec<State> = circuit.iter().rev().map(|(n, o)| (n.clone(), flip(*o))).collect();
        assert_eq!(canonical_circuit(&other_strand), expected);
        assert_eq!(canonical_circuit(&[]), Vec::<State>::new());
    }

    #[test]
    fn test_build_candidate_edges_come_out_sorted() {
        use Orientation::Forward as F;
        let links = vec![
            (b"c".to_vec(), F, b"a".to_vec(), F),
            (b"a".to_vec(), F, b"b".to_vec(), F),
            (b"b".to_vec(), F, b"c".to_vec(), F),
        ];
        let (edges, _) = build_candidate_edges(&HashMap::new(), &HashMap::new(), &build_links_index(&links));
        let keys: Vec<(State, State)> = edges.iter().map(|e| (e.from.clone(), e.to.clone())).collect();
        let mut sorted = keys.clone();
        sorted.sort();
        assert_eq!(keys, sorted);
    }

    #[test]
    fn test_choose_root_is_in_the_largest_circuit() {
        use Orientation::Forward as F;
        // "hub" has the most candidate edges but is in no circuit; the root
        // must come from the larger circuit, its best-connected state
        let circuits = vec![vec![s("x", F), s("y", F)], vec![s("a", F), s("b", F), s("c", F)]];
        let mut edges = vec![
            edge(s("a", F), s("b", F), 1, 9),
            edge(s("b", F), s("c", F), 1, 9),
            edge(s("c", F), s("a", F), 1, 9),
            edge(s("c", F), s("b", F), 2, 1),
            edge(s("x", F), s("y", F), 1, 9),
            edge(s("y", F), s("x", F), 1, 9),
        ];
        for other in ["a", "b", "c", "x", "y"] {
            edges.push(edge(s("hub", F), s(other, F), 2, 1));
        }
        // b and c both have 4 edges; b sorts first
        assert_eq!(choose_root(&circuits, &edges), Some(s("b", F)));
        assert_eq!(choose_root(&[], &edges), None);
    }

    #[test]
    fn test_build_edge_support_pools_both_strands_of_a_junction() {
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        // 3 reads walk a+ -> b+, 4 walk it from the other strand (b- -> a-):
        // 7 for the junction, not two counts each under the threshold of 5
        let mut records = reads(&[s("a", F), s("b", F)], 3);
        records.extend(reads(&[s("b", B), s("a", B)], 4));
        // a junction that is its own reverse complement is counted once
        records.extend(reads(&[s("h", F), s("h", B)], 5));
        // too few reads on both strands together
        records.extend(reads(&[s("x", F), s("y", F)], 2));
        let support = build_edge_support(&records);
        assert_eq!(support[&s("a", F)][&s("b", F)], 7);
        assert_eq!(support[&s("b", B)][&s("a", B)], 7);
        assert_eq!(support[&s("h", F)][&s("h", B)], 5);
        assert!(!support.contains_key(&s("x", F)));
    }

    #[test]
    fn test_detect_bubbles_needs_entry_and_exit_on_one_strand() {
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let circuit_states: HashSet<State> = [s("a", F), s("b", F)].into_iter().collect();
        let excluded: HashSet<Seg> = [b"x".to_vec(), b"y".to_vec()].into_iter().collect();
        let edges = vec![
            // x: in on the forward strand, out on the reverse - not a walk through x
            edge(s("a", F), s("x", F), 1, 9),
            edge(s("x", B), s("b", F), 1, 9),
            // y: walked on its reverse strand from a to b
            edge(s("a", F), s("y", B), 1, 8),
            edge(s("y", B), s("b", F), 1, 6),
        ];
        let bubbles = detect_bubbles(&edges, &excluded, &circuit_states);
        assert_eq!(bubbles.len(), 1);
        assert_eq!(bubbles[0].segment, b"y".to_vec());
        assert_eq!(bubbles[0].orientation, B);
        assert_eq!((bubbles[0].entry.clone(), bubbles[0].exit.clone()), (s("a", F), s("b", F)));
        assert_eq!(bubbles[0].support_reads, 6);
    }

    #[test]
    fn test_try_splice_basic() {
        use Orientation::Forward as F;
        let main = vec![s("a", F), s("b", F), s("c", F)];
        let satellite = vec![s("x", F), s("y", F)];
        let edges = vec![
            edge(s("b", F), s("x", F), 0, 10), // entry
            edge(s("y", F), s("c", F), 0, 10), // exit
        ];
        let merged = try_splice(&main, &satellite, &edges).expect("splice should succeed");
        assert_eq!(
            merged,
            vec![s("a", F), s("b", F), s("x", F), s("y", F), s("c", F)]
        );
    }

    #[test]
    fn test_try_splice_returns_none_without_matching_evidence() {
        use Orientation::Forward as F;
        let main = vec![s("a", F), s("b", F), s("c", F)];
        let satellite = vec![s("x", F), s("y", F)];
        // an entry edge exists, but nothing exits back out -- not spliceable
        let edges = vec![edge(s("b", F), s("x", F), 0, 10)];
        assert!(try_splice(&main, &satellite, &edges).is_none());
    }

    #[test]
    fn test_try_splice_ignores_bare_topology_evidence() {
        use Orientation::Forward as F;
        let main = vec![s("a", F), s("b", F), s("c", F)];
        let satellite = vec![s("x", F), s("y", F)];
        // both edges present, but only as bare graph topology (tier 2) --
        // that's not real evidence, so this must not be treated as spliceable
        let edges = vec![
            edge(s("b", F), s("x", F), 2, 0),
            edge(s("y", F), s("c", F), 2, 0),
        ];
        assert!(try_splice(&main, &satellite, &edges).is_none());
    }

    #[test]
    fn test_splice_circuits_chain_through_growing_main() {
        // C only has real evidence attaching to states that live in B, not
        // in main directly -- it should still get absorbed once B has
        // already been grafted into main and those states are part of it.
        use Orientation::Forward as F;
        let main = vec![s("a", F), s("b", F), s("c", F)];
        let b_sat = vec![s("x", F), s("y", F)];
        let c_sat = vec![s("m", F), s("n", F)];
        let edges = vec![
            edge(s("b", F), s("x", F), 0, 10), // main <-> B entry
            edge(s("y", F), s("c", F), 0, 10), // main <-> B exit
            edge(s("y", F), s("m", F), 0, 10), // B <-> C entry
            edge(s("n", F), s("c", F), 0, 10), // B <-> C exit
        ];
        let (circuits, merges) = splice_circuits(vec![main, b_sat, c_sat], &edges);
        assert_eq!(merges, 2);
        assert_eq!(circuits.len(), 1);
        assert_eq!(circuits[0].len(), 7);
    }

    fn edge(from: State, to: State, tier: u8, reads: u32) -> CandidateEdge {
        CandidateEdge {
            from,
            to,
            tier,
            reads,
        }
    }

    #[test]
    fn test_levenshtein_basic() {
        assert_eq!(levenshtein(b"ACGT", b"ACGT", 1000), Some(0));
        assert_eq!(levenshtein(b"ACGT", b"ACGA", 1000), Some(1)); // substitution
        assert_eq!(levenshtein(b"ACGT", b"ACGTT", 1000), Some(1)); // insertion
        assert_eq!(levenshtein(b"ACGT", b"ACG", 1000), Some(1)); // deletion
        assert_eq!(levenshtein(b"", b"ACGT", 1000), Some(4));
        assert_eq!(levenshtein(b"ACGT", b"", 1000), Some(4));
    }

    #[test]
    fn test_levenshtein_respects_cell_cap() {
        // 100x100 = 10_000 cells, over a cap of 100
        let a = vec![b'A'; 100];
        let b = vec![b'C'; 100];
        assert_eq!(levenshtein(&a, &b, 100), None);
    }

    #[test]
    fn test_find_main_detour_basic() {
        use Orientation::Forward as F;
        let circuit = vec![s("a", F), s("b", F), s("c", F), s("d", F)];
        let detour = find_main_detour(&[circuit], &s("a", F), &s("d", F), 12).unwrap();
        assert_eq!(detour, vec![s("b", F), s("c", F)]);
    }

    #[test]
    fn test_find_main_detour_adjacent_is_empty() {
        // exit immediately follows entry -- main path has no sequence
        // there at all, i.e. the bubble is a pure insertion
        use Orientation::Forward as F;
        let circuit = vec![s("a", F), s("b", F), s("c", F)];
        let detour = find_main_detour(&[circuit], &s("a", F), &s("b", F), 12).unwrap();
        assert!(detour.is_empty());
    }

    #[test]
    fn test_find_main_detour_none_beyond_hop_bound() {
        use Orientation::Forward as F;
        let circuit = vec![s("a", F), s("b", F), s("c", F), s("d", F)];
        assert!(find_main_detour(&[circuit], &s("a", F), &s("d", F), 1).is_none());
    }

    #[test]
    fn test_find_main_detour_picks_shortest_across_repeat_occurrences() {
        // "a" (a repeat) occurs twice; the occurrence right before "d"
        // should win over the one that's much further away
        use Orientation::Forward as F;
        let circuit = vec![
            s("a", F),
            s("x", F),
            s("y", F),
            s("z", F),
            s("a", F),
            s("d", F),
        ];
        let detour = find_main_detour(&[circuit], &s("a", F), &s("d", F), 12).unwrap();
        assert!(detour.is_empty());
    }

    #[test]
    fn test_diff_bubble_arm_reports_comparison_for_similar_sized_sequences() {
        let bubble = b"ACGTACGTAC";
        let detour = b"ACGTACGTAA"; // 1bp different, same length
        match diff_bubble_arm(bubble, detour) {
            BubbleDiff::Compared {
                edit_distance,
                pct_identity,
                compared_bp,
            } => {
                assert_eq!(edit_distance, 1);
                assert_eq!(compared_bp, 10);
                assert!((pct_identity - 90.0).abs() < 1e-9);
            }
            _ => panic!("expected Compared, got a different variant"),
        }
    }

    #[test]
    fn test_diff_bubble_arm_flags_wildly_larger_detour_instead_of_diffing() {
        let bubble = vec![b'A'; 100];
        // 10x the bubble's length -- over the 5x ratio cap
        let detour = vec![b'C'; 1000];
        match diff_bubble_arm(&bubble, &detour) {
            BubbleDiff::DetourTooLarge { detour_bp, ratio } => {
                assert_eq!(detour_bp, 1000);
                assert_eq!(ratio, 10);
            }
            _ => panic!("expected DetourTooLarge"),
        }
    }

    #[test]
    fn test_diff_bubble_arm_within_ratio_cap_still_compares() {
        // exactly at the 5x cap should still be compared, not flagged
        let bubble = vec![b'A'; 100];
        let detour = vec![b'A'; 500];
        assert!(matches!(
            diff_bubble_arm(&bubble, &detour),
            BubbleDiff::Compared { .. }
        ));
    }

    /// Three segments in a ring with 3-base overlaps:
    /// a = AAACCCGGG, b = GGGTTTACG, c = ACGTGCAAA, so a+ -> b+ shares GGG,
    /// b+ -> c+ shares ACG and c+ -> a+ shares AAA. The circle they spell,
    /// each overlap written once, is CCCGGGTTTACGTGCAAA (18 bases).
    fn ring() -> (HashMap<Seg, Vec<u8>>, HashMap<(State, State), usize>) {
        use Orientation::Forward as F;
        let segments: HashMap<Seg, Vec<u8>> = [
            (b"a".to_vec(), b"AAACCCGGG".to_vec()),
            (b"b".to_vec(), b"GGGTTTACG".to_vec()),
            (b"c".to_vec(), b"ACGTGCAAA".to_vec()),
        ]
        .into_iter()
        .collect();
        let links = vec![
            (b"a".to_vec(), F, b"b".to_vec(), F),
            (b"b".to_vec(), F, b"c".to_vec(), F),
            (b"c".to_vec(), F, b"a".to_vec(), F),
        ];
        let cigars = vec![b"3M".to_vec(); 3];
        (segments, build_overlaps(&links, &cigars).unwrap())
    }

    #[test]
    fn test_link_overlap() {
        assert_eq!(link_overlap(b"41M").unwrap(), 41);
        assert_eq!(link_overlap(b"0M").unwrap(), 0);
        assert_eq!(link_overlap(b"*").unwrap(), 0);
        assert!(link_overlap(b"3M1I2M").is_err());
    }

    #[test]
    fn test_path_sequence_circuit_writes_each_overlap_once() {
        use Orientation::Forward as F;
        let (segments, overlaps) = ring();
        let circuit = [s("a", F), s("b", F), s("c", F)];
        let seq = path_sequence(&circuit, &segments, &overlaps, true).unwrap();
        // the closing c -> a overlap (AAA) comes off the start of a
        assert_eq!(seq, b"CCCGGGTTTACGTGCAAA".to_vec());
        // 27 bases of segment, less three 3-base overlaps
        assert_eq!(seq.len(), 27 - 9);
    }

    #[test]
    fn test_path_sequence_circuit_is_the_same_circle_from_any_start() {
        use Orientation::Forward as F;
        let (segments, overlaps) = ring();
        let from_a = path_sequence(&[s("a", F), s("b", F), s("c", F)], &segments, &overlaps, true).unwrap();
        let from_b = path_sequence(&[s("b", F), s("c", F), s("a", F)], &segments, &overlaps, true).unwrap();
        assert_eq!(from_a.len(), from_b.len());
        let doubled = [from_a.clone(), from_a].concat();
        assert!(doubled.windows(from_b.len()).any(|w| w == from_b.as_slice()));
    }

    #[test]
    fn test_path_sequence_reverse_strand_is_the_reverse_complement() {
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let (segments, overlaps) = ring();
        let fwd = path_sequence(&[s("a", F), s("b", F), s("c", F)], &segments, &overlaps, true).unwrap();
        // the same circle walked on the other strand: a- -> c- -> b-
        let rev = path_sequence(&[s("a", B), s("c", B), s("b", B)], &segments, &overlaps, true).unwrap();
        let rc = revcomp(&fwd);
        let doubled = [rc.clone(), rc].concat();
        assert_eq!(rev.len(), fwd.len());
        assert!(doubled.windows(rev.len()).any(|w| w == rev.as_slice()));
    }

    #[test]
    fn test_path_sequence_linear_keeps_the_first_segment_whole() {
        use Orientation::Forward as F;
        let (segments, overlaps) = ring();
        let seq = path_sequence(&[s("a", F), s("b", F)], &segments, &overlaps, false).unwrap();
        assert_eq!(seq, b"AAACCCGGGTTTACG".to_vec());
        assert!(path_sequence(&[], &segments, &overlaps, false).unwrap().is_empty());
    }

    #[test]
    fn test_path_sequence_single_segment_circuit_trims_its_self_link() {
        use Orientation::Forward as F;
        let segments: HashMap<Seg, Vec<u8>> = [(b"x".to_vec(), b"ACGTTTGGCACG".to_vec())].into_iter().collect();
        let overlaps = build_overlaps(&[(b"x".to_vec(), F, b"x".to_vec(), F)], &[b"3M".to_vec()]).unwrap();
        let seq = path_sequence(&[s("x", F)], &segments, &overlaps, true).unwrap();
        assert_eq!(seq, b"TTTGGCACG".to_vec());
    }

    #[test]
    fn test_inexact_junctions_reports_only_overlaps_that_differ() {
        use Orientation::Forward as F;
        let (mut segments, overlaps) = ring();
        let circuit = [s("a", F), s("b", F), s("c", F)];
        assert!(inexact_junctions(&circuit, &segments, &overlaps).is_empty());
        // b now starts GTG where a ends GGG
        segments.insert(b"b".to_vec(), b"GTGTTTACG".to_vec());
        let found = inexact_junctions(&circuit, &segments, &overlaps);
        assert_eq!(found, vec![(s("a", F), s("b", F), 3, 1)]);
    }

    #[test]
    fn test_path_sequence_errors_on_a_transition_with_no_link() {
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let (segments, overlaps) = ring();
        assert!(path_sequence(&[s("a", F), s("b", B)], &segments, &overlaps, false).is_err());
    }

    fn rec(path: &[State]) -> GafRecord {
        GafRecord {
            path: path.to_vec(),
            path_start: 0,
            path_end: 0,
            mapq: 60,
            identity: None,
            counts_for_depth: true,
        }
    }

    fn reads(path: &[State], n: usize) -> Vec<GafRecord> {
        (0..n).map(|_| rec(path)).collect()
    }

    fn is_licensed(links_index: &HashMap<State, HashMap<State, u32>>, e: &CandidateEdge) -> bool {
        links_index
            .get(&e.from)
            .is_some_and(|m| m.contains_key(&e.to))
    }

    /// The u14 shape from vpJunScop1 (Juniperus scopulorum): u72 and u96
    /// both attach only to u14's *minus* end (`L u14 - u72 -`,
    /// `L u14 - u96 -`), and nothing attaches to u14's plus end.
    fn u14_links() -> Vec<(Seg, Orientation, Seg, Orientation)> {
        use Orientation::Backward as B;
        vec![
            (b"u14".to_vec(), B, b"u72".to_vec(), B),
            (b"u14".to_vec(), B, b"u96".to_vec(), B),
            (b"u72".to_vec(), B, b"p".to_vec(), B),
            (b"u96".to_vec(), B, b"q".to_vec(), B),
        ]
    }

    #[test]
    fn test_build_links_index_registers_both_traversals() {
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let idx = build_links_index(&u14_links());
        // as written: u14- -> u72-
        assert!(idx[&s("u14", B)].contains_key(&s("u72", B)));
        // reverse-complement equivalent: u72+ -> u14+
        assert!(idx[&s("u72", F)].contains_key(&s("u14", F)));
        // and never the cross-orientation u72+ -> u14-
        assert!(!idx[&s("u72", F)].contains_key(&s("u14", B)));
    }

    #[test]
    fn test_build_links_index_self_reverse_complement_link_counted_once() {
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let idx = build_links_index(&[(b"x".to_vec(), F, b"x".to_vec(), B)]);
        assert_eq!(idx[&s("x", F)][&s("x", B)], 1);
    }

    #[test]
    fn test_build_triple_lookup_reverse_strand_flips_every_state() {
        // a real read walking u14- -> u72- -> p- must, read from the other
        // strand, register p+ -> u72+ -> u14+ -- not -> u14-
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let records = reads(&[s("u14", B), s("u72", B), s("p", B)], MIN_EDGE_READS as usize);
        let lookup = build_triple_lookup(&records);
        let ctr = &lookup[&(b"p".to_vec(), F, b"u72".to_vec(), F)];
        assert!(ctr.contains_key(&s("u14", F)));
        assert!(!ctr.contains_key(&s("u14", B)));
    }

    #[test]
    fn test_build_candidate_edges_only_emits_gfa_licensed_transitions() {
        // regression: vpJunScop1 resolved `u72+,u14-,u96-`, where
        // u72+ -> u14- is not a valid move in the graph. Real reads in
        // both directions, plus a chimeric read claiming the impossible
        // cross-traversal u72+ -> u14- -> u96- often enough to clear
        // MIN_EDGE_READS.
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let n = MIN_EDGE_READS as usize;
        let mut records = Vec::new();
        records.extend(reads(&[s("u14", B), s("u72", B), s("p", B)], n));
        records.extend(reads(&[s("u14", B), s("u96", B), s("q", B)], n));
        records.extend(reads(&[s("p", F), s("u72", F), s("u14", F)], n));
        records.extend(reads(&[s("u72", F), s("u14", B), s("u96", B)], n));

        let links_index = build_links_index(&u14_links());
        let (edges, n_unlicensed) = build_candidate_edges(
            &build_triple_lookup(&records),
            &build_edge_support(&records),
            &links_index,
        );

        for e in &edges {
            assert!(
                is_licensed(&links_index, e),
                "unlicensed candidate edge {} -> {} (tier {})",
                state_str(&e.from),
                state_str(&e.to),
                e.tier
            );
        }
        assert!(!edges
            .iter()
            .any(|e| e.from == s("u72", F) && e.to == s("u14", B)));
        assert!(n_unlicensed > 0);
        // the real junctions survive, in both traversal directions, with
        // read support rather than falling back to bare topology
        for (u, v) in [
            (s("u14", B), s("u72", B)),
            (s("u72", F), s("u14", F)),
            (s("u14", B), s("u96", B)),
        ] {
            let e = edges.iter().find(|e| e.from == u && e.to == v).unwrap();
            assert!(e.tier < 2, "{} -> {} lost its read support", state_str(&u), state_str(&v));
        }
    }

    #[test]
    fn test_solve_ilp_rejects_orientation_disconnected_repeat() {
        // A budget-2 hub segment used via two *different* orientations in
        // two otherwise-unconnected 2-cycles. The connectivity constraint
        // must not be satisfiable by aggregating flow across both
        // orientations of the same segment -- each orientation's loop has
        // to be genuinely reachable from the root on its own, or the
        // solver should refuse rather than hand back a result that
        // decomposes into two circuits despite claiming success.
        use Orientation::Backward as B;
        use Orientation::Forward as F;
        let edges = vec![
            edge(s("h", F), s("s1", F), 0, 10),
            edge(s("s1", F), s("h", F), 0, 10),
            edge(s("h", B), s("s2", F), 0, 10),
            edge(s("s2", F), s("h", B), 0, 10),
        ];
        let mut budget = HashMap::new();
        budget.insert(b"h".to_vec(), 2u32);
        budget.insert(b"s1".to_vec(), 1u32);
        budget.insert(b"s2".to_vec(), 1u32);
        let mut segments = HashMap::new();
        segments.insert(b"h".to_vec(), b"A".to_vec());
        segments.insert(b"s1".to_vec(), b"A".to_vec());
        segments.insert(b"s2".to_vec(), b"A".to_vec());
        let gene_segs = HashSet::new();

        match solve_ilp(&edges, &budget, &segments, &gene_segs, Some(&s("h", F)), 5.0) {
            Ok((selected, _)) => {
                let circuits = decompose_circuits(&selected, &HashMap::new()).unwrap();
                assert_eq!(
                    circuits.len(),
                    1,
                    "connectivity constraint claimed success but produced {} disjoint circuits: {:?}",
                    circuits.len(),
                    circuits
                );
            }
            Err(_) => {
                // refusing is also an acceptable outcome -- there is no
                // single-circuit solution available at all here, since s1
                // and s2 only ever attach through opposite orientations of
                // h with no edge between them.
            }
        }
    }
}
