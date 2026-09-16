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
    variable, Expression, ProblemVariables, Solution, SolverModel, Variable, WithTimeLimit,
};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::gaf::{read_gaf_records, GafRecord};
use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;

type Seg = Vec<u8>;
type State = (Seg, Orientation);

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
fn compute_depth(records: &[GafRecord], seg_len: &HashMap<Seg, usize>) -> HashMap<Seg, f64> {
    let mut aligned_bases: HashMap<Seg, u64> = HashMap::new();

    for rec in records {
        let mut offset: usize = 0;
        for (seg, _o) in &rec.path {
            let Some(&len) = seg_len.get(seg) else {
                continue;
            };
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

/// The depth of a "typical" single-copy segment (the median of all nonzero
/// depths), and each segment's estimated copy number relative to it,
/// floored at 1 for segments with too little evidence to say otherwise.
fn estimate_budget(depth: &HashMap<Seg, f64>) -> (f64, HashMap<Seg, u32>) {
    let mut values: Vec<f64> = depth.values().copied().filter(|v| *v > 0.0).collect();
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let baseline = if values.is_empty() {
        1.0
    } else if values.len() % 2 == 1 {
        values[values.len() / 2]
    } else {
        let mid = values.len() / 2;
        (values[mid - 1] + values[mid]) / 2.0
    };
    let baseline = if baseline > 0.0 {
        baseline
    } else {
        values.last().copied().unwrap_or(1.0)
    };

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
fn build_triple_lookup(
    records: &[GafRecord],
) -> HashMap<(Seg, Orientation, Seg, Orientation), HashMap<State, u32>> {
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

    let mut lookup: HashMap<(Seg, Orientation, Seg, Orientation), HashMap<State, u32>> =
        HashMap::new();
    for ((a, ao, seg, so, b, bo), cnt) in canon {
        if cnt < MIN_EDGE_READS {
            continue;
        }
        lookup
            .entry((a.clone(), ao, seg.clone(), so))
            .or_default()
            .insert((b.clone(), bo), cnt);
        lookup
            .entry((b, flip(bo), seg, flip(so)))
            .or_default()
            .insert((a, ao), cnt);
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
    edge_support.retain(|_, m| {
        m.retain(|_, &mut cnt| cnt >= MIN_EDGE_READS);
        !m.is_empty()
    });
    edge_support
}

fn build_links_index(
    links: &[(Seg, Orientation, Seg, Orientation)],
) -> HashMap<State, HashMap<State, u32>> {
    let mut idx: HashMap<State, HashMap<State, u32>> = HashMap::new();
    for (fs, fo, ts, to) in links {
        *idx.entry((fs.clone(), *fo))
            .or_default()
            .entry((ts.clone(), *to))
            .or_insert(0) += 1;
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
fn build_candidate_edges(
    triple_lookup: &HashMap<(Seg, Orientation, Seg, Orientation), HashMap<State, u32>>,
    edge_support: &HashMap<State, HashMap<State, u32>>,
    links_index: &HashMap<State, HashMap<State, u32>>,
) -> Vec<CandidateEdge> {
    let mut best: HashMap<(State, State), (u8, u32)> = HashMap::new();

    for ((a, ao, seg, so), ctr) in triple_lookup {
        for (&(ref b, bo), &cnt) in ctr {
            consider(&mut best, (a.clone(), *ao), (seg.clone(), *so), cnt, 0);
            consider(&mut best, (seg.clone(), *so), (b.clone(), bo), cnt, 0);
        }
    }
    for (state, ctr) in edge_support {
        for (next, &cnt) in ctr {
            consider(&mut best, state.clone(), next.clone(), cnt, 1);
        }
    }
    for (state, ctr) in links_index {
        for (next, &cnt) in ctr {
            consider(&mut best, state.clone(), next.clone(), cnt, 2);
        }
    }

    best.into_iter()
        .map(|((from, to), (tier, reads))| CandidateEdge {
            from,
            to,
            tier,
            reads,
        })
        .collect()
}

/// Solve the integer program: pick the subset of candidate edges maximising
/// segment coverage (weighted towards genes), subject to each segment's
/// depth budget and flow conservation at every state. When `enforce_connectivity`
/// is set, an additional single-commodity-flow constraint forces the result
/// into one circuit rather than however many disjoint ones flow conservation
/// alone happens to allow -- see the comment further down for why that's
/// needed. This is the harder problem to solve, and on a few inputs proved
/// too slow for microlp within a reasonable time budget, so the caller
/// retries without it on failure rather than leaving the tool hanging.
fn solve_ilp(
    edges: &[CandidateEdge],
    budget: &HashMap<Seg, u32>,
    segments: &HashMap<Seg, Vec<u8>>,
    gene_segs: &HashSet<Seg>,
    enforce_connectivity: bool,
) -> Result<Vec<(State, State)>> {
    let mut states: HashSet<State> = HashSet::new();
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
    let mut y: HashMap<Seg, Variable> = HashMap::new();
    for seg in segments.keys() {
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
    // existing flow conservation is enough to force one circuit: pick the
    // best-connected segment as a root, and require every other *included*
    // segment be reachable from it via the selected edges. This is a
    // standard single-commodity-flow subtour-elimination technique -- a
    // second set of continuous flow variables piggybacks on the same edges
    // (gated by whether x_e is selected) purely to prove reachability. It
    // reuses the existing per-segment `y` as the consumption indicator
    // rather than introducing a parallel per-state variable with its own
    // big-M linking constraints, which made the solver numerically
    // unstable ("singular matrix") on some inputs.
    let root = if enforce_connectivity {
        let mut seg_degree: HashMap<Seg, usize> = HashMap::new();
        for s in &states {
            let d = incoming.get(s).map_or(0, |v| v.len()) + outgoing.get(s).map_or(0, |v| v.len());
            *seg_degree.entry(s.0.clone()).or_insert(0) += d;
        }
        seg_degree
            .iter()
            .max_by_key(|(_, d)| **d)
            .map(|(s, _)| s.clone())
    } else {
        None
    };

    let mut flow: HashMap<(State, State), Variable> = HashMap::new();
    if root.is_some() {
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
        .with_time_limit(45.0);

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

    if let Some(root) = root {
        // a tighter bound than states.len() would give: total connectivity
        // supply can never exceed the number of segments actually drawing
        // from it, so this is already a valid (if generous) capacity
        let capacity = segments.len() as f64;

        let mut flow_in: HashMap<Seg, Vec<Variable>> = HashMap::new();
        let mut flow_out: HashMap<Seg, Vec<Variable>> = HashMap::new();
        for e in edges {
            let key = (e.from.clone(), e.to.clone());
            let fvar = flow[&key];
            model.add_constraint((1.0 * fvar).leq(capacity * x[&key]));
            flow_out.entry(e.from.0.clone()).or_default().push(fvar);
            flow_in.entry(e.to.0.clone()).or_default().push(fvar);
        }

        model.add_constraint((1.0 * y[&root]).eq(1.0));

        let others_included: Expression = y
            .iter()
            .filter(|(seg, _)| **seg != root)
            .map(|(_, &yvar)| 1.0 * yvar)
            .sum();

        for seg in segments.keys() {
            let inflow: Expression = flow_in
                .get(seg)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .sum();
            let outflow: Expression = flow_out
                .get(seg)
                .cloned()
                .unwrap_or_default()
                .into_iter()
                .sum();
            if *seg == root {
                // the root supplies exactly enough connectivity flow for
                // every other included segment to draw one unit from it
                model.add_constraint((inflow - outflow + others_included.clone()).eq(0.0));
            } else {
                model.add_constraint((inflow - outflow - 1.0 * y[seg]).eq(0.0));
            }
        }
    }

    let solution = model.solve().context("ILP solve failed")?;

    let selected = edges
        .iter()
        .filter(|e| solution.value(x[&(e.from.clone(), e.to.clone())]) > 0.5)
        .map(|e| (e.from.clone(), e.to.clone()))
        .collect();
    Ok(selected)
}

/// Split a flow-conserving edge multiset into disjoint closed circuits
/// (Hierholzer's algorithm). Each returned circuit lists every visited
/// state once; the closing return-to-start edge is not duplicated.
fn decompose_circuits(selected: &[(State, State)]) -> Vec<Vec<State>> {
    let mut adj: HashMap<State, Vec<State>> = HashMap::new();
    let mut remaining: HashMap<(State, State), i64> = HashMap::new();
    for (u, v) in selected {
        adj.entry(u.clone()).or_default().push(v.clone());
        *remaining.entry((u.clone(), v.clone())).or_insert(0) += 1;
    }

    let mut circuits = Vec::new();
    let starts: Vec<State> = adj.keys().cloned().collect();

    for start in starts {
        let has_out = |u: &State, remaining: &HashMap<(State, State), i64>| {
            adj.get(u)
                .map(|vs| {
                    vs.iter()
                        .any(|v| *remaining.get(&(u.clone(), v.clone())).unwrap_or(&0) > 0)
                })
                .unwrap_or(false)
        };
        if !has_out(&start, &remaining) {
            continue;
        }

        let mut stack = vec![start.clone()];
        let mut circuit = Vec::new();
        let mut cur = start;
        loop {
            let mut advanced = false;
            if let Some(vs) = adj.get(&cur).cloned() {
                for v in vs {
                    let key = (cur.clone(), v.clone());
                    if *remaining.get(&key).unwrap_or(&0) > 0 {
                        *remaining.get_mut(&key).unwrap() -= 1;
                        stack.push(v.clone());
                        cur = v;
                        advanced = true;
                        break;
                    }
                }
            }
            if !advanced {
                circuit.push(stack.pop().expect("stack always has the current node"));
                match stack.last() {
                    Some(top) => cur = top.clone(),
                    None => break,
                }
            }
        }
        circuit.reverse();
        if circuit.len() > 1 && circuit.first() == circuit.last() {
            circuit.pop();
        }
        if circuit.len() > 1 {
            circuits.push(circuit);
        }
    }

    circuits
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
    entry: State,
    exit: State,
    support_reads: u32,
}

fn detect_bubbles(
    edges: &[CandidateEdge],
    excluded: &HashSet<Seg>,
    circuit_states: &HashSet<State>,
) -> Vec<BubbleArm> {
    let mut into: HashMap<Seg, (State, State, u32)> = HashMap::new();
    let mut out_of: HashMap<Seg, (State, State, u32)> = HashMap::new();

    for e in edges {
        if e.tier >= 2 {
            continue; // bare topology only, not real read evidence
        }
        if excluded.contains(&e.to.0) && circuit_states.contains(&e.from) {
            let entry = into
                .entry(e.to.0.clone())
                .or_insert((e.from.clone(), e.to.clone(), 0));
            if e.reads > entry.2 {
                *entry = (e.from.clone(), e.to.clone(), e.reads);
            }
        }
        if excluded.contains(&e.from.0) && circuit_states.contains(&e.to) {
            let entry = out_of
                .entry(e.from.0.clone())
                .or_insert((e.from.clone(), e.to.clone(), 0));
            if e.reads > entry.2 {
                *entry = (e.from.clone(), e.to.clone(), e.reads);
            }
        }
    }

    let mut bubbles = Vec::new();
    for seg in excluded {
        if let (Some(in_edge), Some(out_edge)) = (into.get(seg), out_of.get(seg)) {
            bubbles.push(BubbleArm {
                segment: seg.clone(),
                entry: in_edge.0.clone(),
                exit: out_edge.1.clone(),
                support_reads: in_edge.2.min(out_edge.2),
            });
        }
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

fn state_str(s: &State) -> String {
    format!("{}{}", String::from_utf8_lossy(&s.0), orient_char(s.1))
}

fn write_fasta_record(header: &str, seq: &[u8]) {
    println!(">{header}");
    for chunk in seq.chunks(70) {
        println!("{}", String::from_utf8_lossy(chunk));
    }
}

fn circuit_sequence(circuit: &[State], segments: &HashMap<Seg, Vec<u8>>) -> Vec<u8> {
    let mut seq = Vec::new();
    for (seg, o) in circuit {
        let s = &segments[seg];
        match o {
            Orientation::Forward => seq.extend_from_slice(s),
            Orientation::Backward => seq.extend(revcomp(s)),
        }
    }
    seq
}

/// The main entry point for `gfatk resolve`.
pub fn resolve(matches: &clap::ArgMatches) -> Result<()> {
    let gfa_file = matches.get_one::<PathBuf>("GFA");
    let gaf_file = matches
        .get_one::<PathBuf>("gaf")
        .context("`--gaf` is required for `gfatk resolve`")?;
    let gff_file = matches.get_one::<PathBuf>("gff");
    let min_mapq = *matches
        .get_one::<u32>("min-mapq")
        .expect("defaulted by clap");
    let min_identity = *matches
        .get_one::<f64>("min-identity")
        .expect("defaulted by clap");

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
    let records = read_gaf_records(gaf_file, min_mapq, min_identity)?;
    eprintln!(
        "[+]\t{} usable read-path records after quality filtering",
        records.len()
    );

    let depth = compute_depth(&records, &seg_len);
    let (baseline, budget) = estimate_budget(&depth);
    eprintln!("[+]\tBaseline (single-copy) read depth: {baseline:.1}x");

    let triple_lookup = build_triple_lookup(&records);
    let edge_support = build_edge_support(&records);
    let links_index = build_links_index(&links);
    let edges = build_candidate_edges(&triple_lookup, &edge_support, &links_index);
    eprintln!(
        "[+]\t{} candidate edges from read/graph evidence",
        edges.len()
    );

    let gene_segs: HashSet<Seg> = match gff_file {
        Some(path) => {
            let genes = parse_gff_gene_segments(path)?;
            eprintln!(
                "[+]\t{} segments carry a high-confidence annotated gene (from {path:?})",
                genes.len()
            );
            genes.into_keys().collect()
        }
        None => HashSet::new(),
    };

    let selected = match solve_ilp(&edges, &budget, &segments, &gene_segs, true) {
        Ok(selected) => selected,
        Err(e) => {
            eprintln!(
                "[-]\tCould not solve for a single connected circuit within the time budget ({e}); \
                 falling back to the unconstrained solve, which may report multiple circuits."
            );
            solve_ilp(&edges, &budget, &segments, &gene_segs, false)?
        }
    };
    let circuits = {
        let mut c = decompose_circuits(&selected);
        c.sort_by_key(|c| std::cmp::Reverse(c.len()));
        c
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
        let seq = circuit_sequence(circuit, &segments);
        let path_str = circuit.iter().map(state_str).collect::<Vec<_>>().join(",");
        eprintln!(
            "[+]\tCircuit {}: {} segments, {} bp, path={}",
            idx + 1,
            circuit.len(),
            seq.len(),
            path_str
        );
        let species_tag = gfa_file
            .and_then(|p| p.file_stem())
            .map(|s| s.to_string_lossy().to_string())
            .unwrap_or_else(|| "gfatk_resolve".to_string());
        let header = format!(
            "{species_tag}_circuit{}:n_segments={}:bp={}",
            idx + 1,
            circuit.len(),
            seq.len()
        );
        write_fasta_record(&header, &seq);
    }

    if circuits.len() > 1 {
        eprintln!(
            "[-]\t{} separate circuits were resolved -- flow conservation alone guarantees a valid \
             decomposition into closed loops, not that they merge into one, so this can happen even \
             when the underlying GFA graph is a single connected component. Checking for direct \
             evidence between each pair:",
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
                        "[-]\t  Circuit {} <-> Circuit {}: real {} evidence between {} and {} (~{} reads) -- \
                         budget contention likely forced these apart, not a lack of connecting evidence.",
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
    }

    let bubbles = detect_bubbles(&edges, &excluded, &circuit_states);
    if !bubbles.is_empty() {
        eprintln!(
            "[-]\t{} segment(s) excluded from the resolved circuit(s) are recombination-bubble arms \
             (real read evidence connects them to the main circuit on both sides, but the dominant \
             conformation won out) -- reported as separate records below, not treated as missing.",
            bubbles.len()
        );
        for b in &bubbles {
            let genes = gene_segs.contains(&b.segment);
            eprintln!(
                "[-]\t  {} (between {} and {}, ~{} supporting reads{})",
                String::from_utf8_lossy(&b.segment),
                state_str(&b.entry),
                state_str(&b.exit),
                b.support_reads,
                if genes {
                    ", carries an annotated gene"
                } else {
                    ""
                }
            );
            let seq = &segments[&b.segment];
            let header = format!(
                "bubble_arm:{}:between={}-{}:support_reads={}",
                String::from_utf8_lossy(&b.segment),
                state_str(&b.entry),
                state_str(&b.exit),
                b.support_reads
            );
            write_fasta_record(&header, seq);
        }
    }

    let truly_missing: Vec<&Seg> = excluded
        .iter()
        .filter(|s| !bubbles.iter().any(|b| &b.segment == *s))
        .collect();
    if !truly_missing.is_empty() {
        eprintln!(
            "[-]\t{} segment(s) excluded with no bubble evidence either (genuinely unresolved, not just an alternate conformation)",
            truly_missing.len()
        );
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let (baseline, budget) = estimate_budget(&depth);
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
        let (_, budget) = estimate_budget(&depth);
        assert_eq!(budget[&b"weak".to_vec()], 1);
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
        let circuits = decompose_circuits(&selected);
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
        let mut circuits = decompose_circuits(&selected);
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
        let circuits = decompose_circuits(&selected);
        assert_eq!(circuits.len(), 1);
        assert_eq!(circuits[0].len(), 4);
        let repeat_visits = circuits[0].iter().filter(|st| st.0 == b"repeat").count();
        assert_eq!(repeat_visits, 2);
    }
}
