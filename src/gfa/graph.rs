use crate::gfa::gfa::GAFTSVRecord;
use crate::utils::GFAGraphLookups;
use anyhow::{bail, Context, Result};
use gfa::gfa::Orientation;
use gfa::gfa::GFA;
use gfa::optfields::OptFields;
use itertools::Itertools;
use petgraph::{
    dot::{Config, Dot},
    graph::{Graph, IndexType, NodeIndex},
    visit::{EdgeRef, IntoNodeIdentifiers},
    Directed,
    Direction::Outgoing,
    Undirected,
};
use std::collections::HashSet;

pub struct GFAungraph(pub Graph<usize, (), Undirected>);

impl GFAungraph {
    // TEST THIS
    pub fn recursive_search(
        &self,
        sequence_id: usize,
        iterations: i32,
        mut collect_sequence_names: Vec<NodeIndex>,
        graph_indices: GFAGraphLookups,
    ) -> Result<Vec<usize>> {
        let gfa_graph = &self.0;
        eprintln!(
            "[+]\tRecursively searching around node {} at depth {}",
            sequence_id, iterations
        );

        let mut iteration = 0;
        loop {
            // collect all the neighbours of all the current node indices
            for index in collect_sequence_names.clone() {
                for c in gfa_graph.neighbors(index) {
                    collect_sequence_names.push(c);
                }
            }
            // add sorting and deduping here too
            // yes otherwise vectors are enormous.
            collect_sequence_names.sort();
            collect_sequence_names.dedup();

            iteration += 1;
            if iteration == iterations {
                break;
            }
        }

        collect_sequence_names.sort();
        collect_sequence_names.dedup();

        let mut sequences_to_keep = Vec::new();
        // turn node indexes into sequence ID's
        for index in collect_sequence_names {
            let t = graph_indices.node_index_to_seg_id(index)?;
            sequences_to_keep.push(t);
        }

        Ok(sequences_to_keep)
    }
}

// weights are the orientations, used at various points, and an optional
// coverage weight, used in gfatk linear.
// GFA's should always specify Links in a specific direction..?
// so digraphs should be where all the functionality lies.
pub struct GFAdigraph(pub Graph<usize, (Orientation, Orientation, Option<u32>)>);

impl GFAdigraph {
    pub fn debug_with_dot(&self) {
        let gfa_graph = &self.0;
        eprintln!(
            "{:?}",
            Dot::with_config(&gfa_graph, &[Config::GraphContentOnly])
        );
    }
    // we want weakly connected components, as there may only be an edge in one
    // orientation (perhaps unlikely... but still)

    // thanks https://github.com/Qiskit/retworkx/blob/79900cf8da0c0665ac5ce1ccb0f57373434b14b8/src/connectivity/mod.rs
    pub fn weakly_connected_components(
        &self,
        graph_indices: GFAGraphLookups,
    ) -> Result<Vec<Vec<usize>>> {
        let graph = &self.0;
        let mut seen: HashSet<NodeIndex> = HashSet::with_capacity(graph.node_count());
        let mut out_vec: Vec<Vec<usize>> = Vec::new();

        for node in graph.node_indices() {
            if !seen.contains(&node) {
                // BFS node generator

                let mut component_set: std::collections::BTreeSet<NodeIndex> =
                    std::collections::BTreeSet::new();

                let mut bfs_seen: HashSet<NodeIndex> = HashSet::new();

                let mut next_level: HashSet<NodeIndex> = HashSet::new();

                next_level.insert(node);

                while !next_level.is_empty() {
                    let this_level = next_level;

                    next_level = HashSet::new();

                    for bfs_node in this_level {
                        if !bfs_seen.contains(&bfs_node) {
                            component_set.insert(bfs_node);

                            bfs_seen.insert(bfs_node);

                            for neighbor in graph.neighbors_undirected(bfs_node) {
                                next_level.insert(neighbor);
                            }
                        }
                    }
                }
                let set_to_vec: Vec<_> = component_set.iter().cloned().collect();
                // convert node indices to segment ID's
                let x = set_to_vec
                    .iter()
                    .map(|e| {
                        let seg_id = match graph_indices.node_index_to_seg_id(*e) {
                            Ok(s) => s,
                            Err(err) => bail!(
                                "NodeIndex {:?} could not be converted to segment ID.\n{}",
                                e,
                                err
                            ),
                        };
                        Ok(seg_id)
                    })
                    .collect::<Result<Vec<usize>>>();

                out_vec.push(x?);

                seen.extend(bfs_seen);
            }
        }
        Ok(out_vec)
    }

    // mutably insert the coverages
    pub fn add_coverages(
        &mut self,
        coverages: &Option<Result<Vec<GAFTSVRecord>>>,
        graph_indices: &GFAGraphLookups,
    ) -> Result<()> {
        let gfa_graph = &mut self.0;
        // iterate over the edges of the graph
        let mut coverage_vec = Vec::new();
        for edge in gfa_graph.edge_references() {
            // get node ID's
            let from_id = graph_indices.node_index_to_seg_id(edge.source())?;
            let to_id = graph_indices.node_index_to_seg_id(edge.target())?;
            // get the weight of this edge
            let from_orient = edge.weight().0;
            let to_orient = edge.weight().1;

            let coverage = match coverages {
                Some(c) => {
                    // is it okay just to unwrap here?
                    let coverages = c.as_ref().unwrap();
                    let mut res: Option<u32> = None;
                    for cv in coverages {
                        // if everything matches...
                        if from_id == cv.from
                            && to_id == cv.to
                            && from_orient.plus_minus_as_byte() == cv.from_orient as u8
                            && to_orient.plus_minus_as_byte() == cv.to_orient as u8
                        {
                            res = Some(cv.coverage);
                        }
                    }
                    res
                }
                None => None,
            };

            coverage_vec.push((edge.id(), from_orient, to_orient, coverage));
        }
        // had to do this, as cannot immutably borrow, and mutably borrow in the same loop...
        for (edge_id, from_orient, to_orient, coverage) in coverage_vec {
            // get rid of this unwrap.
            *gfa_graph
                .edge_weight_mut(edge_id)
                .with_context(|| format!("Not able to change edge ID: {:?}", edge_id))? =
                (from_orient, to_orient, coverage);
        }
        Ok(())
    }

    // before we search for all the paths, we need to
    // do a check to see if the links indicate both
    // forward and reverse links
    // do this in petgraph by adding nodes and edges.
    pub fn check_both_strands(&mut self) {
        let gfa_graph = &mut self.0;

        let mut edge_vec = Vec::new();
        let mut checks = Vec::new();

        for edge in gfa_graph.edge_references() {
            let from_orient = edge.weight().0;
            let to_orient = edge.weight().1;

            let check = gfa_graph.contains_edge(edge.target(), edge.source());
            checks.push(check);

            edge_vec.push((edge.source(), edge.target(), from_orient, to_orient));
        }

        let return_early = checks.iter().all(|e| *e == true);
        if return_early {
            // don't add loads of extra edges!
            eprintln!("[+]\tGFA contains forward and reverse links.");
            return;
        }

        eprintln!(
            "[+]\tGFA does not contain forward and reverse links. Adding them to internal graph."
        );

        for (source, target, from_orient, to_orient) in edge_vec {
            // if orientations differ, just swap indexes
            if from_orient != to_orient {
                gfa_graph.add_edge(target, source, (from_orient, to_orient, None));
            } else if from_orient == to_orient {
                // reverse the orientation
                let opp_orient = match from_orient.is_reverse() {
                    true => Orientation::Forward,
                    false => Orientation::Backward,
                };
                gfa_graph.add_edge(target, source, (opp_orient, opp_orient, None));
            }
        }
    }

    // core algorithm for gfatk linear
    // iterate over all pairs of nodes to find the paths
    // then search though all these paths to find the
    // longest which doesn't violate sequence orientation

    // need to emit segment ID's that did not make it into the
    // longest path here.

    // return the path, the path ID's, and id's not in the path
    // -> (Vec<NodeIndex>, Vec<usize>, Vec<usize>)

    pub fn all_paths_all_node_pairs(
        &self,
        graph_indices: &GFAGraphLookups,
        coverages: Option<Result<Vec<GAFTSVRecord>>>,
    ) -> Result<(Vec<NodeIndex>, Vec<usize>)> {
        let graph = &self.0;
        let nodes = graph.node_identifiers();
        let node_pairs = nodes.permutations(2).collect::<Vec<_>>();
        // eprintln!("{:?}", node_pairs);

        let all_paths: Vec<_> = node_pairs
            .into_iter()
            .map(|pair| {
                let path = all_paths(&graph, pair[0], pair[1]);
                path
            })
            .collect();

        // make the set of legal paths through the GFA

        let mut valid_paths = Vec::new();
        // iterate over the paths
        for paths in all_paths {
            // iterate over each path
            for path in paths {
                // iterate over adjacent nodes
                let node_pairs = path.windows(2);
                // and the skipped iterator
                let node_pairs_skip = path.windows(2).skip(1);
                // assess whether we should keep a path
                let mut keep = false;
                // so we can compare NodeIndex(0), NodeIndex(1), and NodeIndex(2) directly
                'node_pairs_loop: for (pair1, pair2) in node_pairs.zip(node_pairs_skip) {
                    // first node
                    let from_p1 = pair1[0];
                    // second node
                    let to_p1 = pair1[1];
                    // second node (again)
                    let from_p2 = pair2[0];
                    // third node
                    let to_p2 = pair2[1];

                    // so what we really want is to take the first and second nodes
                    // get all the edges
                    // then get all the edges from the third to the second node
                    // added NodeIndexes here for debugging
                    let a_b_edges: Vec<(
                        NodeIndex,
                        NodeIndex,
                        (Orientation, Orientation, Option<u32>),
                    )> = graph
                        .edges_connecting(from_p1, to_p1)
                        .map(|e| {
                            let s = e.source();
                            let t = e.target();
                            (s, t, *e.weight())
                        })
                        .collect();

                    let c_b_edges: Vec<(
                        NodeIndex,
                        NodeIndex,
                        (Orientation, Orientation, Option<u32>),
                    )> = graph
                        .edges_connecting(to_p2, from_p2)
                        .map(|e| {
                            let s = e.source();
                            let t = e.target();
                            (s, t, *e.weight())
                        })
                        .collect();

                    // we can then compare the orientation of the 'to' Orientation
                    // for a->b and c->b
                    let mut keep_vec = Vec::new();
                    for e in &a_b_edges {
                        for f in &c_b_edges {
                            let (_, a_b_to, _) = e.2;
                            let (_, c_b_to, _) = f.2;
                            // we found a path through!
                            // i.e. the links are not connected to the
                            if a_b_to != c_b_to {
                                // keep = true;
                                keep_vec.push(true);
                            } else if a_b_to == c_b_to {
                                keep_vec.push(false);
                            }
                        }
                    }
                    // keep if any of the elements is true
                    let do_keep = keep_vec.iter().any(|e| *e == true);

                    // debugging...
                    // eprintln!("a->b: {:?}", a_b_edges);
                    // eprintln!("c->b: {:?}", c_b_edges);
                    // eprintln!("Do not keep: {}", do_keep);
                    // eprintln!("");

                    // if we got to here and diff is still false, break out of this path
                    // it's a no-go...
                    if do_keep {
                        keep = true;
                    } else {
                        keep = false;
                        break 'node_pairs_loop;
                    }
                }
                if keep {
                    valid_paths.push(path);
                }
            }
        }
        valid_paths.sort_by(|a, b| b.len().cmp(&a.len()));
        valid_paths.dedup();

        // check if we have coverages
        // if we do then we shall incorporate this information

        let final_path = match coverages.is_none() {
            // if no coverage, take (one of) the longest path(s).
            true => valid_paths.to_vec()[0].clone(),
            false => {
                // now iterate over all the paths again
                let mut index = 0;
                // iterate over all the paths of this length
                // keep track of lengths
                let mut track_coverage = Vec::new();
                let mut final_path = vec![];
                let mut final_coverage = 0;

                for path in &valid_paths {
                    let mut path_coverage = 0;
                    let node_pairs = path.windows(2);
                    for pair in node_pairs {
                        let from = pair[0];
                        let to = pair[1];
                        let connecting = &mut graph.edges_connecting(from, to);
                        let coverage = connecting
                            .next()
                            .with_context(|| {
                                format!("No connecting edges from {:?} to {:?}", from, to)
                            })?
                            .weight()
                            .2;
                        match coverage {
                            Some(c) => path_coverage += c,
                            None => (),
                        }
                    }
                    track_coverage.push(path_coverage);
                    let prev = track_coverage.get(index - 1).unwrap_or(&0);
                    let cur = track_coverage[index];
                    if &cur > prev {
                        final_path = path.to_vec();
                        final_coverage = path_coverage;
                    }

                    index += 1;
                }
                eprintln!("[+]\tHighest cumulative coverage path = {}", final_coverage);
                final_path
            }
        };

        let mut chosen_path_string = Vec::new();
        let final_path_node_pairs = final_path.windows(2);
        for (index, pair) in final_path_node_pairs.enumerate() {
            let from = pair[0];
            let to = pair[1];
            let connecting = &mut graph.edges_connecting(from, to);
            let weight = connecting
                .next()
                .with_context(|| format!("No connecting edges from {:?} to {:?}", from, to))?;
            let from_orient = weight.weight().0;
            let to_orient = weight.weight().1;

            // get segment ID from Node Indices
            let from = graph_indices.node_index_to_seg_id(weight.source())?;
            let to = graph_indices.node_index_to_seg_id(weight.target())?;

            if index == 0 {
                chosen_path_string
                    .push(format!("{} {} -> {} {}", from_orient, from, to_orient, to));
            } else {
                chosen_path_string.push(format!(" -> {} {}", to_orient, to));
            }
        }

        eprintln!(
            "[+]\tChosen path through graph: {}",
            chosen_path_string.join("")
        );

        // sort out the path now
        // we need the strandedness information
        // now sort these overlaps so they are the same order
        // chosen path -> id's
        let chosen_path_ids = final_path
            .iter()
            // .1 needed
            .map(|e| {
                let x = graph_indices
                    .0
                    .iter()
                    .find(|y| y.node_index == *e)
                    .with_context(|| {
                        format!("NodeIndex {:?} could not be converted to segment ID", e)
                    });
                x
            })
            .collect::<Result<Vec<_>>>();

        Ok((
            final_path.to_vec(),
            // error handling a bit annoying here.
            chosen_path_ids?
                .iter()
                .map(|e| e.seg_id)
                .collect::<Vec<usize>>()
                .to_vec(),
        ))
    }

    // simple graph stats
    pub fn node_count(&self) -> usize {
        let gfa_graph = &self.0;

        gfa_graph.node_count()
    }

    pub fn edge_count(&self) -> usize {
        let gfa_graph = &self.0;

        gfa_graph.edge_count()
    }
}

// thanks
// https://github.com/Ninjani/rosalind/blob/e22ecf2c9f0935d970b137684029957c0850d63f/t_ba11b/src/lib.rs
pub fn all_paths<T, U, Ix: IndexType>(
    graph: &Graph<T, U, Directed, Ix>,
    start_node: NodeIndex<Ix>,
    end_node: NodeIndex<Ix>,
) -> Vec<Vec<NodeIndex<Ix>>> {
    let mut visited = HashSet::new();
    visited.insert(start_node);
    all_paths_helper(graph, start_node, end_node, &mut visited)
}

fn all_paths_helper<T, U, Ix: IndexType>(
    graph: &Graph<T, U, Directed, Ix>,
    start_node: NodeIndex<Ix>,
    end_node: NodeIndex<Ix>,
    visited: &mut HashSet<NodeIndex<Ix>>,
) -> Vec<Vec<NodeIndex<Ix>>> {
    if start_node == end_node {
        vec![vec![end_node]]
    } else {
        let mut paths = Vec::new();
        for edge in graph.edges_directed(start_node, Outgoing) {
            let next_node = edge.target();
            if !visited.contains(&next_node) {
                visited.insert(next_node);
                let descendant_paths = all_paths_helper(graph, next_node, end_node, visited);
                visited.remove(&next_node);
                paths.extend(
                    descendant_paths
                        .into_iter()
                        .map(|path| {
                            let mut new_path = vec![start_node];
                            new_path.extend(path);
                            new_path
                        })
                        .collect::<Vec<_>>(),
                )
            }
        }
        paths
    }
}

// Returns a subgraph GFA that only contains elements with the
// provided segment names
// https://github.com/chfi/rs-gfa-utils/blob/master/src/subgraph.rs

pub fn segments_subgraph<T: OptFields + Clone>(
    gfa: &GFA<usize, T>,
    segment_names: Vec<usize>,
) -> GFA<usize, T> {
    let segments = gfa
        .segments
        .iter()
        .filter(|s| segment_names.contains(&s.name))
        .cloned()
        .collect();

    let links = gfa
        .links
        .iter()
        .filter(|l| {
            segment_names.contains(&l.from_segment) && segment_names.contains(&l.to_segment)
        })
        .cloned()
        .collect();

    let containments = gfa
        .containments
        .iter()
        .filter(|l| {
            segment_names.contains(&l.container_name) && segment_names.contains(&l.contained_name)
        })
        .cloned()
        .collect();

    let paths: Vec<_> = gfa
        .paths
        .iter()
        .filter(|p| p.iter().any(|(s, _)| segment_names.contains(&s)))
        .cloned()
        .collect();

    GFA {
        header: gfa.header.clone(),
        segments,
        links,
        paths,
        containments,
    }
}
