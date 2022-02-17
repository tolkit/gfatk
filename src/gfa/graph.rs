use crate::gfa::gfa::GAFTSVRecord;
use gfa::gfa::Orientation;
use itertools::Itertools;
use petgraph::{
    algo::{is_cyclic_directed, kosaraju_scc},
    dot::{Config, Dot},
    graph::{Graph, IndexType, NodeIndex},
    visit::{EdgeRef, IntoNodeIdentifiers, NodeIndexable},
    Directed,
    Direction::Outgoing,
    EdgeDirection, Undirected,
};
use std::collections::HashSet;
use std::error::Error;

pub struct GFAungraph(pub Graph<usize, (), Undirected>);

impl GFAungraph {
    pub fn recursive_search(
        &self,
        sequence_id: usize,
        iterations: i32,
        mut collect_sequence_names: Vec<NodeIndex>,
        graph_indices: Vec<(NodeIndex, usize)>,
    ) -> Vec<usize> {
        let gfa_graph = &self.0;
        eprintln!(
            "[+]\tRecursively searching around node {} at depth {}",
            sequence_id, iterations
        );

        let mut iteration = 0;
        loop {
            for index in collect_sequence_names.clone() {
                for c in gfa_graph.neighbors(index) {
                    let res = graph_indices.iter().find(|x| x.0 == c).unwrap().0;
                    collect_sequence_names.push(res);
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
            let t = graph_indices.iter().find(|x| x.0 == index).unwrap().1;
            sequences_to_keep.push(t);
        }
        sequences_to_keep
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
        graph_indices: Vec<(NodeIndex, usize)>,
    ) -> Vec<Vec<usize>> {
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
                    .map(|e| graph_indices.iter().find(|y| y.0 == *e).unwrap().1)
                    .collect::<Vec<usize>>();

                out_vec.push(x);

                seen.extend(bfs_seen);
            }
        }
        out_vec
    }
    // split the input GFA graph representation into subgraphs
    // where one subgraph might be e.g. the mitochondria.
    // don't know how efficient this function is.
    pub fn split_into_connected_digraphs(&self) -> Vec<Self> {
        let gfa_graph = &self.0;

        // thanks https://github.com/ybyygu/gchemol/blob/f2d139311b71f927b23daca128b8d8d1a240a96d/core/src/molecule/fragment.rs
        let components = kosaraju_scc(gfa_graph);
        println!("Components: {:?}", components);

        println!(
            "Connected component algo: {}",
            petgraph::algo::connected_components(gfa_graph)
        );

        let mut graphs = vec![];
        for nodes in components {
            let g = gfa_graph.filter_map(
                // node closure:
                // keep nodes in the same component
                |i, n| {
                    if nodes.contains(&i) {
                        Some(n.clone())
                    } else {
                        None
                    }
                },
                // edge closure:
                // keep the edge if connected nodes are both in the same component
                |i, e| {
                    let (n1, n2) = gfa_graph.edge_endpoints(i).unwrap();
                    if nodes.contains(&n1) && nodes.contains(&n2) {
                        Some(e.clone())
                    } else {
                        None
                    }
                },
            );
            graphs.push(GFAdigraph(g));
        }

        graphs
    }

    // check if graph is cyclic
    // do we want to check this? Does it need to be cyclic?
    pub fn check_is_cyclic_directed(&self) {
        let gfa_graph = &self.0;
        match is_cyclic_directed(&gfa_graph) {
            true => eprintln!("[+]\tLoaded GFA is a cyclic digraph."),
            false => eprintln!("[+]\tLoaded GFA is not a cyclic digraph."),
        }
    }

    // maybe this should be changed
    // to edge at each end of the segment. how?
    pub fn get_edge_counts(&self) -> Vec<(usize, i32)> {
        let gfa_graph = &self.0;

        // iterate over the nodes
        // just increment until there are no more nodes.
        let mut node_index = 0;
        // collect nodes and their edge counts
        let mut node_edge_counts = Vec::new();
        loop {
            // need to peek ahead so we don't loop forever
            let mut peekable_iter = gfa_graph
                .neighbors_undirected(NodeIndex::new(node_index))
                .peekable();

            // need a different way of breaking here...
            match peekable_iter.peek().is_some() {
                true => {
                    let mut node_count = 0;
                    for _ in gfa_graph.neighbors_undirected(NodeIndex::new(node_index)) {
                        node_count += 1;
                    }
                    node_edge_counts.push((node_index, node_count));
                    node_index += 1;
                }
                false => {
                    if node_index == gfa_graph.node_bound() {
                        break;
                    } else {
                        node_index += 1;
                        continue;
                    }
                }
            }
        }
        node_edge_counts
    }

    // this currently is only reliable for cyclic graphs.
    // how do we get a source target pair from a linear graph
    // if it's not cyclic, we could find paths between all
    // pairs which have only two connections (i.e. only connect to one
    // other segment +/-)
    pub fn get_source_target_pair(
        &self,
        graph_indices: &Vec<(NodeIndex, usize)>,
        node_edge_counts: Vec<(usize, i32)>,
    ) -> Option<(NodeIndex, NodeIndex)> {
        let gfa_graph = &self.0;
        // calculate adjecent nodes, each of which have two neighbours
        // iterate over edges again but stop when we hit a neighbouring
        // pair of nodes which are connected each themselves to at least 2 others
        let mut node_index = 0;
        let mut source_target_pair: Option<(NodeIndex, NodeIndex)> = None;

        loop {
            // need to peek ahead so we don't loop forever
            let mut peekable_iter = gfa_graph.edges(NodeIndex::new(node_index)).peekable();
            match peekable_iter.peek().is_some() {
                true => {
                    for edge in gfa_graph.edges(NodeIndex::new(node_index)) {
                        let source = edge.source();
                        let target = edge.target();

                        let source_id = graph_indices.iter().find(|x| x.0 == source).unwrap().0;
                        let target_id = graph_indices.iter().find(|x| x.0 == target).unwrap().0;

                        let source_edge_count = node_edge_counts
                            .iter()
                            .filter(|e| NodeIndex::new(e.0) == source_id)
                            .collect::<Vec<&(usize, i32)>>()[0]
                            .1;
                        let target_edge_count = node_edge_counts
                            .iter()
                            .filter(|e| NodeIndex::new(e.0) == target_id)
                            .collect::<Vec<&(usize, i32)>>()[0]
                            .1;

                        node_index += 1;

                        if source_edge_count > 1 && target_edge_count > 1 {
                            source_target_pair = Some((source, target));
                            // and break early here
                            break;
                        }
                    }
                }
                false => {
                    if node_index == gfa_graph.node_bound() {
                        break;
                    } else {
                        node_index += 1;
                        continue;
                    }
                }
            }
        }
        source_target_pair
    }

    // mutably insert the coverages
    pub fn add_coverages(
        &mut self,
        coverages: &Option<Result<Vec<GAFTSVRecord>, Box<dyn Error>>>,
        graph_indices: &Vec<(NodeIndex, usize)>,
    ) {
        let gfa_graph = &mut self.0;
        // iterate over the edges of the graph
        let mut coverage_vec = Vec::new();
        for edge in gfa_graph.edge_references() {
            // get node ID's
            let from_id = graph_indices
                .iter()
                .find(|x| x.0 == edge.source())
                .unwrap()
                .1;
            let to_id = graph_indices
                .iter()
                .find(|x| x.0 == edge.target())
                .unwrap()
                .1;
            // get the weight of this edge
            let from_orient = edge.weight().0;
            let to_orient = edge.weight().1;

            let coverage = match coverages {
                Some(ref c) => {
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
            *gfa_graph.edge_weight_mut(edge_id).unwrap() = (from_orient, to_orient, coverage);
        }
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

    // get the path out as both a vec of node indices
    // and as actual node names

    // this is wrong when there are edges coming only
    // from a single orientation from a segment.

    // pub fn find_hamiltonian_path(
    //     &self,
    //     source_target_pair: Option<(NodeIndex, NodeIndex)>,
    //     graph_indices: &Vec<(NodeIndex, usize)>,
    //     coverages: Option<Result<Vec<GAFTSVRecord>, Box<dyn Error>>>,
    // ) -> (Vec<NodeIndex>, Vec<usize>) {
    //     let gfa_graph = &self.0;

    //     // debugging
    //     eprintln!(
    //         "[+]\tSearching between nodes {} and {}",
    //         graph_indices
    //             .iter()
    //             .find(|y| y.0 == source_target_pair.unwrap().0)
    //             .unwrap()
    //             .1,
    //         graph_indices
    //             .iter()
    //             .find(|y| y.0 == source_target_pair.unwrap().1)
    //             .unwrap()
    //             .1
    //     );

    //     // for cyclic digraphs this below works
    //     // for more linear graphs, choosing the
    //     // start and end node is a bit more complex.
    //     // maybe user can input start & end nodes.

    //     let mut paths = all_paths(
    //         &gfa_graph,
    //         source_target_pair.unwrap().1,
    //         source_target_pair.unwrap().0,
    //     );
    //     // debugging
    //     // eprintln!("{:?}", source_target_pair);
    //     for path in &paths {
    //         eprintln!("{:?}", path);
    //     }
    //     // TODO: in all the paths, check that there are links at either end.

    //     // need a conditional here, if we have the data...
    //     // else...
    //     let final_path = match coverages.is_none() {
    //         true => {
    //             // just take the longest path (slight randomness to this)
    //             // sort & dedup paths and pick longest one
    //             paths.sort_by(|a, b| b.len().cmp(&a.len()));
    //             // paths.iter().for_each(|e| eprintln!("{}", e.len()));
    //             let final_path = paths[0].to_vec();
    //             final_path
    //         }
    //         false => {
    //             // or by taking into account coverages
    //             // length of the longest path
    //             let longest_path_length_vec = paths.iter().map(|e| e.len()).collect::<Vec<usize>>();
    //             let longest_path_length = longest_path_length_vec.iter().max().unwrap();
    //             // iterate over all the paths of this length
    //             // keep track of lengths
    //             let mut track_lens = Vec::new();
    //             let mut index = 0;
    //             let mut final_path = vec![];
    //             let mut final_coverage = 0;
    //             for path in &paths {
    //                 // if the length of the path is maximal
    //                 if &path.len() == longest_path_length {
    //                     let mut path_len = 0;
    //                     // for each node in the path
    //                     // find the corresponding
    //                     let node_pairs = path.windows(2);
    //                     for pair in node_pairs {
    //                         let from = pair[0];
    //                         let to = pair[1];
    //                         let connecting = &mut gfa_graph.edges_connecting(from, to);
    //                         let coverage = connecting.next().unwrap().weight().2;
    //                         match coverage {
    //                             Some(c) => path_len += c,
    //                             None => (),
    //                         }
    //                     }
    //                     track_lens.push(path_len);
    //                     let prev = track_lens.get(index - 1).unwrap_or(&0);
    //                     let cur = track_lens[index];
    //                     if &cur > prev {
    //                         final_path = path.to_vec();
    //                         final_coverage = path_len;
    //                     }

    //                     index += 1;
    //                 }
    //             }
    //             eprintln!("[+]\tTotal coverage for chosen path is {}", final_coverage);
    //             final_path
    //         }
    //     };

    //     // eprintln!("Final path is {:?}", final_path);
    //     // calculate the total coverage for each path

    //     // if paths length is zero, bail out here.
    //     if paths.is_empty() {
    //         panic!("No simple paths were found.")
    //     }

    //     let mut chosen_path_string2 = Vec::new();
    //     let final_path_node_pairs = final_path.windows(2);
    //     for (index, pair) in final_path_node_pairs.enumerate() {
    //         let from = pair[0];
    //         let to = pair[1];
    //         let connecting = &mut gfa_graph.edges_connecting(from, to);
    //         let weight = connecting.next().unwrap();
    //         let from_orient = weight.weight().0;
    //         let to_orient = weight.weight().1;
    //         let from = graph_indices
    //             .iter()
    //             .find(|y| y.0 == weight.source())
    //             .unwrap()
    //             .1;
    //         let to = graph_indices
    //             .iter()
    //             .find(|y| y.0 == weight.target())
    //             .unwrap()
    //             .1;

    //         if index == 0 {
    //             chosen_path_string2
    //                 .push(format!("{} {} -> {} {}", from_orient, from, to_orient, to));
    //         } else {
    //             chosen_path_string2.push(format!(" -> {} {}", to_orient, to));
    //         }
    //     }

    //     eprintln!(
    //         "[+]\tChosen path through graph: {}",
    //         chosen_path_string2.join("")
    //     );

    //     // sort out the path now
    //     // we need the strandedness information
    //     // now sort these overlaps so they are the same order
    //     // chosen path -> id's
    //     let chosen_path_ids = final_path
    //         .iter()
    //         .map(|e| graph_indices.iter().find(|y| y.0 == *e).unwrap().1)
    //         .collect::<Vec<_>>();

    //     (final_path.to_vec(), chosen_path_ids.to_vec())
    // }

    // this might be crazy
    pub fn all_paths_all_node_pairs(
        &self,
        graph_indices: &Vec<(NodeIndex, usize)>,
        coverages: Option<Result<Vec<GAFTSVRecord>, Box<dyn Error>>>,
    ) -> (Vec<NodeIndex>, Vec<usize>) {
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

        // check if we have coverages
        // if we do then we shall incorporate this information

        let final_path = match coverages.is_none() {
            true => {
                // let mut valid_paths = Vec::new();
                // iterate over the paths
                for paths in all_paths {
                    // iterate over each path
                    for path in paths {
                        // iterate over adjacent nodes
                        let node_pairs = path.windows(2);
                        // and the skipped iterator
                        let node_pairs_skip = path.windows(2).skip(1);
                        // so we can compare NodeIndex(0), NodeIndex(1), and NodeIndex(2) directly
                        for (pair1, pair2) in node_pairs.zip(node_pairs_skip) {
                            let from_p1 = pair1[0];
                            let to_p1 = pair1[1];
                            let from_p2 = pair2[0];
                            let to_p2 = pair2[1];

                            // get the edge outgoing NodeIndex(0) -> NodeIndex(1)
                            let edge_outgoing_p1 = graph.find_edge(from_p1, to_p1).unwrap();
                            // get the edge incoming NodeIndex(2) -> NodeIndex(1) 
                            let edge_incoming_p2 = graph.find_edge(to_p2, from_p2).unwrap();
                            // get the weights for the focal node
                            let w_out = graph.edge_weight(edge_outgoing_p1).unwrap().1;
                            let w_in = graph.edge_weight(edge_incoming_p2).unwrap().1;
                            // debug
                            // YEEEEEEEEEEEEEEEES
                            if w_out == w_in {
                                eprintln!(
                                    "Between nodes: {:?}, {:?}, {:?}, {:?} :: {}/{}",
                                    from_p1, to_p1, from_p2, to_p2, w_out, w_in
                                );
                            }
                        }
                    }
                }
                vec![petgraph::graph::NodeIndex::new(0)]
            }
            false => {
                // or by taking into account coverages
                // length of the longest path

                // keep track of longest path length
                let mut longest_path_length = 0usize;
                // search through all_paths to find length of longest path
                for paths in &all_paths {
                    let longest_path_length_vec =
                        paths.iter().map(|e| e.len()).collect::<Vec<usize>>();
                    let longest_path_length_temp = longest_path_length_vec.iter().max().unwrap();

                    if longest_path_length_temp > &longest_path_length {
                        longest_path_length = *longest_path_length_temp;
                    }
                }

                // now iterate over all the paths again
                let mut final_paths_and_coverages = Vec::new();

                for paths in &all_paths {
                    // iterate over all the paths of this length
                    // keep track of lengths
                    let mut track_lens = Vec::new();
                    let mut index = 0;
                    let mut final_path = vec![];
                    let mut final_coverage = 0;
                    for path in paths.to_vec() {
                        // if the length of the path is maximal
                        if path.len() == longest_path_length {
                            let mut path_len = 0;
                            // for each node in the path
                            // find the corresponding
                            let node_pairs = path.windows(2);
                            for pair in node_pairs {
                                let from = pair[0];
                                let to = pair[1];
                                let connecting = &mut graph.edges_connecting(from, to);
                                let coverage = connecting.next().unwrap().weight().2;
                                match coverage {
                                    Some(c) => path_len += c,
                                    None => (),
                                }
                            }
                            track_lens.push(path_len);
                            let prev = track_lens.get(index - 1).unwrap_or(&0);
                            let cur = track_lens[index];
                            if &cur > prev {
                                final_path = path.to_vec();
                                final_coverage = path_len;
                            }

                            index += 1;
                        }
                    }
                    final_paths_and_coverages.push((final_coverage, final_path));
                }
                final_paths_and_coverages.sort_by(|(a, _), (b, _)| b.cmp(&a));
                final_paths_and_coverages[0].1.to_vec()
            }
        };

        let mut chosen_path_string2 = Vec::new();
        let final_path_node_pairs = final_path.windows(2);
        for (index, pair) in final_path_node_pairs.enumerate() {
            let from = pair[0];
            let to = pair[1];
            let connecting = &mut graph.edges_connecting(from, to);
            let weight = connecting.next().unwrap();
            let from_orient = weight.weight().0;
            let to_orient = weight.weight().1;
            let from = graph_indices
                .iter()
                .find(|y| y.0 == weight.source())
                .unwrap()
                .1;
            let to = graph_indices
                .iter()
                .find(|y| y.0 == weight.target())
                .unwrap()
                .1;

            if index == 0 {
                chosen_path_string2
                    .push(format!("{} {} -> {} {}", from_orient, from, to_orient, to));
            } else {
                chosen_path_string2.push(format!(" -> {} {}", to_orient, to));
            }
        }

        eprintln!(
            "[+]\tChosen path through graph: {}",
            chosen_path_string2.join("")
        );

        // sort out the path now
        // we need the strandedness information
        // now sort these overlaps so they are the same order
        // chosen path -> id's
        let chosen_path_ids = final_path
            .iter()
            .map(|e| graph_indices.iter().find(|y| y.0 == *e).unwrap().1)
            .collect::<Vec<_>>();

        (final_path.to_vec(), chosen_path_ids.to_vec())
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

use gfa::gfa::GFA;
use gfa::optfields::OptFields;

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
