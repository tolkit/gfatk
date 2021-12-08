use crate::gfa::gfa::GAFTSVRecord;
use gfa::gfa::Orientation;
use petgraph::{
    algo::is_cyclic_directed,
    graph::{Graph, IndexType, NodeIndex},
    visit::{EdgeRef, NodeIndexable},
    Directed,
    Direction::Outgoing,
    Undirected,
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
pub struct GFAdigraph(pub Graph<usize, (Orientation, Orientation, Option<u32>)>);

impl GFAdigraph {
    // check if graph is cyclic
    // do we want to check this? Does it need to be cyclic?
    pub fn check_is_cyclic_directed(&self) {
        let gfa_graph = &self.0;
        match is_cyclic_directed(&gfa_graph) {
            true => eprintln!("[+]\tLoaded GFA is a cyclic digraph."),
            false => eprintln!("[+]\tLoaded GFA is not a cyclic digraph."),
        }
    }

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
    pub fn find_hamiltonian_path(
        &self,
        source_target_pair: Option<(NodeIndex, NodeIndex)>,
        graph_indices: &Vec<(NodeIndex, usize)>,
        coverages: Option<Result<Vec<GAFTSVRecord>, Box<dyn Error>>>,
    ) -> (Vec<NodeIndex>, Vec<usize>) {
        let gfa_graph = &self.0;
        // debugging
        // use petgraph::dot::{Config, Dot};
        // eprintln!(
        //     "{:?}",
        //     Dot::with_config(&gfa_graph, &[Config::GraphContentOnly])
        // );

        // debugging
        eprintln!(
            "[+]\tSearching between nodes {} and {}",
            graph_indices
                .iter()
                .find(|y| y.0 == source_target_pair.unwrap().0)
                .unwrap()
                .1,
            graph_indices
                .iter()
                .find(|y| y.0 == source_target_pair.unwrap().1)
                .unwrap()
                .1
        );

        // for cyclic digraphs this below works
        // for more linear graphs, choosing the
        // start and end node is a bit more complex.
        // maybe user can input start & end nodes.

        let mut paths = all_paths(
            &gfa_graph,
            source_target_pair.unwrap().1,
            source_target_pair.unwrap().0,
        );
        // debugging
        // eprintln!("{:?}", source_target_pair);
        // eprintln!("{:?}", paths);

        // need a conditional here, if we have the data...
        // else...
        let final_path = match coverages.is_none() {
            true => {
                // just take the longest path (slight randomness to this)
                // sort & dedup paths and pick longest one
                paths.sort_by(|a, b| b.len().cmp(&a.len()));
                let final_path = paths[0].to_vec();
                final_path
            }
            false => {
                // or by taking into account coverages
                // length of the longest path
                let longest_path_length_vec = paths.iter().map(|e| e.len()).collect::<Vec<usize>>();
                let longest_path_length = longest_path_length_vec.iter().max().unwrap();
                // iterate over all the paths of this length
                // keep track of lengths
                let mut track_lens = Vec::new();
                let mut index = 0;
                let mut final_path = vec![];
                let mut final_coverage = 0;
                for path in &paths {
                    // if the length of the path is maximal
                    if &path.len() == longest_path_length {
                        let mut path_len = 0;
                        // for each node in the path
                        // find the corresponding
                        let node_pairs = path.windows(2);
                        for pair in node_pairs {
                            let from = pair[0];
                            let to = pair[1];
                            let connecting = &mut gfa_graph.edges_connecting(from, to);
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
                eprintln!("[+]\tTotal coverage for chosen path is {}", final_coverage);
                final_path
            }
        };

        // eprintln!("Final path is {:?}", final_path);
        // calculate the total coverage for each path

        // if paths length is zero, bail out here.
        if paths.is_empty() {
            panic!("No simple paths were found.")
        }
        // sort & dedup paths and pick longest one
        // paths.sort_by(|a, b| b.len().cmp(&a.len()));

        // format string just for my (and user) interest
        let chosen_path_string = final_path
            .iter()
            .map(|e| format!("{}", graph_indices.iter().find(|y| y.0 == *e).unwrap().1))
            .collect::<Vec<String>>();

        eprintln!(
            "[+]\tChosen path through graph: {}",
            chosen_path_string.join(" -> ")
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
