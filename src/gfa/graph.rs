use crate::gfa::gfa::GFAtk;
use crate::utils::{format_usize_to_kb, GFAGraphLookups};
use anyhow::{bail, Context, Result};
use gfa::gfa::Orientation;
use gfa::gfa::GFA;
use gfa::optfields::OptFields;
use itertools::Itertools;
use petgraph::{
    graph::{Graph, IndexType, NodeIndex},
    visit::{EdgeRef, IntoNodeIdentifiers, IntoNodeReferences, NodeIndexable, NodeRef},
    Directed,
    Direction::Outgoing,
    Undirected,
};
use std::collections::HashMap;
use std::collections::HashSet;

/// A wrapper of petgraph's undirected `Graph` struct, applied to a GFA. No weights.
pub struct GFAungraph(pub Graph<usize, (), Undirected>);

impl GFAungraph {
    /// The algorithm called in `gfatk extract`.
    ///
    /// The number of iterations of searching for neighbouring nodes can be modified.
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

/// A wrapper of petgraph's directed `Graph` struct, applied to a GFA. The edge weights included are the `Orientation`'s of the adjacent segments, and the coverage of this edge.
pub struct GFAdigraph(pub Graph<usize, (Orientation, Orientation, Option<i64>)>);

impl GFAdigraph {
    /// The main function called from `gfatk dot`.
    ///
    /// It is a somewhat modified, simplified version of this:
    /// <https://docs.rs/petgraph/latest/src/petgraph/dot.rs.html#1-349>
    ///
    /// Generating a DOT language output of a GFA file.
    pub fn dot(&self, gfa: GFAtk) -> Result<()> {
        let gfa_graph = &self.0;
        static INDENT: &str = "    ";

        println!("digraph GFA {{");
        // print nodes
        for node in gfa_graph.node_references() {
            let e = gfa_graph.to_index(node.id());
            let w = node.weight();
            let meta = gfa.node_seq_len_and_cov(*w)?;
            println!(
                // see https://stackoverflow.com/questions/20516143/graphviz-dot-different-fontsizes-in-same-label
                "{}{} [ label = <<FONT POINT-SIZE=\'20\'>{}</FONT><br/><FONT POINT-SIZE=\'10\'>L: {}</FONT><br/><FONT POINT-SIZE=\'10\'>C: {}</FONT>> ];",
                INDENT, e, w, format_usize_to_kb(meta.0), meta.1
            );
        }
        // print edges
        for edge in gfa_graph.edge_references() {
            let from = gfa_graph.to_index(edge.source());
            let to = gfa_graph.to_index(edge.target());
            let from_o = edge.weight().0;
            let to_o = edge.weight().1;

            let arrowhead_shape = match to_o {
                Orientation::Forward => "ornormal",
                Orientation::Backward => "olnormal",
            };

            let ec = edge
                .weight()
                .2
                .context(format!("No edge weight for edge {:?}", edge))?;

            println!("{}{} -> {} [ label = \"  {}  \" taillabel = \"  {}  \" headlabel = \"  {}  \" arrowhead = \"{}\" ];", 
                INDENT,
                from,
                to,
                ec,
                from_o,
                to_o,
                arrowhead_shape
            );
        }

        println!("}}");

        Ok(())
    }
    // we want weakly connected components, as there may only be an edge in one
    // orientation (perhaps unlikely... but still)

    /// Split the GFA digraph into subgraphs which are the weakly connected components of the graph.
    ///
    /// Taken from <https://github.com/Qiskit/retworkx/blob/79900cf8da0c0665ac5ce1ccb0f57373434b14b8/src/connectivity/mod.rs>
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

    // core algorithm for gfatk linear
    // iterate over all pairs of nodes to find the paths
    // then search though all these paths to find the
    // longest which doesn't violate sequence orientation

    // need to emit segment ID's that did not make it into the
    // longest path here.

    // return the path, the path ID's, id's not in the path, and the
    // `String` of a fasta header.

    // -> (Vec<NodeIndex>, Vec<usize>, Vec<usize>)

    /// The main function called from `gfatk linear`.
    ///
    /// This function will generate the longest path through the GFA, by
    /// filtering the output of `all_paths`, and choosing the path with
    /// the highest cumulative edge coverage.
    pub fn all_paths_all_node_pairs(
        &self,
        graph_indices: &GFAGraphLookups,
        rel_coverage_map: Option<&HashMap<NodeIndex, usize>>,
    ) -> Result<(Vec<NodeIndex>, Vec<usize>, Vec<usize>, String)> {
        let graph = &self.0;
        let nodes = graph.node_identifiers();
        let node_pairs = nodes.permutations(2).collect::<Vec<_>>();
        // eprintln!("{:?}", node_pairs);

        let all_paths: Result<Vec<_>> = node_pairs
            .into_iter()
            .map(|pair| {
                let path = all_paths(&graph, pair[0], pair[1], rel_coverage_map);
                path
            })
            .collect();

        // make the set of legal paths through the GFA

        let mut valid_paths = Vec::new();
        // iterate over the paths
        for paths in all_paths? {
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
                        (Orientation, Orientation, Option<i64>),
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
                        (Orientation, Orientation, Option<i64>),
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
        // maybe add an expected number of segments?

        let final_path = {
            // now iterate over all the paths again
            // let mut index = 0;
            // iterate over all the paths of this length
            // keep track of lengths

            // let mut final_path = Vec::new();
            // let mut final_coverage = 0;
            // push path and coverage into map
            // don't care about memory allocations for the moment.
            let mut map = HashMap::new();

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

                map.insert(path, path_coverage);
            }

            let highest_coverage_path_op =
                map.iter().max_by(|a, b| a.1.cmp(&b.1)).map(|(k, v)| (k, v));

            // explicit error out here
            let highest_coverage_path = match highest_coverage_path_op {
                Some(p) => p,
                None => bail!("There was no highest coverage path."),
            };

            eprintln!(
                "[+]\tHighest cumulative coverage path = {}",
                highest_coverage_path.1
            );
            (highest_coverage_path.0.to_vec(), *highest_coverage_path.1)
        };

        let mut chosen_path_string = Vec::new();
        let final_path_node_pairs = final_path.0.windows(2);
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

            // no spaces between the formatted strings
            if index == 0 {
                chosen_path_string.push(format!("{}{},{}{}", from, from_orient, to, to_orient));
            } else {
                chosen_path_string.push(format!("->{}{}", to, to_orient));
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
            .0
            .iter()
            // .1 needed
            .map(|e| graph_indices.node_index_to_seg_id(*e))
            .collect::<Result<Vec<_>>>();

        // make a vector of nodes not in the final path
        // these will be passed later and printed to a fasta.
        let final_path_set: HashSet<_> = final_path.0.iter().collect();
        let nodes: Vec<_> = graph.node_identifiers().collect();
        let difference: Vec<_> = nodes
            .into_iter()
            .filter(|item| !final_path_set.contains(item))
            .collect();

        let difference_ids: Result<Vec<usize>> = difference
            .iter()
            .map(|e| {
                let seg_id = graph_indices.node_index_to_seg_id(*e);
                seg_id
            })
            .collect();

        // the fasta header should contain the tool, path information, and coverage
        let fasta_header = format!(
            "gfatk_linear:path={}:coverage={}",
            chosen_path_string.join(""),
            final_path.1
        );

        Ok((
            final_path.0.to_vec(),
            // error handling a bit annoying here.
            chosen_path_ids?,
            difference_ids?,
            fasta_header,
        ))
    }

    /// Simple wrapper of `Graph.node_count()` in petgraph.
    pub fn node_count(&self) -> usize {
        let gfa_graph = &self.0;

        gfa_graph.node_count()
    }

    /// Simple wrapper of `Graph.edge_count()` in petgraph.
    pub fn edge_count(&self) -> usize {
        let gfa_graph = &self.0;

        gfa_graph.edge_count()
    }

    /// Trim a graph to include only nodes connected to two or more other nodes.
    ///
    /// This algorithm will loop for as long as the longest branch in the GFA yields a segment connected to only a single node.
    pub fn trim(&self, graph_indices: GFAGraphLookups) -> Vec<usize> {
        let gfa_graph = &self.0;

        let mut all_nodes = HashSet::new();
        // get all node indices into a hashset
        for (node_index, _) in gfa_graph.node_references() {
            all_nodes.insert(node_index);
        }

        // initiate new hashset for the nodes we remove
        let mut removed_nodes = HashSet::new();
        // keep track of the number of removed nodes in a vector
        let mut track_removed_nodes = Vec::new();
        // index for the above vector, keeping track of iterations
        let mut index = 0;
        loop {
            // iterate over the nodes
            for (node_index, _) in gfa_graph.node_references() {
                // how many neighbours for this particular node?
                let neighbours = gfa_graph.neighbors(node_index).collect::<HashSet<_>>();

                // if there are fewer than two neighbours
                // OR the difference between neighbours & removed nodes == 1
                if neighbours.len() < 2
                    || neighbours
                        .difference(&removed_nodes)
                        .collect::<HashSet<_>>()
                        .len()
                        == 1
                {
                    removed_nodes.insert(node_index);
                }
            }
            // push the length
            // if the previous length is the same as the current,
            // there are no more nodes to delete.
            track_removed_nodes.push(removed_nodes.len());
            if track_removed_nodes.get(index).unwrap()
                == track_removed_nodes.get(index - 1).unwrap_or(&0)
            {
                break;
            }
            index += 1;
        }
        // print for user info
        for el in &removed_nodes {
            let seg_id = graph_indices.node_index_to_seg_id(*el).unwrap();
            eprintln!("[+]\tRemoved segment {} from GFA.", seg_id);
        }

        all_nodes
            .difference(&removed_nodes)
            .map(|e| graph_indices.node_index_to_seg_id(*e).unwrap())
            .collect::<Vec<_>>()
    }
}

/// A function generic over certain types of `Directed` petgraph `Graph`s.
///
/// Given a graph, a start node, an end node, and optionally a map of the coverage of each node, compute all simple paths between these nodes.
///
/// Modified from: <https://github.com/Ninjani/rosalind/blob/e22ecf2c9f0935d970b137684029957c0850d63f/t_ba11b/src/lib.rs>
pub fn all_paths<T, U, Ix: IndexType>(
    graph: &Graph<T, U, Directed, Ix>,
    start_node: NodeIndex<Ix>,
    end_node: NodeIndex<Ix>,
    rel_coverage_map: Option<&HashMap<NodeIndex<Ix>, usize>>,
) -> Result<Vec<Vec<NodeIndex<Ix>>>> {
    match rel_coverage_map {
        Some(cov_map) => {
            // for the set of visited nodes
            let mut visited = HashMap::new();
            visited.insert(start_node, 1);
            recursive_path_finder_incl_coverage(graph, start_node, end_node, &mut visited, cov_map)
        }
        None => {
            let mut visited = HashSet::new();
            visited.insert(start_node);
            recursive_path_finder_no_coverage(graph, start_node, end_node, &mut visited)
        }
    }
}

/// Function called by `all_paths` where a `HashMap` is supplied instead of a
/// `HashSet` in order to keep track of how many times a segment/node has been
/// passed in a path.
///
/// Causes a stack overflow if the variance in node coverages
/// is too high.
fn recursive_path_finder_incl_coverage<T, U, Ix: IndexType>(
    graph: &Graph<T, U, Directed, Ix>,
    start_node: NodeIndex<Ix>,
    end_node: NodeIndex<Ix>,
    visited: &mut HashMap<NodeIndex<Ix>, usize>,
    rel_coverage_map: &HashMap<NodeIndex<Ix>, usize>,
) -> Result<Vec<Vec<NodeIndex<Ix>>>> {
    // if the start node is the same as the end
    // the path is just to the end node
    if start_node == end_node {
        Ok(vec![vec![end_node]])
    } else {
        let mut paths = Vec::new();
        for edge in graph.edges_directed(start_node, Outgoing) {
            let next_node = edge.target();

            let test = *rel_coverage_map.get(&next_node).unwrap();

            if !visited.contains_key(&next_node) || !(*visited.get(&next_node).unwrap() == test) {
                *visited.entry(next_node).or_insert(0) += 1;
                let descendant_paths = recursive_path_finder_incl_coverage(
                    graph,
                    next_node,
                    end_node,
                    visited,
                    rel_coverage_map,
                )?;
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
        Ok(paths)
    }
}

/// A safer and more reliable alternative to `recursive_path_finder_incl_coverage` where
/// a `HashSet` determines whether a segment/node has already been seen or not.
fn recursive_path_finder_no_coverage<T, U, Ix: IndexType>(
    graph: &Graph<T, U, Directed, Ix>,
    start_node: NodeIndex<Ix>,
    end_node: NodeIndex<Ix>,
    visited: &mut HashSet<NodeIndex<Ix>>,
) -> Result<Vec<Vec<NodeIndex<Ix>>>> {
    if start_node == end_node {
        Ok(vec![vec![end_node]])
    } else {
        let mut paths = Vec::new();
        for edge in graph.edges_directed(start_node, Outgoing) {
            let next_node = edge.target();
            if !visited.contains(&next_node) {
                visited.insert(next_node);
                let descendant_paths =
                    recursive_path_finder_no_coverage(graph, next_node, end_node, visited)?;
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
        Ok(paths)
    }
}

/// Returns a subgraph GFA that only contains elements with the provided segment names.
///
/// Taken from <https://github.com/chfi/rs-gfa-utils/blob/master/src/subgraph.rs>
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

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
