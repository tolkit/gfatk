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
pub struct GFAungraph(pub Graph<Vec<u8>, (), Undirected>);

impl GFAungraph {
    /// The algorithm called in `gfatk extract`.
    ///
    /// The number of iterations of searching for neighbouring nodes can be modified.
    ///
    /// It's a naive algorithm, but it's fast enough for our purposes.
    pub fn recursive_search(
        &self,
        sequence_id: Vec<Vec<u8>>,
        iterations: i32,
        collect_sequence_names: Vec<NodeIndex>,
        graph_indices: GFAGraphLookups,
    ) -> Result<Vec<Vec<u8>>> {
        let gfa_graph = &self.0;

        let sequence_id_d = sequence_id
            .iter()
            .map(|e| String::from_utf8_lossy(e).to_string())
            .join(", ");

        eprintln!(
            "[+]\tRecursively searching around nodes {} at depth {}",
            sequence_id_d, iterations
        );

        let mut collect_sequence_set: HashSet<_> = collect_sequence_names.iter().copied().collect();

        for _ in 0..iterations {
            // collect all the neighbours of all the current node indices
            for index in collect_sequence_set.clone() {
                for c in gfa_graph.neighbors(index) {
                    // could possibly add a conditional in here.
                    collect_sequence_set.insert(c);
                }
            }
        }

        // turn node indexes into sequence ID's
        collect_sequence_set
            .into_iter()
            .map(|index| graph_indices.node_index_to_seg_id(index))
            .collect()
    }
}

// weights are the orientations, used at various points, and an optional
// coverage weight, used in gfatk linear.
// GFA's should always specify Links in a specific direction..?
// so digraphs should be where all the functionality lies.

/// A wrapper of petgraph's directed `Graph` struct, applied to a GFA. The edge weights included are the `Orientation`'s of the adjacent segments, and the coverage of this edge.
pub struct GFAdigraph(pub Graph<Vec<u8>, (Orientation, Orientation, Option<i64>)>);

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
            let w_d = std::str::from_utf8(&w)?;
            let meta = gfa.node_seq_len_and_cov(w.to_vec())?;
            println!(
                // see https://stackoverflow.com/questions/20516143/graphviz-dot-different-fontsizes-in-same-label
                "{}{} [ label = <<FONT POINT-SIZE=\'20\'>{}</FONT><br/><FONT POINT-SIZE=\'10\'>L: {}</FONT><br/><FONT POINT-SIZE=\'10\'>C: {}</FONT>> ];",
                INDENT, e, w_d, format_usize_to_kb(meta.0), meta.1
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
    ) -> Result<Vec<Vec<Vec<u8>>>> {
        let graph = &self.0;
        let mut seen: HashSet<NodeIndex> = HashSet::with_capacity(graph.node_count());
        let mut out_vec: Vec<Vec<Vec<u8>>> = Vec::new();

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
                    .collect::<Result<Vec<Vec<u8>>>>();

                out_vec.push(x?);

                seen.extend(bfs_seen);
            }
        }
        Ok(out_vec)
    }

    /// The main function called from `gfatk linear`.
    ///
    /// This function will generate the longest path through the GFA, by
    /// filtering the output of `all_paths`, and choosing the path with
    /// the highest cumulative edge coverage.
    pub fn all_paths_all_node_pairs(
        &self,
        graph_indices: &GFAGraphLookups,
        rel_coverage_map: Option<&HashMap<NodeIndex, usize>>,
    ) -> Result<(Vec<(NodeIndex, Orientation)>, Vec<Vec<u8>>, String)> {
        let graph = &self.0;
        let nodes = graph.node_identifiers();

        let all_paths: Result<Vec<_>> = nodes
            .permutations(2)
            .enumerate()
            .map(|(index, pair)| all_paths(graph, pair[0], pair[1], rel_coverage_map, index))
            .collect();

        // this is kind of annoying to add this, but I could not think
        // of another way to overcome the eprint!() in `all_paths()`
        eprintln!();
        // make the set of legal paths through the GFA

        let mut valid_paths = Vec::new();
        // iterate over the paths
        for paths in all_paths? {
            // iterate over each path
            for path in paths {
                // I think easiest to just append all paths length = 2
                // does this make sense? if there are paths of longer length
                // than two, these will *always* be the highest coverage
                // so no need to filter later.
                if path.len() == 2 {
                    valid_paths.push(path.clone());
                }
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

                    // 1. we can then compare the orientation of the 'to' Orientation
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
                    let do_keep = keep_vec.iter().any(|e| *e);

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
        valid_paths.sort_by_key(|b| std::cmp::Reverse(b.len()));
        valid_paths.dedup();

        // now make the final path
        let final_path = {
            // push path and coverage into map
            // don't care about memory allocations for the moment.
            let mut map = HashMap::new();

            // test this please.
            'outer: for path in &valid_paths {
                let mut path_coverage = 0;

                // test this
                let mut path_orientations = Vec::new();

                // try a different method
                let mut fi = 0;
                let mut se = 1;

                let path_len = path.len();
                for _ in 0..path_len {
                    let node_1 = match path.get(fi) {
                        Some(p) => *p,
                        None => continue,
                    };
                    let node_2 = match path.get(se) {
                        Some(p) => *p,
                        None => continue,
                    };

                    let pair_connecting = &mut graph.edges_connecting(node_1, node_2);
                    // from node 1 to node 2, we just choose the first edge
                    // as I think it doesn't matter which edge is chosen (they will have the same coverage in MBG)
                    if fi == 0 {
                        let pair_weight = pair_connecting
                            .next()
                            .with_context(|| {
                                format!("No connecting edges from {:?} to {:?}", node_1, node_2)
                            })?
                            .weight();

                        let node_1_orientation = pair_weight.0;
                        let node_2_orientation = pair_weight.1;
                        // this will be path_orientations[0]
                        path_orientations.push(node_1_orientation);
                        // this will be path_orientations[1]
                        path_orientations.push(node_2_orientation);

                        let coverage = pair_weight.2;
                        if let Some(c) = coverage {
                            path_coverage += c;
                        }
                    } else {
                        // this might be wrong...
                        let prev_orientation = path_orientations[fi];
                        let pair_weight =
                            pair_connecting.find(|e| e.weight().0 == prev_orientation);

                        // if:
                        // A -> B (orientation)
                        // is not equal to
                        // B (orientation) -> C
                        // we skip this path.
                        if pair_weight.is_none() {
                            continue 'outer;
                        }
                        let node_2_orientation = pair_weight.unwrap().weight().1;
                        path_orientations.push(node_2_orientation);

                        let coverage = pair_weight.unwrap().weight().2;
                        if let Some(c) = coverage {
                            path_coverage += c;
                        }
                    }

                    // increment the indices.
                    fi += 1;
                    se += 1;
                }

                let path_orientation_tuple = path
                    .iter()
                    .zip(path_orientations.iter())
                    .map(|(e, f)| (*e, *f))
                    .collect::<Vec<_>>();

                map.insert(path_orientation_tuple, path_coverage);
            }

            let highest_coverage_path_op =
                map.iter().max_by(|a, b| a.1.cmp(b.1)).map(|(k, v)| (k, v));

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
            let from = pair[0].0;
            let from_orient = pair[0].1;
            let to = pair[1].0;
            let to_orient = pair[1].1;

            // get segment ID from Node Indices
            let from_inner = graph_indices.node_index_to_seg_id(from)?;
            let to_inner = graph_indices.node_index_to_seg_id(to)?;
            let from = std::str::from_utf8(&from_inner).unwrap();
            let to = std::str::from_utf8(&to_inner).unwrap();

            // no spaces between the formatted strings
            if index == 0 {
                chosen_path_string.push(format!("{}{},{}{}", from, from_orient, to, to_orient));
            } else {
                chosen_path_string.push(format!(",{}{}", to, to_orient));
            }
        }

        eprintln!(
            "[+]\tChosen path through graph: {}",
            chosen_path_string.join("")
        );

        // make a vector of nodes not in the final path
        // these will be passed later and printed to a fasta.
        let final_path_set: HashSet<_> = final_path.0.iter().map(|(e, _f)| *e).collect();

        let difference: Vec<_> = graph
            .node_identifiers()
            .filter(|item| !final_path_set.contains(item))
            .collect();

        let difference_ids: Result<Vec<Vec<u8>>> = difference
            .iter()
            .map(|e| graph_indices.node_index_to_seg_id(*e))
            .collect();

        // the fasta header should contain the tool, path information, and coverage
        let fasta_header = format!(
            "gfatk_linear:path={}:coverage={}",
            chosen_path_string.join(""),
            final_path.1
        );

        Ok((final_path.0.to_vec(), difference_ids?, fasta_header))
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
    pub fn trim(&self, graph_indices: GFAGraphLookups) -> Vec<Vec<u8>> {
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
            eprintln!(
                "[+]\tRemoved segment {} from GFA.",
                std::str::from_utf8(&seg_id).unwrap()
            );
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
    index: usize,
) -> Result<Vec<Vec<NodeIndex<Ix>>>> {
    match rel_coverage_map {
        Some(cov_map) => {
            // for the set of visited nodes
            let mut visited = HashMap::new();
            visited.insert(start_node, 1);
            // a counter for recursion depth.
            let depth = 0;
            match recursive_path_finder_incl_coverage(
                graph,
                start_node,
                end_node,
                &mut visited,
                cov_map,
                depth,
                index,
            ) {
                Some(p) => p,
                None => {
                    // copy of the chunk below!
                    // so if we go past the self imposed stack limit
                    // we default to our other (not including coverage) method.
                    let mut visited = HashSet::new();
                    visited.insert(start_node);
                    recursive_path_finder_no_coverage(graph, start_node, end_node, &mut visited)
                }
            }
        }
        None => {
            let mut visited = HashSet::new();
            visited.insert(start_node);
            recursive_path_finder_no_coverage(graph, start_node, end_node, &mut visited)
        }
    }
}

/// A recursion depth limit, so we don't hit a stack overflow
/// and instead, abort and call another function.
///
/// Why is it 1000? Seemed sensible, and that's what python's is.
const MAX_RECURSION_DEPTH: usize = 1000;

/// Function called by `all_paths` where a `HashMap` is supplied instead of a
/// `HashSet` in order to keep track of how many times a segment/node has been
/// passed in a path.
///
/// Should be no more stack overflows.
fn recursive_path_finder_incl_coverage<T, U, Ix: IndexType>(
    graph: &Graph<T, U, Directed, Ix>,
    start_node: NodeIndex<Ix>,
    end_node: NodeIndex<Ix>,
    visited: &mut HashMap<NodeIndex<Ix>, usize>,
    rel_coverage_map: &HashMap<NodeIndex<Ix>, usize>,
    depth: usize,
    index: usize,
) -> Option<Result<Vec<Vec<NodeIndex<Ix>>>>> {
    if depth > MAX_RECURSION_DEPTH {
        eprint!(
            "\r[-]\tRecursion depth limit ({}) exceeded in permutation {}. Switching to default path finder.",
            MAX_RECURSION_DEPTH,
            index + 1
        );
        return None;
    }
    // if the start node is the same as the end
    // the path is just to the end node
    if start_node == end_node {
        Some(Ok(vec![vec![end_node]]))
    } else {
        let mut paths = Vec::new();
        for edge in graph.edges_directed(start_node, Outgoing) {
            let next_node = edge.target();

            let test = *rel_coverage_map.get(&next_node).unwrap();

            if !visited.contains_key(&next_node) || *visited.get(&next_node).unwrap() != test {
                *visited.entry(next_node).or_insert(0) += 1;
                let descendant_paths = match recursive_path_finder_incl_coverage(
                    graph,
                    next_node,
                    end_node,
                    visited,
                    rel_coverage_map,
                    depth + 1,
                    index,
                ) {
                    // can I get rid of this unwrap? is it safe?
                    Some(p) => p.unwrap(),
                    None => return None,
                };
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
        Some(Ok(paths))
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
    gfa: &GFA<Vec<u8>, T>,
    segment_names: Vec<Vec<u8>>,
) -> GFA<Vec<u8>, T> {
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
        .filter(|p| p.iter().any(|(s, _)| segment_names.contains(&s.to_vec())))
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

    use super::*;

    // we want to make a test graph to play with
    // make the inner graph representation of:
    // ./examples/mito_NC_037304.1.MZ323108.1.fasta.BOTH.HiFiMapped.bam.filtered.1k.gfa

    fn make_graph() -> GFAdigraph {
        let mut graph = Graph::<Vec<u8>, (Orientation, Orientation, Option<i64>)>::new();

        // node weights are usize
        let node0 = graph.add_node("0".as_bytes().to_vec());
        let node1 = graph.add_node("1".as_bytes().to_vec());
        let node2 = graph.add_node("2".as_bytes().to_vec());
        let node3 = graph.add_node("3".as_bytes().to_vec());
        let node4 = graph.add_node("4".as_bytes().to_vec());
        let node5 = graph.add_node("5".as_bytes().to_vec());

        // we create the following graph
        //
        //  0 <-----> 3 <-----> 1
        //    \     / ˄ \     /
        //      \ /   |   \ /
        //      / \   |   / \
        //    /     \ ˅ /     \
        //  5 <-----> 2 <-----> 4
        //

        graph.extend_with_edges(&[
            (
                node0,
                node3,
                (Orientation::Backward, Orientation::Backward, Some(379)),
            ),
            (
                node0,
                node2,
                (Orientation::Forward, Orientation::Backward, Some(338)),
            ),
            (
                node1,
                node3,
                (Orientation::Backward, Orientation::Backward, Some(380)),
            ),
            (
                node1,
                node2,
                (Orientation::Forward, Orientation::Backward, Some(374)),
            ),
            (
                node2,
                node4,
                (Orientation::Backward, Orientation::Forward, Some(347)),
            ),
            (
                node2,
                node5,
                (Orientation::Backward, Orientation::Forward, Some(399)),
            ),
            (
                node2,
                node1,
                (Orientation::Forward, Orientation::Backward, Some(374)),
            ),
            (
                node2,
                node0,
                (Orientation::Forward, Orientation::Backward, Some(338)),
            ),
            (
                node3,
                node5,
                (Orientation::Backward, Orientation::Backward, Some(397)),
            ),
            (
                node3,
                node4,
                (Orientation::Backward, Orientation::Backward, Some(349)),
            ),
            (
                node3,
                node1,
                (Orientation::Forward, Orientation::Forward, Some(380)),
            ),
            (
                node3,
                node0,
                (Orientation::Forward, Orientation::Forward, Some(379)),
            ),
            (
                node4,
                node2,
                (Orientation::Backward, Orientation::Forward, Some(347)),
            ),
            (
                node4,
                node3,
                (Orientation::Forward, Orientation::Forward, Some(349)),
            ),
            (
                node5,
                node2,
                (Orientation::Backward, Orientation::Forward, Some(399)),
            ),
            (
                node5,
                node3,
                (Orientation::Forward, Orientation::Forward, Some(397)),
            ),
        ]);

        GFAdigraph(graph)
    }

    // there are 6 nodes in this graph
    #[test]
    fn test_node_count() {
        let graph = make_graph();

        assert_eq!(graph.node_count(), 6);
    }

    // there are 16 edges in this graph (incl +/- orientations)
    #[test]
    fn test_edge_count() {
        let graph = make_graph();

        assert_eq!(graph.edge_count(), 16);
    }

    // there are two possible paths between node indexes 0 and 2
    #[test]
    fn test_path_generation() {
        let graph = make_graph();

        let paths = all_paths(&graph.0, NodeIndex::new(0), NodeIndex::new(2), None, 0).unwrap();

        // there should be two paths
        let path1: Vec<NodeIndex> = vec![NodeIndex::new(0), NodeIndex::new(2)];
        let path2: Vec<NodeIndex> = vec![
            NodeIndex::new(0),
            NodeIndex::new(3),
            NodeIndex::new(1),
            NodeIndex::new(2),
        ];
        let path3: Vec<NodeIndex> = vec![
            NodeIndex::new(0),
            NodeIndex::new(3),
            NodeIndex::new(4),
            NodeIndex::new(2),
        ];
        let path4: Vec<NodeIndex> = vec![
            NodeIndex::new(0),
            NodeIndex::new(3),
            NodeIndex::new(5),
            NodeIndex::new(2),
        ];

        assert!(paths.contains(&path1));
        assert!(paths.contains(&path2));
        assert!(paths.contains(&path3));
        assert!(paths.contains(&path4));
    }

    //
    #[test]
    fn test_path_generation_incl_node_cov() {
        let graph = make_graph();

        let mut map: HashMap<NodeIndex, usize> = HashMap::new();

        // we can provide a map to say we want to visit certain nodes twice
        map.insert(NodeIndex::new(0), 1);
        map.insert(NodeIndex::new(1), 1);
        map.insert(NodeIndex::new(2), 2);
        map.insert(NodeIndex::new(3), 2);
        map.insert(NodeIndex::new(4), 1);
        map.insert(NodeIndex::new(5), 1);

        let lookup = GFAGraphLookups(vec![
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(0),
                seg_id: "4".as_bytes().to_vec(),
            },
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(2),
                seg_id: "6".as_bytes().to_vec(),
            },
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(5),
                seg_id: "9".as_bytes().to_vec(),
            },
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(3),
                seg_id: "7".as_bytes().to_vec(),
            },
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(1),
                seg_id: "5".as_bytes().to_vec(),
            },
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(4),
                seg_id: "8".as_bytes().to_vec(),
            },
        ]);

        // generate the paths
        let paths = graph.all_paths_all_node_pairs(&lookup, Some(&map));

        // either this path
        let longest_path1: Vec<NodeIndex> = vec![
            NodeIndex::new(2),
            NodeIndex::new(5),
            NodeIndex::new(3),
            NodeIndex::new(1),
            NodeIndex::new(2),
            NodeIndex::new(4),
            NodeIndex::new(3),
            NodeIndex::new(0),
        ];

        // or this path
        let longest_path2: Vec<NodeIndex> = vec![
            NodeIndex::new(2),
            NodeIndex::new(4),
            NodeIndex::new(3),
            NodeIndex::new(1),
            NodeIndex::new(2),
            NodeIndex::new(5),
            NodeIndex::new(3),
            NodeIndex::new(0),
        ];

        // will be chosen
        let both = vec![longest_path1, longest_path2];

        let path = &paths.unwrap().0.iter().map(|(a, _)| *a).collect::<Vec<_>>();

        assert!(both.contains(path));
    }
}
