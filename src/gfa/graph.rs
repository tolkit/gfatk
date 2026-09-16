use crate::gfa::gfa::GFAtk;
use crate::utils::{format_usize_to_kb, GFAGraphLookups};
use anyhow::{bail, Context, Result};
use gfa::gfa::Orientation;
use gfa::gfa::GFA;
use gfa::optfields::OptFields;
use itertools::Itertools;
use petgraph::{
    graph::{Graph, NodeIndex},
    visit::{EdgeRef, IntoNodeIdentifiers, IntoNodeReferences, NodeIndexable, NodeRef},
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
    /// Builds a linear representation of the GFA with a greedy, budget-bounded
    /// walk of the graph, rather than exhaustively enumerating every simple
    /// path between every pair of nodes (which is combinatorially explosive
    /// and prone to stack overflows on complex, repeat-rich graphs).
    ///
    /// Each node is given a "visit budget": how many times the walk may pass
    /// through it. Without coverage information every node gets a budget of
    /// 1. With `-i`, the budget is the node's coverage relative to the
    /// lowest-coverage node in the graph (`rel_coverage_map`), so repeat
    /// segments can legitimately be walked through more than once.
    ///
    /// The walk is attempted from every node, in both orientations, and the
    /// walk with the highest cumulative edge coverage is kept. Because each
    /// walk is bounded by the total visit budget and touches each edge a
    /// bounded number of times, this scales to large graphs without
    /// recursion.
    pub fn all_paths_all_node_pairs(
        &self,
        graph_indices: &GFAGraphLookups,
        rel_coverage_map: Option<&HashMap<NodeIndex, usize>>,
    ) -> Result<(Vec<(NodeIndex, Orientation)>, Vec<Vec<u8>>, String)> {
        let graph = &self.0;

        let base_budget: HashMap<NodeIndex, usize> = graph
            .node_identifiers()
            .map(|node| {
                let budget = rel_coverage_map
                    .and_then(|map| map.get(&node).copied())
                    .unwrap_or(1)
                    .max(1);
                (node, budget)
            })
            .collect();

        let mut best: Option<(Vec<(NodeIndex, Orientation)>, i64)> = None;

        for start in graph.node_identifiers() {
            for start_orientation in [Orientation::Forward, Orientation::Backward] {
                let (path, coverage) =
                    greedy_walk(graph, start, start_orientation, base_budget.clone());

                let is_better = match &best {
                    None => true,
                    Some((best_path, best_coverage)) => {
                        coverage > *best_coverage
                            || (coverage == *best_coverage && path.len() > best_path.len())
                    }
                };
                if is_better {
                    best = Some((path, coverage));
                }
            }
        }

        let final_path = best.context("There was no highest coverage path.")?;

        eprintln!("[+]\tHighest cumulative coverage path = {}", final_path.1);

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

/// Perform a single greedy, budget-bounded walk through the graph, starting
/// at `start_node` with orientation `start_orientation`.
///
/// At each step, among the outgoing edges whose orientation continues the
/// current strand (`edge.weight().0 == current_orientation`, the same rule
/// used to stitch the final path together) and whose target node still has
/// visit budget remaining, the edge with the highest coverage is followed.
/// The walk stops at a dead end. Because every node's budget is finite and
/// strictly decreases on each visit, the walk is guaranteed to terminate in
/// a bounded number of steps with no recursion, so it scales to large or
/// highly cyclic/repeat-rich graphs that would otherwise force an
/// exhaustive path search to blow the stack.
fn greedy_walk(
    graph: &Graph<Vec<u8>, (Orientation, Orientation, Option<i64>)>,
    start_node: NodeIndex,
    start_orientation: Orientation,
    mut budget: HashMap<NodeIndex, usize>,
) -> (Vec<(NodeIndex, Orientation)>, i64) {
    let start_budget = budget.entry(start_node).or_insert(1);
    if *start_budget == 0 {
        return (Vec::new(), 0);
    }
    *start_budget -= 1;

    let mut path = vec![(start_node, start_orientation)];
    let mut total_coverage: i64 = 0;
    let mut current = start_node;
    let mut current_orientation = start_orientation;

    loop {
        let next_edge = graph
            .edges_directed(current, Outgoing)
            .filter(|edge| {
                edge.weight().0 == current_orientation
                    && budget.get(&edge.target()).copied().unwrap_or(0) > 0
            })
            .max_by_key(|edge| {
                (
                    edge.weight().2.unwrap_or(0),
                    std::cmp::Reverse(edge.target().index()),
                )
            });

        let edge = match next_edge {
            Some(edge) => edge,
            None => break,
        };

        let next_node = edge.target();
        let next_orientation = edge.weight().1;
        total_coverage += edge.weight().2.unwrap_or(0);

        *budget.get_mut(&next_node).unwrap() -= 1;
        path.push((next_node, next_orientation));

        current = next_node;
        current_orientation = next_orientation;
    }

    (path, total_coverage)
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

    // a single greedy walk should actually move through the graph, following
    // real, orientation-consistent edges, rather than getting stuck immediately.
    #[test]
    fn test_greedy_walk() {
        let graph = make_graph();

        let budget: HashMap<NodeIndex, usize> = graph.0.node_indices().map(|n| (n, 1)).collect();

        let (path, _coverage) =
            greedy_walk(&graph.0, NodeIndex::new(0), Orientation::Forward, budget);

        assert_eq!(path[0], (NodeIndex::new(0), Orientation::Forward));
        assert!(path.len() > 1);
    }

    // this is a regression test for the bug this greedy walk replaced: exhaustively
    // enumerating every simple path between every pair of nodes via recursion blew
    // the stack on large, cyclic/repeat-rich graphs (like complex mitochondria
    // assemblies). this builds a large graph with many cycles and checks that
    // linearising it completes without overflowing the stack.
    #[test]
    fn test_large_cyclic_graph_does_not_overflow() {
        let mut graph = Graph::<Vec<u8>, (Orientation, Orientation, Option<i64>)>::new();

        let n = 2000;
        let nodes: Vec<NodeIndex> = (0..n)
            .map(|i| graph.add_node(i.to_string().into_bytes()))
            .collect();

        for i in 0..n - 1 {
            graph.add_edge(
                nodes[i],
                nodes[i + 1],
                (Orientation::Forward, Orientation::Forward, Some(1)),
            );
            // add a few edges back to earlier nodes to create lots of cycles,
            // mimicking the repeat structure of a complex mitochondrial graph.
            if i >= 3 {
                graph.add_edge(
                    nodes[i],
                    nodes[i - 3],
                    (Orientation::Forward, Orientation::Forward, Some(1)),
                );
            }
        }

        let gfa_graph = GFAdigraph(graph);

        let lookup = GFAGraphLookups(
            nodes
                .iter()
                .enumerate()
                .map(|(i, &node_index)| crate::utils::GFAGraphPair {
                    node_index,
                    seg_id: i.to_string().into_bytes(),
                })
                .collect(),
        );

        let result = gfa_graph.all_paths_all_node_pairs(&lookup, None);
        assert!(result.is_ok());
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

        // generate the path
        let (path, _not_in_path, fasta_header) =
            graph.all_paths_all_node_pairs(&lookup, Some(&map)).unwrap();

        // exhaustive search (see git history of this test) established 2625 as the
        // true maximum cumulative coverage achievable under this coverage map, via
        // one of two nodewise-equivalent optimal routes. the greedy walk explores a
        // superset of starting points/orientations, so it may land on either of
        // those two routes, or another route of equal coverage - what matters is
        // that it actually finds the true optimum and respects every node's visit
        // budget (node2 and node3 twice, everything else once), not the exact node
        // order.
        assert!(fasta_header.contains("coverage=2625"));

        let mut visit_counts: HashMap<NodeIndex, usize> = HashMap::new();
        for (node, _orientation) in &path {
            *visit_counts.entry(*node).or_insert(0) += 1;
        }
        for (node, expected_visits) in &map {
            assert_eq!(
                visit_counts.get(node).copied().unwrap_or(0),
                *expected_visits
            );
        }
    }
}
