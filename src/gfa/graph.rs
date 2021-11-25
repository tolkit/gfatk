use gfa::gfa::Orientation;
// common functions for graph structures
use petgraph::{
    algo::{all_simple_paths, is_cyclic_directed},
    graph::{Graph, NodeIndex},
    visit::{EdgeRef, NodeIndexable},
    Undirected,
};
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

pub struct GFAdigraph(pub Graph<usize, (Orientation, Orientation)>);

impl GFAdigraph {
    // check if graph is cyclic
    // do we want to check this? Does it need to be cyclic?
    pub fn check_is_cyclic_directed(&self) {
        let gfa_graph = &self.0;
        match is_cyclic_directed(&gfa_graph) {
            true => eprintln!("[+]\tLoaded GFA is a cyclic digraph."),
            false => panic!("Loaded GFA is not a cyclic digraph"),
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

    // get the path out as both a vec of node indices
    // and as actual node names
    pub fn find_hamiltonian_path(
        &self,
        source_target_pair: Option<(NodeIndex, NodeIndex)>,
        min_intermediate_nodes: usize,
        max_intermediate_nodes: Option<usize>,
        graph_indices: &Vec<(NodeIndex, usize)>,
    ) -> (Vec<NodeIndex>, Vec<usize>) {
        let gfa_graph = &self.0;
        // enumerating all simple paths is computationally expensive
        let mut paths = all_simple_paths::<Vec<_>, _>(
            &gfa_graph,
            // i.e. from target -> source in a digraph
            source_target_pair.unwrap().1,
            source_target_pair.unwrap().0,
            min_intermediate_nodes,
            max_intermediate_nodes,
        )
        .collect::<Vec<_>>();

        // if paths length is zero, bail out here.
        if paths.is_empty() {
            panic!("No simple paths were found.")
        }
        // sort & dedup paths and pick longest one
        paths.sort_by(|a, b| b.len().cmp(&a.len()));

        // format string just for my (and user) interest
        let chosen_path_string = &paths[0]
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
        let chosen_path_ids = &paths[0]
            .iter()
            .map(|e| graph_indices.iter().find(|y| y.0 == *e).unwrap().1)
            .collect::<Vec<_>>();

        (paths[0].to_vec(), chosen_path_ids.to_vec())
    }
}
