use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;
use indexmap::IndexMap;

pub fn force_linear(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").unwrap();
    let fasta_header = matches.value_of("fasta-header").unwrap();

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    // load gfa into graph structure
    let (graph_indices, gfa_graph) = gfa.into_digraph();

    // check if it's cyclic + directed
    gfa_graph.check_is_cyclic_directed();

    // get the counts of number of links per node
    let node_edge_counts = gfa_graph.get_edge_counts();

    // if there are any nodes with fewer than 2 connections
    // get rid of these, as they are not part of the cycle.
    let filtered_nodes = node_edge_counts
        .iter()
        .filter(|e| e.1 > 1)
        .collect::<Vec<&(usize, i32)>>();

    // for starting node and ending node
    // if the graph is cyclic but not filtered
    // might get stuck in the wrong subgraph
    let min_intermediate_nodes = filtered_nodes.len() - 2;

    // no max intermediate nodes
    let max_intermediate_nodes = None;

    let source_target_pair = gfa_graph.get_source_target_pair(&graph_indices, node_edge_counts);

    if source_target_pair.is_none() {
        panic!("Did not find a pair of adjacent segments, each themselves with connections.");
    }

    let (chosen_path, chosen_path_ids) = gfa_graph.find_hamiltonian_path(
        source_target_pair,
        min_intermediate_nodes,
        max_intermediate_nodes,
        &graph_indices,
    );

    let sorted_chosen_path_overlaps =
        gfa.determine_path_overlaps(&chosen_path, &graph_indices, &chosen_path_ids);

    // merge constrand on the ids
    // as order is critical, we use an IndexMap
    let mut merged_sorted_chosen_path_overlaps = IndexMap::new();
    for (id, orientation, overlap, side) in sorted_chosen_path_overlaps {
        merged_sorted_chosen_path_overlaps
            .entry(id)
            .or_insert(Vec::new())
            .push((orientation, overlap, side));
    }

    gfa.print_path_to_fasta(merged_sorted_chosen_path_overlaps, fasta_header);

    Ok(())
}
