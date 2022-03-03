use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;
use anyhow::{Context, Result};
use indexmap::IndexMap;

// what happens if the graph is not circular?

pub fn force_linear(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").context("No gfa file specified")?;
    let fasta_header = matches
        .value_of("fasta-header")
        .context("No fasta header specified")?;
    let include_node_coverage = matches.is_present("include-node-coverage");

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    // load gfa into graph structure
    let (graph_indices, gfa_graph) = gfa.into_digraph()?;

    // don't evaluate the coverage if we don't care about it
    let rel_coverage_map = match include_node_coverage {
        true => Some(gfa.gen_cov_hash(&graph_indices)?),
        false => None,
    };

    let (chosen_path, chosen_path_ids, segments_not_in_path) =
        gfa_graph.all_paths_all_node_pairs(&graph_indices, rel_coverage_map.as_ref())?;

    let sorted_chosen_path_overlaps =
        gfa.determine_path_overlaps(&chosen_path, graph_indices, &chosen_path_ids)?;

    // merge constrand on the ids
    // as order is critical, we use an IndexMap
    let mut merged_sorted_chosen_path_overlaps = IndexMap::new();
    for (id, orientation, overlap, side) in sorted_chosen_path_overlaps {
        merged_sorted_chosen_path_overlaps
            .entry(id)
            .or_insert(Vec::new())
            .push((orientation, overlap, side));
    }

    gfa.print_path_to_fasta(
        merged_sorted_chosen_path_overlaps,
        fasta_header,
        segments_not_in_path,
    )?;

    Ok(())
}
