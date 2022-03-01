use crate::gfa::gfa::{read_gaf_to_records, GFAtk};
use crate::load::load_gfa;
use anyhow::Result;
use indexmap::IndexMap;

// what happens if the graph is not circular?

pub fn force_linear(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").unwrap();
    let fasta_header = matches.value_of("fasta-header").unwrap();
    let coverage_file = matches.value_of("coverage-file");
    let include_node_coverage = matches.is_present("include-node-coverage");

    let coverages = match coverage_file {
        Some(f) => Some(read_gaf_to_records(f)),
        None => None,
    };

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    // load gfa into graph structure
    let (graph_indices, mut gfa_graph) = gfa.into_digraph()?;

    // don't evaluate the coverage if we don't care about it
    let rel_coverage_map = match include_node_coverage {
        true => Some(gfa.gen_cov_hash(&graph_indices)?),
        false => None,
    };

    // add in the coverages from the file, if there is one.
    gfa_graph.add_coverages(&coverages, &graph_indices)?;

    let (chosen_path, chosen_path_ids, segments_not_in_path) =
        gfa_graph.all_paths_all_node_pairs(&graph_indices, coverages, rel_coverage_map.as_ref())?;

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
