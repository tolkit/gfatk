use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Context, Result};
use indexmap::IndexMap;

/// Force a linear representation of the GFA.
///
/// This function finds *all* legal paths through a GFA, and returns the longest path, with the highest cumulative edge coverage.
///
/// If the `-i` option is included, node coverages are taken into account, and paths are created with nodes appearing in the final path the number of times they relatively occur according to coverage information.
///
/// For example:
/// ```bash
/// # simple
/// gfatk linear in.gfa > out.fasta
/// # account for node coverage
/// gfatk -i linear in.gfa > out.fasta
/// ```
pub fn force_linear(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("GFA");
    let fasta_header = matches
        .value_of("fasta-header")
        .context("No fasta header specified")?;
    let include_node_coverage = matches.is_present("include-node-coverage");

    let gfa: GFAtk = match gfa_file {
        Some(f) => {
            if !f.ends_with(".gfa") {
                bail!("Input file is not a GFA.")
            }
            GFAtk(load_gfa(f)?)
        }
        None => match utils::is_stdin() {
            true => GFAtk(load_gfa_stdin(std::io::stdin().lock())?),
            false => bail!("No input from STDIN. Run `gfatk extract -h` for help."),
        },
    };

    // load gfa into graph structure
    let (graph_indices, gfa_graph) = gfa.into_digraph()?;

    // check how many subgraphs there are
    let no_subgraphs = gfa_graph
        .weakly_connected_components(graph_indices.clone())?
        .len();

    if no_subgraphs > 1 {
        eprintln!(
            "[-]\tThe input GFA has multiple subgraphs ({}).",
            no_subgraphs
        )
    }

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
