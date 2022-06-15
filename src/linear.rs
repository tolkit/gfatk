use crate::gfa::gfa::GFAtk;
use crate::gfa::graph::{segments_subgraph, GFAdigraph};
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils::{self, GFAGraphLookups};
use anyhow::{bail, Result};
use indexmap::IndexMap;
use petgraph::algo::is_cyclic_directed;

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
pub fn linear(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("GFA");
    let include_node_coverage = matches.is_present("include-node-coverage");
    let evaluate_subgraphs = matches.is_present("evaluate-subgraphs");
    // this unwrap_or is redundant.
    let node_threshold: usize = matches.value_of_t("node-threshold").unwrap_or(60);

    let gfa: GFAtk = match gfa_file {
        Some(f) => {
            if !f.ends_with(".gfa") {
                bail!("Input file is not a GFA.")
            }
            GFAtk(load_gfa(f)?)
        }
        None => match utils::is_stdin() {
            true => GFAtk(load_gfa_stdin(std::io::stdin().lock())?),
            false => bail!("No input from STDIN. Run `gfatk linear -h` for help."),
        },
    };

    // load gfa into graph structure
    let (graph_indices, gfa_graph) = gfa.into_digraph()?;

    // if we have only one node (segment) then all we can do
    // is print the sequence
    // otherwise we hit this error: `Error: There was no highest coverage path.`
    // makes sense as you can't have a path of length 1.
    if gfa_graph.node_count() == 1 {
        // as we would in `gfatk fasta`
        eprintln!("[+]\tOnly a single segment detected. Printing sequence and exiting.");
        gfa.print_sequences(None)?;
        return Ok(());
    }

    // check how many subgraphs there are
    let subgraphs = gfa_graph.weakly_connected_components(graph_indices.clone())?;

    // Warn user if there is more than one subgraph
    if subgraphs.len() > 1 {
        eprintln!(
            "[-]\tThe input GFA has multiple subgraphs ({}).",
            subgraphs.len()
        );
        if !evaluate_subgraphs {
            eprintln!("[-]\tYou did not specify the `-e` option, so only the first subgraph will be linearised.");
        }
    }

    match evaluate_subgraphs {
        true => {
            for (mut index, id_set) in subgraphs.iter().enumerate() {
                // so we don't zero index on the fasta headers
                index += 1;
                // have to make the extra allocation here.
                let gfa = gfa.clone();
                // make the new GFA
                let subgraph_gfa = GFAtk(segments_subgraph(&gfa.0, id_set.to_vec()));
                let (graph_indices_subgraph, subgraph) = subgraph_gfa.into_digraph()?;

                // it can be useful to see here whether the subgraph is cyclic
                let is_circular = is_cyclic_directed(&subgraph.0);
                // check the node count here. If there's one segment, then we can just print the sequence.
                // otherwise we go ahead and linearise the subgraph.

                if subgraph.node_count() == 1 {
                    let subgraph_index_header =
                        Some(format!(" subgraph-{}:is_circular-{}", index, is_circular));
                    subgraph_gfa.print_sequences(subgraph_index_header)?;
                } else if subgraph.node_count() > node_threshold {
                    eprintln!(
                        "[-]\tDetected {} nodes in a subgraph. Skipping.",
                        subgraph.node_count()
                    );
                    continue;
                } else {
                    // add a subgraph index to the fasta header
                    let subgraph_index_header =
                        Some(format!(" subgraph-{}:is_circular-{}", index, is_circular));
                    linear_inner(
                        subgraph_gfa,
                        include_node_coverage,
                        graph_indices_subgraph,
                        subgraph,
                        subgraph_index_header,
                    )?;
                }
            }
        }
        false => {
            if gfa_graph.node_count() > node_threshold {
                bail!(
                    "Detected {} nodes in this graph! It may be possible to linearise some subgraphs (if present) with the `-e` flag. Exiting.",
                    gfa_graph.node_count()
                );
            }
            linear_inner(gfa, include_node_coverage, graph_indices, gfa_graph, None)?;
        }
    }

    Ok(())
}

/// Reusable function to call on subgraphs in a GFA if necessary.
fn linear_inner(
    gfa: GFAtk,
    include_node_coverage: bool,
    graph_indices: GFAGraphLookups,
    gfa_graph: GFAdigraph,
    subgraph_index_header: Option<String>,
) -> Result<()> {
    // don't evaluate the coverage if we don't care about it
    let rel_coverage_map = match include_node_coverage {
        true => Some(gfa.gen_cov_hash(&graph_indices)?),
        false => None,
    };

    let (chosen_path, chosen_path_ids, segments_not_in_path, fasta_header) =
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
        &fasta_header,
        segments_not_in_path,
        subgraph_index_header,
    )?;
    Ok(())
}
