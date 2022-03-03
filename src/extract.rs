use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;
use anyhow::{Context, Result};

pub fn extract(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").context("No gfa file specified")?;
    let sequence_id: usize = matches.value_of_t("sequence-id")?;
    let iterations: i32 = matches.value_of_t("iterations")?;

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    let (graph_indices, gfa_graph) = gfa.into_ungraph()?;

    // get the node index of the target sequence ID.
    let target_index = graph_indices
        .seg_id_to_node_index(sequence_id)
        .context("Sequence ID not found.")?;
    let mut collect_sequence_names = Vec::new();
    // and this is the first item in this vector
    collect_sequence_names.push(target_index);

    let sequences_to_keep = gfa_graph.recursive_search(
        sequence_id,
        iterations,
        collect_sequence_names,
        graph_indices,
    )?;

    gfa.print_extract(sequences_to_keep);

    Ok(())
}
