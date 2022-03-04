use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Context, Result};

pub fn extract(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("GFA");
    let sequence_id: usize = matches.value_of_t("sequence-id")?;
    let iterations: i32 = matches.value_of_t("iterations")?;

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
