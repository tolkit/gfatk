use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Result};
use petgraph::graph::NodeIndex;

/// Supply a sequence/segment ID from the GFA, and extract the GFA with all nodes connected to the input node.
///
/// For example:
/// ```bash
/// gfatk extract in.gfa -s 1 > out.gfa
/// ```
pub fn extract(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("GFA");
    let sequence_ids: Vec<usize> = matches.values_of_t("sequence-ids")?.into_iter().collect();
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
    let target_indices = sequence_ids
        .iter()
        .map(|e| graph_indices.seg_id_to_node_index(*e).unwrap())
        .collect::<Vec<NodeIndex>>();

    let sequences_to_keep =
        gfa_graph.recursive_search(sequence_ids, iterations, target_indices, graph_indices)?;

    gfa.print_extract(sequences_to_keep);

    Ok(())
}
