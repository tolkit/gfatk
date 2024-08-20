use std::path::PathBuf;

use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Context, Result};
use petgraph::graph::NodeIndex;

/// Supply a sequence/segment ID from the GFA, and extract the GFA with all nodes connected to the input node.
///
/// For example:
/// ```bash
/// gfatk extract in.gfa -s 1 > out.gfa
/// ```
pub fn extract(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.get_one::<PathBuf>("GFA");
    let sequence_ids = matches
        .get_many::<String>("sequence-ids")
        .expect("errored by clap")
        .collect::<Vec<_>>();
    let iterations = *matches
        .get_one::<i32>("iterations")
        .expect("defaulted by clap");

    let gfa: GFAtk = match gfa_file {
        Some(f) => {
            let ext = f.extension();
            match ext {
                Some(e) => {
                    if e == "gfa" {
                        GFAtk(load_gfa(f)?)
                    } else {
                        bail!("Input is not a GFA.")
                    }
                }
                None => bail!("Could not read file."),
            }
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
        .map(|e| graph_indices.seg_id_to_node_index(e.as_bytes().to_vec()))
        .collect::<Result<Vec<NodeIndex>>>();

    let target_indices =
        target_indices.context("One of your input segment ID's does not exist in the graph.")?;

    let sequences_to_keep = gfa_graph.recursive_search(
        sequence_ids.iter().map(|e| e.as_bytes().to_vec()).collect(),
        iterations,
        target_indices,
        graph_indices,
    )?;

    gfa.print_extract(sequences_to_keep);

    Ok(())
}
