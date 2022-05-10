use crate::gfa::gfa::GFAtk;
use crate::gfa::gfa_string;
use crate::gfa::graph::segments_subgraph;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Result};

/// Trim a GFA file of segments which are connected only to one other segment.
///
/// For example:
/// ```bash
/// gfatk trim in.gfa > out.gfa
/// ```
pub fn trim(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("GFA");

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

    let (graph_indices, gfa_graph) = gfa.into_digraph()?;

    let trimmed = gfa_graph.trim(graph_indices);

    let subgraph = segments_subgraph(&gfa.0, trimmed);

    println!("{}", gfa_string(&subgraph));

    Ok(())
}
