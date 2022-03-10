// make a dot language representation
// of the GFA

use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Result};

/// Make a DOT (<https://graphviz.org/doc/info/lang.html>) language representation of a GFA.
///
/// For example:
/// ```bash
/// gfatk dot in.gfa | dot -Tsvg out.svg
/// ```
pub fn dot(matches: &clap::ArgMatches) -> Result<()> {
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
            false => bail!("No input from STDIN. Run `gfatk dot -h` for help."),
        },
    };

    let (_, gfa_graph) = gfa.into_digraph()?;

    gfa_graph.dot(gfa)?;

    Ok(())
}
