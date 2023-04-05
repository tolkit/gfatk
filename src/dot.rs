// make a dot language representation
// of the GFA

use std::path::PathBuf;

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
    let gfa_file = matches.get_one::<PathBuf>("GFA");

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
            false => bail!("No input from STDIN. Run `gfatk dot -h` for help."),
        },
    };

    let (_, gfa_graph) = gfa.into_digraph()?;

    gfa_graph.dot(gfa)?;

    Ok(())
}
