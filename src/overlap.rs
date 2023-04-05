use std::path::PathBuf;

use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Result};

/// Generate overlaps between segments, with an optional parameter of how large to make these overlaps.
/// For example:
/// ```bash
/// gfatk overlap in.gfa > out.fasta
/// ```
pub fn overlap(matches: &clap::ArgMatches) -> Result<()> {
    // required so unwrap safely
    let gfa_file = matches.get_one::<PathBuf>("GFA");
    let extend_length = *matches.get_one::<usize>("size").expect("defaulted by clap");

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
            false => bail!("No input from STDIN. Run `gfatk path -h` for help."),
        },
    };

    let overlaps = gfa.make_overlaps(extend_length)?;

    overlaps.print(extend_length);

    Ok(())
}
