use crate::gfa::gfa::GFAtk;
use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Result};

/// Print a fasta representation of the sequences in a GFA.
///
/// For example:
/// ```bash
/// gfatk fasta in.gfa > out.fasta
/// ```
pub fn fasta(matches: &clap::ArgMatches) -> Result<()> {
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

    // None here, as we aren't lookiing/care about
    // subgraphs.
    gfa.print_sequences(None)?;

    Ok(())
}
