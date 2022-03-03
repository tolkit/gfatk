// extract sequences from a GFA
// easy peasy

use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;
use anyhow::{Context, Result};

pub fn fasta(matches: &clap::ArgMatches) -> Result<()> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").context("No gfa file specified")?;

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    gfa.print_sequences()?;

    Ok(())
}
