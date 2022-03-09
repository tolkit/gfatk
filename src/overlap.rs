use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;
use anyhow::{Context, Result};

/// Generate overlaps between segments, with an optional parameter of how large to make these overlaps.
pub fn overlap(matches: &clap::ArgMatches) -> Result<()> {
    // required so unwrap safely
    let gfa_file = matches.value_of("gfa").context("No gfa file specified")?;
    let extend_length: usize = matches.value_of_t("size")?;

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    let overlaps = gfa.make_overlaps(extend_length)?;

    overlaps.print(extend_length);

    Ok(())
}
