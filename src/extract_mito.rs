use crate::gfa::gfa::GFAtk;
use crate::gfa::graph::segments_subgraph;
use crate::gfa::writer::gfa_string;
use crate::load::load_gfa;
use crate::stats;
use anyhow::{Context, Result};

pub fn extract_mito(matches: &clap::ArgMatches) -> Result<()> {
    let segments =
        stats::stats(matches, true)?.context("Should never reach here with further = true")?;

    let gfa_file = matches.value_of("gfa").context("No gfa file specified")?;

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    let subgraph = segments_subgraph(&gfa.0, segments);

    println!("{}", gfa_string(&subgraph));

    Ok(())
}
