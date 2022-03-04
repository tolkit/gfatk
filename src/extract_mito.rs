use crate::gfa::graph::segments_subgraph;
use crate::gfa::writer::gfa_string;
use crate::stats;
use anyhow::{Context, Result};

pub fn extract_mito(matches: &clap::ArgMatches) -> Result<()> {
    let result =
        stats::stats(matches, true)?.context("Should never reach here with further = true")?;

    let subgraph = segments_subgraph(&result.0 .0, result.1);

    println!("{}", gfa_string(&subgraph));

    Ok(())
}
