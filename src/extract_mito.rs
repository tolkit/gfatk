use crate::gfa::graph::segments_subgraph;
use crate::gfa::writer::gfa_string;
use crate::stats;
use anyhow::{Context, Result};

/// Using a combination of GC% of the segments, relative coverage of the segments, and expected length of the mitochondrial genome, extract the putative mitochondrial subgraph.
///
/// For example:
/// ```bash
/// gfatk extract-mito in.gfa > out.gfa
/// ```
pub fn extract_mito(matches: &clap::ArgMatches) -> Result<()> {
    let result =
        stats::stats(matches, true)?.context("Should never reach here with further = true")?;

    let subgraph = segments_subgraph(&result.0 .0, result.1);

    println!("{}", gfa_string(&subgraph));

    Ok(())
}
