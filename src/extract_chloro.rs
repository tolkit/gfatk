use crate::gfa::gfa_string;
use crate::{gfa::graph::segments_subgraph, stats, stats::GenomeType};
use anyhow::{Context, Result};

/// Using a combination of GC% of the segments, relative coverage of the
/// segments, and these parameters:
///
/// For example:
/// ```bash
/// gfatk extract-chloro in.gfa > out.gfa
/// ```
pub fn extract_chloro(matches: &clap::ArgMatches, genome_type: GenomeType) -> Result<()> {
    let result = stats::stats(matches, genome_type)?
        .context("Should never reach here with `stats::GenomeType::Chloroplast`")?;

    let subgraph = segments_subgraph(&result.0 .0, result.1);

    println!("{}", gfa_string(&subgraph));

    Ok(())
}
