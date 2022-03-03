// make a dot language representation
// of the GFA

use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;
use anyhow::{Context, Result};

pub fn dot(matches: &clap::ArgMatches) -> Result<()> {
    let gfa_file = matches.value_of("gfa").context("No gfa file specified")?;

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    // load gfa into graph structure
    let (_graph_indices, gfa_graph) = gfa.into_digraph()?;

    gfa_graph.dot()?;

    Ok(())
}
