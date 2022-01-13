// extract sequences from a GFA
// easy peasy

use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;

pub fn fasta(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").unwrap();

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    gfa.print_sequences();

    Ok(())
}
