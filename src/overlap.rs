use clap::value_t;

use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;

pub fn overlap(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // required so unwrap safely
    let gfa_file = matches.value_of("gfa").unwrap();
    let extend_length = value_t!(matches.value_of("size"), usize).unwrap_or_else(|e| e.exit());

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    let overlaps = gfa.make_overlaps(extend_length);

    overlaps.print(extend_length);

    Ok(())
}
