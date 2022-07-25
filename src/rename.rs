use crate::load::{load_gfa, load_gfa_stdin};
use crate::utils;
use anyhow::{bail, Result};
use gfa::{gfa::name_conversion::NameMap, gfa::GFA, optfields::OptionalFields};

pub fn rename_gfa(matches: &clap::ArgMatches) -> Result<()> {
    let gfa_file = matches.value_of("GFA");

    let gfa: GFA<Vec<u8>, OptionalFields> = match gfa_file {
        Some(f) => {
            if !f.ends_with(".gfa") {
                bail!("Input file is not a GFA.")
            }
            load_gfa(f)?
        }
        None => match utils::is_stdin() {
            true => load_gfa_stdin(std::io::stdin().lock())?,
            false => bail!("No input from STDIN. Run `gfatk path -h` for help."),
        },
    };

    let name_map = NameMap::build_from_gfa(&gfa);

    if let Some(new_gfa) = name_map.gfa_bytestring_to_usize(&gfa, false) {
        let gfa_string = crate::gfa::gfa_string(&new_gfa);
        println!("{}", gfa_string);
    } else {
        bail!("Could not convert segment ID's to usize.")
    }

    Ok(())
}
