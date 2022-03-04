use anyhow::Result;
use clap::{crate_version, Arg, Command};
use gfatk::dot;
use gfatk::extract;
use gfatk::extract_mito;
use gfatk::fasta;
use gfatk::linear;
use gfatk::overlap;
use gfatk::stats;

fn main() -> Result<()> {
    let matches = Command::new("gfatk")
        .version(crate_version!())
        .propagate_version(true)
        .arg_required_else_help(true)
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Explore and linearise (mitochondrial) GFA files.")
        .subcommand(
            Command::new("overlap")
                .about("Extract overlaps from a GFA.")
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("size")
                        .short('s')
                        .long("size")
                        .default_value("1000")
                        .help("Region around overlap to extract."),
                ),
        )
        .subcommand(
            Command::new("extract")
                .about("Extract subgraph from a GFA, given a segment name.")
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("sequence-id")
                        .short('s')
                        .long("sequence-id")
                        .takes_value(true)
                        .required(true)
                        .help("Extract subgraph of which this sequence is part of."),
                )
                .arg(
                    Arg::new("iterations")
                        .short('i')
                        .long("iterations")
                        .default_value("3")
                        .help("Number of iterations to recursively search for connecting nodes."),
                ),
        )
        .subcommand(
            Command::new("linear")
                .about("Force a linear representation of the graph.")
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("fasta-header")
                        .short('f')
                        .long("fasta-header")
                        .takes_value(true)
                        .required(true)
                        .default_value("gfatk-linear")
                        .help("Name of the fasta header in the output file."),
                )
                .arg(
                    Arg::new("include-node-coverage")
                        .short('i')
                        .long("include-node-coverage")
                        .help("Should the coverage information of the segments be incorporated into linearisation?")
                ),
        )
        .subcommand(
            Command::new("fasta")
                .about(
                    "Extract a fasta file.\nAlmost as simple as: awk \'/^S/{print \">\"$2\"\\n\"$3}\'.",
                )
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
                ),
        )
        .subcommand(
            Command::new("stats")
                .about(
                    "Some stats about the input GFA.",
                )
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                ),
        )
        .subcommand(
            Command::new("extract-mito")
                .about(
                    "Extract the mitochondria from a GFA.",
                )
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                ),
        )
        .subcommand(
            Command::new("dot")
                .about("Return the dot representation of a GFA.")
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                ),
        )
        .get_matches();

    match matches.subcommand() {
        Some(("overlap", matches)) => {
            overlap::overlap(matches)?;
        }
        Some(("extract", matches)) => {
            extract::extract(matches)?;
        }
        Some(("linear", matches)) => {
            linear::force_linear(matches)?;
        }
        Some(("fasta", matches)) => {
            fasta::fasta(matches)?;
        }
        Some(("stats", matches)) => {
            stats::stats(matches, false)?;
        }
        Some(("extract-mito", matches)) => {
            extract_mito::extract_mito(matches)?;
        }
        Some(("dot", matches)) => {
            dot::dot(matches)?;
        }
        _ => {
            println!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            std::process::exit(1);
        }
    }

    Ok(())
}
