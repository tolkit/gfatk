use clap::{App, AppSettings, Arg};
use gfatk::extract;
use gfatk::extract_mito;
use gfatk::fasta;
use gfatk::gaf;
use gfatk::linear;
use gfatk::overlap;
use gfatk::stats;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = App::new("gfatk")
        .global_setting(AppSettings::PropagateVersion)
        .global_setting(AppSettings::UseLongFormatForHelpSubcommand)
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Some functions to process GFA files.")
        .subcommand(
            App::new("overlap")
                .about("Extract overlaps from a GFA.")
                // output file name
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
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
            App::new("extract")
                .about("Extract subgraph from a GFA, given a segment name.")
                // output file name
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
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
            App::new("gaf")
                .about("Extract coverages from GAF file.")
                // output file name
                .arg(
                    Arg::new("gaf")
                        .short('g')
                        .long("gaf")
                        .takes_value(true)
                        .required(true)
                        .help("Input GAF file."),
                ),
        )
        .subcommand(
            App::new("linear")
                .about(
                    "Force a linear representation of the graph.\nEach node is included once only.",
                )
                // output file name
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
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
                    Arg::new("coverage-file")
                        .short('c')
                        .long("coverage-file")
                        .takes_value(true)
                        .help("Name of the text file indicating the oriented coverage of links in a GFA."),
                ),
        )
        .subcommand(
            App::new("fasta")
                .about(
                    "Extract a fasta file.\nAlmost as simple as: awk \'/^S/{print \">\"$2\"\\n\"$3}\'.",
                )
                // output file name
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
            App::new("stats")
                .about(
                    "Some stats about the input GFA.",
                )
                // output file name
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
            App::new("extract-mito")
                .about(
                    "Extract the mitochondria from a GFA.",
                )
                // output file name
                .arg(
                    Arg::new("gfa")
                        .short('g')
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
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
        Some(("gaf", matches)) => {
            gaf::gaf_to_graph(matches)?;
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
        _ => {
            println!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            std::process::exit(1);
        }
    }

    Ok(())
}
