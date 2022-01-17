use clap::{App, Arg};
use gfatk::extract;
use gfatk::fasta;
use gfatk::gaf;
use gfatk::linear;
use gfatk::overlap;
use gfatk::stats;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let matches = App::new("gfatk")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Some functions to process GFA files.")
        .subcommand(
            clap::SubCommand::with_name("overlap")
                .about("Extract overlaps from a GFA.")
                // output file name
                .arg(
                    Arg::with_name("gfa")
                        .short("g")
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
                )
                .arg(
                    Arg::with_name("size")
                        .short("s")
                        .long("size")
                        .default_value("1000")
                        .help("Region around overlap to extract."),
                ),
        )
        .subcommand(
            clap::SubCommand::with_name("extract")
                .about("Extract subgraph from a GFA, given a segment name.")
                // output file name
                .arg(
                    Arg::with_name("gfa")
                        .short("g")
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
                )
                .arg(
                    Arg::with_name("sequence-id")
                        .short("s")
                        .long("sequence-id")
                        .takes_value(true)
                        .required(true)
                        .help("Extract subgraph of which this sequence is part of."),
                )
                .arg(
                    Arg::with_name("iterations")
                        .short("i")
                        .long("iterations")
                        .default_value("3")
                        .help("Number of iterations to recursively search for connecting nodes."),
                ),
        )
        .subcommand(
            clap::SubCommand::with_name("gaf")
                .about("gaf")
                // output file name
                .arg(
                    Arg::with_name("gaf")
                        .short("g")
                        .long("gaf")
                        .takes_value(true)
                        .required(true)
                        .help("Input GAF file."),
                ),
        )
        .subcommand(
            clap::SubCommand::with_name("linear")
                .about(
                    "Force a linear representation of the graph.\nEach node is included once only.",
                )
                // output file name
                .arg(
                    Arg::with_name("gfa")
                        .short("g")
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
                )
                .arg(
                    Arg::with_name("fasta-header")
                        .short("f")
                        .long("fasta-header")
                        .takes_value(true)
                        .required(true)
                        .default_value("gfatk-linear")
                        .help("Name of the fasta header in the output file."),
                )
                .arg(
                    Arg::with_name("coverage-file")
                        .short("c")
                        .long("coverage-file")
                        .takes_value(true)
                        .help("Name of the text file indicating the oriented coverage of links in a GFA."),
                ),
        )
        .subcommand(
            clap::SubCommand::with_name("fasta")
                .about(
                    "Extract a fasta file.\nAlmost as simple as: awk \'/^S/{print \">\"$2\"\\n\"$3}\'.",
                )
                // output file name
                .arg(
                    Arg::with_name("gfa")
                        .short("g")
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
                ),
        )
        .subcommand(
            clap::SubCommand::with_name("stats")
                .about(
                    "Some stats about the input GFA.",
                )
                // output file name
                .arg(
                    Arg::with_name("gfa")
                        .short("g")
                        .long("gfa")
                        .takes_value(true)
                        .required(true)
                        .help("Input GFA file."),
                ),
        )
        .get_matches();

    let subcommand = matches.subcommand();
    match subcommand.0 {
        "overlap" => {
            let matches = subcommand.1.unwrap();
            overlap::overlap(matches)?;
        }
        "extract" => {
            let matches = subcommand.1.unwrap();
            extract::extract(matches)?;
        }
        "gaf" => {
            let matches = subcommand.1.unwrap();
            gaf::gaf_to_graph(matches)?;
        }
        "linear" => {
            let matches = subcommand.1.unwrap();
            linear::force_linear(matches)?;
        }
        "fasta" => {
            let matches = subcommand.1.unwrap();
            fasta::fasta(matches)?;
        }
        "stats" => {
            let matches = subcommand.1.unwrap();
            stats::stats(matches)?;
        }
        _ => {
            println!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            std::process::exit(1);
        }
    }

    Ok(())
}
