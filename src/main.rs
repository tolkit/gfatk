// Max Brown
// Wellcome Sanger Institute 2022

use anyhow::Result;
use clap::{crate_version, Arg, Command};
use gfatk::{
    dot, extract, extract_chloro, extract_mito, fasta, linear, overlap, path,
    stats::{self, GenomeType},
    trim,
};

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
                    Arg::new("sequence-ids")
                        .short('s')
                        .long("sequence-ids")
                        .takes_value(true)
                        .required(true)
                        .multiple_values(true)
                        .use_value_delimiter(true)
                        .require_value_delimiter(true)
                        .help("Extract subgraph of which this sequence is part of. Specifying multiple segments requires a delimiter, e.g. 1,2,3"),
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
                    Arg::new("include-node-coverage")
                        .short('i')
                        .long("include-node-coverage")
                        .help("Should the coverage information of the segments be incorporated into linearisation?")
                )
                .arg(
                    Arg::new("evaluate-subgraphs")
                        .short('e')
                        .long("evaluate-subgraphs")
                        .help("If there are multiple subgraphs within a GFA, evaluate linear on each of these.")
                )
        )
        .subcommand(
            Command::new("fasta")
                .about(
                    "Extract a fasta file.\nAlmost as simple as: awk \'/^S/{print \">\"$2\"\\n\"$3}\'.",
                )
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
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
                )
                .arg(
                    Arg::new("tabular")
                        .short('t')
                        .long("tabular")
                        .help("Output tabular stats.")
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
                )
                .arg(
                    Arg::new("size-lower")
                        .long("size-lower")
                        // 200,000 default
                        .default_value("200000")
                        .help("Minimum size (bp) of expected mitochondria."),
                )
                .arg(
                    Arg::new("size-upper")
                        .long("size-upper")
                        // 1 million is high enough for most species?
                        .default_value("1000000")
                        .help("Maximum size (bp) of expected mitochondria."),
                )
                .arg(
                    Arg::new("gc-lower")
                        .long("gc-lower")
                        .default_value("0.42")
                        .help("Minimum GC% of expected mitochondria."),
                )
                .arg(
                    Arg::new("gc-upper")
                        .long("gc-upper")
                        .default_value("0.50")
                        .help("Maximum GC% of expected mitochondria."),
                ),
        )
        .subcommand(
            Command::new("extract-chloro")
                .about(
                    "Extract the plastid from a GFA.",
                )
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("size-lower")
                        .long("size-lower")
                        .default_value("126000")
                        .help("Minimum size (bp) of expected plastid."),
                )
                .arg(
                    Arg::new("size-upper")
                        .long("size-upper")
                        .default_value("180000")
                        .help("Maximum size (bp) of expected plastid."),
                )
                .arg(
                    Arg::new("gc-lower")
                        .long("gc-lower")
                        .default_value("0.35")
                        .help("Minimum GC% of expected plastid."),
                )
                .arg(
                    Arg::new("gc-upper")
                        .long("gc-upper")
                        .default_value("0.39")
                        .help("Maximum GC% of expected plastid."),
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
        .subcommand(
            Command::new("trim")
                .about("Trim a GFA to remove nodes of degree < 4 (i.e. only has one neighbour).")
                .arg(
                    Arg::new("GFA")
                        .help("Input GFA file.")
                ),
        )
        .subcommand(
            Command::new("path")
                .about("Supply an input path to evaluate a linear representation of. Input must be a text file of a single comma separated line with node ID's and orientations. E.g.:\n\t1+,2-,3+")
                .arg(
                    Arg::new("GFA")
                        .index(1)
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("path_cli")
                        .index(2)
                        .help("Input path from CLI.")
                )
                .arg(
                    Arg::new("path_file")
                        .short('p')
                        .long("path")
                        .takes_value(true)
                        .help("Input path from file.")
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
            linear::linear(matches)?;
        }
        Some(("fasta", matches)) => {
            fasta::fasta(matches)?;
        }
        Some(("stats", matches)) => {
            stats::stats(matches, GenomeType::None)?;
        }
        Some(("extract-mito", matches)) => {
            extract_mito::extract_mito(matches, GenomeType::Mitochondria)?;
        }
        Some(("extract-chloro", matches)) => {
            extract_chloro::extract_chloro(matches, GenomeType::Chloroplast)?;
        }
        Some(("dot", matches)) => {
            dot::dot(matches)?;
        }
        Some(("trim", matches)) => {
            trim::trim(matches)?;
        }
        Some(("path", matches)) => {
            path::path(matches)?;
        }
        _ => {
            eprintln!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            std::process::exit(1);
        }
    }

    Ok(())
}
