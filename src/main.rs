// Max Brown
// Wellcome Sanger Institute 2023

use std::path::PathBuf;

use anyhow::Result;
use clap::{crate_version, value_parser, Arg, ArgAction, Command};
use gfatk::{
    dot, extract, extract_chloro, extract_mito, fasta, linear, overlap, path, rename, resolve,
    stats::{self, GenomeType},
    trim,
};

fn main() -> Result<()> {
    let matches = Command::new("gfatk")
        .version(crate_version!())
        .propagate_version(true)
        .arg_required_else_help(true)
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("Explore and linearise (plant organellar) GFA files.")
        .subcommand(
            Command::new("overlap")
                .about("Extract overlaps from a GFA.")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("size")
                        .short('s')
                        .long("size")
                        .default_value("1000")
                        .value_parser(value_parser!(usize))
                        .help("Region around overlap to extract."),
                ),
        )
        .subcommand(
            Command::new("extract")
                .about("Extract subgraph from a GFA, given a segment name.")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("sequence-ids")
                        .short('s')
                        .long("sequence-ids")
                        .required(true)
                        .value_delimiter(',')
                        .value_parser(value_parser!(String))
                        .help("Extract subgraph of which this sequence is part of. Specifying multiple segments requires a delimiter, e.g. 1,2,3 - note there should be no spaces between delimited segments."),
                )
                .arg(
                    Arg::new("iterations")
                        .short('i')
                        .long("iterations")
                        .default_value("3")
                        .value_parser(value_parser!(i32))
                        .help("Number of iterations to recursively search for connecting nodes."),
                ),
        )
        .subcommand(
            Command::new("linear")
                .about("Force a linear representation of the graph.")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("include-node-coverage")
                        .short('i')
                        .long("include-node-coverage")
                        .action(ArgAction::SetTrue)
                        .help("Should the coverage information of the segments be incorporated into linearisation?")
                )
                .arg(
                    Arg::new("evaluate-subgraphs")
                        .short('e')
                        .long("evaluate-subgraphs")
                        .action(ArgAction::SetTrue)
                        .help("If there are multiple subgraphs within a GFA, evaluate linear on each of these.")
                )
                .arg(
                    Arg::new("node-threshold")
                        .short('n')
                        .long("node-threshold")
                        .default_value("10000")
                        .value_parser(value_parser!(usize))
                        .help("Skip (sub)graphs with more nodes than this, as a safety cap.")
                )
        )
        .subcommand(
            Command::new("fasta")
                .about(
                    "Extract a fasta file.\nAlmost as simple as: awk \'/^S/{print \">\"$2\"\\n\"$3}\'.",
                )
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
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
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("tabular")
                        .short('t')
                        .long("tabular")
                        .action(ArgAction::SetTrue)
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
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("size-lower")
                        .long("size-lower")
                        // 200,000 default
                        .default_value("200000")
                        .value_parser(value_parser!(usize))
                        .help("Minimum size (bp) of expected mitochondria."),
                )
                .arg(
                    Arg::new("size-upper")
                        .long("size-upper")
                        // 1 million is high enough for most species?
                        .default_value("1000000")
                        .value_parser(value_parser!(usize))
                        .help("Maximum size (bp) of expected mitochondria."),
                )
                .arg(
                    Arg::new("gc-lower")
                        .long("gc-lower")
                        .default_value("0.42")
                        .value_parser(value_parser!(f32))
                        .help("Minimum GC% of expected mitochondria."),
                )
                .arg(
                    Arg::new("gc-upper")
                        .long("gc-upper")
                        .default_value("0.50")
                        .value_parser(value_parser!(f32))
                        .help("Maximum GC% of expected mitochondria."),
                )
                .arg(
                    Arg::new("tabular")
                        .short('t')
                        .long("tabular")
                        .action(ArgAction::SetTrue)
                        .help("Output tabular stats.")
                )
                ,
        )
        .subcommand(
            Command::new("extract-chloro")
                .about(
                    "Extract the plastid from a GFA.",
                )
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("size-lower")
                        .long("size-lower")
                        .default_value("126000")
                        .value_parser(value_parser!(usize))
                        .help("Minimum size (bp) of expected plastid."),
                )
                .arg(
                    Arg::new("size-upper")
                        .long("size-upper")
                        .default_value("180000")
                        .value_parser(value_parser!(usize))
                        .help("Maximum size (bp) of expected plastid."),
                )
                .arg(
                    Arg::new("gc-lower")
                        .long("gc-lower")
                        .default_value("0.35")
                        .value_parser(value_parser!(f32))
                        .help("Minimum GC% of expected plastid."),
                )
                .arg(
                    Arg::new("gc-upper")
                        .long("gc-upper")
                        .default_value("0.39")
                        .value_parser(value_parser!(f32))
                        .help("Maximum GC% of expected plastid."),
                )
                .arg(
                    Arg::new("tabular")
                        .short('t')
                        .long("tabular")
                        .action(ArgAction::SetTrue)
                        .help("Output tabular stats.")
                ),
        )
        .subcommand(
            Command::new("dot")
                .about("Return the dot representation of a GFA.")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                ),
        )
        .subcommand(
            Command::new("trim")
                .about("Trim a GFA to remove nodes of degree < 4 (i.e. only has one neighbour).")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                ),
        )
        .subcommand(
            Command::new("path")
                .about("Supply an input path to evaluate a linear representation of.
\r\t\t  Input must be a text file of a single comma separated line with node ID's and orientations. E.g. 1+,2-,3+")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .index(1)
                        .help("Input GFA file.")
                )
                // TODO: this is broke
                .arg(
                    Arg::new("path_cli")
                        .index(2)
                        .value_parser(value_parser!(String))
                        .help("Input path from CLI.")
                )
                .arg(
                    Arg::new("path_file")
                        .short('p')
                        .long("path")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input path from file.")
                )
                .arg(
                    Arg::new("all_paths")
                        .short('a')
                        .long("all")
                        .action(ArgAction::SetTrue)
                        .help("If there are path (P) lines in the input, output all paths in fasta format.")
                ),
        )
        .subcommand(
            Command::new("resolve")
                .about("Resolve a repeat-rich assembly graph into a circular genome using PacBio HiFi read-path (GAF) evidence.")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
                .arg(
                    Arg::new("gaf")
                        .long("gaf")
                        .required(true)
                        .value_parser(value_parser!(PathBuf))
                        .help("GAF file of long reads aligned to the graph (e.g. from GraphAligner).")
                )
                .arg(
                    Arg::new("gff")
                        .long("gff")
                        .value_parser(value_parser!(PathBuf))
                        .help("Optional GFF gene annotation, used to weight gene-bearing segments during resolution.")
                )
                .arg(
                    Arg::new("bed")
                        .long("bed")
                        .value_parser(value_parser!(PathBuf))
                        .help("Optional oatk-style BED gene annotation (seq_name, align_from, align_to, gene_name, score_capped_at_1000, strand), used the same way as --gff. Can be given alongside --gff; segments from both are combined.")
                )
                .arg(
                    Arg::new("min-mapq")
                        .long("min-mapq")
                        .default_value("1")
                        .value_parser(value_parser!(u32))
                        .help("Minimum GAF mapping quality to keep a read.")
                )
                .arg(
                    Arg::new("min-identity")
                        .long("min-identity")
                        .default_value("0.90")
                        .value_parser(value_parser!(f64))
                        .help("Minimum alignment identity (id:f: tag) to keep a read.")
                )
                .arg(
                    Arg::new("time-limit")
                        .long("time-limit")
                        .default_value("45.0")
                        .value_parser(value_parser!(f64))
                        .help("Solver time budget in seconds for the single-circuit connectivity constraint, per attempt. Larger/more tangled graphs may need more.")
                )
                .arg(
                    Arg::new("no-bubble-fasta")
                        .long("no-bubble-fasta")
                        .action(ArgAction::SetTrue)
                        .help("Don't write bubble-arm sequences to the FASTA output. They are still reported (counts and evidence) on stderr.")
                ),
        )
        .subcommand(
            Command::new("rename")
                .about("Rename the segment ID's of a GFA.")
                .arg(
                    Arg::new("GFA")
                        .value_parser(value_parser!(PathBuf))
                        .help("Input GFA file.")
                )
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
        Some(("rename", matches)) => {
            rename::rename_gfa(matches)?;
        }
        Some(("resolve", matches)) => {
            resolve::resolve(matches)?;
        }
        _ => {
            eprintln!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            std::process::exit(1);
        }
    }

    Ok(())
}
