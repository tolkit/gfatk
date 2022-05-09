//! `gfatk` is a tool for Graphical Fragment Assembly (GFA) manipulation.
//!
//! GFA's are at their heart, simple, directed graphs. As such, all internal
//! representations of GFA's are [`petgraph::Graph`'s](https://docs.rs/petgraph/latest/petgraph/graph/struct.Graph.html).
//!
//! `gfatk` is designed mainly for reasonably small GFA files, and was designed for
//! de-tangling and linearising plant mitochondrial genomes on the command line.
//!
//! Almost all of the core functionality of `gfatk` resides in the [`gfatk::gfa`](./gfa/index.html)
//! module. The other modules are entry points for the command line application.

/// Make a DOT language representation of a GFA.
pub mod dot;
/// Extract a subgraph from a GFA.
pub mod extract;
/// Extract the putative plastid subgraph in a GFA.
pub mod extract_chloro;
/// Extract the putative mitochondrial subgraph in a GFA.
pub mod extract_mito;
/// Print all the sequences in a GFA to fasta format.
pub mod fasta;
/// A module with all the methods to manipulate GFA's in.
pub mod gfa;
/// Coerce a GFA into a fasta, finding the longest path through the graph.
pub mod linear;
/// Helper functions to load a GFA from a file, or read from STDIN.
/// Modified from <https://github.com/chfi/rs-gfa-utils/blob/2065b001d107ee9f5d7abe04d65ab82193fc5904/src/commands.rs>
pub mod load;
/// Generate overlapping sequences between segments in a GFA.
pub mod overlap;
/// Extract a fasta given a path.
pub mod path;
/// Generate statistics about the input GFA file.
pub mod stats;
/// Utility to trim a GFA of isolated nodes.
pub mod trim;
/// Utility functions used throughout.
pub mod utils;
