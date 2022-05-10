use ::gfa::writer::write_gfa;
use ::gfa::{gfa::GFA, optfields::OptionalFields};

/// A module where all the methods of GFA manipulations are defined.
pub mod gfa;
/// A module where a GFA is coerced into a petgraph `Graph` structure, with associated methods.
pub mod graph;

/// Writes a GFA to a string.
///
/// Modified from function of the same name in gfa crate.
pub fn gfa_string(gfa: &GFA<usize, OptionalFields>) -> String {
    let mut result = String::new();
    write_gfa(gfa, &mut result);
    result
}
