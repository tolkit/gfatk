use gfa::writer::write_gfa;
use gfa::{gfa::GFA, optfields::OptionalFields};

// modified from function of the same name in gfa crate.
pub fn gfa_string(gfa: &GFA<usize, OptionalFields>) -> String {
    let mut result = String::new();
    write_gfa(gfa, &mut result);
    result
}
