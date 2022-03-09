// why re-invent every wheel.
// taken from:
// https://github.com/chfi/rs-gfa-utils/blob/2065b001d107ee9f5d7abe04d65ab82193fc5904/src/commands.rs

use anyhow::{Context, Result};
use bstr::io::*;
use gfa::{
    gfa::{SegmentId, GFA},
    optfields::OptFields,
    parser::{GFAParser, ParseError},
};
use std::io::{BufReader, Read, StdinLock};

/// Iterate over the byte lines of a file.
pub fn byte_lines_iter<'a, R: Read + 'a>(reader: R) -> Box<dyn Iterator<Item = Vec<u8>> + 'a> {
    Box::new(BufReader::new(reader).byte_lines().map(|l| l.unwrap()))
}

/// Given a path, load the GFA into a `GFA` struct.
pub fn load_gfa<N, T, P>(path: P) -> Result<GFA<N, T>>
where
    N: SegmentId,
    T: OptFields,
    P: AsRef<std::path::Path>,
{
    let parser = GFAParser::new();
    let gfa = parser.parse_file(path.as_ref()).with_context(|| {
        format!(
            "Failed to parse GFA from path: {:?}",
            path.as_ref().as_os_str()
        )
    })?;
    Ok(gfa)
}

// take input from stdin, instead of a file.
// we'll lock on to it, saves a bit of code repitition

/// If the file is coming from STDIN, this function reads a GFA in.
pub fn load_gfa_stdin<N, T>(stdin: StdinLock) -> Result<GFA<N, T>, ParseError>
where
    N: SegmentId,
    T: OptFields,
{
    let parser = GFAParser::new();
    let lines = BufReader::new(stdin).byte_lines();

    let mut gfa = GFA::new();

    for line in lines {
        let line = line?;
        // if this not added then
        if line.is_empty() {
            continue;
        }
        match parser.parse_gfa_line(line.as_ref()) {
            Ok(parsed) => gfa.insert_line(parsed),
            // I don't have access to the .tolerance field...
            // Err(err) if err.can_safely_continue(&parser.tolerance) => (),
            Err(err) => return Err(err),
        };
    }

    Ok(gfa)
}
