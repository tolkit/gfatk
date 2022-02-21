// why re-invent every wheel.
// taken from:
// https://github.com/chfi/rs-gfa-utils/blob/2065b001d107ee9f5d7abe04d65ab82193fc5904/src/commands.rs

use anyhow::{Context, Result};
use bstr::io::*;
use gfa::{
    gfa::{SegmentId, GFA},
    optfields::OptFields,
    parser::GFAParser,
};
use std::io::{BufReader, Read};

pub fn byte_lines_iter<'a, R: Read + 'a>(reader: R) -> Box<dyn Iterator<Item = Vec<u8>> + 'a> {
    Box::new(BufReader::new(reader).byte_lines().map(|l| l.unwrap()))
}

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
