//! Parsing of GAF (Graph Alignment Format) records: how a long read (e.g.
//! PacBio HiFi) actually walked through a GFA assembly graph, as produced by
//! tools like `GraphAligner`.
//!
//! A GAF path field such as `>u66<u67>u28` is the read's literal walk: `u66`
//! forward, `u67` reverse, `u28` forward. That's much stronger evidence for
//! resolving repeats and branch points than the graph topology alone, since
//! it comes from real sequenced molecules rather than an assumption about
//! coverage.

use anyhow::{bail, Context, Result};
use gfa::gfa::Orientation;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A single parsed and quality-filtered GAF record.
#[derive(Debug, Clone)]
pub struct GafRecord {
    /// The segments the read traversed, each with the orientation it was
    /// walked in.
    pub path: Vec<(Vec<u8>, Orientation)>,
    /// Start of the aligned region, in path-length coordinates (i.e. summed
    /// segment lengths along `path`, ignoring link overlaps).
    pub path_start: usize,
    /// End of the aligned region, in path-length coordinates.
    pub path_end: usize,
    /// Mapping quality (GAF column 12).
    pub mapq: u32,
    /// Alignment identity, from the `id:f:` optional tag, when present.
    pub identity: Option<f64>,
}

/// Parse a GAF path field, e.g. `>u66<u67>u28`, into segment/orientation
/// pairs.
pub fn parse_gaf_path(path_str: &str) -> Vec<(Vec<u8>, Orientation)> {
    let mut out = Vec::new();
    let mut current_orient: Option<Orientation> = None;
    let mut current_name = String::new();

    for c in path_str.chars() {
        match c {
            '>' | '<' => {
                if let Some(o) = current_orient.take() {
                    if !current_name.is_empty() {
                        out.push((current_name.as_bytes().to_vec(), o));
                    }
                }
                current_orient = Some(if c == '>' {
                    Orientation::Forward
                } else {
                    Orientation::Backward
                });
                current_name = String::new();
            }
            _ => current_name.push(c),
        }
    }
    if let Some(o) = current_orient {
        if !current_name.is_empty() {
            out.push((current_name.as_bytes().to_vec(), o));
        }
    }
    out
}

/// Parse one tab-delimited GAF line, applying the mapping-quality and
/// identity thresholds. Returns `Ok(None)` for a record that fails the
/// quality filter or has too short a path to be useful (fewer than 2
/// segments), rather than an error -- most GAF lines will legitimately be
/// filtered this way.
pub fn parse_gaf_line(line: &str, min_mapq: u32, min_identity: f64) -> Result<Option<GafRecord>> {
    let fields: Vec<&str> = line.trim_end().split('\t').collect();
    if fields.len() < 12 {
        bail!("Malformed GAF line, fewer than 12 columns: {line}");
    }

    let mapq: u32 = fields[11].parse().context("GAF mapq not an integer")?;

    let mut identity = None;
    for tag in &fields[12..] {
        if let Some(rest) = tag.strip_prefix("id:f:") {
            identity = rest.parse::<f64>().ok();
            break;
        }
    }

    if mapq < min_mapq || identity.map(|i| i < min_identity).unwrap_or(false) {
        return Ok(None);
    }

    let path = parse_gaf_path(fields[5]);
    if path.len() < 2 {
        return Ok(None);
    }

    let path_start: usize = fields[7].parse().context("GAF path start not an integer")?;
    let path_end: usize = fields[8].parse().context("GAF path end not an integer")?;

    Ok(Some(GafRecord {
        path,
        path_start,
        path_end,
        mapq,
        identity,
    }))
}

/// Read and quality-filter every usable record from a GAF file.
pub fn read_gaf_records<P: AsRef<Path>>(
    path: P,
    min_mapq: u32,
    min_identity: f64,
) -> Result<Vec<GafRecord>> {
    let file = File::open(path.as_ref())
        .with_context(|| format!("Failed to open GAF file: {:?}", path.as_ref()))?;
    let reader = BufReader::new(file);

    let mut records = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        if let Some(rec) = parse_gaf_line(&line, min_mapq, min_identity)? {
            records.push(rec);
        }
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_gaf_path() {
        let path = parse_gaf_path(">u66<u67>u28>u42<u43");
        assert_eq!(
            path,
            vec![
                (b"u66".to_vec(), Orientation::Forward),
                (b"u67".to_vec(), Orientation::Backward),
                (b"u28".to_vec(), Orientation::Forward),
                (b"u42".to_vec(), Orientation::Forward),
                (b"u43".to_vec(), Orientation::Backward),
            ]
        );
    }

    #[test]
    fn test_parse_gaf_line_filters_low_quality() {
        let line = "read1\t1000\t0\t998\t+\t>u1>u2\t2000\t100\t1098\t900\t998\t0\tid:f:0.5";
        let rec = parse_gaf_line(line, 1, 0.9).unwrap();
        assert!(rec.is_none());
    }

    #[test]
    fn test_parse_gaf_line_keeps_good_record() {
        let line = "read1\t1000\t0\t998\t+\t>u1>u2\t2000\t100\t1098\t900\t998\t60\tid:f:0.99";
        let rec = parse_gaf_line(line, 1, 0.9).unwrap().unwrap();
        assert_eq!(rec.path.len(), 2);
        assert_eq!(rec.path_start, 100);
        assert_eq!(rec.path_end, 1098);
        assert_eq!(rec.mapq, 60);
        assert!((rec.identity.unwrap() - 0.99).abs() < 1e-9);
    }

    #[test]
    fn test_parse_gaf_line_drops_short_paths() {
        let line = "read1\t1000\t0\t998\t+\t>u1\t2000\t100\t1098\t900\t998\t60\tid:f:0.99";
        let rec = parse_gaf_line(line, 1, 0.9).unwrap();
        assert!(rec.is_none());
    }
}
