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
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A single parsed and quality-filtered GAF record.
#[derive(Debug, Clone)]
pub struct GafRecord {
    /// The segments the read traversed, each with the orientation it was
    /// walked in.
    pub path: Vec<(Vec<u8>, Orientation)>,
    /// Start of the aligned region, in path coordinates: positions along
    /// `path` with the bases two linked segments share counted once (the
    /// path is as long as its segments less the overlaps between them).
    pub path_start: usize,
    /// End of the aligned region, in path-length coordinates.
    pub path_end: usize,
    /// Mapping quality (GAF column 12).
    pub mapq: u32,
    /// Alignment identity, from the `id:f:` optional tag, when present.
    pub identity: Option<f64>,
    /// Whether the read is aligned over enough of its length for this
    /// alignment to count towards read depth (see [`read_gaf_evidence`]).
    /// Adjacency evidence uses every record regardless.
    pub counts_for_depth: bool,
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
/// quality filter or has no path, rather than an error. A read lying
/// within one segment is kept: it says nothing about adjacency, but it is
/// read depth like any other, and most of it - on long segments nearly all
/// reads start and end inside the segment.
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
    if path.is_empty() {
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
        counts_for_depth: true,
    }))
}

/// What a GAF says about a graph: the reads that cross at least one
/// junction, kept whole, and - for read depth - how many bases the reads
/// lying within a single segment align to each segment. Those are most of
/// the reads (over 90% on an organelle graph) and say nothing about
/// adjacency, so they are tallied as they are read, not stored.
#[derive(Debug, Default)]
pub struct GafEvidence {
    /// Reads whose path has two or more segments.
    pub records: Vec<GafRecord>,
    /// Aligned bases of single-segment alignments that count for depth,
    /// per segment.
    pub single_segment_bases: HashMap<Vec<u8>, u64>,
    /// Reads with at least one alignment passing the quality filter.
    pub n_reads: usize,
    /// Of those, the reads aligned over enough of their length to count
    /// for depth.
    pub n_reads_for_depth: usize,
}

/// Read and quality-filter every usable record from a GAF file.
///
/// A read counts for depth only if its alignments, together, cover at
/// least `min_read_aligned` of its length. Reads recruited to an organelle
/// graph are often only partly organellar - a nuclear read carrying an
/// organelle-derived insert, say - and the piece that aligns is not depth
/// of the organelle genome: on mitochondrial graphs two reads in three can
/// be like this, and counting them puts depth out by segment (most of all
/// against short segments, which few of them reach). Alignments of one
/// read are taken to be consecutive lines, as GraphAligner writes them.
pub fn read_gaf_evidence<P: AsRef<Path>>(
    path: P,
    min_mapq: u32,
    min_identity: f64,
    min_read_aligned: f64,
) -> Result<GafEvidence> {
    let file = File::open(path.as_ref())
        .with_context(|| format!("Failed to open GAF file: {:?}", path.as_ref()))?;
    let reader = BufReader::new(file);

    let mut evidence = GafEvidence::default();
    let mut read_name = String::new();
    let mut read_len: usize = 0;
    // this read's alignments: (query start, query end, record)
    let mut pending: Vec<(usize, usize, GafRecord)> = Vec::new();

    fn flush(
        pending: &mut Vec<(usize, usize, GafRecord)>,
        read_len: usize,
        min_read_aligned: f64,
        evidence: &mut GafEvidence,
    ) {
        if pending.is_empty() {
            return;
        }
        // bases of the read covered by any of its alignments
        let mut spans: Vec<(usize, usize)> = pending.iter().map(|(s, e, _)| (*s, *e)).collect();
        spans.sort_unstable();
        let (mut covered, mut reached) = (0usize, 0usize);
        for (s, e) in spans {
            if e > reached {
                covered += e - s.max(reached);
                reached = e;
            }
        }
        let counts = covered as f64 >= min_read_aligned * read_len as f64;
        evidence.n_reads += 1;
        if counts {
            evidence.n_reads_for_depth += 1;
        }
        for (_, _, mut rec) in pending.drain(..) {
            if rec.path.len() == 1 {
                if counts {
                    let (seg, _) = rec.path.pop().expect("one segment");
                    *evidence.single_segment_bases.entry(seg).or_insert(0) +=
                        rec.path_end.saturating_sub(rec.path_start) as u64;
                }
            } else {
                rec.counts_for_depth = counts;
                evidence.records.push(rec);
            }
        }
    }

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut fields = line.splitn(5, '\t');
        let name = fields.next().unwrap_or("");
        let len: usize = fields.next().and_then(|f| f.parse().ok()).unwrap_or(0);
        let start: usize = fields.next().and_then(|f| f.parse().ok()).unwrap_or(0);
        let end: usize = fields.next().and_then(|f| f.parse().ok()).unwrap_or(0);
        if name != read_name {
            flush(&mut pending, read_len, min_read_aligned, &mut evidence);
            read_name.clear();
            read_name.push_str(name);
            read_len = len;
        }
        if let Some(rec) = parse_gaf_line(&line, min_mapq, min_identity)? {
            pending.push((start, end, rec));
        }
    }
    flush(&mut pending, read_len, min_read_aligned, &mut evidence);
    Ok(evidence)
}

/// Read and quality-filter the records of a GAF file whose path crosses at
/// least one junction (two or more segments).
pub fn read_gaf_records<P: AsRef<Path>>(
    path: P,
    min_mapq: u32,
    min_identity: f64,
) -> Result<Vec<GafRecord>> {
    Ok(read_gaf_evidence(path, min_mapq, min_identity, 0.0)?.records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_gaf_evidence_counts_only_well_aligned_reads_for_depth() {
        // read1 (1000 bp): two alignments covering 0-600 and 500-1000 = all of it
        // read2 (1000 bp): one 300 bp alignment = 30% of it
        // read3 (1000 bp): crosses a junction but only 40% aligned
        let gaf = "\
read1\t1000\t0\t600\t+\t>u1\t5000\t100\t700\t600\t600\t60\tid:f:0.99
read1\t1000\t500\t1000\t+\t>u2>u3\t9000\t0\t500\t500\t500\t60\tid:f:0.99
read2\t1000\t0\t300\t+\t>u1\t5000\t0\t300\t300\t300\t60\tid:f:0.99
read3\t1000\t0\t400\t+\t>u2>u3\t9000\t0\t400\t400\t400\t60\tid:f:0.99
";
        let path = std::env::temp_dir().join(format!("gfatk_test_{}.gaf", std::process::id()));
        std::fs::write(&path, gaf).unwrap();
        let ev = read_gaf_evidence(&path, 1, 0.9, 0.9).unwrap();
        let all = read_gaf_evidence(&path, 1, 0.9, 0.0).unwrap();
        std::fs::remove_file(&path).unwrap();
        assert_eq!((ev.n_reads, ev.n_reads_for_depth), (3, 1));
        assert_eq!(ev.single_segment_bases[&b"u1".to_vec()], 600);
        // both junction-crossing reads are kept as adjacency evidence
        assert_eq!(ev.records.len(), 2);
        assert_eq!(ev.records.iter().filter(|r| r.counts_for_depth).count(), 1);
        // with no threshold every alignment counts
        assert_eq!(all.single_segment_bases[&b"u1".to_vec()], 900);
        assert!(all.records.iter().all(|r| r.counts_for_depth));
    }

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
    fn test_parse_gaf_line_keeps_single_segment_reads() {
        let line = "read1\t1000\t0\t998\t+\t>u1\t2000\t100\t1098\t900\t998\t60\tid:f:0.99";
        let rec = parse_gaf_line(line, 1, 0.9).unwrap().unwrap();
        assert_eq!(rec.path.len(), 1);
        assert_eq!((rec.path_start, rec.path_end), (100, 1098));
    }

    #[test]
    fn test_parse_gaf_line_drops_empty_paths() {
        let line = "read1\t1000\t0\t998\t+\t*\t2000\t100\t1098\t900\t998\t60\tid:f:0.99";
        assert!(parse_gaf_line(line, 1, 0.9).unwrap().is_none());
    }
}
