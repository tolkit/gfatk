use clap::value_t;
use gfa::gfa::Orientation;
use gfa::gfa::GFA;

use crate::load::load_gfa;
use crate::utils::{parse_cigar, reverse_complement};

pub fn overlap(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // required so unwrap safely
    let gfa_file = matches.value_of("gfa").unwrap();
    let extend_length = value_t!(matches.value_of("size"), usize).unwrap_or_else(|e| e.exit());

    let gfa: GFA<usize, ()> = load_gfa(gfa_file)?;

    // tuple of (from: overlap - length (incl. overlap), to: overlap + length)
    let mut from_to = Vec::new();
    // outer loop over links
    for link in &gfa.links {
        // get all the info out of each link
        let from_segment = link.from_segment;
        let from_orient = link.from_orient;
        let to_segment = link.to_segment;
        let to_orient = link.to_orient;
        let overlap = parse_cigar(&link.overlap);

        eprintln!(
            "From segment {} ({}) to segment {} ({})\nOverlap: {}",
            from_segment, from_orient, to_segment, to_orient, overlap
        );

        let mut from_seq: &[u8] = &[];
        let mut to_seq: &[u8] = &[];

        // get the from and to sequences.
        for line in gfa.lines_iter() {
            // if we meet a segment, let's do something
            match line.some_segment() {
                Some(s) => {
                    if s.name == from_segment && s.name == to_segment {
                        from_seq = &s.sequence;
                        to_seq = &s.sequence;
                    } else if s.name == from_segment {
                        from_seq = &s.sequence;
                    } else if s.name == to_segment {
                        to_seq = &s.sequence;
                    }
                }
                None => (),
            }
        }

        // initiate so we can append to vec
        let mut overlap_str_from_f: Option<String> = None;
        let mut overlap_str_from_r: Option<String> = None;
        let mut overlap_str_to_f: Option<String> = None;
        let mut overlap_str_to_r: Option<String> = None;

        // deal with the from's
        match from_orient {
            Orientation::Forward => {
                // do nothing
                // let overlap_seq = &from_seq[from_seq.len() - overlap - extend_length..];
                let overlap_seq = &from_seq.get(from_seq.len() - overlap - extend_length..);
                // if the extend length is too long, it means that
                // we hit the start of the sequence, so take full slice.
                let overlap_str = match overlap_seq {
                    Some(sl) => std::str::from_utf8(sl).unwrap(),
                    None => std::str::from_utf8(&from_seq[..]).unwrap(),
                };
                overlap_str_from_f = Some(overlap_str.to_string());
            }
            // if the relative negative strand matches
            // revcomp and take the end.
            Orientation::Backward => {
                let revcomp = reverse_complement(&from_seq);
                // let overlap_revcomp = revcomp[revcomp.len() - overlap - extend_length..].to_vec();
                let overlap_revcomp = revcomp.get(revcomp.len() - overlap - extend_length..);

                let overlap_str = match overlap_revcomp {
                    Some(sl) => String::from_utf8(sl.to_vec()).unwrap(),
                    // take the whole thing.
                    None => String::from_utf8(revcomp).unwrap(),
                };

                overlap_str_from_r = Some(overlap_str);
            }
        }
        // deal with the to's
        // here we ignore the overlap, as that is
        // captured above.
        match to_orient {
            Orientation::Forward => {
                // do nothing
                // let overlap = &to_seq[overlap..overlap + extend_length];
                let overlap_seq = &to_seq.get(overlap..overlap + extend_length);

                let overlap_str = match overlap_seq {
                    Some(sl) => std::str::from_utf8(sl).unwrap(),
                    // from end of overlap to the end of the sequence
                    None => std::str::from_utf8(&to_seq[overlap..]).unwrap(),
                };

                overlap_str_to_f = Some(overlap_str.to_string());
            }
            // if the relative negative strand matches
            // revcomp and take the start.
            Orientation::Backward => {
                let revcomp = reverse_complement(&to_seq);
                // let overlap_revcomp = revcomp[overlap..overlap + extend_length].to_vec();
                let overlap_revcomp = revcomp.get(overlap..overlap + extend_length);

                let overlap_str = match overlap_revcomp {
                    Some(sl) => String::from_utf8(sl.to_vec()).unwrap(),
                    None => String::from_utf8(revcomp[overlap..].to_vec()).unwrap(),
                };

                overlap_str_to_r = Some(overlap_str);
            }
        }

        from_to.push((
            overlap_str_from_f,
            overlap_str_from_r,
            overlap_str_to_f,
            overlap_str_to_r,
            from_segment,
            to_segment,
            from_orient,
            to_orient,
        ));
    }

    // long winded...
    for tup in from_to {
        let ff = tup.0.unwrap_or("".to_string());
        let fr = tup.1.unwrap_or("".to_string());
        let tf = tup.2.unwrap_or("".to_string());
        let tr = tup.3.unwrap_or("".to_string());
        let from_seg = tup.4;
        let to_seg = tup.5;
        let from_orient = tup.6;
        let to_orient = tup.7;

        println!(
            ">{}({}) -> {}({}): extend = {}\n{}{}{}{}",
            from_seg, from_orient, to_seg, to_orient, extend_length, ff, fr, tf, tr
        );
    }

    Ok(())
}
