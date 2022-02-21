// common functions for gfa structures

use crate::gfa::graph::{segments_subgraph, GFAdigraph, GFAungraph};
use crate::gfa::writer;
use crate::utils;
use crate::utils::{parse_cigar, reverse_complement, GFAGraphLookups, GFAGraphPair};
use anyhow::{bail, Context, Result};
use csv::ReaderBuilder;
use gfa::gfa::{Orientation, GFA};
use gfa::optfields::{OptFieldVal, OptionalFields};
use indexmap::IndexMap;
use petgraph::graph::{Graph, NodeIndex, UnGraph};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct GAFTSVRecord {
    pub from_orient: char,
    pub from: usize,
    pub to_orient: char,
    pub to: usize,
    pub coverage: u32,
}

pub fn read_gaf_to_records(file: &str) -> Result<Vec<GAFTSVRecord>> {
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(file)
        .with_context(|| format!("File: {} not found.", file))?;
    let mut res = Vec::new();
    for result in rdr.deserialize() {
        let record: GAFTSVRecord = result.context("Could not parse record.")?;
        res.push(record);
    }
    Ok(res)
}

// the GFA type used throughout
pub struct GFAtk(pub GFA<usize, OptionalFields>);

impl GFAtk {
    // return a tuple of graph indices
    // and the graph itself
    pub fn into_ungraph(&self) -> Result<(GFAGraphLookups, GFAungraph)> {
        // alias to get GFA out
        let gfa = &self.0;
        // we're reading in now
        eprintln!("[+]\tReading GFA into an undirected graph.");
        let mut gfa_graph: UnGraph<usize, ()> = Graph::new_undirected();

        let mut graph_indices = GFAGraphLookups::new();
        // read the segments into graph nodes
        // save the indexes for populating the edges
        for node in &gfa.segments {
            let index = gfa_graph.add_node(node.name);
            graph_indices.push(GFAGraphPair {
                node_index: index,
                seg_id: node.name,
            });
        }

        // populate the edges
        for edge in &gfa.links {
            let from = edge.from_segment;
            let to = edge.to_segment;

            // get the node index for a given edge (like a map)
            let from_index = graph_indices.seg_id_to_node_index(from)?;
            let to_index = graph_indices.seg_id_to_node_index(to)?;

            // add the edges
            gfa_graph.add_edge(from_index, to_index, ());
        }

        Ok((graph_indices, GFAungraph(gfa_graph)))
    }

    // the Option<u32> here is for the coverages of the edges
    // which we optionally annotate later
    // maybe all these functions should be ported to ungraphs with weights?
    // not sure what the digraph does for GFA?

    pub fn into_digraph(&self) -> Result<(GFAGraphLookups, GFAdigraph)> {
        let gfa = &self.0;
        // eprintln!("[+]\tReading GFA into a directed graph.");
        let mut gfa_graph: Graph<usize, (Orientation, Orientation, Option<u32>)> = Graph::new();

        let mut graph_indices = GFAGraphLookups::new();
        // read the segments into graph nodes
        // save the indexes for populating the edges
        for node in &gfa.segments {
            let index = gfa_graph.add_node(node.name);
            graph_indices.push(GFAGraphPair {
                node_index: index,
                seg_id: node.name,
            });
        }

        // populate the edges
        for edge in &gfa.links {
            let from = edge.from_segment;
            let to = edge.to_segment;
            let from_orient = edge.from_orient;
            let to_orient = edge.to_orient;

            // get the node index for a given edge
            let from_index = graph_indices.seg_id_to_node_index(from)?;
            let to_index = graph_indices.seg_id_to_node_index(to)?;

            // add the edges
            gfa_graph.add_edge(from_index, to_index, (from_orient, to_orient, None));
        }

        Ok((graph_indices, GFAdigraph(gfa_graph)))
    }

    // print the GFA with only the sequences to keep
    pub fn print_extract(&self, sequences_to_keep: Vec<usize>) {
        let gfa = &self.0;
        let subgraph_gfa = GFAtk(segments_subgraph(&gfa, sequences_to_keep));

        print!("{}", writer::gfa_string(&subgraph_gfa.0));
    }

    // detect all the overlaps between joins in a GFA.
    pub fn make_overlaps(&self, extend_length: usize) -> Result<Overlaps> {
        let gfa = &self.0;
        // tuple of (from: overlap - length (incl. overlap), to: overlap + length)
        let mut from_to = Overlaps::new();
        // outer loop over links
        for link in &gfa.links {
            // get all the info out of each link
            let from_segment = link.from_segment;
            let from_orient = link.from_orient;
            let to_segment = link.to_segment;
            let to_orient = link.to_orient;
            let overlap = parse_cigar(&link.overlap)?;

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
                    // length - overlap - extend length at the end of the sequence.
                    let overlap_seq = &from_seq.get(from_seq.len() - overlap - extend_length..);
                    // if the extend length is too long, it means that
                    // we hit the start of the sequence, so take full slice.
                    let overlap_str = match overlap_seq {
                        Some(sl) => std::str::from_utf8(sl)
                            .with_context(|| format!("Malformed UTF8: {:?}", sl))?,
                        None => std::str::from_utf8(&from_seq[..])
                            .with_context(|| format!("Malformed UTF8: {:?}", &from_seq[..]))?,
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
                        Some(sl) => String::from_utf8(sl.to_vec())
                            .with_context(|| format!("Malformed UTF8: {:?}", sl))?,
                        // take the whole thing.
                        None => String::from_utf8(revcomp).context("Malformed UTF8.")?,
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
                        Some(sl) => std::str::from_utf8(sl)
                            .with_context(|| format!("Malformed UTF8: {:?}", sl))?,
                        // from end of overlap to the end of the sequence
                        None => std::str::from_utf8(&to_seq[overlap..])
                            .with_context(|| format!("Malformed UTF8: {:?}", &to_seq[overlap..]))?,
                    };

                    overlap_str_to_f = Some(overlap_str.to_string());
                }
                // if the relative negative strand matches
                // revcomp and take the start.
                Orientation::Backward => {
                    let revcomp = reverse_complement(&to_seq);
                    // let overlap_revcomp = revcomp[overlap..overlap + extend_length].to_vec();
                    let overlap_revcomp = revcomp.get(overlap..overlap + extend_length);

                    let overlap_str =
                        match overlap_revcomp {
                            Some(sl) => String::from_utf8(sl.to_vec())
                                .with_context(|| format!("Malformed UTF8: {:?}", sl))?,
                            None => String::from_utf8(revcomp[overlap..].to_vec()).with_context(
                                || format!("Malformed UTF8: {:?}", revcomp[overlap..].to_vec()),
                            )?,
                        };

                    overlap_str_to_r = Some(overlap_str);
                }
            }

            from_to.push(Overlap {
                overlap_str_from_f,
                overlap_str_from_r,
                overlap_str_to_f,
                overlap_str_to_r,
                from_segment,
                to_segment,
                from_orient,
                to_orient,
            });
        }
        Ok(from_to)
    }

    pub fn determine_path_overlaps(
        &self,
        chosen_path: &Vec<NodeIndex>,
        graph_indices: GFAGraphLookups,
        chosen_path_ids: &Vec<usize>,
    ) -> Result<Vec<(usize, Orientation, usize, &str)>> {
        let gfa = &self.0;
        // another new vec of (id, overlap, start/end)
        // to sort out for fasta generation
        let mut chosen_path_overlaps = Vec::new();

        for edge in &gfa.links {
            let from = edge.from_segment;
            let to = edge.to_segment;
            let from_orient = edge.from_orient;
            let to_orient = edge.to_orient;

            // overlap
            let overlap = parse_cigar(&edge.overlap)?;
            // here look at the strandedness between pairs of nodes
            let path_pairs = chosen_path.windows(2);
            for path in path_pairs {
                let from_path = graph_indices.node_index_to_seg_id(path[0])?;
                let to_path = graph_indices.node_index_to_seg_id(path[1])?;

                // okay this appears to work.
                if from == from_path && to == to_path {
                    chosen_path_overlaps.push((from, from_orient, overlap, "end"));
                    chosen_path_overlaps.push((to, to_orient, overlap, "start"));
                }
            }
        }

        chosen_path_overlaps.sort();
        chosen_path_overlaps.dedup();

        let mut sorted_chosen_path_overlaps = Vec::new();
        // sorting (allocate to new vec)
        for path_id in chosen_path_ids {
            for path_id2 in &chosen_path_overlaps {
                if &path_id2.0 == path_id {
                    sorted_chosen_path_overlaps.push(*path_id2)
                }
            }
        }
        Ok(sorted_chosen_path_overlaps)
    }

    // logic follows Marcela's python script now (thanks!)
    // first segment added, then segment[n+1][overlap..] added
    pub fn print_path_to_fasta(
        &self,
        merged_sorted_chosen_path_overlaps: IndexMap<usize, Vec<(Orientation, usize, &str)>>,
        fasta_header: &str,
    ) -> Result<()> {
        let gfa = &self.0;

        // debugging
        // for (k, v) in &merged_sorted_chosen_path_overlaps {
        //     println!("Index {}: {:?}", k, v);
        // }

        println!(">{}", fasta_header);
        // create iterator over the paths
        let mut path_iter = merged_sorted_chosen_path_overlaps.iter();
        // get the first element of the iterator.
        let (id_init, vector_init) = path_iter
            .next()
            .context("First element of the path not found.")?;

        for line in gfa.lines_iter() {
            match line.some_segment() {
                Some(s) => {
                    if s.name == *id_init {
                        let orientation = vector_init[0].0;
                        match orientation {
                            Orientation::Forward => {
                                // do nothing except print
                                print!(
                                    "{}",
                                    std::str::from_utf8(&s.sequence).with_context(|| format!(
                                        "Malformed UTF8: {:?}",
                                        &s.sequence
                                    ))?
                                );
                            }
                            Orientation::Backward => {
                                let revcomp_seq = reverse_complement(&s.sequence);
                                print!(
                                    "{}",
                                    std::str::from_utf8(&revcomp_seq).with_context(|| format!(
                                        "Malformed UTF8: {:?}",
                                        revcomp_seq
                                    ))?
                                );
                            }
                        }
                    }
                }
                None => (),
            }
        }

        // now do the rest of the paths
        for (id, vector) in path_iter {
            // get the from and to sequences.
            for line in gfa.lines_iter() {
                // if we meet a segment, let's do something
                match line.some_segment() {
                    Some(s) => {
                        // we've got the segment we wanted
                        if s.name == *id {
                            // remove start overlap
                            // but we need to watch out for orientation
                            // only need the starting overlaps
                            let start_overlap = vector
                                .iter()
                                .find(|(_or, _ov, side)| side == &"start")
                                .context("Start overlap could not be found.")?;

                            // start and end orientation should be the same
                            // hence just matching on start here.
                            let seq_minus_overlap = match start_overlap.0 {
                                Orientation::Forward => {
                                    // do nothing
                                    let seq_minus_overlap = s
                                        .sequence
                                        .get(start_overlap.1 - 1..)
                                        .with_context(|| {
                                        format!(
                                            "{} is outside the bounds of the sequence.",
                                            start_overlap.1 - 1
                                        )
                                    })?;
                                    seq_minus_overlap.to_vec()
                                }
                                Orientation::Backward => {
                                    let revcomp_seq = reverse_complement(&s.sequence);
                                    let seq_minus_overlap = revcomp_seq
                                        .get(start_overlap.1 - 1..)
                                        .with_context(|| {
                                            format!(
                                                "{} is outside the bounds of the sequence.",
                                                start_overlap.1 - 1
                                            )
                                        })?;
                                    seq_minus_overlap.to_vec()
                                }
                            };
                            print!(
                                "{}",
                                std::str::from_utf8(&seq_minus_overlap).with_context(|| {
                                    format!("Malformed UTF8: {:?}", &seq_minus_overlap)
                                })?
                            );
                        }
                    }
                    None => (),
                }
            }
        }
        Ok(())
    }

    pub fn print_sequences(&self) -> Result<()> {
        let gfa = &self.0;

        for line in gfa.lines_iter() {
            match line.some_segment() {
                Some(s) => {
                    let seq = std::str::from_utf8(&s.sequence)
                        .with_context(|| format!("Malformed UTF8: {:?}", &s.sequence))?;
                    let id = s.name;
                    println!(">{}\n{}", id, seq);
                }
                None => (),
            }
        }
        Ok(())
    }

    // get the average of the coverages from a GFA.
    fn parse_coverage_opt(opt: &OptFieldVal) -> Result<&f32> {
        let ll = match opt {
            OptFieldVal::Float(f) => f,
            _ => bail!("ll: coverage should be Float()"),
        };
        Ok(ll)
    }

    fn get_coverage(&self) -> Result<f32> {
        let gfa = &self.0;

        let ll_tag: [u8; 2] = [108, 108];
        let mut ll_tag_vec = Vec::new();

        for seg in &gfa.segments {
            let opts = &seg.optional;
            for opt in opts {
                if opt.tag == ll_tag {
                    ll_tag_vec.push(Self::parse_coverage_opt(&opt.value)?);
                }
            }
        }
        let len = ll_tag_vec.len() as f32;
        let sum: f32 = ll_tag_vec.iter().fold(0.0, |a, b| a + **b);

        Ok(sum / len)
    }

    // we actually return only a float here
    // but could change this later...

    pub fn sequence_stats(&self, further: bool) -> Result<(f32, f32)> {
        let gfa = &self.0;

        let cov = Self::get_coverage(&self)?;

        let mut total_overlap_length = 0;
        for link in &gfa.links {
            total_overlap_length += parse_cigar(&link.overlap)?;
        }

        let mut total_sequence_length = 0;
        let mut gc_vec = Vec::new();

        for segment in &gfa.segments {
            let seq = &segment.sequence;
            total_sequence_length += seq.len();

            let gc = utils::gc_content(&seq);
            gc_vec.push(gc);
        }

        let avg_gc = gc_vec.iter().sum::<f32>() / gc_vec.len() as f32;

        if !further {
            println!("\tTotal sequence length:\t{}", total_sequence_length);
            println!("\tTotal sequence overlap length:\t{}", total_overlap_length);
            println!(
                "\tSequence length minus overlaps:\t{}",
                total_sequence_length as i32 - total_overlap_length as i32
            );
            println!("\tGC content of total sequence:\t{}", avg_gc);
            println!("\tAverage coverage of total segments:\t{}", cov);
        }

        Ok((avg_gc, cov))
    }
}

pub struct Overlap {
    pub overlap_str_from_f: Option<String>,
    pub overlap_str_from_r: Option<String>,
    pub overlap_str_to_f: Option<String>,
    pub overlap_str_to_r: Option<String>,
    pub from_segment: usize,
    pub to_segment: usize,
    pub from_orient: Orientation,
    pub to_orient: Orientation,
}

pub struct Overlaps(Vec<Overlap>);

impl Overlaps {
    fn new() -> Self {
        Self(Vec::new())
    }
    fn push(&mut self, add: Overlap) {
        self.0.push(add)
    }
    pub fn print(self, extend_length: usize) {
        // long winded...
        for o in self.0 {
            // unwrap None -> zero length string.
            let ff = o.overlap_str_from_f.unwrap_or("".to_string());
            let fr = o.overlap_str_from_r.unwrap_or("".to_string());
            let tf = o.overlap_str_to_f.unwrap_or("".to_string());
            let tr = o.overlap_str_to_r.unwrap_or("".to_string());
            let from_seg = o.from_segment;
            let to_seg = o.to_segment;
            let from_orient = o.from_orient;
            let to_orient = o.to_orient;

            println!(
                ">{}({})->{}({}): extend = {}\n{}{}{}{}",
                from_seg, from_orient, to_seg, to_orient, extend_length, ff, fr, tf, tr
            );
        }
    }
}
