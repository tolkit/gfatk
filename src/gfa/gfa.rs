use crate::gfa::graph::{segments_subgraph, GFAdigraph, GFAungraph};
use crate::gfa::writer;
use crate::path::GFAPath;
use crate::stats::GenomeType;
use crate::utils::{
    self, get_edge_coverage, parse_cigar, reverse_complement, GFAGraphLookups, GFAGraphPair,
};
use anyhow::{bail, Context, Result};
use gfa::gfa::{Orientation, GFA};
use gfa::optfields::{OptFieldVal, OptionalFields};
use indexmap::IndexMap;
use petgraph::graph::{Graph, NodeIndex, UnGraph};
use std::collections::HashMap;

/// A wrapper around GFA from the gfa crate
#[derive(Clone)]
pub struct GFAtk(pub GFA<usize, OptionalFields>);

impl GFAtk {
    /// Returns a tuple of GFAGraphLookups (a struct of indices/node names)
    /// and an undirected GFA graph structure.
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

    /// Returns a tuple of GFAGraphLookups (a struct of indices/node names) and an directed GFA graph structure.
    ///
    /// Most functionality of this binary is on directed graph structures
    pub fn into_digraph(&self) -> Result<(GFAGraphLookups, GFAdigraph)> {
        let gfa = &self.0;
        // eprintln!("[+]\tReading GFA into a directed graph.");
        let mut gfa_graph: Graph<usize, (Orientation, Orientation, Option<i64>)> = Graph::new();

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

            let ec = get_edge_coverage(&edge.optional)?;

            // get the node index for a given edge
            let from_index = graph_indices.seg_id_to_node_index(from)?;
            let to_index = graph_indices.seg_id_to_node_index(to)?;

            // add the edges
            gfa_graph.add_edge(from_index, to_index, (from_orient, to_orient, Some(ec)));
        }

        Ok((graph_indices, GFAdigraph(gfa_graph)))
    }

    /// A method to print a GFA to STDOUT, given a vector of sequence ID's to keep.
    pub fn print_extract(&self, sequences_to_keep: Vec<usize>) {
        let gfa = &self.0;
        let subgraph_gfa = GFAtk(segments_subgraph(&gfa, sequences_to_keep));

        print!("{}", writer::gfa_string(&subgraph_gfa.0));
    }

    /// Returns the overlaps between all the segments in a GFA.
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

    /// Computes the overlaps between segments as a vector of:
    /// `(usize, Orientation, usize, &str)`
    /// Which is a tuple of from/to, the orientation of the overlap, the extent of the overlap, and whether the overlap is at the start or end of a segment.
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

    /// A method to print to STDOUT a fasta, given a path through the GFA.
    ///
    /// The first segment is added first, then all subsequent segments
    /// (-overlap with previous segment).
    pub fn print_path_to_fasta(
        &self,
        merged_sorted_chosen_path_overlaps: IndexMap<usize, Vec<(Orientation, usize, &str)>>,
        fasta_header: &str,
        segments_not_in_path: Vec<usize>,
    ) -> Result<()> {
        let gfa = &self.0;

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
                                    let seq_minus_overlap =
                                        s.sequence.get(start_overlap.1..).with_context(|| {
                                            format!(
                                                "{} is outside the bounds of the sequence.",
                                                start_overlap.1
                                            )
                                        })?;
                                    seq_minus_overlap.to_vec()
                                }
                                Orientation::Backward => {
                                    let revcomp_seq = reverse_complement(&s.sequence);
                                    let seq_minus_overlap =
                                        revcomp_seq.get(start_overlap.1..).with_context(|| {
                                            format!(
                                                "{} is outside the bounds of the sequence.",
                                                start_overlap.1
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

        // just a blank new line.
        // as the last statements were just print!.
        println!("");

        // print the rest of the segments
        if !segments_not_in_path.is_empty() {
            for segment in segments_not_in_path {
                for line in gfa.lines_iter() {
                    match line.some_segment() {
                        Some(seg) => {
                            if seg.name == segment {
                                println!(
                                    ">{}\n{}",
                                    segment,
                                    std::str::from_utf8(&seg.sequence).with_context(|| format!(
                                        "Malformed UTF8: {:?}",
                                        &seg.sequence
                                    ))?
                                )
                            }
                        }
                        None => {}
                    }
                }
            }
        }
        Ok(())
    }

    /// The internal function called when `gfatk fasta` is called.
    ///
    /// Prints all segments of the GFA as-is.
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

    /// Two internal functions below to parse coverage of a GFA segment.
    ///
    /// Used in `gfatk stats`.
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

    /// Return the coverage and sequence length for a segment, given a segment name.
    ///
    /// Note segment names are always `usize`.
    pub fn node_seq_len_and_cov(&self, node: usize) -> Result<(usize, f32)> {
        let gfa = &self.0;

        let ll_tag: [u8; 2] = [108, 108];
        let mut seq_len = None;
        let mut cov = None;

        for segment in &gfa.segments {
            if segment.name == node {
                seq_len = Some(segment.sequence.len());
                let opt = &segment.optional;
                for c in opt {
                    if c.tag == ll_tag {
                        cov = Some(*Self::parse_coverage_opt(&c.value)?);
                        break;
                    }
                }
            }
        }

        Ok((
            seq_len.context("No sequence length for each segment in GFA.")?,
            cov.context("No segment coverage for each segment in GFA.")?,
        ))
    }

    /// The internal function called in `gfatk stats`.
    ///
    /// Returns average GC%, average coverage, and total sequence length for a GFA (sub)graph.
    pub fn sequence_stats(&self, genome_type: GenomeType) -> Result<(f32, f32, usize)> {
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

        if genome_type == GenomeType::None {
            println!("\tTotal sequence length:\t{}", total_sequence_length);
            println!("\tTotal sequence overlap length:\t{}", total_overlap_length);
            println!(
                "\tSequence length minus overlaps:\t{}",
                total_sequence_length as i32 - total_overlap_length as i32
            );
            println!("\tGC content of total sequence:\t{}", avg_gc);
            println!("\tAverage coverage of total segments:\t{}", cov);
        }

        Ok((avg_gc, cov, total_sequence_length))
    }

    /// Returns a `HashMap` of relative coverage of each node (segment) in the GFA.
    ///
    /// Relative here indicates that each segment coverage is divided by the lowest coverage node, and rounded.
    pub fn gen_cov_hash(
        &self,
        graph_lookup: &GFAGraphLookups,
    ) -> Result<HashMap<NodeIndex, usize>> {
        let gfa = &self.0;

        let ll_tag: [u8; 2] = [108, 108];
        let mut node_cov_map = HashMap::new();

        // the initial map contains node index and coverage
        for seg in &gfa.segments {
            // get the coverage
            let opts = &seg.optional;
            let node_index = graph_lookup.seg_id_to_node_index(seg.name)?;

            for opt in opts {
                if opt.tag == ll_tag {
                    let cov = Self::parse_coverage_opt(&opt.value)?;
                    node_cov_map.insert(node_index, *cov);
                }
            }
        }

        // we want to convert the node index and coverage
        // to node index and *relative* coverage
        let mut lowest_cov_iter = node_cov_map.iter().map(|(_k, v)| v).enumerate();
        let init = lowest_cov_iter
            .next()
            .context("No coverage information found in this GFA.")?;
        // we process the rest
        let result = lowest_cov_iter.try_fold(init, |acc, x| {
            // return None if x is NaN
            let cmp = x.1.partial_cmp(&acc.1)?;
            // if x is less than the acc
            let min = if let std::cmp::Ordering::Less = cmp {
                x
            } else {
                acc
            };
            Some(min)
        });

        // allocate to a new map as we want u32's
        let mut rel_cov_map = HashMap::new();

        for (k, v) in &node_cov_map {
            rel_cov_map.insert(*k, (v / result.unwrap().1).round() as usize);
        }

        Ok(rel_cov_map)
    }

    /// Take a [`GFAPath`] and print out the path
    /// from a GFA.
    ///
    /// Currently implemented requires two loops of the GFA, and storage
    /// of the sequences in a [`HashMap`].
    pub fn from_path_cli(&self, path: GFAPath, link_map: HashMap<String, usize>) -> Result<()> {
        let gfa = &self.0;

        // put all the segments in memory - easiest way for now.
        let mut seg_map = HashMap::new();

        for seg in &gfa.segments {
            let id = seg.name;
            let seq = seg.sequence.clone();

            seg_map.insert(id, seq);
        }

        println!(">{}", path.to_string());

        // now iterate over the path itself
        for path_el in path.inner.windows(2) {
            // from
            let seg_id_from = path_el[0].segment_id;
            let orientation_from = path_el[0].orientation;
            // to
            let seg_id_to = path_el[1].segment_id;
            let orientation_to = path_el[1].orientation;

            // format so we can match on the links map
            let cigar_match = format!(
                "{}{}|{}{}",
                seg_id_from, orientation_from, seg_id_to, orientation_to
            );

            let overlap = *link_map.get(&cigar_match).context(format!(
                "This link: {} - does not occur in the input GFA. Perhaps re-consider the input path?",
                cigar_match
            ))?;

            // for the first element in the path
            // we print the entire sequence
            if path_el[0].index == 0 {
                let sequence = match orientation_from {
                    Orientation::Forward => seg_map.get(&seg_id_from).unwrap().to_vec(),
                    Orientation::Backward => {
                        let seq = seg_map.get(&seg_id_from).unwrap();
                        utils::reverse_complement(seq)
                    }
                };

                // and we must print out this second sequence, otherwise it's skipped!
                let sequence2 = match orientation_to {
                    Orientation::Forward => seg_map.get(&seg_id_to).unwrap()[overlap..].to_vec(),
                    Orientation::Backward => {
                        let seq = seg_map.get(&seg_id_to).unwrap();
                        utils::reverse_complement(seq)[overlap..].to_vec()
                    }
                };

                print!("{}", std::str::from_utf8(&sequence)?);
                print!("{}", std::str::from_utf8(&sequence2)?);
            } else {
                // otherwise we print the second element of the windows :)
                // and these are all dealt with in the same way
                let sequence = match orientation_to {
                    Orientation::Forward => seg_map.get(&seg_id_to).unwrap()[overlap..].to_vec(),
                    Orientation::Backward => {
                        let seq = seg_map.get(&seg_id_to).unwrap();
                        utils::reverse_complement(seq)[overlap..].to_vec()
                    }
                };

                print!("{}", std::str::from_utf8(&sequence)?);
            }
        }

        // newline
        println!("");

        Ok(())
    }
}

/// Overlap from one segment to another.
pub struct Overlap {
    /// From segment forward.
    pub overlap_str_from_f: Option<String>,
    /// From segment reverse.
    pub overlap_str_from_r: Option<String>,
    /// To segment forward.
    pub overlap_str_to_f: Option<String>,
    /// To segment reverse.
    pub overlap_str_to_r: Option<String>,
    /// ID of from segment.
    pub from_segment: usize,
    /// ID of to segment.
    pub to_segment: usize,
    /// Orientation of from segment.
    pub from_orient: Orientation,
    /// Orientation of to segment.
    pub to_orient: Orientation,
}

/// A vector of `Overlap` structs.
pub struct Overlaps(Vec<Overlap>);

impl Overlaps {
    /// Create a new instance of `Overlaps`.
    fn new() -> Self {
        Self(Vec::new())
    }
    /// Append to `Overlaps`, adding another `Overlap`.
    fn push(&mut self, add: Overlap) {
        self.0.push(add)
    }
    /// Print overlaps to STDOUT.
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

#[cfg(test)]
mod tests {

    use super::*;
    use crate::load::load_gfa;
    use crate::stats::GenomeType;

    // the GFA -> GFAtk structure used in tests below.
    fn make_gfa(path: &str) -> GFAtk {
        GFAtk(load_gfa(path).unwrap())
    }

    #[test]
    fn test_gfa_sequence_stats() {
        let gfa = make_gfa("./tests/test_linear.gfa");

        // could be mitochondria/chloroplast
        let (gc, cov, len) = gfa.sequence_stats(GenomeType::Mitochondria).unwrap();

        assert!(cov == 40.0);
        assert!(gc == 0.39523807);
        assert!(len == 18);
    }

    #[test]
    fn test_gen_cov_hash() {
        let gfa = make_gfa("./tests/test_linear.gfa");

        let lookup = GFAGraphLookups(vec![
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(0),
                seg_id: 11,
            },
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(1),
                seg_id: 12,
            },
            crate::utils::GFAGraphPair {
                node_index: NodeIndex::new(2),
                seg_id: 13,
            },
        ]);

        let cov_hash = gfa.gen_cov_hash(&lookup).unwrap();

        assert_eq!(cov_hash.get(&NodeIndex::new(0)).unwrap(), &1);
        assert_eq!(cov_hash.get(&NodeIndex::new(1)).unwrap(), &2);
        assert_eq!(cov_hash.get(&NodeIndex::new(2)).unwrap(), &1);
    }
}
