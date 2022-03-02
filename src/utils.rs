use anyhow::{bail, Context, Result};
use gfa::optfields::{OptField, OptFieldVal::*};
use petgraph::graph::NodeIndex;
use std::collections::HashMap;

pub fn get_edge_coverage(options: &Vec<OptField>) -> Result<i64> {
    for op in options {
        match op.tag {
            // ec
            [101, 99] => match op.value {
                Int(i) => return Ok(i),
                _ => bail!("Could not find integer ec:i:<i64> tag."),
            },
            _ => bail!("Could not find ec (edge coverage) tag."),
        };
    }
    bail!("Edge coverage not found.")
}

pub fn get_option_string(options: Vec<OptField>) -> Result<String> {
    let mut tag_val = String::new();
    for op in options {
        let tag = std::str::from_utf8(&op.tag)
            .with_context(|| format!("Malformed UTF8: {:?}", op.tag))?;
        let value = match op.value {
            Float(f) => format!(":f:{:.3}", f),
            A(a) => format!(":A:{}", a.to_string()),
            Int(i) => format!(":i:{}", i.to_string()),
            Z(z) => format!(
                ":Z:{}",
                std::str::from_utf8(&z).with_context(|| format!("Malformed UTF8: {:?}", z))?
            ),
            // J(j) => ???,
            // a hexadecimal array
            H(h) => format!(":H:{}", h.iter().map(|x| x.to_string()).collect::<String>()),
            // B is a general array
            // is it capital B?
            BInt(bi) => format!(
                ":B:{}",
                bi.iter().map(|x| x.to_string()).collect::<String>()
            ),
            BFloat(bf) => format!(
                ":B:{}",
                bf.iter().map(|x| format!("{:.3}", x)).collect::<String>()
            ),
            _ => "".to_string(),
        };
        tag_val += &format!("{}{}\t", tag, value);
    }
    // should always end in \t ^
    let tag_val_op_un = tag_val
        .strip_suffix("\t")
        .context("Could not strip a tab from the suffix.")?;
    Ok(tag_val_op_un.to_string())
}

// not a very safe function, but works
pub fn parse_cigar(cigar: &[u8]) -> Result<usize> {
    // check it ends with an M
    if !cigar.ends_with(&[77]) {
        panic!("CIGAR did not end with M.");
    }
    let stripped = cigar.strip_suffix(&[77]);
    match stripped {
        Some(s) => {
            let string_rep =
                std::str::from_utf8(s).with_context(|| format!("Malformed UTF8: {:?}", s))?;
            Ok(string_rep
                .parse::<usize>()
                .with_context(|| format!("{} could not be parsed to <usize>", string_rep))?)
        }
        None => bail!("Could not strip suffix (M) of the CIGAR string."),
    }
}

pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    let dna_vec = dna.to_vec();
    let mut revcomp = Vec::new();

    for base in dna_vec.iter() {
        revcomp.push(switch_base(*base))
    }
    revcomp.as_mut_slice().reverse();
    revcomp
}

fn switch_base(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'T' => b'A',
        b't' => b'a',
        b'G' => b'C',
        b'g' => b'c',
        b'N' => b'N',
        b'n' => b'n',
        _ => b'N',
    }
}

// pinched from past Max
// https://github.com/tolkit/fasta_windows/blob/master/src/seq_statsu8.rs

fn nucleotide_counts(dna: &[u8]) -> HashMap<&u8, i32> {
    let mut map = HashMap::new();
    for nucleotide in dna {
        let count = map.entry(nucleotide).or_insert(0);
        *count += 1;
    }
    map
}

pub fn gc_content(dna: &[u8]) -> f32 {
    // G/C/A/T counts
    let counts = nucleotide_counts(dna);
    let g_counts = counts.get(&71).unwrap_or(&0) + counts.get(&103).unwrap_or(&0);
    let c_counts = counts.get(&67).unwrap_or(&0) + counts.get(&99).unwrap_or(&0);
    let a_counts = counts.get(&65).unwrap_or(&0) + counts.get(&97).unwrap_or(&0);
    let t_counts = counts.get(&84).unwrap_or(&0) + counts.get(&116).unwrap_or(&0);

    (g_counts + c_counts) as f32 / (g_counts + c_counts + a_counts + t_counts) as f32
}

// convert Node Index to segment ID and vice versa
// I rely a lot on this tuple:
// (NodeIndex, usize)
// which stores the node index and it's corresponding segment ID
// I just realise this should 100000% be a hashmap... change that later.
#[derive(Clone, Copy)]
pub struct GFAGraphPair {
    pub node_index: NodeIndex,
    pub seg_id: usize,
}
#[derive(Clone)]
pub struct GFAGraphLookups(pub Vec<GFAGraphPair>);

impl GFAGraphLookups {
    // simple new
    pub fn new() -> Self {
        Self(Vec::new())
    }
    // simple push
    pub fn push(&mut self, other: GFAGraphPair) {
        self.0.push(other);
    }

    // return seg_id from node index
    pub fn node_index_to_seg_id(&self, node_index: NodeIndex) -> Result<usize> {
        let seg_id = &self
            .0
            .iter()
            .find(|e| e.node_index == node_index)
            .with_context(|| {
                format!(
                    "Node index {:?} could not be converted to segment ID",
                    node_index
                )
            })?
            .seg_id;

        Ok(*seg_id)
    }
    // and the reverse operation
    pub fn seg_id_to_node_index(&self, seg_id: usize) -> Result<NodeIndex> {
        let node_index = &self
            .0
            .iter()
            .find(|e| e.seg_id == seg_id)
            .with_context(|| {
                format!(
                    "Segment ID {:?} could not be converted to NodeIndex",
                    seg_id
                )
            })?
            .node_index;

        Ok(*node_index)
    }
}

use std::fmt;

impl fmt::Display for GFAGraphLookups {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut output = String::new();
        output += "\n\tSegment ID's:\n\t";

        let mut seg_ids: String = self
            .0
            .iter()
            .map(|pair| pair.seg_id.to_string() + ", ")
            .collect();
        seg_ids.drain(seg_ids.len() - 2..);

        output += &seg_ids;

        write!(f, "{}\n", output)
    }
}
