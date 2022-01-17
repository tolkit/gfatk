use gfa::optfields::{OptField, OptFieldVal::*};
use std::collections::HashMap;

pub fn get_option_string(options: Vec<OptField>) -> String {
    let mut tag_val = String::new();
    for op in options {
        let tag = std::str::from_utf8(&op.tag).unwrap();
        let value = match op.value {
            Float(f) => format!(":f:{:.3}", f),
            A(a) => format!(":A:{}", a.to_string()),
            Int(i) => format!(":i:{}", i.to_string()),
            Z(z) => format!(":Z:{}", std::str::from_utf8(&z).unwrap()),
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
    let tag_val_op_un = tag_val.strip_suffix("\t").unwrap();
    tag_val_op_un.to_string()
}

// not a very safe function, but works
pub fn parse_cigar(cigar: &[u8]) -> usize {
    // check it ends with an M
    if !cigar.ends_with(&[77]) {
        panic!("CIGAR did not end with M.");
    }
    let stripped = cigar.strip_suffix(&[77]);
    match stripped {
        Some(s) => {
            let string_rep = std::str::from_utf8(s).unwrap();
            string_rep.parse::<usize>().unwrap()
        }
        None => panic!(""),
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
