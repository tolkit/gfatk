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
