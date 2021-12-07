use bstr::ByteSlice;
use gfa::gafpaf::GAFStep;
use gfa::{gafpaf::parse_gaf, gafpaf::GAFPath, optfields::OptionalFields};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

// can we get it to print out in a way like this:
// id   (+/-)   id2 (+/-)   coverage
// I think we can.
// this can then be passed to gfatk linear for
// finding the most covered path through the gfa.

// > is forward, < is reverse.

type GAF = gfa::gafpaf::GAF<OptionalFields>;

pub fn gaf_to_graph(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    let gaf_filename = matches.value_of("gaf").unwrap();
    let directed = true;

    let gaf_file = File::open(gaf_filename).unwrap();
    let mut lines = BufReader::new(gaf_file);
    let mut line: Vec<u8> = Vec::new();

    let mut hash = HashMap::new();

    loop {
        // clear the line
        line.clear();
        // read the bytes
        let bytes_read = lines.read_until(0xA, &mut line);
        // break if there's an error
        if bytes_read.is_err() || bytes_read.unwrap() == 0 {
            break;
        }
        // split the fields
        let fields: bstr::Split = line[0..line.len()].split_str(b"\t");
        // now we get to the juicy bit
        if let Some::<GAF>(gaf) = parse_gaf(fields) {
            let path = gaf.path;
            match path {
                // the ID of a stable rGFA identifier
                // what is this??
                GAFPath::StableId(e) => {
                    eprintln!("{}", std::str::from_utf8(&e).unwrap())
                }
                // a list of oriented steps (what I think we are interested in)
                GAFPath::OrientIntv(e) => {
                    let len = e.len();
                    let mut vector_to_tuple = Vec::new();
                    // only interested in overlaps between two segments
                    if len > 1 && len < 3 {
                        // sort the tuple
                        // why are we just unwrapping this .get()?
                        let x = match e.get(0).unwrap() {
                            GAFStep::SegId(o, id) => {
                                let y = std::str::from_utf8(id).unwrap();
                                let orientation = match o.plus_minus_as_byte() {
                                    43 => "+",
                                    45 => "-",
                                    _ => "",
                                };
                                (orientation, y.parse::<usize>().unwrap())
                            }
                            _ => panic!("No segment ID"),
                        };
                        let y = match e.get(1).unwrap() {
                            GAFStep::SegId(o, id) => {
                                let y = std::str::from_utf8(id).unwrap();
                                let orientation = match o.plus_minus_as_byte() {
                                    43 => "+",
                                    45 => "-",
                                    _ => "",
                                };
                                (orientation, y.parse::<usize>().unwrap())
                            }
                            _ => panic!("No segment ID"),
                        };

                        vector_to_tuple.push(x);
                        vector_to_tuple.push(y);
                        // undirected might be desirable in future?
                        if !directed {
                            vector_to_tuple.sort();
                        }
                        // add to map
                        *hash
                            .entry((vector_to_tuple[0], vector_to_tuple[1]))
                            .or_insert(0) += 1;
                    } else if len > 3 {
                        // if we want to do something with
                        // the paths which go through multiple segments...
                    }
                }
            }
        }
    }

    // headers
    println!("from_orient\tfrom\tto_orient\tto\tcoverage");
    // body
    for ((k1, k2), v) in hash {
        println!("{}\t{}\t{}\t{}\t{}", k1.0, k1.1, k2.0, k2.1, v);
    }

    Ok(())
}
