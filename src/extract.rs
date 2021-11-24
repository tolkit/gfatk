use crate::load::load_gfa;
use crate::utils::get_option_string;
use clap::value_t;
use gfa::gfa::GFA;
use gfa::optfields::OptionalFields;
use petgraph::graph::{Graph, UnGraph};

// parsing, editing, and making V1 GFA's
const HEADER: &str = "H\tVN:Z:1.0";

// TODO: I haven't added in paths/containments etc...
// just segements and links.
// also the optional extra bits are useful - e.g. sequence length
// if they are present in a gfa. Work out the best bits to keep.

pub fn extract(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").unwrap();
    let sequence_id = value_t!(matches.value_of("sequence-id"), usize).unwrap_or_else(|e| e.exit());
    let iterations = value_t!(matches.value_of("iterations"), usize).unwrap_or_else(|e| e.exit());

    let gfa: GFA<usize, OptionalFields> = load_gfa(gfa_file)?;

    eprintln!("[+]\tReading GFA into an undirected graph.");
    let mut gfa_graph: UnGraph<usize, ()> = Graph::new_undirected();

    let mut graph_indices = Vec::new();
    // read the segments into graph nodes
    // save the indexes for populating the edges
    for node in &gfa.segments {
        let index = gfa_graph.add_node(node.name);
        graph_indices.push((index, node.name));
    }

    // populate the edges
    for edge in &gfa.links {
        let from = edge.from_segment;
        let to = edge.to_segment;

        // get the node index for a given edge (like a map)
        let from_index = graph_indices.iter().find(|x| x.1 == from).unwrap().0;
        let to_index = graph_indices.iter().find(|x| x.1 == to).unwrap().0;

        // add the edges
        gfa_graph.add_edge(from_index, to_index, ());
    }

    // get the node index of the target sequence ID.
    let target_index = graph_indices.iter().find(|x| x.1 == sequence_id).unwrap().0;
    let mut collect_sequence_names = Vec::new();
    // and this is the first item in this vector
    collect_sequence_names.push(target_index);

    eprintln!(
        "[+]\tRecursively searching around node {} at depth {}",
        sequence_id, iterations
    );

    let mut iteration = 0;
    loop {
        for index in collect_sequence_names.clone() {
            for c in gfa_graph.neighbors(index) {
                let res = graph_indices.iter().find(|x| x.0 == c).unwrap().0;
                collect_sequence_names.push(res);
            }
        }
        // add sorting and deduping here too
        // yes otherwise vectors are enormous.
        collect_sequence_names.sort();
        collect_sequence_names.dedup();

        iteration += 1;
        if iteration == iterations {
            break;
        }
    }
    collect_sequence_names.sort();
    collect_sequence_names.dedup();

    let mut sequences_to_keep = Vec::new();
    // turn node indexes into sequence ID's
    for index in collect_sequence_names {
        let t = graph_indices.iter().find(|x| x.0 == index).unwrap().1;
        sequences_to_keep.push(t);
    }

    eprintln!("[+]\tGenerating GFA subgraph");

    println!("{}", HEADER);
    // segments are easy to print.
    // if they match our sequence ID's we print them
    for segment in &gfa.segments {
        let name = segment.name;
        // hacky options parsing.
        // it works but needs to me made much better.
        let options = segment.optional.clone();
        let tag_val = get_option_string(options);

        if sequences_to_keep.contains(&name) {
            println!(
                "S\t{}\t{}\t{}",
                segment.name,
                std::str::from_utf8(&segment.sequence).unwrap(),
                tag_val,
            )
        }
    }

    // keep track of the segment pairs to avoid unneccesary
    // printing of duplicate pairs in opposite orientations.
    let mut keep_track_pairs = Vec::new();

    for link in &gfa.links {
        let from = link.from_segment;
        let to = link.to_segment;
        let options = link.optional.clone();
        let tag_val = get_option_string(options);

        if !keep_track_pairs.contains(&(from, to)) || !keep_track_pairs.contains(&(to, from)) {
            if sequences_to_keep.contains(&from) || sequences_to_keep.contains(&to) {
                println!(
                    "L\t{}\t{}\t{}\t{}\t{}\t{}",
                    from,
                    link.from_orient,
                    to,
                    link.to_orient,
                    std::str::from_utf8(&link.overlap).unwrap(),
                    tag_val
                )
            }
        }
        // keep only unique pairs
        keep_track_pairs.push((from, to));
        keep_track_pairs.push((to, from));
    }

    Ok(())
}
