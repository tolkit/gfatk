use crate::gfa::gfa::GFAtk;
use crate::load::load_gfa;
use clap::value_t;

// TODO: I haven't added in paths/containments etc...
// just segements and links.
// also the optional extra bits are useful - e.g. sequence length
// if they are present in a gfa. Work out the best bits to keep.

pub fn extract(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").unwrap();
    let sequence_id = value_t!(matches.value_of("sequence-id"), usize).unwrap_or_else(|e| e.exit());
    let iterations = value_t!(matches.value_of("iterations"), i32).unwrap_or_else(|e| e.exit());

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    let (graph_indices, gfa_graph) = gfa.into_ungraph();

    // get the node index of the target sequence ID.
    let target_index = graph_indices.iter().find(|x| x.1 == sequence_id).unwrap().0;
    let mut collect_sequence_names = Vec::new();
    // and this is the first item in this vector
    collect_sequence_names.push(target_index);

    let sequences_to_keep = gfa_graph.recursive_search(
        sequence_id,
        iterations,
        collect_sequence_names,
        graph_indices,
    );

    gfa.print_extract(sequences_to_keep);

    Ok(())
}
