use crate::gfa::gfa::GFAtk;
use crate::gfa::graph::segments_subgraph;
// use crate::gfa::writer;
use crate::load::load_gfa;

pub fn stats(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // required so unwrap safely
    let gfa_file = matches.value_of("gfa").unwrap();

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    // load gfa into graph structure
    let (graph_indices, gfa_graph) = gfa.into_digraph();

    let subgraphs = gfa_graph.weakly_connected_components(graph_indices);

    let mut no_subgraphs = 0;
    for id_set in subgraphs {
        let subgraph_gfa = GFAtk(segments_subgraph(&gfa.0, id_set));

        let (_, subgraph) = subgraph_gfa.into_digraph();

        // print stats
        println!("Subgraph {}:", no_subgraphs + 1);
        println!("\tNumber of nodes: {}", subgraph.node_count());
        println!("\tNumber of edges: {}", subgraph.edge_count());
        subgraph_gfa.sequence_stats();

        // that appears to work
        // println!("{}\n", writer::gfa_string(&subgraph_gfa.0));
        no_subgraphs += 1;
    }

    println!("Total number of subgraphs: {}", no_subgraphs);

    Ok(())
}
