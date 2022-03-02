use crate::gfa::gfa::GFAtk;
use crate::gfa::graph::segments_subgraph;
use crate::load::load_gfa;
use anyhow::Result;

// so we can extract the segments with highest GC content.
#[derive(Clone)]
pub struct Stat {
    pub index: usize,
    pub gc: f32,
    pub cov: f32,
    pub segments: Vec<usize>,
    pub total_sequence_length: usize,
}

pub struct Stats(pub Vec<Stat>);

impl Stats {
    pub fn push(&mut self, stat: Stat) {
        let stats = &mut self.0;
        stats.push(stat);
    }
    // we extract the mito by ordering our stats
    // by gc content & coverage.
    pub fn extract_mito(&mut self) -> Vec<usize> {
        let stat_vec = &mut self.0;
        // reverse the cov..
        stat_vec.sort_by(|a, b| (a.cov, b.gc).partial_cmp(&(b.cov, a.gc)).unwrap());

        let stat_vec_len = stat_vec.len();

        if stat_vec_len > 1 {
            // now check that the length is sufficiently high
            let z: Vec<&Stat> = stat_vec
                .iter()
                .filter(|e| e.total_sequence_length > 100_000)
                .collect();
            z[0].segments.clone()
        } else {
            stat_vec[0].segments.clone()
        }
    }
}

// I've handled 'further' here really badly...
// I want node indices & segment names printed too (maybe optionally.)
pub fn stats(matches: &clap::ArgMatches, further: bool) -> Result<Option<Vec<usize>>> {
    // required so unwrap safely
    let gfa_file = matches.value_of("gfa").unwrap();

    let gfa: GFAtk = GFAtk(load_gfa(gfa_file)?);

    // load gfa into graph structure
    let (graph_indices, gfa_graph) = gfa.into_digraph()?;

    let subgraphs = gfa_graph.weakly_connected_components(graph_indices)?;

    let mut no_subgraphs = 0;
    let mut store_stats = Stats(Vec::new());

    for id_set in &subgraphs {
        let subgraph_gfa = GFAtk(segments_subgraph(&gfa.0, id_set.to_vec()));

        let (graph_indices_subgraph, subgraph) = subgraph_gfa.into_digraph()?;

        // print stats
        if !further {
            println!("Subgraph {}:", no_subgraphs + 1);
            println!("\tNumber of nodes/segments: {}", subgraph.node_count());
            println!("\tNumber of edges/links: {}", subgraph.edge_count());
            // equivalent to id_set
            println!("{}", graph_indices_subgraph);
        }
        let (avg_gc, cov, total_sequence_length) = subgraph_gfa.sequence_stats(further)?;

        store_stats.push(Stat {
            index: no_subgraphs,
            gc: avg_gc,
            cov,
            segments: id_set.clone(),
            total_sequence_length,
        });
        no_subgraphs += 1;
    }

    // if we want to do more stat things
    if further {
        return Ok(Some(store_stats.extract_mito()));
    } else {
        println!("Total number of subgraphs: {}", no_subgraphs);
    }

    Ok(None)
}
