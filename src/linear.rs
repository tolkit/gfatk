// the idea of this is to take
// a path through the gfa which includes all
// segments only once

// I think load into an DiGraph (as strand matters)?
// how?

// need a digraph
// weight is strandedness

// maybe I can filter the links
// so only one strand is left
// then apply all simple paths

use crate::load::load_gfa;
use crate::utils::{parse_cigar, reverse_complement};
use gfa::gfa::{Orientation, GFA};
use gfa::optfields::OptionalFields;
use indexmap::IndexMap;
use petgraph::algo::all_simple_paths;
use petgraph::algo::is_cyclic_directed;
use petgraph::graph::{Graph, NodeIndex};

// debugging
use petgraph::visit::{EdgeRef, NodeIndexable};

pub fn force_linear(matches: &clap::ArgMatches) -> Result<(), Box<dyn std::error::Error>> {
    // read in path and parse gfa
    let gfa_file = matches.value_of("gfa").unwrap();

    let gfa: GFA<usize, OptionalFields> = load_gfa(gfa_file)?;

    eprintln!("[+]\tReading GFA into a directed graph.");
    let mut gfa_graph: Graph<usize, (Orientation, Orientation)> = Graph::new();

    let mut graph_indices = Vec::new();
    // read the segments into graph nodes
    // save the indexes for populating the edges
    for node in &gfa.segments {
        let index = gfa_graph.add_node(node.name);
        graph_indices.push((index, node.name));
    }

    // populate the edges
    let mut graph_edge_weights = Vec::new();
    for edge in &gfa.links {
        let from = edge.from_segment;
        let to = edge.to_segment;
        let from_orient = edge.from_orient;
        let to_orient = edge.to_orient;

        // overlap
        let overlap = parse_cigar(&edge.overlap);

        // here filter out one strand...

        // get the node index for a given edge (like a map)
        let from_index = graph_indices.iter().find(|x| x.1 == from).unwrap().0;
        let to_index = graph_indices.iter().find(|x| x.1 == to).unwrap().0;

        // add the edges
        // conditional on strandedness
        if from_orient == to_orient {
            gfa_graph.add_edge(from_index, to_index, (from_orient, to_orient));
            // populate weights too
            graph_edge_weights.push((
                (from_index, from_orient),
                (to_index, to_orient),
                overlap,
                "start",
            ));
        } else {
            gfa_graph.add_edge(to_index, from_index, (to_orient, from_orient));
            // populate weights too
            graph_edge_weights.push((
                (to_index, from_orient),
                (from_index, to_orient),
                overlap,
                "end",
            ));
        }
    }

    // check if graph is cyclic
    // do we want to check this? Does it need to be cyclic?
    match is_cyclic_directed(&gfa_graph) {
        true => eprintln!("[+]\tLoaded GFA is a cyclic digraph."),
        false => panic!("Loaded GFA is not a cyclic digraph"),
    }

    // iterate over the nodes
    // just increment until there are no more nodes.
    let mut node_index = 0;
    // collect nodes and their edge counts
    let mut node_edge_counts = Vec::new();
    loop {
        // need to peek ahead so we don't loop forever
        let mut peekable_iter = gfa_graph
            .neighbors_undirected(NodeIndex::new(node_index))
            .peekable();

        // need a different way of breaking here...
        match peekable_iter.peek().is_some() {
            true => {
                let mut node_count = 0;
                for _ in gfa_graph.neighbors_undirected(NodeIndex::new(node_index)) {
                    node_count += 1;
                }
                node_edge_counts.push((node_index, node_count));
                node_index += 1;
            }
            false => {
                if node_index == gfa_graph.node_bound() {
                    break;
                } else {
                    node_index += 1;
                    continue;
                }
            }
        }
    }

    // if there are any nodes with fewer than 2 connections
    // get rid of these, as they are not part of the cycle.
    let filtered_nodes = node_edge_counts
        .iter()
        .filter(|e| e.1 > 1)
        .collect::<Vec<&(usize, i32)>>();

    // for starting node and ending node
    // if the graph is cyclic but not filtered
    // might get stuck in the wrong subgraph
    let min_intermediate_nodes = filtered_nodes.len() - 2;

    // no max intermediate nodes
    let max_intermediate_nodes = None;
    // calculate adjecent nodes, each of which have two neighbours
    // iterate over edges again but stop when we hit a neighbouring
    // pair of nodes which are connected each themselves to at least 2 others
    node_index = 0;
    let mut source_target_pair: Option<(NodeIndex, NodeIndex)> = None;

    loop {
        // need to peek ahead so we don't loop forever
        let mut peekable_iter = gfa_graph.edges(NodeIndex::new(node_index)).peekable();
        match peekable_iter.peek().is_some() {
            true => {
                for edge in gfa_graph.edges(NodeIndex::new(node_index)) {
                    let source = edge.source();
                    let target = edge.target();

                    let source_id = graph_indices.iter().find(|x| x.0 == source).unwrap().0;
                    let target_id = graph_indices.iter().find(|x| x.0 == target).unwrap().0;

                    let source_edge_count = node_edge_counts
                        .iter()
                        .filter(|e| NodeIndex::new(e.0) == source_id)
                        .collect::<Vec<&(usize, i32)>>()[0]
                        .1;
                    let target_edge_count = node_edge_counts
                        .iter()
                        .filter(|e| NodeIndex::new(e.0) == target_id)
                        .collect::<Vec<&(usize, i32)>>()[0]
                        .1;

                    node_index += 1;

                    if source_edge_count > 1 && target_edge_count > 1 {
                        source_target_pair = Some((source, target));
                        // and break early here
                        break;
                    }
                }
            }
            false => {
                if node_index == gfa_graph.node_bound() {
                    break;
                } else {
                    node_index += 1;
                    continue;
                }
            }
        }
    }

    if source_target_pair.is_none() {
        panic!("Did not find a pair of adjacent segments, each themselves with connections.");
    }

    // enumerating all simple paths is computationally expensive
    let mut paths = all_simple_paths::<Vec<_>, _>(
        &gfa_graph,
        // i.e. from target -> source in a digraph
        source_target_pair.unwrap().1,
        source_target_pair.unwrap().0,
        min_intermediate_nodes,
        max_intermediate_nodes,
    )
    .collect::<Vec<_>>();

    // if paths length is zero, bail out here.
    if paths.is_empty() {
        panic!("No simple paths were found.")
    }
    // sort & dedup paths and pick longest one
    paths.sort_by(|a, b| b.len().cmp(&a.len()));

    // sort out the path now
    // we need the strandedness information
    let chosen_path = &paths[0];

    // another new vec of (id, overlap, start/end)
    // to sort out for fasta generation
    let mut chosen_path_overlaps = Vec::new();

    for edge in &gfa.links {
        let from = edge.from_segment;
        let to = edge.to_segment;
        let from_orient = edge.from_orient;
        let to_orient = edge.to_orient;

        // println!("{} {} -> {} {}", from, from_orient, to, to_orient);

        // overlap
        let overlap = parse_cigar(&edge.overlap);
        // here look at the strandedness between pairs of nodes
        let path_pairs = chosen_path.windows(2);
        for path in path_pairs {
            let from_path = graph_indices.iter().find(|y| y.0 == path[0]).unwrap().1;
            let to_path = graph_indices.iter().find(|y| y.0 == path[1]).unwrap().1;

            // okay this appears to work.
            if from == from_path && to == to_path {
                chosen_path_overlaps.push((from, from_orient, overlap, "end"));
                chosen_path_overlaps.push((to, to_orient, overlap, "start"));
            } else if from == to_path && to == from_path {
                chosen_path_overlaps.push((from, from_orient, overlap, "start"));
                chosen_path_overlaps.push((to, to_orient, overlap, "end"));
            }
        }
    }

    chosen_path_overlaps.sort();
    chosen_path_overlaps.dedup();
    // now sort these overlaps so they are the same order
    // chosen path -> id's
    let chosen_path_ids = chosen_path
        .iter()
        .map(|e| graph_indices.iter().find(|y| y.0 == *e).unwrap().1)
        .collect::<Vec<_>>();

    let mut sorted_chosen_path_overlaps = Vec::new();
    // sorting (allocate to new vec)
    for path_id in chosen_path_ids {
        for path_id2 in &chosen_path_overlaps {
            if path_id2.0 == path_id {
                sorted_chosen_path_overlaps.push(path_id2)
            }
        }
    }
    // println!("{:?}", sorted_chosen_path_overlaps);

    // format string just for my (and user) interest
    let chosen_path_string = chosen_path
        .iter()
        .map(|e| format!("{}", graph_indices.iter().find(|y| y.0 == *e).unwrap().1))
        .collect::<Vec<String>>();
    eprintln!(
        "[+]\tChosen path through graph: {}",
        chosen_path_string.join(" -> ")
    );

    // merge constrand on the ids
    // as order is critical, we use an IndexMap
    let mut merged_constrand = IndexMap::new();
    for (id, orientation, overlap, side) in sorted_chosen_path_overlaps {
        merged_constrand
            .entry(id)
            .or_insert(Vec::new())
            .push((orientation, overlap, side));
    }

    // println!("{:?}", merged_constrand);

    // finally make the fasta.
    // iterate over the segments

    println!(">TEST");

    for (id, vector) in merged_constrand {
        // get the from and to sequences.
        for line in gfa.lines_iter() {
            // if we meet a segment, let's do something
            match line.some_segment() {
                Some(s) => {
                    // we've got the segment we wanted
                    if s.name == *id {
                        // remove overlap
                        // but we need to watch out for orientation

                        match vector.len() {
                            1 => {
                                // this is either the start or the end.
                                // as the start and end of each inner segment
                                // is stripped, this can be printed as-is
                                print!("{}", std::str::from_utf8(&s.sequence).unwrap());
                            }
                            2 => {
                                // this is an inner segment
                                let start_overlap = vector
                                    .iter()
                                    .find(|(_or, _ov, side)| side == &&"start")
                                    .unwrap();
                                let end_overlap = vector
                                    .iter()
                                    .find(|(_or, _ov, side)| side == &&"end")
                                    .unwrap();

                                // start and end orientation should be the same
                                // hence just matching on start here.
                                let seq_minus_overlap = match start_overlap.0 {
                                    Orientation::Forward => {
                                        // do nothing
                                        let seq_minus_overlap = s
                                            .sequence
                                            .get(*start_overlap.1..s.sequence.len() - end_overlap.1)
                                            .unwrap();
                                        seq_minus_overlap.to_vec()
                                    }
                                    Orientation::Backward => {
                                        let revcomp_seq = reverse_complement(&s.sequence);
                                        let seq_minus_overlap = revcomp_seq
                                            .get(*start_overlap.1..s.sequence.len() - end_overlap.1)
                                            .unwrap();
                                        seq_minus_overlap.to_vec()
                                    }
                                };
                                print!("{}", std::str::from_utf8(&seq_minus_overlap).unwrap());
                            }
                            _ => {
                                // anything else is should not happen
                            }
                        }
                    }
                }
                None => (),
            }
        }
    }

    Ok(())
}
