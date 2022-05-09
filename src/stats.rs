use crate::gfa::graph::segments_subgraph;
use crate::load::load_gfa;
use crate::utils;
use crate::{gfa::gfa::GFAtk, load::load_gfa_stdin};
use anyhow::{bail, Result};

/// Enumeration of the genomes we are interested in.
#[derive(PartialEq, Clone, Copy)]
pub enum GenomeType {
    /// The mitochondrial genome
    Mitochondria,
    /// The chloroplast/plastid genome
    Chloroplast,
    /// This will process the stats and return nothing
    None,
}

/// The statistics associated with a subgraph in a GFA.
#[derive(Clone)]
pub struct Stat {
    /// Arbitrary index of the subgraph(s).
    pub index: usize,
    /// The average GC% across a subgraph.
    pub gc: f32,
    /// The average coverage across a subgraph.
    pub cov: f32,
    /// Names of the segments.
    pub segments: Vec<usize>,
    /// Total sequence length of all the segments.
    pub total_sequence_length: usize,
}

/// A vector of `Stat`.
pub struct Stats(pub Vec<Stat>);

impl Stats {
    /// Add a new `Stat` to `Stats`.
    pub fn push(&mut self, stat: Stat) {
        let stats = &mut self.0;
        stats.push(stat);
    }
    // we extract the mito by ordering our stats
    // by gc content & coverage.

    /// Extract the putative chloroplast genome from a GFA
    /// file.
    ///
    /// The upper and lower limits of genome size and GC content are supplied through the
    /// CLI. Defaults should probably be

    #[allow(unused_variables)]
    pub fn extract_chloro(
        &mut self,
        size_lower: usize,
        mut size_upper: usize,
        gc_lower: f32,
        gc_upper: f32,
    ) -> Result<Vec<usize>> {
        // just going to hard code these for the moment
        // these values are taken from GoaT
        // these values are within 2 stddevs of the mean,
        // so most chloroplasts should pop out

        // adjust because of overlaps between segments
        // kind of arbitrary...
        let seq_len_adj = 20000;
        size_upper += seq_len_adj;

        let stat_vec = &mut self.0;
        let stat_vec_len = stat_vec.len();
        // filter this vector to have stats in line with the span/gc

        if stat_vec_len > 1 {
            // apply the filter
            let mut z: Vec<&Stat> = stat_vec
                .iter()
                .filter(
                    |Stat {
                         index,
                         gc,
                         cov,
                         segments,
                         total_sequence_length,
                     }| {
                        (gc > &gc_lower && gc < &gc_upper)
                            && (total_sequence_length > &size_lower
                                && total_sequence_length < &size_upper)
                    },
                )
                .collect();
            // We could just return the complete set of deduped segments
            // that were found
            z.sort_by(|a, b| (b.cov, a.gc).partial_cmp(&(a.cov, b.gc)).unwrap());

            let res = match z.get(0) {
                Some(stat) => stat,
                None => bail!(
                    "No subgraphs within the bounds:\nsize_upper: {size_upper}\nsize_lower: {size_lower}\ngc_upper: {gc_upper}\ngc_lower: {gc_lower}\nTry changing limits?"
                ),
            };
            Ok(res.segments.clone())
        } else {
            // we only have 1 or no segments
            let extracted_segments_op = stat_vec.get(0);
            let extracted_segments = match extracted_segments_op {
                Some(s) => s,
                None => bail!("There were no segments to be extracted. Check input GFA file."),
            };
            Ok(extracted_segments.segments.clone())
        }
    }

    /// The function called from `gfatk extract-mito`.
    ///
    /// Extracts the putative mitochondrial subgraph from a GFA.
    pub fn extract_mito(&mut self, size: usize) -> Result<Vec<usize>> {
        let stat_vec = &mut self.0;
        // reverse the cov..
        stat_vec.sort_by(|a, b| (a.cov, b.gc).partial_cmp(&(b.cov, a.gc)).unwrap());

        let stat_vec_len = stat_vec.len();

        if stat_vec_len > 1 {
            // now check that the length is sufficiently high
            let z: Vec<&Stat> = stat_vec
                .iter()
                // hardcoded for now, but does not have to be.
                .filter(|e| e.total_sequence_length > size)
                .collect();
            let res = match z.get(0) {
                Some(stat) => stat,
                None => bail!(
                    "No subgraphs with minimal sequence length of {}. Maybe reduce <size>?",
                    size
                ),
            };
            Ok(res.segments.clone())
        } else {
            // TODO: make a better error message here?
            let extracted_segments_op = stat_vec.get(0);
            let extracted_segments = match extracted_segments_op {
                Some(s) => s,
                None => bail!("There were no segments to be extracted. Check input GFA file."),
            };
            Ok(extracted_segments.segments.clone())
        }
    }
}

// I've handled 'further' here really badly...
// I want node indices & segment names printed too (maybe optionally.)

/// Internal function called in `gfatk stats`.
///
/// Used in `gfatk stats`, `gfatk extract-mito`, and `gfatk extract-chloro`.
///
/// For example:
/// ```bash
/// gfatk stats in.gfa
/// ```
pub fn stats(
    matches: &clap::ArgMatches,
    genome_type: GenomeType,
) -> Result<Option<(GFAtk, Vec<usize>)>> {
    // required so unwrap safely
    let gfa_file = matches.value_of("GFA");
    // only passed through extract_mito
    let mito_args = if matches!(genome_type, GenomeType::Mitochondria) {
        let size: usize = matches.value_of_t("size")?;
        Some(size)
    } else {
        None
    };
    // only required for extract_chloro
    let chloro_args = if matches!(genome_type, GenomeType::Chloroplast) {
        let size_lower: usize = matches.value_of_t("size-lower")?;
        let size_upper: usize = matches.value_of_t("size-upper")?;
        let gc_lower: f32 = matches.value_of_t("gc-lower")?;
        let gc_upper: f32 = matches.value_of_t("gc-upper")?;
        Some((size_lower, size_upper, gc_lower, gc_upper))
    } else {
        None
    };

    let error_string = match genome_type {
        GenomeType::Chloroplast => "extract-chloro",
        GenomeType::Mitochondria => "extract-mito",
        GenomeType::None => "stats",
    };

    let gfa = match gfa_file {
        Some(f) => {
            if !f.ends_with(".gfa") {
                bail!("Input file is not a GFA.")
            }

            GFAtk(load_gfa(f)?)
        }
        None => match utils::is_stdin() {
            true => GFAtk(load_gfa_stdin(std::io::stdin().lock())?),
            false => bail!(
                "No input from STDIN. Run `gfatk {} -h` for help.",
                error_string
            ),
        },
    };

    // load gfa into graph structure
    let (graph_indices, gfa_graph) = gfa.into_digraph()?;

    let subgraphs = gfa_graph.weakly_connected_components(graph_indices)?;

    let mut no_subgraphs = 0;
    let mut store_stats = Stats(Vec::new());

    for id_set in &subgraphs {
        let subgraph_gfa = GFAtk(segments_subgraph(&gfa.0, id_set.to_vec()));

        let (graph_indices_subgraph, subgraph) = subgraph_gfa.into_digraph()?;

        // print stats
        if genome_type == GenomeType::None {
            println!("Subgraph {}:", no_subgraphs + 1);
            println!("\tNumber of nodes/segments: {}", subgraph.node_count());
            println!("\tNumber of edges/links: {}", subgraph.edge_count());
            // equivalent to id_set
            println!("{}", graph_indices_subgraph);
        }
        let (avg_gc, cov, total_sequence_length) = subgraph_gfa.sequence_stats(genome_type)?;

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
    match genome_type {
        GenomeType::Mitochondria => {
            return Ok(Some((gfa, store_stats.extract_mito(mito_args.unwrap())?)))
        }
        GenomeType::Chloroplast => {
            let chloro_args = chloro_args.unwrap();
            return Ok(Some((
                gfa,
                store_stats.extract_chloro(
                    chloro_args.0,
                    chloro_args.1,
                    chloro_args.2,
                    chloro_args.3,
                )?,
            )));
        }
        GenomeType::None => println!("Total number of subgraphs: {}", no_subgraphs),
    }

    Ok(None)
}
