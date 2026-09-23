# gfatk

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22870219.svg)](https://doi.org/10.5281/zenodo.22870219)

<p align="center">
    <img width="300" height="132" src="https://www.darwintreeoflife.org/wp-content/themes/dtol/dist/assets/gfx/dtol-logo-w.png">
</p>

> [!IMPORTANT]  
> `oatk` (<a href="https://github.com/c-zhou/oatk">here</a>) is a complete de-novo assembler which built upon the ideas generated from this work, as well as <a href="https://github.com/tolkit/fpma">fpma</a> and <a href="https://github.com/tolkit/fppa">fppa</a>.
> Work still continues here in the space of GFA linearisation, which is an open problem in complex cases, and general maintenance.

A command line utility to explore, extract, and linearise plant mitochondrial assemblies. The Graphical Fragment Assembly files (GFA's) used to refine the code in this repository are almost exclusively generated from the assembly program <a href="https://github.com/maickrau/MBG">`MBG`</a> or `oatk`. See the testing section below for caveats.

## Install

Grab from the releases (Mac & Linux only):

```bash
# for mac (Intel or Apple Silicon, via Rosetta 2)
curl -L "https://github.com/tolkit/gfatk/releases/download/0.5.3/gfatk_0.5.3_x86_64-apple-darwin.zip" -o gfatk.zip && unzip gfatk.zip && chmod +x gfatk
# and linux (musl, static binary)
curl -L "https://github.com/tolkit/gfatk/releases/download/0.5.3/gfatk_0.5.3_x86_64-unknown-linux-musl.tar.zst" | tar --zstd -xf - && chmod +x gfatk
```

Or build from source.

```bash
# e.g. get rustup!
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# get directly from crates.io
cargo install gfatk

# or clone this repo!
git clone https://github.com/tolkit/gfatk
# cd!
cd gfatk
# build!
cargo build --release
# or install into your path!
cargo install --path .
```

## Features

The features of the toolkit reflect only their usefulness in debugging, visualising, and linearising GFA's from (especially) plant mitochondrial genome assemblies output from <a href="https://github.com/maickrau/MBG">`MBG`</a> or `oatk`. These genomes are usually pretty small (up to 2Mb), and in many cases have circular or branching paths.

Current help:

```
Explore and linearise (plant organellar) GFA files.

Usage: gfatk [COMMAND]

Commands:
  overlap         Extract overlaps from a GFA.
  extract         Extract subgraph from a GFA, given a segment name.
  linear          Force a linear representation of the graph.
  fasta           Extract a fasta file.
                      Almost as simple as: awk '/^S/{print ">"$2"\n"$3}'.
  stats           Some stats about the input GFA.
  extract-mito    Extract the mitochondria from a GFA.
  extract-chloro  Extract the plastid from a GFA.
  dot             Return the dot representation of a GFA.
  trim            Trim a GFA to remove nodes of degree < 4 (i.e. only has one neighbour).
  path            Supply an input path to evaluate a linear representation of.
                  Input must be a text file of a single comma separated line with node ID's and orientations. E.g. 1+,2-,3+
  resolve         Resolve a repeat-rich assembly graph into a circular genome using PacBio HiFi read-path (GAF) evidence.
  rename          Rename the segment ID's of a GFA.
  help            Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

To explain each of these briefly:

- `gfatk dot <GFA>` - generates a `<a href="https://graphviz.org/doc/info/lang.html">`DOT language`</a>` representation of the GFA.
- `gfatk extract <GFA> -s <segment-ids> -i <iterations>` - extracts the subgraph from the GFA, given a segment name, or multiple (if multiple, these must be comma separated without space). Number of iterations may need to be increased for large graphs.
- `gfatk extract-chloro <GFA>` - extracts the plastid from the GFA. It has default parameters which seem to work okay.
- `gfatk extract-mito <GFA>` - extracts the mitochondria from the GFA. It has default parameters which seem to work okay.
- `gfatk fasta <GFA>` - extracts a fasta file from the GFA. This simply prints each of the segments from the GFA. I say it's almost as simple as the `awk` version, but the toolkit does some checks to see if we are actually dealing with a GFA or not.
- `gfatk linear <GFA> -e -i -n <node-threshold>` - forces the longest linear legal representation of the graph, via a budget-bounded greedy walk. You can evaluate within subgraphs (`-e`), or include node coverage information (`-i`).
- `gfatk overlap <GFA> -s <size>` - extracts the overlaps from the GFA. These are taken from the CIGAR string from each of the links, and optionally extended (e.g. `-s 1000` to 1000bp either side of the overlap).
- `gfatk path <GFA> <path> (-p path/to/path.txt)` - evaluates a linear representation of the graph, given an input path. The input path can be on the command line, or a file. Simply, it must be an comma separated list of node ID's and orientations (1+,2-,3+ ... ).
- `gfatk resolve <GFA> --gaf <reads.gaf> --gff <annot.gff> --bed <annot.bed>` - resolves a repeat-rich assembly graph into a circular genome using PacBio HiFi read-path (GAF) evidence, rather than assembly-graph coverage alone. Segment copy number is estimated directly from real read depth, and an integer program selects the best-supported, single-circuit set of graph edges. Segments left out of the main circuit, but with direct read evidence connecting them to it, are reported separately as candidate recombination-mediated structural variants ("bubble arms") rather than discarded -- pass `--no-bubble-fasta` to keep them out of the FASTA output (they're still reported, with their supporting evidence, on stderr). Each bubble arm is also compared, where a genuinely local, size-comparable alternative can be found, to whatever the main path does between the same two attachment points, reporting an edit distance and percent identity -- useful for telling a near-identical minor variant apart from a real, structurally distinct alternative. If the graph is too tangled to prove a single circuit within `--time-limit` seconds (default 45), any remaining separate circuits are automatically grafted back together wherever real read evidence supports it; stderr explains why any that are left over couldn't be.
- `gfatk stats <GFA> -t` - some stats about the input GFA. Can be quite verbose for large, unconnected graphs. `-t` outputs tabular data (TSV).
- `gfatk trim <GFA>` - removes segments if they have only a single neighbour. Useful for trimming GFA's which have segments attached at low coverage.

These are not all the options for each subcommand. Run:

`gfatk help <subcommand>` for more information on each subcommand.

Many of these commands can be chained in a pipeline, e.g. `gfatk extract-chloro in.gfa | gfatk linear > out.fa`.

## Examples and docs

A couple of more detailed examples can be seen in the `examples` directory, where there is a `README.md` file. To view the auto-generated documentation of the binary itself, including details of all underlying functions, see:

<p align="center">
    <b>
        <a href="https://tolkit.github.io/gfatk/">API documentation</a>
    </b>
</p>

## Requirements and testing

Some unit tests are now provided in the `tests` directory. To run these (you'll need Rust):

```bash
cargo test --release
```

For full functionality of the toolkit, two tags are required, node coverage and edge coverage. Other functionality will fail if the CIGAR string is not purely an overlap; i.e. in the format `<integer>M`. Only GFA version 1 supported. Only header (`H`), segment (`S`), and link (`L`) lines are required. `P` lines are used in `gfatk path --all <GFA>`.

`gfatk resolve` additionally requires a GAF file of long reads aligned back to the assembly graph (e.g. from <a href="https://github.com/maickrau/GraphAligner">`GraphAligner`</a>), and optionally a gene annotation to weight gene-bearing segments during resolution -- either a GFF (`--gff`) or, directly, an oatk-style BED file (`--bed`, columns `seq_name align_from align_to gene_name score_capped_at_1000 strand`, `seq_name` referencing the GFA's own segment IDs, no conversion step needed). A BED hit only counts as a confident gene call once its score reaches the 1000 cap, the same standard `--gff` applies via `coverage=full`. Both can be given together; the resulting gene-bearing segments are combined.

```
H	VN:Z:1.0
S	11	ACCTT	ll:f:30.0 <- this tag indicates node/segment coverage (here it's 30.0)
S	12	TCAAGG	ll:f:60.0
S	13	CTTGATT	ll:f:30.0
L	11	+	12	-	4M	ec:i:1 <- this tag indicates edge coverage (here it's 1; EC:i is also accepted)
L	12	-	13	+	5M	ec:i:1
L	11	+	13	+	3M	ec:i:1
L	12	+	11	-	4M	ec:i:1
L	13	-	12	+	5M	ec:i:1
L	13	-	11	-	3M <- simple overlap on the CIGAR string (overlap == 3)	ec:i:1

```

## Thanks

Many thanks to the developers of MBG, and partners in the Tree of Life program, and beyond:

- Marcela Uliano-Silva
- Sergey Nurk
- Alex Twyford
- Lucia Campos
- Chenxi Zhou
- Mark Blaxter
