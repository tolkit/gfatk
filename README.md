# gfatk

<p align="center">
    <img width="300" height="132" src="https://www.darwintreeoflife.org/wp-content/themes/dtol/dist/assets/gfx/dtol-logo-round.png">
</p>

A command line utility to explore, extract, and linearise plant mitochondrial assemblies. The Graphical Fragment Assembly files (GFA's) used to refine the code in this repository are almost exclusively generated from the assembly program <a href="https://github.com/maickrau/MBG">MBG</a>. See the testing section below for caveats.

## Install

Grab from the releases (Mac & Linux only):

```bash
# for mac
curl -L "https://github.com/tolkit/gfatk/releases/download/0.2.16/gfatk_mac_0.2.2" > gfatk && chmod +x gfatk
# and linux (ubuntu)
curl -L "https://github.com/tolkit/gfatk/releases/download/0.2.16/gfatk_ubuntu_0.2.2" > gfatk && chmod +x gfatk
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

The features of the toolkit reflect only their usefulness in debugging, visualising, and linearising GFA's from (especially) plant mitochondrial genome assemblies output from <a href="https://github.com/maickrau/MBG">MBG</a>. These genomes are usually pretty small (up to 2Mb), and in many cases have circular or branching paths.

Current help:

```
gfatk 0.2.2
Max Brown <mb39@sanger.ac.uk>
Explore and linearise (mitochondrial) GFA files.

USAGE:
    gfatk [SUBCOMMAND]

OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

SUBCOMMANDS:
    dot               Return the dot representation of a GFA.
    extract           Extract subgraph from a GFA, given a segment name.
    extract-chloro    Extract the plastid from a GFA.
    extract-mito      Extract the mitochondria from a GFA.
    fasta             Extract a fasta file.
                          Almost as simple as: awk '/^S/{print ">"$2"\n"$3}'.
    help              Print this message or the help of the given subcommand(s)
    linear            Force a linear representation of the graph.
    overlap           Extract overlaps from a GFA.
    path              Supply an input path to evaluate a linear representation of. Input must be
                          a text file of a single comma separated line with node ID's and
                          orientations. E.g.:
                                1+,2-,3+
    stats             Some stats about the input GFA.
    trim              Trim a GFA to remove nodes of degree < 4 (i.e. only has one neighbour).
```

To explain each of these briefly:
- `gfatk dot <GFA>` - generates a <a href="https://graphviz.org/doc/info/lang.html">DOT language</a> representation of the GFA.
- `gfatk extract <GFA> -s <segment-ids> -i <iterations>` - extracts the subgraph from the GFA, given a segment name, or multiple (if multiple, these must be comma separated without space). Number of iterations may need to be increased for large graphs.
- `gfatk extract-chloro <GFA>` - extracts the plastid from the GFA. It has default parameters which seem to work okay.
- `gfatk extract-mito <GFA>` - extracts the mitochondria from the GFA. It has default parameters which seem to work okay.
- `gfatk fasta <GFA>` - extracts a fasta file from the GFA. This simply prints each of the segments from the GFA. I say it's almost as simple as the `awk` version, but the toolkit does some checks to see if we are actually dealing with a GFA or not.
- `gfatk linear <GFA> -e -i -n <node-threshold>` - forces the longest linear legal representation of the graph. You can evaluate within subgraphs (`-e`), or include node coverage information (`-i`).
- `gfatk overlap <GFA> -s <size>` - extracts the overlaps from the GFA. These are taken from the CIGAR string from each of the links, and optionally extended (e.g. `-s 1000` to 1000bp either side of the overlap).
- `gfatk path <GFA> <path> (-p path/to/path.txt)` - evaluates a linear representation of the graph, given an input path. The input path can be on the command line, or a file. Simply, it must be an comma separated list of node ID's and orientations (1+,2-,3+ ... ).
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

For full functionality of the toolkit, two tags are required, node coverage and edge coverage. Other functionality will fail if the CIGAR string is not purely an overlap; i.e. in the format `<integer>M`. Only GFA version 1 supported. Only header (`H`), segment (`S`), and link (`L`) lines are required. Other lines (e.g. path, `P`), will I think be ignored.

```
H	VN:Z:1.0
S	11	ACCTT	ll:f:30.0 <- this tag indicates node/segment coverage (here it's 30.0)
S	12	TCAAGG	ll:f:60.0
S	13	CTTGATT	ll:f:30.0
L	11	+	12	-	4M	ec:i:1 <- this tag indicates edge coverage (here it's 1)
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
