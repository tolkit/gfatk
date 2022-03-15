# gfatk

<p align="center">
    <img width="300" height="132" src="https://www.darwintreeoflife.org/wp-content/themes/dtol/dist/assets/gfx/dtol-logo-round.png">
</p>

A command line utility to explore, extract, and linearise plant mitochondrial assemblies. The Graphical Fragment Assembly files (GFAs) used to refine the code in this repository are almost exclusively generated from the assembly program <a href="https://github.com/maickrau/MBG">MBG</a>. See the testing section below for caveats.

## Install

Grab from the releases (Mac & Linux only):

```bash
# for mac
curl -L "https://github.com/tolkit/gfatk/releases/download/0.1.4/gfatk_mac_0.1.4" > goat && chmod +x goat
# and linux (ubuntu)
curl -L "https://github.com/tolkit/gfatk/releases/download/0.1.4/gfatk_ubuntu_0.1.4" > goat && chmod +x goat
```

Or build from source.

```bash
# e.g. get rustup!
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
# clone!
git clone https://github.com/tolkit/gfatk
# cd!
cd gfatk
# build!
cargo build --release
# or install into your path!
cargo install --path .
```

## Pipeline

The pipeline is still being ironed out, but it progresses like this:

- Assemble plant mitome/plastome using MBG
- Extract mitochondrial subgraph from the MBG output GFA
- Linearise the GFA by traversing the graph to find the longest path (simple or not)

In terms of `gfatk` commands, it looks like this:

```bash
# install rust and `cargo install --path .` etc.
# MBG output = output.gfa

# optionally check what your GFA structure is
gfatk stats output.gfa
# or check what the graph looks like visually
# using `dot` from `graphviz`
gfatk dot output.gfa | dot -Tsvg > output.svg

# extract the mitochondrial subgraph
# this is done using coverage, GC content info
# and expected minimum genome size (100,000bp)
# then check stats
gfatk extract-mito output.gfa | gfatk stats

# if this isn't the right subgraph, manually intervene
# select one segment name from the subgraph of interest
gfatk extract output.gfa -s 2 | gfatk stats

# now we linearise
gfatk extract-mito output.gfa | gfatk linear > linearised.fasta
# alternatively pass the `-i` flag which will include
# coverage information from the GFA itself
gfatk extract-mito output.gfa | gfatk linear -i > linearised.fasta
# if the graph is even slightly complex, this last command
# will fail with a stack overflow.
```

## Examples and docs

A couple of more detailed examples can be seen in the `examples` directory, where there is a `README.md` file. To view the auto-generated documentation of the binary itself, including details of all underlying functions, see:

https://tolkit.github.io/gfatk/

## Testing

Some unit tests are now provided in the `tests` directory. To run these (you'll need Rust):

```bash
cargo test --release
```

For full functionality of the toolkit, two tags are required, node coverage and edge coverage. Other functionality will fail if the CIGAR string is not purely an overlap; i.e. in the format `<integer>M`. Only GFA version 1 supported.

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

More tests to come.

## Thanks

Many thanks to the developers of MBG, and partners in the Tree of Life program, and beyond:
- Marcela Uliano-Silva
- Sergey Nurk
- Alex Twyford
- Lucia Campos
- Chenxi Zhou