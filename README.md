# gfatk

<p align="center">
    <img width="300" height="132" src="https://www.darwintreeoflife.org/wp-content/themes/dtol/dist/assets/gfx/dtol-logo-round.png">
</p>

A command line utility to explore, extract, and linearise plant mitochondrial assemblies. The Graphical Fragment Assembly files (GFAs) used to refine the code in this repository are almost exclusively generated from the assembly program <a href="https://github.com/maickrau/MBG">MBG</a>. Little testing has been done beyond the output of MBG.

## Install

For now, I haven't put any precompiled binaries out as the code has been changing so much. Nearly there but in the meantime, install rust and build yourself.

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

## Thanks

Many thanks to the developers of MBG, and partners in the Tree of Life program and beyond:
- Marcela Uliano-Silva
- Sergey Nurk
- Alex Twyford
- Lucia Campos
- Chenxi Zhou