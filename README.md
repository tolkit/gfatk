# gfatk

<p align="center">
    <img width="300" height="132" src="https://www.darwintreeoflife.org/wp-content/themes/dtol/dist/assets/gfx/dtol-logo-round.png">
</p>

A command line utility to explore, extract, and linearise plant mitochondrial assemblies. The Graphical Fragment Assembly files (GFAs) used to refine the code in this repository are almost exclusively generated from the assembly program <a href="https://github.com/maickrau/MBG">MBG</a>. Little testing has been done beyond the output of MBG.

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
gfatk stats -g output.gfa

# extract the mitochondrial subgraph
# this is done using coverage, GC content info
# and expected minimum genome size (100,000bp)
gfatk extract-mito -g output.gfa > mito.gfa

# now we linearise
gfatk linear -g mito.gfa > linearised.fasta
# alternatively pass the `-i` flag which will include
# coverage information from the GFA itself
gfatk linear -ig mito.gfa > linearised.fasta
# if the graph is even slightly complex, this last command
# will fail with a stack overflow.
```

## Caveats

Still in active development, so I expect there will be bugs & API changes.

## Thanks

Many thanks to the developers of MBG and GraphAligner, and partners in the Tree of Life program and beyond:
- Marcela Uliano-Silva
- Sergey Nurk
- Alex Twyford
- Lucia Campos
- Chenxi Zhou