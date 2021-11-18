# gfatk

A GFA toolkit. 

## Usage

### gfatk overlap

Extract overlapping regions in a GFA between two segments, and extend the sequence either side of the overlap.

```
gfatk-overlap 
Extract overlaps from a GFA.

USAGE:
    gfatk overlap [OPTIONS] --gfa <gfa>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -g, --gfa <gfa>      Input GFA file.
    -s, --size <size>    Region around overlap to extract. [default: 1000]
```

E.g. 

`gfatk overlap -g test.gfa > test.fa`

### gfatk extract

Extract a connected subgraph given a node. Supply a gfa, a sequence ID connected to the subgraph you want to extract, and optionally a number of iterations to search for neighbours for.

```
gfatk-extract 
Extract subgraph from a GFA, given a segment name.

USAGE:
    gfatk extract [OPTIONS] --gfa <gfa> --sequence-id <sequence-id>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -g, --gfa <gfa>                    Input GFA file.
    -i, --iterations <iterations>      Number of iterations to recursively search for connecting nodes. [default: 3]
    -s, --sequence-id <sequence-id>    Extract subgraph of which this sequence is part of.
```