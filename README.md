# gfatk

Exploratory functions to manipulate Graphical Fragment Assembly Format (GFA), and Graph Alignment Format (GAF). Not very useful yet.

## Usage

### gfatk linear

Turn a GFA into a linear sequence by traversing the graph, using each segment only once. This is an NP hard problem (Hamiltonian Path). The code currently is very much prototype and untested.

```
gfatk-linear 
Force a linear representation of the graph.
Each node is included once only.

USAGE:
    gfatk linear --gfa <gfa>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -g, --gfa <gfa>    Input GFA file.
```

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

I realise now this is almost identical to `gfatools view -l [node] -r [some number] gfa.gfa`.

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

### gfatk gaf

Merge forward and reverse mappings from a .gaf.

```
gfatk-gaf 
gaf

USAGE:
    gfatk gaf --gaf <gaf>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -g, --gaf <gaf>    Input GAF file.
```

### TODO's

In general tidy up code into struct implementations. E.g.
- Loading gfa -> petgraph.
- Methods on the petgraph.