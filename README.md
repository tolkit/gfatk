# gfatk

A GFA toolkit. 

## Usage

Only functionality currently is to extract overlapping regions in a GFA between two segments, and extend the sequence either side of the overlap.

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

## TODO's

- Extract a connected subgraph given a node?