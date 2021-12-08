# gfatk

Exploratory functions to manipulate Graphical Fragment Assembly Format (GFA), and Graph Alignment Format (GAF) files.

All code should be considered a prototype, with minimal testing.

## Usage

### gfatk linear

Turn a GFA into a linear sequence by traversing the graph, using each segment only once.

```
gfatk-linear 
Force a linear representation of the graph.
Each node is included once only.

USAGE:
    gfatk linear [OPTIONS] --fasta-header <fasta-header> --gfa <gfa>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --coverage-file <coverage-file>    Name of the text file indicating the oriented coverage of links in a GFA.
    -f, --fasta-header <fasta-header>      Name of the fasta header in the output file. [default: gfatk-linear]
    -g, --gfa <gfa>                        Input GFA file.
```

Can take the output from `gfatk gaf`, but using different starting points in the graph gives different coverages. Not currently sure how to solve this, apart from iterating through every longest path.

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

Count the coverage of each junction/overlap from a GAF file.

Outputs a TSV of:

```
from_orient	from	to_orient	to	coverage
+       232     -       228     122
-       229     -       231     97
+       231     +       228     148
-       227     -       228     184
+       232     -       229     116
+       229     +       227     189
-       233     +       231     83
-       232     +       233     106
+       230     -       227     583
-       230     +       231     108
-       230     +       232     178
+       228     -       232     133
-       228     -       231     115
+       230     -       2       3
-       228     -       233     168
-       233     +       232     165
-       231     +       230     112
+       229     -       232     145
+       231     +       229     104
+       233     +       229     143
+       2       -       230     21
-       232     +       230     132
-       227     -       229     181
+       228     +       227     250
-       229     -       233     160
+       233     +       228     154
-       231     +       233     121
+       227     -       230     535
```

Without the headers though. The idea is to passs these as the edge weights in the GFA graph.

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