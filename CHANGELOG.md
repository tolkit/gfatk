# Changelog

All notable changes to gfatk are listed here, newest first. Versions follow
[semantic versioning](https://semver.org/); release binaries are on the
[GitHub releases page](https://github.com/tolkit/gfatk/releases).

## [0.7.0] - 2026-11-27

### Changed

- `gfatk resolve` may walk a junction more than once. The solver chose each
  transition between two segments at most once, so a repeat spanning several
  segments could be walked only once, however many copies its depth showed,
  and everything a circuit could reach only by walking it again was left
  out. A transition may now be walked as many times as the smaller budget of
  its two segments allows. On fern mitochondrial graphs, this was most of
  what resolve left out.
- `gfatk resolve` keeps as much of the sequence read depth shows as it can,
  without overstating it, rather than including as many segments as it can.
  Each segment counted once whatever its length or depth, so the solver
  would give up a long segment to keep two short ones, or a single-copy
  segment to keep one at a few percent of the baseline depth. A segment's
  weight is now its length (times 20 for a gene-bearing segment, as before).
  A segment at or above the single-copy baseline depth earns its weight when
  included. One below it, at copy ratio r, earns r of its weight but
  overstates its depth by 1 - r, so it nets 2r - 1: sequence at 5% of the
  baseline (typically from the other organelle, nuclear copies or an
  assembly error) is walked only where that connects more sequence than it
  overstates, and sequence at half the baseline breaks even. Segments left
  out are still written as unplaced records.
- `gfatk resolve` forces its single circuit through the circuit, found
  without that constraint, that scores best in those terms. It used to take
  the one with the most segment visits, which can be one walking a few short
  repeats many times rather than the one holding the genome.
- Read support now only decides between solutions to within 10 bp of
  sequence. If the solver fails on a graph with it, resolve solves again
  without it and says so.
- Together, on 54 plant organelle graphs (31 fern and lycophyte
  mitochondria), against 0.6.3:
  - 43 graphs now place at least 70% (mitochondria) or 90% (plastids) of
    their graph sequence in the circuit, up from 25.
  - The median share placed rises from 0.70 to 0.84.
  - How far each segment's walks are from its copy ratio, weighted by
    length, falls from a median of 0.43 to 0.24 of the depth-implied genome
    size.
  - Sequence placed at under 0.2x the baseline depth falls from 188 to 40
    kb. Because of that, Geranium phaeum's plastid now places 88% of its
    graph rather than 95%: it no longer walks 46 kb of such sequence.
  - All 54 resolve to a single circuit, output is identical between runs,
    and the 54 together take 686 s rather than 1201 s.
  - 94 of 810 repeat runs, in 32 graphs, are flagged as not crossed by
    reads.

### Added

- `gfatk resolve` flags repeat runs whose arrangement the reads don't show. A
  repeat run is a stretch of consecutive repeat segments (budget above 1)
  between two single-copy segments. Which flank each copy joins is shown by
  reads that cross the whole run, flank to flank; where fewer than 5 do, the
  circuit's arrangement there is one the graph allows, not one the reads
  show. Each circuit's header now carries `repeat_runs=`,
  `unspanned_repeat_runs=` and `unspanned_at=` (1-based path positions of the
  repeat segments of each such run, or `none`), before `path=`. stderr lists
  each such run with its segments, length, flanks and spanning reads.

### Fixed

- When the single-circuit solve fails, the warning now gives the solver's own
  error instead of blaming the time limit.

## [0.6.3] - 2026-11-19

### Fixed

- `gfatk resolve` wrote each circuit as its segments joined end to end, so
  the bases two linked segments share (the link's overlap, `41M` say) came
  out twice at every junction. Each segment after the first is now written
  without its overlap with the one before, as `gfatk path` and `gfatk linear`
  already do, and the closing overlap, from the last segment back to the
  first, is trimmed off the start of the first segment. That is where oatk
  puts the origin of a circular contig, so a `path=` gives the same
  coordinates in either tool's output. `bp=` in the header and the length on
  stderr are those of the trimmed sequence. The path chosen is unchanged.
- The main-path sequence a bubble arm is compared with is trimmed the same
  way.
- `gfatk resolve` gave a different answer each time it was run on the same
  input: the order of candidate edges, of the solver's variables and of the
  walk through the chosen edges all followed hash-map order. Everything is
  now in sorted order, and each circuit is written in one fixed form (of all
  its rotations on either strand, the one that sorts first, so it starts at
  its first segment by name). The same GFA and GAF always give the same
  FASTA, byte for byte.
- Where a circuit passes through a segment more than once, which way it
  leaves after each entry was arbitrary. It is now chosen by the reads: at
  each state, entries and exits are paired to maximise the support of the
  `(previous, segment, next)` triples walked, and closed walks sharing a
  state are joined where that costs the least support.
- A circuit of a single segment (a circular unitig linked to itself) was
  dropped, and the segment written as `unresolved_segment` with its overlap
  untrimmed. It is now a circuit like any other.
- The single-circuit constraint was rooted at the best-connected state of
  the whole graph, which forced that state into the circuit whether or not
  the best circuit passes through it, and could cost segments (or leave no
  solution, and several circuits). The root is now the best-connected state
  of the largest circuit found without the constraint; if that is already a
  single circuit, the constrained solve is skipped.

- `gfatk linear`, `gfatk path`, `gfatk trim` and `gfatk extract` only
  followed a link in the direction it is written: `L a + b +` let a walk go
  from `a+` to `b+` but not from `b-` to `a-`, which is the same junction on
  the other strand. On a graph that writes each link once (oatk's, for one)
  `linear` therefore stopped wherever the next link was written the other
  way round, `path` refused such a step (including any `path=` that
  `resolve` or oatk reports), `trim` counted only the segments a segment
  links *to* as its neighbours, and `extract` only searched outwards along
  links as written. All four now treat a link as usable from either strand.
  On 300 oatk graphs, the longest `linear` walk gets longer in 157 and
  shorter in 1.
- `gfatk resolve` counted the reads supporting a junction separately for
  each strand, and needed 5 on one of them; reads from both strands are now
  pooled before the threshold, and both traversals of a junction carry the
  same support.
- `gfatk resolve` could report a segment as a bubble arm from an entry onto
  one of its strands and an exit from the other, which is not a walk
  through it. Entry and exit must now be on the same strand; such segments
  are reported as unresolved. The arm is written, and compared with the
  main path, on the strand it is walked on (before, always forward, so the
  identity reported for a reverse-strand arm was meaningless), and its
  header gains `:strand=+` or `:strand=-` after `support_reads`.

- `gfatk resolve` measured read depth, and from it each segment's copy
  number, wrongly in three ways, now fixed together. Segments were placed
  along a read's path end to end, though GAF path coordinates count the
  bases two linked segments share once, so bases were credited to the
  wrong segment. Reads lying within a single segment were ignored, which
  left long segments (few junction-crossing reads for their length) at
  0.4x the depth of the rest. And single-copy depth was the median over
  segments, which the many short ones set; it is now the median weighted
  by segment length. Only reads aligned over at least 90% of their length
  count for depth (`--min-read-aligned`): on mitochondrial graphs two
  recruited reads in three are only partly organellar, and counting their
  aligned pieces puts depth out again. Against the graphs' k-mer coverage
  on 54 oatk graphs, the copy number agrees for 90% of segments (77%
  before), and a single-copy segment is called two-copy in 3% of cases
  (17% before), a two-copy one single-copy in 8% (16%).

### Added

- `gfatk resolve --min-read-aligned` (default 0.90; 0 counts every
  alignment), and a line on stderr saying how many reads count for depth.
- `gfatk resolve` reports on stderr each junction of a circuit where the two
  segments differ within their link's overlap (an overlap of `41M` is 41
  aligned bases, not necessarily 41 identical ones). The earlier segment's
  bases are the ones written.
- `gfatk resolve` says so on stderr when the solver stopped at `--time-limit`
  without proving its answer the best one. Before, such an answer was used
  without comment.

### Changed

- `gfatk resolve` now needs every link's overlap to be `<integer>M` or `*`
  (read as no overlap), and stops with an error naming the link otherwise.
- Tests: `linear`, `trim` and `stats` on a GFA without edge coverage are
  expected to succeed, as they have since edge coverage became optional.

## [0.6.2] - 2026-10-15

### Fixed

- `gfatk resolve` no longer produces circuits containing transitions that no
  GFA link allows. When registering the reverse-strand reading of a read
  triple `a seg b`, the orientation of `a` was not flipped, so every valid
  triple also produced an impossible tier-0 edge (a real `<u14<u72` read
  yielded `u72+ -> u14-`). On vpJunScop1 (*Juniperus scopulorum*) this put 31
  impossible transitions into a 130-step circuit; with the fix there are none,
  and the dead-end segments involved (u14, u72, u96) are reported as unplaced
  instead of being bridged by an invented junction.
- The bare-topology evidence tier now registers both traversals of each GFA
  link (as written and its reverse complement), not only the literal one.

### Added

- Read-derived transitions (tiers 0 and 1) are checked against the GFA's
  links before they reach the solver. Any that are not allowed are discarded,
  and their number is reported on stderr, as a guard against a GAF aligned to
  a different graph.

## [0.6.1] - 2026-09-24

### Changed

- Each resolved circuit's FASTA header now includes its segment path
  (`path=u13-,u6-,...`), matching oatk's convention.

## [0.6.0] - 2026-09-23

### Added

- `gfatk resolve --bed` reads oatk-style gene annotation directly, with no
  GFF conversion step. Every BED row counts as a gene call, whatever its
  score.
- Segments with no evidence connecting them to any resolved circuit are
  written to the FASTA as `unresolved_segment` records instead of only being
  counted.

## [0.5.3] - 2026-09-21

### Changed

- CITATION.cff lists a single author; the README uses a DOI-based Zenodo
  badge.

## [0.5.2] - 2026-09-21

### Fixed

- Corrected the DOI and affiliations in CITATION.cff and pointed the README
  install URLs at the current release.

## [0.5.1] - 2026-09-21

### Added

- CITATION.cff with a Zenodo DOI, and a DOI badge in the README.

## [0.5.0] - 2026-09-18

### Added

- `gfatk resolve` compares each bubble arm with the sequence the main circuit
  uses between the same attachment points and reports the edit distance and
  percent identity on stderr and in the FASTA header. The search covers the
  whole circuit, and comparisons against a detour more than 5× the arm's size
  are reported as not meaningful rather than diffed.

## [0.4.0] - 2026-09-17

### Fixed

- `gfatk resolve`'s single-circuit constraint now tracks flow per
  (segment, orientation) state, not per segment. Previously a repeat segment
  used in both orientations could make two disconnected loops look connected.

### Added

- Evidence-based circuit splicing: leftover circuits are grafted into the
  main one where real read evidence supports a join, and stderr explains why
  any that cannot be joined were left apart.
- `--time-limit` for the solver (previously fixed at 45 s).
- `--no-bubble-fasta` to keep bubble-arm sequences out of the FASTA output
  (they are still reported on stderr).

### Changed

- GFA node/segment lookups are backed by hash maps, so loading very large
  graphs no longer takes hours (#16).

## [0.3.0] - 2026-09-16

### Added

- `gfatk resolve`: resolves a repeat-rich assembly graph into a circular
  genome from HiFi read paths (GAF), with copy number taken from read depth
  and edge selection posed as an integer program.

### Fixed

- `gfatk linear` no longer overflows the stack on complex, repeat-rich
  graphs: exhaustive path enumeration was replaced with a coverage-budgeted
  greedy walk, and the `--node-threshold` default was raised.
- Release CI builds again (current Rust toolchain, macOS deployment target).

## Earlier releases

Versions 0.1.x–0.2.x (2022–2025) predate this changelog; see the git history.

[0.7.0]: https://github.com/tolkit/gfatk/compare/0.6.3...0.7.0
[0.6.3]: https://github.com/tolkit/gfatk/compare/0.6.2...0.6.3
[0.6.2]: https://github.com/tolkit/gfatk/compare/0.6.1...0.6.2
[0.6.1]: https://github.com/tolkit/gfatk/compare/0.6.0...0.6.1
[0.6.0]: https://github.com/tolkit/gfatk/compare/0.5.3...0.6.0
[0.5.3]: https://github.com/tolkit/gfatk/compare/0.5.2...0.5.3
[0.5.2]: https://github.com/tolkit/gfatk/compare/0.5.1...0.5.2
[0.5.1]: https://github.com/tolkit/gfatk/compare/0.5.0...0.5.1
[0.5.0]: https://github.com/tolkit/gfatk/compare/0.4.0...0.5.0
[0.4.0]: https://github.com/tolkit/gfatk/compare/0.3.0...0.4.0
[0.3.0]: https://github.com/tolkit/gfatk/compare/0.2.16...0.3.0
