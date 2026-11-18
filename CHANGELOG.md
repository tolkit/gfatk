# Changelog

All notable changes to gfatk are listed here, newest first. Versions follow
[semantic versioning](https://semver.org/); release binaries are on the
[GitHub releases page](https://github.com/tolkit/gfatk/releases).

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

[0.6.2]: https://github.com/tolkit/gfatk/compare/0.6.1...0.6.2
[0.6.1]: https://github.com/tolkit/gfatk/compare/0.6.0...0.6.1
[0.6.0]: https://github.com/tolkit/gfatk/compare/0.5.3...0.6.0
[0.5.3]: https://github.com/tolkit/gfatk/compare/0.5.2...0.5.3
[0.5.2]: https://github.com/tolkit/gfatk/compare/0.5.1...0.5.2
[0.5.1]: https://github.com/tolkit/gfatk/compare/0.5.0...0.5.1
[0.5.0]: https://github.com/tolkit/gfatk/compare/0.4.0...0.5.0
[0.4.0]: https://github.com/tolkit/gfatk/compare/0.3.0...0.4.0
[0.3.0]: https://github.com/tolkit/gfatk/compare/0.2.16...0.3.0
