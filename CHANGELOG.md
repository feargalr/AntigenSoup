# Changelog

The version lives in the `VERSION` file at the repository root. `antigensoup` and
both R scripts read it at runtime. `CITATION.cff` is a static file and carries its
own copy, so a release bumps those two files and nothing else.

## v0.5.1
- Aligned the version reported by `scripts/build_iedb_db.R` and
  `scripts/generate_spikein_benchmark.R` with the release version. Both declared
  1.0.0 independently, which appeared in `--help` and, more importantly, in the
  `script_version` field of the build report and benchmark manifest. A provenance
  record should not name a version that does not exist.

## v0.5.0
- Removed the bundled `iedb.fasta.gz`. The database is now built from the current
  IEDB release with `scripts/build_iedb_db.R`, which records the release it used
  and every filter applied. The bundled file had no such record and, as a result
  of being assembled from two separate length queries, contained no 12mers at all
  despite IEDB holding 267,718 of them. It remains retrievable from the git history
  (`git show 396d94b:iedb.fasta.gz`) for anyone needing to reproduce older results.
- Install instructions now build the database instead of unzipping the bundled one.

## v0.4.0
- Renamed the pipeline executable from `daedalus` to `antigensoup`, along with its
  conda environment, install paths and user-facing messages. This also fixes the
  install instructions, which already documented the new names while the wrapper
  still looked for the old ones.
- Fixed `CITATION.cff`, which still named the project Daedalus and pointed at the
  old repository.
- Added `scripts/build_iedb_db.R`, which builds the epitope database directly from
  the current IEDB export: linear 8-25mer peptides with positive experimental
  T-cell or B-cell evidence, filtered for amino-acid alphabet and sequence
  complexity, deduplicated, with a metadata TSV and a build report. Replaces the
  manually assembled `iedb.fasta.gz`, which had no record of its own filtering and
  contained no 12mers at all. Documented in `docs/iedb_database_build.md`.
- Added `scripts/generate_spikein_benchmark.R`, which builds artificial epitope
  spike-in benchmarks for sensitivity and specificity testing. Artificial peptides
  are composition-matched shuffles of retained IEDB epitopes, proven absent from
  both the database and the target proteins before insertion at recorded positions.
  Documented in `docs/spikein_benchmark.md`, including the recorded validation run
  against the E. coli K-12 reference proteome: 1000 spikes, 13/13 internal checks,
  100% recovery by `ac_match.py` at the correct protein and coordinate, and zero
  hits against the unspiked proteome.

## v0.3.3
- Added gen variants script, scaffolds as input. 

## v0.3.2
- Added skip steps if output files already present 

## v0.3.1
- Updated handling of additional acmatch script and conda env

## v0.3.0
- Changed to AC Match for string matching epitopes

## v0.2.0
- Expanded functionality, added filtering of human sequences

## v0.1.0
- Initial build
