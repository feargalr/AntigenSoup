# Changelog

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
