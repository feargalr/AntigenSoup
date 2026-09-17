# Building an AntigenSoup epitope database from IEDB

`scripts/build_iedb_db.R` turns the public [IEDB](https://www.iedb.org) database
export into an AntigenSoup epitope FASTA. It replaces the manual process that
produced the bundled `iedb.fasta.gz`, which carried no record of how it was
filtered.

```bash
Rscript scripts/build_iedb_db.R --outdir databases/iedb
```

That is the whole thing. It downloads the export, filters it, and writes the
database. Requires R with `data.table`, plus the `curl` and `unzip` commands.
A build takes about 5 minutes and needs ~6 GB of memory.

## Outputs

For a release stamped `2026-09-08`:

| File | Contents |
|---|---|
| `iedb_antigensoup_2026-09-08.fasta` | The database. One record per unique peptide. |
| `iedb_antigensoup_2026-09-08.tsv` | 20 metadata columns, one row per FASTA record, same order. |
| `iedb_antigensoup_2026-09-08.report.txt` | Provenance, settings, per-stage retention, output distributions. |

Source archives are cached in `<outdir>/source` (~530 MB) and reused on the next
build. Use the same `--outdir` to avoid re-downloading.

The FASTA is deliberately minimal, because `ac_match.py` echoes the header into
its `epitope_ids` column and long headers make the hit table unreadable:

```
>AS_0000001
AAAACTTMK
```

Join a hit table back to the metadata on that identifier.

## Options

```
--outdir DIR              Output directory (default: databases/iedb)
--min-length N            Minimum peptide length (default: 8)
--max-length N            Maximum peptide length (default: 25)
--include-mhc-only        Also include peptides whose only positive evidence is
                          MHC binding or elution
--include-low-complexity  Skip the sequence complexity filter
--natural-only            Also drop neo-epitopes
--force-download          Re-download even if cached
--threads N               data.table threads
-h, --help
```

```bash
# Recommended default database
Rscript scripts/build_iedb_db.R --outdir databases/iedb

# Broader: include eluted MHC ligands and longer linear B-cell epitopes
Rscript scripts/build_iedb_db.R --outdir databases/iedb_broad \
    --include-mhc-only --max-length 50
```

## Which IEDB files are used

Four archives from `https://www.iedb.org/downloader.php`: `epitope_full_v3.zip`
as the master epitope table, and `tcell_full_v3.zip`,
`bcell_full_v3_multi_file.zip` and `mhc_ligand_full_multi_file.zip` for assay
evidence.

Two practical notes, both learned the hard way and both handled in the code.

The CSVs carry **two header rows**, a group row and a name row, so column
identity is the pair, for example `Assay :: Qualitative Measurement`. The builder
re-reads those headers on every run and aborts with an explicit message if a name
has moved, rather than silently producing a wrong database.

`data.table` **cannot read the single-file B-cell and MHC ligand exports**. They
are 3.2 GB and 9.2 GB uncompressed, and `fread` fails with a 2 GB string limit
and then a bus error. The multi-file variants split them into sub-2 GB members,
which read cleanly. Only member `00` carries the header rows, and `fread`
mis-detects the column count on the headerless members, so the expected width is
passed explicitly and every read is checked for malformed identifiers.

The export has no version number, so the release stamp is the latest build date
across the archive members.

## Sequence handling

Sequences come from `Epitope :: Name`. IEDB writes a chemically modified epitope
as `SEQUENCE + MOD(position)`, for example `AAIAAVKEEAF + METH(A10)`. The builder
keeps the residue string before the ` + `, which is the parent amino-acid
sequence a microbial protein could actually match, and flags the record as
modified. Residues are never rewritten. A sequence is dropped only if *every*
record for it is modified.

Sequences are upper-cased and stripped of whitespace. Only the 20 standard amino
acids are accepted; anything else is rejected rather than repaired.

**On U and O.** Selenocysteine and pyrrolysine do not occur in IEDB linear
peptide sequences. Apparent hits are inside modification names, not sequences.
The only non-standard residue that does occur is X, in 550 records, which are
rejected because an unknown residue cannot participate in exact matching.

### Length

Default 8 to 25 amino acids.

The minimum matters most for AntigenSoup. Exact matches of very short peptides
occur by chance in large metagenomic protein spaces: expected hits scale roughly
as `N / 20^L` for a search space of `N` residues, so a 6-mer is expected many
times over in a single assembly while an 8-mer is not. Eight is also the floor of
the MHC class I peptide range.

The maximum keeps the database to peptides that plausibly represent a single
presented or antibody-bound determinant. Raising it admits tiled peptide-array
constructs: the current release contains 228,233 47-mers from one *Trypanosoma
cruzi* array study alone. Extend to 50 with `--max-length 50` for broader linear
B-cell epitope work.

The range is applied as one continuous test to a single master table. The old
bundled FASTA contained no 12-mers at all, the signature of two separate queries
covering 8-11 and 13-25 being concatenated. That cannot happen here.

## Evidence

A peptide being present in IEDB says nothing on its own about whether it was ever
shown to do anything. Each assay row carries a curated `Qualitative Measurement`
(`Qualitative Measure` in the B-cell table, a naming inconsistency the builder
handles). Values are `Positive`, `Positive-Low`, `Positive-Intermediate`,
`Positive-High` and `Negative`. Anything beginning `Positive` counts as positive.

| Category | Meaning | In the FASTA |
|---|---|---|
| `IMMUNE_SUPPORTED` | At least one positive T-cell or B-cell assay. | Yes |
| `MHC_ONLY` | At least one positive MHC binding or elution assay, no positive T-cell or B-cell assay. | Only with `--include-mhc-only` |
| `NEGATIVE_ONLY` | Has experimental assays, none positive. | No |
| `NO_ASSAY` | No assay row in any table. | No |

**There is no `PREDICTION_ONLY` category, because the export cannot populate
one.** Every row in the three assay tables is a curated experimental
measurement. The only prediction-related field is `MHC Restriction :: Evidence
Code`, whose values include `MHC binding prediction`, but that describes how the
restricting allele was assigned for an otherwise experimental assay. It is not
evidence about the epitope.

The split between `IMMUNE_SUPPORTED` and `MHC_ONLY` is the most consequential
default. Mass-spectrometry immunopeptidomics has made eluted ligands by far the
largest class in IEDB, and an eluted ligand shows that a peptide was presented,
not that any lymphocyte responded to it. Including them roughly doubles the
database without adding demonstrated immunogenicity.

`evidence_strength` is reported separately from the inclusion decision, so a
single well-supported study is still eligible: `HIGH` means immune-supported with
at least 2 independent IEDB references, `SUPPORTED` means one. The underlying
counts are all columns in the TSV, so you can impose a stricter rule by filtering
it.

## Natural, synthetic and neo-epitopes

IEDB curates this in `Related Object :: Epitope Relation`. An empty value means
the epitope has no related object, that is, it is recorded as a straight fragment
of a natural molecule.

| Field value | `epitope_nature` | Default |
|---|---|---|
| *(empty)* | `NATURAL` | Included |
| `in-frame`, `frameshift`, `fusion`, `unspecified neo-epitope` | `NEOEPITOPE` | Included |
| `analog`, `mimotope` | `ANALOG_OR_MIMOTOPE` | Excluded |
| anything unrecognised | `UNDETERMINED` | Included, with a warning |

Neo-epitopes are kept and labelled separately. They are not synthetic constructs;
they are real sequences arising from somatic variation, and microbial proteins
matching a tumour neo-epitope is a biologically interesting result.
`--natural-only` drops them if you want the strictest set.

Anything the export does not classify is retained and flagged, never discarded. A
peptide attested anywhere as a natural fragment counts as natural even if another
study used the same string as an analogue backbone.

## Sequence complexity

Exact matching makes low-complexity peptides actively harmful. A poly-alanine or
poly-glutamate run will match many unrelated proteins, and those matches carry no
information about shared antigenicity.

Four fixed thresholds, applied together:

| Rule | Threshold | Candidates removed |
|---|---|---|
| distinct residues | >= 4 | 0.33% |
| fraction from the commonest residue | <= 0.50 | 0.67% |
| longest single-residue run | <= 4 | 0.31% |
| normalised Shannon entropy | >= 0.55 | 1.36% |
| **any of the four** | | **1.52%** |

Entropy is normalised by `log2(min(length, 20))`, the maximum attainable at that
length rather than the absolute maximum, which is what makes an 8-mer and a
25-mer comparable. In the current release the fraction flagged stays between 0.6%
and 5% across every length from 8 to 25, so the filter is not acting as a
disguised length filter.

These values were checked against the real distribution before being fixed, not
assumed. The entropy threshold sits at roughly the 1.5th percentile of the
candidate distribution, a real tail rather than a haircut; tightening it to 0.70
would remove 8.6%, which starts discarding ordinary peptides. They are constants
in the script rather than options because loosening them individually is rarely
what anyone wants; `--include-low-complexity` turns the filter off entirely.

**On SEG.** A peptide of 8 to 25 residues is shorter than SEG's default
12-residue window is designed to slide within, and SEG's complexity measure is
itself the entropy of the residue composition of a window. Computing that entropy
over the whole peptide is the same measure applied at the right scale, and needs
no dependency.

Sequences removed are the intended ones: `AVRRRRRRRV`, the His-tag-like
`PHHHHHHHRHPQPAT`, the gliadin repeat `PQPQQPQQSFPQQQ`, and
`GQGPGAPQGPGAPQGPGAPQGPGAP`, the Epstein-Barr EBNA-1 Gly-Ala repeat. That last
one is worth pausing on: some low-complexity sequences are genuinely
immunologically important. They are excluded because they behave badly under
exact matching, not because they are uninteresting.

## Source organism

Source organism is **not** used as a filter. Bacterial, viral, fungal, allergen,
human, tumour and autoimmune epitopes are all retained, because microbial
proteins matching human or tumour epitopes are exactly the molecular mimicry
signal AntigenSoup is built to find. Organism information is aggregated into
`source_organisms` so organism-specific subsets can be cut from the TSV later.

## Deduplication and identifiers

Identical sequences are collapsed to one FASTA record, **after** all evidence is
aggregated, so a peptide appearing in many records keeps the union of its IEDB
identifiers, assay counts, references, organisms, antigens, hosts and MHC
restriction. List-valued columns hold at most 5 entries, ordered by frequency
then alphabetically, with a `...(+N more)` suffix.

Identifiers are `AS_0000001` upward over the retained sequences sorted
lexicographically. Given the same release and the same options this is exactly
reproducible, verified by building twice with different thread counts and
comparing checksums.

They are **not stable across IEDB releases**: adding one peptide near the start of
the alphabet shifts every identifier after it. The peptide sequence is the stable
key, and it is present in both the FASTA and the metadata, so join on `sequence`
when comparing results across two builds.

## Metadata columns

| Column | Meaning |
|---|---|
| `antigensoup_id` | FASTA identifier |
| `sequence`, `length` | Peptide and its length |
| `iedb_epitope_ids`, `n_iedb_records` | Contributing IEDB epitope IDs and how many there were |
| `n_positive_tcell`, `n_negative_tcell` | T-cell assay counts |
| `n_positive_bcell`, `n_negative_bcell` | B-cell assay counts |
| `n_positive_mhc` | Positive MHC binding/elution assays |
| `n_references` | Distinct IEDB references reporting on this peptide |
| `source_organisms`, `source_antigens` | Curated source species and molecules |
| `hosts` | Assay host organisms |
| `mhc_classes`, `mhc_alleles` | MHC restriction |
| `evidence_category`, `evidence_strength` | See above |
| `epitope_nature` | `NATURAL`, `NEOEPITOPE`, `ANALOG_OR_MIMOTOPE` or `UNDETERMINED` |
| `complexity_status` | `PASS` or `LOW_COMPLEXITY` |

## Reproducibility

The report records the IEDB release stamp, each source archive with its size and
md5, the script version, the R and `data.table` versions, the build timestamp and
every setting. To reproduce a build, use the archives its report names and the
same options.
