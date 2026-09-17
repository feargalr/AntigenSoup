# Artificial epitope spike-in benchmark

`scripts/generate_spikein_benchmark.R` builds a controlled benchmark for
measuring AntigenSoup's sensitivity and specificity. It inserts artificial
peptides at known positions in real proteins, having first proven those peptides
occur nowhere in the inputs, so every later detection is attributable to the
spike-in alone and every miss is a true false negative.

```bash
Rscript scripts/generate_spikein_benchmark.R \
    --iedb databases/iedb/iedb_antigensoup_2026-09-08.tsv \
    --fasta proteins.faa \
    --n-spikes 1000 --seed 12345 \
    --outdir benchmark_spikein
```

Amino-acid sequences only. No nucleotide, CDS, translation, codon or
reading-frame handling, by design.

## How it works

1. Sample a real epitope from the metadata TSV written by `build_iedb_db.R`.
   That file contains only peptides the standard database retains, so the
   artificial set is modelled on epitopes AntigenSoup would actually search for,
   not on IEDB records the builder rejected.
2. Shuffle its residues. This preserves length and amino-acid composition
   exactly, and therefore also unique residue count, maximum residue fraction and
   composition-derived Shannon entropy. Only order-dependent properties, chiefly
   the longest homopolymer, can change.
3. Reject the candidate if it equals its source, already exists in the epitope
   database, already occurs anywhere in the input proteins, or duplicates an
   earlier artificial peptide. After 50 failed shuffles, draw a different source
   peptide rather than getting stuck.
4. Insert accepted peptides at random positions in randomly chosen proteins.
5. Validate the finished benchmark, then write it.

The generator deliberately stays this simple. A peptide shuffle is transparent,
auditable and exactly composition-matched; a generative model would be neither,
and there is no evidence shuffling fails. Each artificial peptide keeps its
source identifier and sequence in the truth table.

## Outputs

| File | Contents |
|---|---|
| `benchmark_spiked.fasta` | Every input protein, with artificial epitopes inserted. Identifiers unchanged; unspiked records byte-identical to the input. |
| `benchmark_epitopes.fasta` | The artificial epitopes, ready to pass to AntigenSoup as `-e`. |
| `benchmark_truth.tsv` | One row per spike: epitope, source, target, coordinates, absence checks, sequence metrics. |
| `benchmark_manifest.tsv` | Provenance, settings, checksums, rejection counts, validation results. |

The manifest is a two-column `field`/`value` TSV rather than JSON, to avoid a
dependency for a flat key-value record.

```
>SPIKE_000001
QRLALPFPGIDLFERGDNE
```

## Coordinates

Both coordinates are 1-based, and the distinction matters when a protein receives
more than one insertion.

- **`insertion_position_original`** is the position in the **original** protein
  *before which* the epitope was inserted. Original residues `1 .. pos-1` precede
  the epitope and `pos .. L` follow it.
- **`insertion_position_final`** is the position of the epitope's **first
  residue** in the **final** spiked protein. It equals
  `insertion_position_original` plus the combined length of every epitope
  inserted into the same protein at a strictly smaller coordinate.

With one insertion per protein the two are equal. Insertions into one protein are
applied from the highest coordinate downwards, so lower coordinates are never
displaced while the sequence is being built.

Coordinates within a protein are always distinct, so two artificial epitopes can
never abut or concatenate; at least one original residue separates them.

## Insertion placement

Positions are drawn uniformly from `margin+1 .. L-margin+1`, where `L` excludes
any trailing stop character and `margin` is `--terminal-margin` (default 5). A
protein must therefore be at least `2 x margin` residues long to be eligible.
Placement is uniform along the protein: in a 1000-spike run against the E. coli
K-12 proteome the median insertion sat at 50.0% of protein length, with 7.3% in
the first tenth and 7.7% in the last.

Targets are filled in passes, so every eligible protein is used once before any
is used twice. `--max-spikes-per-sequence` defaults to 1 and is raised
automatically, with a log message, only when the requested spike count cannot
otherwise be placed. There is no bias toward long proteins.

## Input validation

The protein FASTA is checked before anything else, and the utility fails with an
explicit message rather than guessing:

- unreadable or malformed FASTA
- missing, blank or duplicated identifiers (uniqueness is enforced on the first
  whitespace-delimited token, which is what AntigenSoup reports)
- empty sequences
- characters that are not amino acids, including alignment gaps and digits
- internal `*` stop characters, which indicate a raw six-frame translation

**Nucleotide input is detected and rejected.** A, C, G and T are all valid
amino-acid letters, so `Biostrings` would accept a DNA FASTA silently. If more
than 90% of records of at least 50 characters consist solely of A/C/G/T/U/N, the
run aborts.

**Prodigal-GV trailing stops are handled.** Prodigal writes a `*` at the end of
each protein. Those bytes are preserved exactly, so unspiked records stay
identical, but the trailing stop is excluded from the insertable region, so no
epitope is ever placed after it.

## Absence checking

Before acceptance a candidate must be absent from all three of:

- the epitope database, by exact set membership on the `sequence` column
- the input proteins, as a **substring** of the complete sequences, not merely
  compared against whole records
- the artificial peptides already accepted

The protein check concatenates all records into one `AAString` joined by `*`.
Because `*` cannot occur in a candidate, no match can straddle a record boundary,
and one `Biostrings::countPattern` call then scans the whole dataset.

One scope note worth stating plainly: the database check covers exactly the
peptides in the TSV you supply. Given the builder's output that is the retained
set, which is what AntigenSoup searches, so a spike-in can never collide with a
real database entry. It does **not** prove the peptide is absent from IEDB
records the builder rejected.

## Final validation

The run fails, loudly, unless all 13 checks pass:

every epitope is present in the spiked FASTA; occurs there exactly once; lies in
its recorded target protein; begins at its recorded final coordinate; was absent
from the original proteins; was absent from the epitope database; all epitopes
are unique; unspiked proteins are unchanged; spiked lengths match the truth
table; the spiked FASTA reads back with identical sequences and headers; the
epitope FASTA matches the truth table exactly; the truth row count equals the
requested spike count; spike identifiers are unique.

The "occurs exactly once" check is what catches an insertion junction
accidentally recreating another artificial peptide.

## Reproducibility

Identical inputs, parameters and `--seed` produce byte-identical outputs. Verified
by building twice with seed 12345 and comparing checksums of all three files, and
by confirming a different seed produces a different epitope set. The seed, both
input checksums and the full command line are recorded in the manifest.

## Options

```
--iedb FILE                     Epitope metadata TSV from build_iedb_db.R (required)
--fasta FILE                    Amino-acid protein FASTA to spike into (required)
--outdir DIR                    Output directory (default: benchmark_spikein)
--n-spikes N                    Artificial epitopes to insert (default: 100)
--seed N                        Random seed (default: 1)
--terminal-margin N             Residues kept free at each terminus (default: 5)
--max-spikes-per-sequence N     Cap per protein (default: 1, raised only if needed)
--max-generation-attempts N     Cap on total candidate shuffles
                                (default: larger of 50000 and 100 x --n-spikes)
-h, --help
```

## Measuring recovery

```bash
python scripts/ac_match.py \
    --epitopes benchmark_spikein/benchmark_epitopes.fasta \
    --proteins benchmark_spikein/benchmark_spiked.fasta \
    --out hits.tsv
```

Join `hits.tsv` on `epitope_ids` to `spike_id` in the truth table. Recovery should
be 100%, because these are exact insertions. `start` in the hit table should equal
`insertion_position_final`, which is an independent check of the coordinate
arithmetic. If anything is missing, investigate the cause rather than adjusting
the benchmark to make the test pass.

## Validation record

This is the reference run the generator was validated against. It is recorded so
the numbers below can be reproduced and so a future change that breaks something
is visible as a difference from a known result.

### Inputs

**Proteins.** The UniProt reference proteome for *Escherichia coli* K-12, strain
MG1655 / ATCC 47076, accession `UP000000625`, derived from genome assembly
`GCA_000005845.2`. Fetched as amino acids:

```bash
curl -sSL "https://rest.uniprot.org/uniprotkb/stream?query=proteome:UP000000625&format=fasta&compressed=true" \
  -o databases/test_data/ecoli_k12.faa.gz && gunzip databases/test_data/ecoli_k12.faa.gz
```

These are **gene products**, one canonical protein per protein-coding gene. Not
the genome, not gene nucleotide sequences.

| Property | Value |
|---|---|
| Records | 4,403 |
| Total residues | 1,354,442 |
| Median protein length | 271 aa |
| Records matching only A/C/G/T/N | 0, confirming these are not nucleotides |
| Header prefixes | all `sp|`, every entry manually reviewed Swiss-Prot |
| Trailing stop characters | 0 (UniProt omits them; Prodigal-GV adds them) |
| Selenocysteine (U) | 3 residues across 3 proteins |
| Unknown residues (X) | 8 residues across 7 proteins |

The U and X residues are real, not corruption: E. coli uses selenocysteine in its
formate dehydrogenases. They exercise the validator's tolerance of ambiguity codes
in input proteins, which are never generated in artificial epitopes.

**Epitopes.** `iedb_antigensoup_2026-09-08.tsv`, 240,971 retained peptides from
the IEDB release of 2026-09-08.

### Command

```bash
Rscript scripts/generate_spikein_benchmark.R \
    --iedb databases/iedb/iedb_antigensoup_2026-09-08.tsv \
    --fasta databases/test_data/ecoli_k12.faa \
    --n-spikes 1000 --seed 12345 --outdir benchmark_spikein
```

Runs in about 4 seconds.

### Results

| Measure | Result |
|---|---|
| Artificial epitopes requested / accepted | 1000 / 1000 |
| Candidate shuffles needed | 1000, i.e. no candidate was rejected |
| Built-in validation checks | 13 of 13 passed |
| Distinct target proteins used | 1000, one spike each |

**Recovery by the matching stage.** Running `ac_match.py` on the spiked proteome
with the artificial epitope FASTA:

| Measure | Result |
|---|---|
| Epitopes recovered | 1000 of 1000 (100%) |
| Found in the protein the truth table names | 1000 |
| Start coordinate equal to `insertion_position_final` | 1000 |

The coordinate agreement is an independent check of the insertion arithmetic,
because `ac_match.py` computes position without any knowledge of the truth table.

**Controls.**

| Search | Hits | Interpretation |
|---|---|---|
| Artificial epitopes vs the **spiked** proteome | 1000 | Every spike found, no extras |
| Artificial epitopes vs the **original** proteome | 0 | No false positives; the peptides really were novel |
| Real IEDB database vs the original proteome | 290 | Genuine background matches, for scale |

The middle row matters most. It confirms the generator's novelty proof using the
actual matching tool rather than the generator checking its own work.

**Property fidelity.** Comparing each artificial peptide with the source peptide
it was shuffled from:

| Property | Outcome |
|---|---|
| Peptide length | preserved exactly |
| Amino-acid composition | preserved exactly, maximum difference 0.0000 percentage points across all 20 residues |
| Shannon entropy | preserved exactly |
| Longest homopolymer | changed in 483 of 1000, as expected for an order-dependent property |

**Insertion placement.** Median insertion sat at 50.0% of protein length, with
7.3% in the first tenth and 7.7% in the last, consistent with uniform placement
inside the terminal margin.

**Determinism.** Two runs at seed 12345 produced byte-identical FASTA, epitope and
truth files. A run at seed 999 produced a different epitope set.

### What this run does and does not establish

It establishes that the matching stage, `ac_match.py`, reliably finds peptides
present in the proteins it is given, and does not report peptides that are absent.

It does **not** test the rest of the pipeline. AntigenSoup runs six stages: SRA
download, fastp, nohuman, MetaSPAdes assembly, Prodigal-GV gene prediction, then
epitope matching. Only the last was exercised, because the benchmark inserts into
proteins rather than into reads or contigs. Epitope loss during assembly or gene
calling is invisible to this benchmark by design, since nucleotide handling is out
of scope.

A curated UniProt proteome is also tidier than real input. A Prodigal-GV `.faa`
from a metagenomic assembly has contig-derived identifiers, a trailing `*` on every
protein, gene fragments truncated at contig ends, and is far larger. For results
you intend to report, point `--fasta` at one of those instead.

### Rejection-path testing

The reference run rejected no candidates, which is expected: a random shuffle of a
9mer colliding with 240,971 database entries or 1.35 million residues of proteome
has a probability of roughly 1 in a million. Zero counts therefore prove nothing,
so each path was forced separately with deliberately hostile fixtures.

| Fixture | Paths exercised |
|---|---|
| Epitope pool seeded with 900 permutations of one 8mer, proteins containing 800 more | `present_in_iedb` 45, `present_in_protein_fasta` 30, `duplicate_of_spike` 2 |
| Pool dominated by single-residue peptides | `identical_to_source` 100, `shuffle_attempts_exhausted` 2 |
| Pool containing only single-residue peptides | Fails cleanly at the attempt cap with an actionable message, rather than hanging |

All 13 validation checks passed in the first two cases, including with multiple
insertions per protein.

### Input rejection testing

Each of these fails with a specific message rather than being silently accepted:
nucleotide FASTA, duplicated identifiers, empty sequences, alignment gaps,
internal stop characters. Prodigal-style trailing stops are accepted, preserved
byte-for-byte, and excluded from the insertable region so no epitope is placed
after a stop.

## Dependencies

`Biostrings` (Bioconductor) for FASTA parsing and writing, amino-acid alphabet
validation, and exact substring search. `data.table` (CRAN) for the metadata TSV
and fast `%chin%` set membership. Argument parsing, checksums and shuffling use
base R.
