# AntigenSoup
**Version:** 0.5.1

This pipeline was built for metagenomic assembly and the identification of cross-reactive epitopes in the human microbiome. AntigenSoup is a wrapper pipeline that orchestrates several established bioinformatics tools, and so we *strongly encourage* users cite the underlying software components appropriately in any resulting publications.

AntigenSoup searches for short, exact peptide matches between predicted microbial proteins and known epitopes, typically in the 8–15 amino acid range. This reflects the biology of antigen presentation, as MHC class I epitopes are usually 8–11 aa, while MHC class II epitopes contain shorter core motifs. T cell cross-reactivity can arise from identical short peptides embedded within otherwise unrelated proteins. For this reason, epitope detection is treated as a string-matching problem, rather than a homology search.
To efficiently detect these matches at scale, AntigenSoup uses the Aho–Corasick algorithm, which enables simultaneous, exact matching of millions of epitope sequences against large protein databases in a single pass.

## **Database**
The epitope database is built from the current IEDB release with
`scripts/build_iedb_db.R`, rather than shipped as a file in this repository. A
build takes about five minutes and records exactly which IEDB release it came
from and how it was filtered, so results stay traceable.

Earlier versions bundled a pre-filtered `iedb.fasta.gz`. It was removed in v0.5.0
because it carried no record of how it had been filtered, and was missing whole
peptide lengths as a result. It remains in the git history if you need it for
comparison; the last commit that shipped it is `396d94b`:

```bash
git show 396d94b:iedb.fasta.gz > iedb.fasta.gz
```

### Building a database from the current IEDB release

`scripts/build_iedb_db.R` rebuilds the epitope database from the current IEDB
database export. One command, no manual downloads:

```bash
Rscript scripts/build_iedb_db.R --outdir databases/iedb
```

It fetches the IEDB export, keeps linear peptides with positive experimental
T cell or B cell evidence, filters on length, amino-acid alphabet and sequence
complexity, collapses duplicate sequences while aggregating their evidence, and
writes three files:

- `iedb_antigensoup_<RELEASE>.fasta` — the database, ready for `-e`
- `iedb_antigensoup_<RELEASE>.tsv` — one metadata row per FASTA record
- `iedb_antigensoup_<RELEASE>.report.txt` — provenance, settings and retention counts

The FASTA drops straight into the pipeline:

```bash
antigensoup --scaffolds my_assembly.fasta -e databases/iedb/iedb_antigensoup_<RELEASE>.fasta -n 16
```

Headers are short, stable identifiers, so `ac_match.py` output stays readable and
joins back to the metadata TSV on `epitope_ids`:

```
>AS_0000001
AAAACTTMK
```

Defaults build the recommended database with no extra arguments. To vary it:

```bash
--include-mhc-only        # also include eluted / binding-only MHC ligands
--max-length 50           # broader database for longer linear B cell epitopes
--natural-only            # also drop neo-epitopes
--include-low-complexity  # skip the complexity filter
--force-download          # refresh the cached IEDB export
```

**Requirements:** R with `data.table` (`install.packages("data.table")`), plus the
`curl` and `unzip` commands. Nothing else. A build takes about 5 minutes, uses
~6GB of memory, and caches ~530MB of IEDB archives in `<outdir>/source` for
reuse, so rerunning with the same `--outdir` does not re-download.

See [docs/iedb_database_build.md](docs/iedb_database_build.md) for the biological
rationale behind each default filter, what the evidence categories mean, and a
description of every output column.

### Benchmarking detection with artificial spike-ins

`scripts/generate_spikein_benchmark.R` builds a controlled benchmark for
measuring sensitivity and specificity. It generates artificial peptides by
shuffling real retained IEDB epitopes, proves each one is absent from both the
epitope database and the input proteins, then inserts them at known positions:

```bash
Rscript scripts/generate_spikein_benchmark.R \
    --iedb databases/iedb/iedb_antigensoup_<RELEASE>.tsv \
    --fasta proteins.faa \
    --n-spikes 1000 --seed 12345 \
    --outdir benchmark_spikein
```

This writes a spiked protein FASTA, an epitope FASTA usable directly as `-e`, a
ground-truth table with 1-based coordinates, and a manifest. Because every
insertion is exact and provably novel, recovery should be 100% and any miss is a
true false negative. Amino-acid sequences only.

See [docs/spikein_benchmark.md](docs/spikein_benchmark.md) for coordinate
semantics, options, and the recorded validation run against the E. coli K-12
reference proteome (1000 spikes, 100% recovery, zero false positives).

## **Epitope Matching Strategy**

### Exact matching at scale
At the scale AntigenSoup operates — hundreds of thousands of epitopes searched across large metagenomic datasets — the epitope search uses **exact string matching only**. This is a deliberate design choice. Short peptides (8–16 aa) are particularly prone to false positives under mismatch-tolerant search: allowing even 1–2 substitutions in an 8-mer permits 12–25% sequence divergence, which at metagenomic scale produces an unacceptable number of spurious hits. Exact matching is unambiguous, interpretable, and fast.

### Variant search for specific epitopes
We recognise that users may wish to search for near-identical variants of specific epitopes of interest — for example, to identify microbial mimics of a known T cell epitope that differ by one or two residues. For these targeted use cases, AntigenSoup provides `gen_variants.py`, a utility that generates a FASTA file of all sequences within a given Hamming distance of one or more query epitopes. This variant FASTA can then be used directly as input to `ac_match.py` for exact matching, with full provenance encoded in the sequence headers.

This approach is intentionally kept as a separate, opt-in step rather than a built-in mismatch mode. Variant expansion should only be applied to specific epitopes a user has prior reason to care about — not applied globally — in order to keep false positive rates under control.

**Generating variants for a single epitope:**
```bash
python gen_variants.py SIINFEKL --mismatches 1
# Output: SIINFEKL_variants_d1.fasta
```

**Generating variants from a FASTA of epitopes of interest:**
```bash
python gen_variants.py my_epitopes.fasta --mismatches 1
# Output: my_epitopes_variants_d1.fasta
```

**Then search with ac_match as normal:**
```bash
python ac_match.py --epitopes SIINFEKL_variants_d1.fasta --proteins my_proteins.fasta --out hits.tsv
```

Each matched variant is fully traceable via the header format `>original_epitope|d<N>|substitutions|original_sequence`, for example:
```
>SIINFEKL|d1|S1C|SIINFEKL
CIINFEKL
```
This encodes the original epitope name, the Hamming distance, the exact substitutions (position and amino acid change), and the original sequence — so downstream analysis can easily filter by distance or identify which positions tolerate substitution.

## **Overview**
**AntigenSoup** assembles metagenomes, predicts genes, and identifies cross-reactive epitopes from metagenomic data. It integrates:

- **SRA-Tools** - for downloading sequenceing data from the SRA (optional)
- **fastp** - for read quality control and filtering
- **nohuman**- for remove of human reads
- **MetaSPAdes** – for assembly  
- **Prodigal-GV** – for gene prediction  
- **Aho-Corasick algorithm** – for searching epitope sequences  
- **Pigz** – for compressing output  

Everything is wrapped into a single executable: `antigensoup`. You can simply give AntigenSoup an SRA ID and it will return a list of epitopes identified in that sample post-assembly.

## **Pipeline Workflow**

1. **Input acquisition**
   - Download reads from the SRA using `fasterq-dump`, or
   - Use local paired-end FASTQ files.

2. **Read quality control**
   - Quality filtering and adapter trimming using **fastp**.

3. **Host read removal**
   - Removal of human reads using **nohuman**. Please see the [nohuman GitHub repository](https://github.com/mbhall88/nohuman) for the most up-to-date recommendations on which database to use.
4. **Metagenomic assembly**
   - Assembly of filtered reads using **MetaSPAdes**.

5. **Gene prediction**
   - Prediction of protein-coding genes using **Prodigal-GV**.

6. **Epitope search**
   - Identification of epitope matches in predicted proteins using **Aho–Corasick string-matching algorithm implemented in acmatch in Python**.

7. **Output compression**
   - Compression of large intermediate and final outputs using **pigz**.

## **Installation**

### **1. Clone the repository**
```bash
git clone https://github.com/feargalr/AntigenSoup.git
cd AntigenSoup
```

### **2. Install dependencies**
```bash

#First create conda envs
conda env create -f conda_ymls/antigensoup_env.yml
conda env create -f conda_ymls/acmatch_env.yml

#Second ensure AntigenSoup scripts are in paths for individual envs
conda activate antigensoup
cp antigensoup "$CONDA_PREFIX/bin/antigensoup"

conda activate acmatch
mkdir -p "$CONDA_PREFIX/share/antigensoup"
cp scripts/ac_match.py "$CONDA_PREFIX/share/antigensoup/ac_match.py"
cp scripts/gen_variants.py "$CONDA_PREFIX/share/antigensoup/gen_variants.py"


#Third.For multi-threaded gene prediction we use the parallel-prodigal-gv.py
git clone https://github.com/apcamargo/prodigal-gv
conda activate antigensoup
mkdir -p "$CONDA_PREFIX/share/antigensoup"
cp prodigal-gv/parallel-prodigal-gv.py "$CONDA_PREFIX/share/antigensoup/parallel-prodigal-gv.py"

# Fourth. Download the nohuman db
conda activate antigensoup
nohuman --download --db /example_directory/nohuman_db

# Add this to your ~/.bashrc (or ~/.zshrc)
export NOHUMAN_DB="/example_directory/nohuman_db"

# Fifth. Build your epitope database (needs R with data.table)
Rscript scripts/build_iedb_db.R --outdir databases/iedb

#Alternatively prepare your own database of sequences. 

```


## **Usage**
```bash
Usage:
  antigensoup -e <epitope_fasta> -n <num_cores> -m <memory_gb> [--sra <SRA_ID>] [--read1 <read1.fastq.gz> --read2 <read2.fastq.gz>]
  antigensoup -e <epitope_fasta> -n <num_cores> --scaffolds <scaffolds.fasta>

Flags:
  -e, --epitopes        Path to epitope FASTA file (required)
  -n, --num-cores       Number of cores (default: 1)
  -m, --memory          Memory in GB (default: 8)
  --sra <SRA_ID>        SRA accession ID (optional)
  --read1 <file>        Path to local R1 FASTQ file (required if no SRA)
  --read2 <file>        Path to local R2 FASTQ file (required if no SRA)
  --scaffolds <file>    Path to pre-assembled scaffolds FASTA (skips QC, host removal, and assembly)
  --nohuman-db <path>   Path to nohuman database (or set $NOHUMAN_DB; not required with --scaffolds)
  -h, --help            Show this help message
  -V, --version         Print version and exit

Notes:
- If --sra is provided, local --read1 and --read2 are ignored.
- If no SRA is provided, both --read1 and --read2 must be specified.
- If --scaffolds is provided, all upstream steps (fastp, nohuman, SPAdes) are skipped.
  Gene prediction and epitope search run directly on the provided scaffolds.
  Output is written to <scaffolds_basename>_output/.

Examples:
  antigensoup --sra SRR123456 -e epitopes.fasta -n 32 -m 64
  antigensoup --read1 sample_1.fastq.gz --read2 sample_2.fastq.gz -e epitopes.fasta -n 16 -m 32
  antigensoup --scaffolds my_assembly.fasta -e epitopes.fasta -n 16
```

## **Inputs**
- **Epitope FASTA**: Required – a list of peptide sequences to search for.
- **SRA accession** *or* paired FASTQ files: Provide one or the other.

## **Output**
- `sra_fastq/`: SRA downloaded files.
- `fastp_output/`: Quality filtered reads and nohuman removed reads.
- `spades_output/`: Contains assembled scaffolds and protein predictions.
- `ac_matches.txt`: List of matched epitopes and corresponding genes


## **Notes**
- SRA downloads are handled via `fasterq-dump`.
- Intermediate and final files are compressed with `pigz` for efficiency.
- Designed for paired-end Illumina shotgun metagenomic data
- Assumes human host for read removal
- Assembly-based approach may miss low-abundance epitopes
