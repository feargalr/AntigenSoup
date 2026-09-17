#!/usr/bin/env Rscript
#
# generate_spikein_benchmark.R -- build an artificial epitope spike-in benchmark
# for measuring AntigenSoup sensitivity and specificity.
#
#   Rscript scripts/generate_spikein_benchmark.R \
#       --iedb databases/iedb/iedb_antigensoup_2026-09-08.tsv \
#       --fasta proteins.faa --n-spikes 1000 --seed 12345 \
#       --outdir benchmark_spikein
#
# Artificial epitopes are made by shuffling real retained IEDB epitopes, which
# preserves length and amino-acid composition exactly. Each candidate is proven
# absent from the IEDB database and from the input proteins before it is
# inserted, so any later detection is attributable to the spike-in alone.
#
# Amino-acid sequences only. No nucleotide, CDS, translation or frame handling.
#
# Requires: R with Biostrings (Bioconductor) and data.table (CRAN).
# See docs/spikein_benchmark.md.

## Version comes from the VERSION file at the repository root, the single source
## of truth. Reported in --help and written to the provenance record, so it falls
## back to "unknown" rather than to a stale hard-coded string.
read_version <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  here <- if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else NA
  cand <- c(if (!is.na(here)) file.path(dirname(here), "VERSION"),
            if (!is.na(here)) file.path(here, "VERSION"),
            file.path(getwd(), "VERSION"))
  for (f in cand) if (file.exists(f)) {
    v <- trimws(readLines(f, warn = FALSE)[1])
    if (nzchar(v)) return(v)
  }
  "unknown"
}
VERSION <- read_version()

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
})

## Shuffles tried per sampled source peptide before giving up on it and drawing
## a different one. Short or repetitive peptides have few distinct permutations.
MAX_SHUFFLE_ATTEMPTS <- 50L

## Residues an artificial epitope may be built from and inserted proteins may
## contain in their insertable region.
STD_AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

## Ambiguity codes tolerated in the input proteins but never generated.
AMBIG_AA <- c("X", "U", "O", "B", "Z", "J")


## ---------------------------------------------------------------------------
## Helpers
## ---------------------------------------------------------------------------

.start <- Sys.time()

log_msg <- function(...) cat(sprintf("[%5.0fs] %s\n",
  as.numeric(difftime(Sys.time(), .start, units = "secs")), paste0(...)),
  file = stderr())

die <- function(...) {
  cat(paste0("\nERROR: ", paste0(...), "\n"), file = stderr())
  quit(status = 1, save = "no")
}

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)


## ---------------------------------------------------------------------------
## Command line
## ---------------------------------------------------------------------------

USAGE <- paste0(
"generate_spikein_benchmark.R ", VERSION, " -- artificial epitope spike-in benchmark.

Usage:
  Rscript scripts/generate_spikein_benchmark.R --iedb <tsv> --fasta <faa> [options]

Required:
  --iedb FILE               Epitope metadata TSV from build_iedb_db.R
  --fasta FILE              Amino-acid protein FASTA to spike into

Options:
  --outdir DIR              Output directory (default: benchmark_spikein)
  --n-spikes N              Artificial epitopes to insert (default: 100)
  --seed N                  Random seed (default: 1)
  --terminal-margin N       Residues kept free at each protein terminus (default: 5)
  --max-spikes-per-sequence N
                            Cap per protein (default: 1; raised automatically only
                            if the requested spike count cannot otherwise be placed)
  --max-generation-attempts N
                            Cap on total candidate shuffles
                            (default: the larger of 50000 and 100 x --n-spikes)
  -h, --help                Show this message

Outputs, written to --outdir:
  benchmark_spiked.fasta    Input proteins with artificial epitopes inserted
  benchmark_epitopes.fasta  The artificial epitopes, for use as the AntigenSoup -e database
  benchmark_truth.tsv       One row per spike: epitope, target, coordinates, checks
  benchmark_manifest.tsv    Provenance, settings, checksums, rejection counts

Coordinates in the truth table are 1-based. insertion_position_original is the
position in the ORIGINAL protein before which the epitope was inserted;
insertion_position_final is the position of the epitope's FIRST residue in the
final spiked protein.
")

parse_args <- function(argv) {
  o <- list(iedb = NULL, fasta = NULL, outdir = "benchmark_spikein",
            n_spikes = 100L, seed = 1L, terminal_margin = 5L,
            max_per_seq = 1L, max_attempts = NA_integer_)
  i <- 1L
  need <- function(a) { if (i == length(argv)) die(sprintf("%s requires a value", a)); argv[i + 1L] }
  as_int <- function(a, v) {
    n <- suppressWarnings(as.integer(v))
    if (is.na(n)) die(sprintf("%s expects an integer, got '%s'", a, v))
    n
  }
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("-h", "--help")) { cat(USAGE); quit(status = 0, save = "no") }
    switch(a,
      "--iedb"                    = { o$iedb <- need(a); i <- i + 2L },
      "--fasta"                   = { o$fasta <- need(a); i <- i + 2L },
      "--outdir"                  = { o$outdir <- need(a); i <- i + 2L },
      "--n-spikes"                = { o$n_spikes <- as_int(a, need(a)); i <- i + 2L },
      "--seed"                    = { o$seed <- as_int(a, need(a)); i <- i + 2L },
      "--terminal-margin"         = { o$terminal_margin <- as_int(a, need(a)); i <- i + 2L },
      "--max-spikes-per-sequence" = { o$max_per_seq <- as_int(a, need(a)); i <- i + 2L },
      "--max-generation-attempts" = { o$max_attempts <- as_int(a, need(a)); i <- i + 2L },
      die(sprintf("unknown option '%s' (see --help)", a))
    )
  }
  if (is.null(o$iedb))  die("--iedb is required (see --help)")
  if (is.null(o$fasta)) die("--fasta is required (see --help)")
  for (f in c(o$iedb, o$fasta)) if (!file.exists(f)) die("file not found: ", f)
  if (o$n_spikes < 1L) die("--n-spikes must be at least 1")
  if (o$terminal_margin < 0L) die("--terminal-margin must be 0 or more")
  if (o$max_per_seq < 1L) die("--max-spikes-per-sequence must be at least 1")
  if (is.na(o$max_attempts)) o$max_attempts <- max(50000L, 100L * o$n_spikes)
  o
}


## ---------------------------------------------------------------------------
## Peptide metrics (same definitions as build_iedb_db.R, so the two are
## directly comparable)
## ---------------------------------------------------------------------------

peptide_metrics <- function(seqs) {
  n <- nchar(seqs)
  cnt <- vapply(STD_AA, function(a) n - nchar(gsub(a, "", seqs, fixed = TRUE)),
                integer(length(seqs)))
  if (length(seqs) == 1L) cnt <- matrix(cnt, nrow = 1L)
  p <- cnt / n
  H <- -rowSums(p * ifelse(p > 0, log2(p), 0))
  homo <- rep(1L, length(seqs)); k <- 2L; idx <- seq_along(seqs)
  while (length(idx) && k <= max(n)) {
    idx <- idx[grepl(sprintf("(.)\\1{%d}", k - 1L), seqs[idx])]
    if (length(idx)) homo[idx] <- k
    k <- k + 1L
  }
  data.table(peptide_length = n,
             unique_amino_acids = as.integer(rowSums(cnt > 0L)),
             max_aa_fraction = round(apply(cnt, 1L, max) / n, 4),
             shannon_entropy = round(H, 4),
             normalised_shannon_entropy = round(
               ifelse(log2(pmin(n, 20L)) > 0, H / log2(pmin(n, 20L)), 0), 4),
             longest_homopolymer = homo)
}


## ---------------------------------------------------------------------------
## Inputs
## ---------------------------------------------------------------------------

load_epitope_db <- function(path) {
  db <- fread(path, sep = "\t", showProgress = FALSE, colClasses = list(character = "sequence"))
  if (!"sequence" %in% names(db))
    die("the epitope table has no 'sequence' column.\n",
        "Expected the metadata TSV written by scripts/build_iedb_db.R.")
  ## Tolerate an all-sequences table by keeping only what the standard database
  ## would contain, so artificial peptides are modelled on retained epitopes.
  if ("included_in_default_db" %in% names(db)) {
    before <- nrow(db)
    db <- db[as.logical(included_in_default_db)]
    log_msg(sprintf("epitope table carries inclusion flags: kept %s of %s retained peptides",
                    fmt_n(nrow(db)), fmt_n(before)))
  }
  db[, sequence := toupper(trimws(sequence))]
  db <- db[nzchar(sequence)]
  if (!nrow(db)) die("the epitope table contains no usable sequences")
  if (!"antigensoup_id" %in% names(db)) db[, antigensoup_id := paste0("row", .I)]
  ## Only peptides over the standard 20 residues can be shuffled into a valid
  ## artificial peptide; the builder emits nothing else, so this should be a no-op.
  bad <- grepl("[^ACDEFGHIKLMNPQRSTVWY]", db$sequence)
  if (any(bad)) {
    log_msg(sprintf("WARNING: ignoring %s epitope(s) with non-standard residues",
                    fmt_n(sum(bad))))
    db <- db[!bad]
  }
  unique(db, by = "sequence")[]
}

read_protein_fasta <- function(path) {
  aa <- tryCatch(readAAStringSet(path),
    error = function(e) die("could not read '", path, "' as FASTA:\n  ", conditionMessage(e)))
  if (length(aa) == 0L) die("'", path, "' contains no FASTA records")

  ids <- names(aa)
  if (is.null(ids) || any(is.na(ids)) || any(!nzchar(trimws(ids))))
    die("every FASTA record must have a non-empty identifier; ",
        sum(is.null(ids) | is.na(ids) | !nzchar(trimws(ids))), " do not")
  ## AntigenSoup reports the first whitespace-delimited token, so that is what
  ## must be unique for the truth table to address a record.
  short <- sub("\\s.*$", "", ids)
  if (anyDuplicated(short)) {
    d <- unique(short[duplicated(short)])
    die("FASTA identifiers must be unique; ", length(d), " are repeated, e.g. ",
        paste(head(d, 3), collapse = ", "))
  }
  if (any(width(aa) == 0L))
    die(sum(width(aa) == 0L), " FASTA record(s) have an empty sequence")

  seqs <- as.character(aa)

  ## Guard against a nucleotide FASTA: A, C, G and T are all valid amino-acid
  ## letters, so Biostrings would accept one silently.
  nuc <- grepl("^[ACGTUN]+$", seqs) & nchar(seqs) >= 50L
  if (mean(nuc) > 0.9)
    die("'", path, "' looks like a NUCLEOTIDE FASTA: ",
        round(100 * mean(nuc)), "% of records contain only A/C/G/T/U/N.\n",
        "This utility operates on amino-acid sequences only. Translate first, ",
        "or supply the protein FASTA (for AntigenSoup, the Prodigal-GV .faa).")

  ## Prodigal-GV appends a stop character. Keep the bytes exactly as given so
  ## unspiked records stay identical, but exclude a trailing stop from the
  ## insertable region. An internal stop means a bad translation.
  core <- sub("\\*$", "", seqs)
  if (any(grepl("*", core, fixed = TRUE)))
    die(sum(grepl("*", core, fixed = TRUE)), " sequence(s) contain an internal '*' stop ",
        "character. Supply mature protein sequences, not raw six-frame translations.")
  n_stop <- sum(seqs != core)

  allowed <- paste0("[^", paste(c(STD_AA, AMBIG_AA), collapse = ""), "]")
  bad <- grepl(allowed, core)
  if (any(bad)) {
    ex <- unique(unlist(regmatches(core[bad], gregexpr(allowed, core[bad]))))
    die(sum(bad), " sequence(s) contain characters that are not amino acids: ",
        paste(sprintf("'%s'", head(ex, 8)), collapse = ", "), ".\n",
        "Gaps, digits and nucleotide ambiguity codes are not accepted.")
  }
  n_ambig <- sum(grepl(paste0("[", paste(AMBIG_AA, collapse = ""), "]"), core))

  list(ids = short, headers = ids, seqs = seqs, core_len = nchar(core),
       n_stop = n_stop, n_ambig = n_ambig)
}


## ---------------------------------------------------------------------------
## Artificial epitope generation
## ---------------------------------------------------------------------------

## Candidates are checked against the concatenated proteins rather than record by
## record. Records are joined with '*', which cannot occur in a candidate, so no
## match can straddle a boundary.
generate_epitopes <- function(db, prot_subject, n_spikes, max_attempts) {
  n_src <- nrow(db)
  src_seq <- db$sequence
  src_id  <- db$antigensoup_id
  iedb_set <- src_seq

  acc_seq <- character(n_spikes); acc_src <- integer(n_spikes)
  n_acc <- 0L; attempts <- 0L
  reject <- c(identical_to_source = 0L, present_in_iedb = 0L,
              present_in_protein_fasta = 0L, duplicate_of_spike = 0L,
              shuffle_attempts_exhausted = 0L, invalid_sequence = 0L)

  while (n_acc < n_spikes) {
    if (attempts >= max_attempts)
      die("gave up after ", fmt_n(attempts), " candidate shuffles with only ",
          fmt_n(n_acc), " of ", fmt_n(n_spikes), " epitopes accepted.\n",
          "Raise --max-generation-attempts, or lower --n-spikes.")
    s <- sample.int(n_src, 1L)
    src <- src_seq[s]
    chars <- strsplit(src, "", fixed = TRUE)[[1]]
    placed <- FALSE
    for (k in seq_len(MAX_SHUFFLE_ATTEMPTS)) {
      attempts <- attempts + 1L
      cand <- paste(sample(chars), collapse = "")
      if (identical(cand, src))                { reject["identical_to_source"] <- reject["identical_to_source"] + 1L; next }
      if (!grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cand)) { reject["invalid_sequence"] <- reject["invalid_sequence"] + 1L; next }
      if (cand %chin% iedb_set)                { reject["present_in_iedb"] <- reject["present_in_iedb"] + 1L; next }
      if (n_acc > 0L && cand %chin% acc_seq[seq_len(n_acc)]) {
                                                 reject["duplicate_of_spike"] <- reject["duplicate_of_spike"] + 1L; next }
      if (countPattern(cand, prot_subject) > 0L) {
                                                 reject["present_in_protein_fasta"] <- reject["present_in_protein_fasta"] + 1L; next }
      n_acc <- n_acc + 1L
      acc_seq[n_acc] <- cand; acc_src[n_acc] <- s
      placed <- TRUE
      break
    }
    if (!placed) reject["shuffle_attempts_exhausted"] <- reject["shuffle_attempts_exhausted"] + 1L
    if (n_acc %% 200L == 0L && placed) log_msg(sprintf("  generated %s / %s", fmt_n(n_acc), fmt_n(n_spikes)))
  }
  list(sequence = acc_seq, source_row = acc_src,
       source_id = src_id[acc_src], source_sequence = src_seq[acc_src],
       attempts = attempts, reject = reject)
}


## ---------------------------------------------------------------------------
## Target and position assignment
## ---------------------------------------------------------------------------

assign_targets <- function(core_len, n_spikes, margin, max_per_seq) {
  ## A protein can host an insertion only if at least `margin` residues can sit
  ## on each side of it, i.e. positions margin+1 .. core_len-margin+1 exist.
  n_pos <- core_len - 2L * margin + 1L
  eligible <- which(n_pos >= 1L)
  if (!length(eligible))
    die("no protein is long enough to accept an insertion with --terminal-margin ",
        margin, " (the longest is ", max(core_len), " aa; at least ", 2L * margin,
        " aa is required). Lower --terminal-margin.")

  ## Capacity is limited by the per-protein cap and by how many distinct
  ## insertion coordinates each protein actually offers.
  cap <- pmin(max_per_seq, n_pos[eligible])
  if (sum(cap) < n_spikes)
    die("cannot place ", fmt_n(n_spikes), " spikes: ", fmt_n(length(eligible)),
        " eligible protein(s) can hold at most ", fmt_n(sum(cap)),
        " with --max-spikes-per-sequence ", max_per_seq, ".\n",
        "Raise --max-spikes-per-sequence, lower --n-spikes, or lower --terminal-margin.")

  ## Fill by passes so every protein is used once before any is used twice.
  picked <- integer(0); used <- setNames(integer(length(eligible)), eligible)
  pass <- 0L
  while (length(picked) < n_spikes) {
    pass <- pass + 1L
    avail <- eligible[used[as.character(eligible)] < cap]
    if (!length(avail)) break
    take <- sample(avail)
    take <- head(take, n_spikes - length(picked))
    picked <- c(picked, take)
    used[as.character(take)] <- used[as.character(take)] + 1L
  }
  list(target = picked, passes = pass)
}


## ---------------------------------------------------------------------------
## Build
## ---------------------------------------------------------------------------

main <- function(argv) {
  opt <- parse_args(argv)
  cmd <- paste(c("Rscript scripts/generate_spikein_benchmark.R", argv), collapse = " ")
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$seed)

  log_msg("reading epitope database: ", opt$iedb)
  db <- load_epitope_db(opt$iedb)
  log_msg(sprintf("  %s unique retained epitopes", fmt_n(nrow(db))))

  log_msg("validating protein FASTA: ", opt$fasta)
  fa <- read_protein_fasta(opt$fasta)
  log_msg(sprintf("  %s records, %s aa total%s", fmt_n(length(fa$ids)),
                  fmt_n(sum(fa$core_len)),
                  if (fa$n_stop) sprintf(", %s with a trailing stop", fmt_n(fa$n_stop)) else ""))

  ## One subject for all absence checks.
  prot_subject <- AAString(paste(fa$seqs, collapse = "*"))

  log_msg("generating ", fmt_n(opt$n_spikes), " artificial epitopes by shuffling")
  gen <- generate_epitopes(db, prot_subject, opt$n_spikes, opt$max_attempts)
  log_msg(sprintf("  accepted %s from %s candidate shuffles",
                  fmt_n(opt$n_spikes), fmt_n(gen$attempts)))

  ## Targets. Raise the per-sequence cap only if the request cannot be met.
  max_per_seq <- opt$max_per_seq
  n_pos_all <- fa$core_len - 2L * opt$terminal_margin + 1L
  n_elig <- sum(n_pos_all >= 1L)
  if (opt$n_spikes > n_elig * max_per_seq && n_elig > 0L) {
    max_per_seq <- as.integer(ceiling(opt$n_spikes / n_elig))
    log_msg("  raising max spikes per sequence to ", max_per_seq,
            " (", fmt_n(opt$n_spikes), " spikes, ", fmt_n(n_elig), " eligible proteins)")
  }
  asg <- assign_targets(fa$core_len, opt$n_spikes, opt$terminal_margin, max_per_seq)
  target <- asg$target

  ## Distinct coordinates per protein, so two insertions can never concatenate.
  pos <- integer(opt$n_spikes)
  for (t in unique(target)) {
    idx <- which(target == t)
    lo <- opt$terminal_margin + 1L; hi <- fa$core_len[t] - opt$terminal_margin + 1L
    pos[idx] <- sort(sample(seq.int(lo, hi), length(idx)))
  }

  log_msg("inserting epitopes")
  spike_id <- sprintf("SPIKE_%06d", seq_len(opt$n_spikes))
  spiked <- fa$seqs
  final_pos <- integer(opt$n_spikes)
  for (t in unique(target)) {
    idx <- which(target == t)
    idx <- idx[order(pos[idx])]
    ## Offset of each insertion in the final sequence: everything inserted at a
    ## strictly smaller coordinate shifts it right.
    lens <- nchar(gen$sequence[idx])
    final_pos[idx] <- pos[idx] + c(0L, cumsum(lens)[-length(lens)])
    ## Apply from the highest coordinate down so lower ones are not displaced.
    s <- spiked[t]
    for (k in rev(idx))
      s <- paste0(substr(s, 1L, pos[k] - 1L), gen$sequence[k],
                  substr(s, pos[k], nchar(s)))
    spiked[t] <- s
  }

  ## Outputs ----------------------------------------------------------------
  log_msg("writing outputs")
  spiked_path <- file.path(opt$outdir, "benchmark_spiked.fasta")
  epi_path    <- file.path(opt$outdir, "benchmark_epitopes.fasta")
  truth_path  <- file.path(opt$outdir, "benchmark_truth.tsv")
  man_path    <- file.path(opt$outdir, "benchmark_manifest.tsv")

  out <- AAStringSet(spiked); names(out) <- fa$headers
  writeXStringSet(out, spiked_path)
  epi <- AAStringSet(gen$sequence); names(epi) <- spike_id
  writeXStringSet(epi, epi_path)

  m <- peptide_metrics(gen$sequence)
  truth <- data.table(
    spike_id = spike_id,
    synthetic_epitope = gen$sequence,
    peptide_length = m$peptide_length,
    source_iedb_id = gen$source_id,
    source_iedb_sequence = gen$source_sequence,
    target_fasta_id = fa$ids[target],
    original_protein_length = nchar(fa$seqs[target]),
    final_protein_length = nchar(spiked[target]),
    insertion_position_original = pos,
    insertion_position_final = final_pos,
    terminal_margin = opt$terminal_margin,
    seed = opt$seed,
    confirmed_absent_from_original_fasta = TRUE,
    confirmed_absent_from_iedb = TRUE,
    confirmed_unique_among_spikes = TRUE,
    unique_amino_acids = m$unique_amino_acids,
    max_aa_fraction = m$max_aa_fraction,
    shannon_entropy = m$shannon_entropy,
    normalised_shannon_entropy = m$normalised_shannon_entropy,
    longest_homopolymer = m$longest_homopolymer)
  setorder(truth, spike_id)
  fwrite(truth, truth_path, sep = "\t", quote = FALSE)

  ## Validation -------------------------------------------------------------
  log_msg("validating benchmark")
  checks <- validate(truth, spiked, fa, gen, prot_subject, db, spiked_path, epi_path)

  man <- data.table(field = c(
    "script", "script_version", "date_generated", "command",
    "protein_fasta", "protein_fasta_md5", "n_input_sequences", "total_input_aa",
    "iedb_database", "iedb_database_md5", "n_iedb_epitopes_available",
    "n_spikes_requested", "n_spikes_generated", "generation_method",
    "seed", "terminal_margin", "max_spikes_per_sequence",
    "n_target_proteins_used", "total_generation_attempts", "total_rejected",
    paste0("rejected_", names(gen$reject)),
    "r_version", "Biostrings_version", "data.table_version",
    paste0("validation_", names(checks))),
    value = c(
    "generate_spikein_benchmark.R", VERSION,
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), cmd,
    normalizePath(opt$fasta), unname(tools::md5sum(opt$fasta)),
    length(fa$ids), sum(fa$core_len),
    normalizePath(opt$iedb), unname(tools::md5sum(opt$iedb)), nrow(db),
    opt$n_spikes, nrow(truth), "shuffle_of_retained_iedb_epitope",
    opt$seed, opt$terminal_margin, max_per_seq,
    uniqueN(target), gen$attempts, sum(gen$reject),
    unname(gen$reject),
    as.character(getRversion()), as.character(packageVersion("Biostrings")),
    as.character(packageVersion("data.table")),
    ifelse(unlist(checks), "PASS", "FAIL")))
  fwrite(man, man_path, sep = "\t", quote = FALSE)

  cat("\n")
  for (n in names(checks))
    cat(sprintf("  [%s] %s\n", ifelse(checks[[n]], "PASS", "FAIL"), gsub("_", " ", n)))
  cat("\n  candidate rejections:\n")
  for (n in names(gen$reject)) cat(sprintf("    %-28s %8s\n", n, fmt_n(gen$reject[[n]])))
  cat(sprintf("\n  outputs in %s\n", normalizePath(opt$outdir)))

  if (!all(unlist(checks)))
    die("benchmark validation failed; see the checks above. The outputs were ",
        "written so the failure can be inspected, but must not be used as a benchmark.")
  log_msg("done.")
}


## ---------------------------------------------------------------------------
## Validation of the finished benchmark
## ---------------------------------------------------------------------------

validate <- function(truth, spiked, fa, gen, prot_subject, db, spiked_path, epi_path) {
  final_subject <- AAString(paste(spiked, collapse = "*"))
  tgt <- match(truth$target_fasta_id, fa$ids)

  ## Each epitope appears in the final data exactly once, in its recorded target,
  ## starting at its recorded final coordinate.
  occ <- vapply(truth$synthetic_epitope, function(p) countPattern(p, final_subject),
                integer(1), USE.NAMES = FALSE)
  in_target <- vapply(seq_len(nrow(truth)), function(i)
    grepl(truth$synthetic_epitope[i], spiked[tgt[i]], fixed = TRUE), logical(1))
  at_coord <- substr(spiked[tgt], truth$insertion_position_final,
                     truth$insertion_position_final + truth$peptide_length - 1L) ==
              truth$synthetic_epitope

  ## Absence from the untouched inputs.
  absent_orig <- vapply(truth$synthetic_epitope, function(p) countPattern(p, prot_subject),
                        integer(1), USE.NAMES = FALSE) == 0L
  absent_iedb <- !(truth$synthetic_epitope %chin% db$sequence)

  ## Untouched proteins are byte-identical to the input.
  untouched <- setdiff(seq_along(fa$seqs), unique(tgt))
  reread <- readAAStringSet(spiked_path)
  epi_re <- readAAStringSet(epi_path)

  list(
    every_epitope_present_in_spiked_fasta   = all(occ >= 1L),
    every_epitope_occurs_exactly_once       = all(occ == 1L),
    every_epitope_in_recorded_target        = all(in_target),
    every_epitope_at_recorded_final_coord   = all(at_coord),
    every_epitope_absent_from_original      = all(absent_orig),
    every_epitope_absent_from_iedb          = all(absent_iedb),
    all_epitopes_unique                     = !anyDuplicated(truth$synthetic_epitope),
    unspiked_proteins_unchanged             = all(spiked[untouched] == fa$seqs[untouched]),
    spiked_length_matches_truth             = all(nchar(spiked[tgt]) == truth$final_protein_length),
    spiked_fasta_reads_back_as_protein      = length(reread) == length(fa$seqs) &&
                                              identical(unname(as.character(reread)), unname(spiked)) &&
                                              identical(names(reread), fa$headers),
    epitope_fasta_valid_and_matches_truth   = length(epi_re) == nrow(truth) &&
                                              identical(names(epi_re), truth$spike_id) &&
                                              identical(unname(as.character(epi_re)), truth$synthetic_epitope),
    truth_rows_equal_requested_spikes       = nrow(truth) == length(gen$sequence),
    spike_ids_unique                        = !anyDuplicated(truth$spike_id)
  )
}

if (sys.nframe() == 0L && !interactive()) main(commandArgs(trailingOnly = TRUE))
