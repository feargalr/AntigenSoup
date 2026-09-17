#!/usr/bin/env Rscript
#
# build_iedb_db.R -- build an AntigenSoup epitope FASTA from the IEDB export.
#
#   Rscript scripts/build_iedb_db.R --outdir databases/iedb
#
# Downloads the current IEDB database export, keeps linear peptides with
# positive experimental T-cell or B-cell evidence, collapses duplicate
# sequences, and writes three files:
#
#   iedb_antigensoup_<release>.fasta       the database, ready for -e
#   iedb_antigensoup_<release>.tsv         one metadata row per FASTA record
#   iedb_antigensoup_<release>.report.txt  provenance, settings, retention
#
# Requires: R with data.table, plus the `curl` and `unzip` commands.
# Rationale for every default is in docs/iedb_database_build.md.

VERSION <- "1.0.0"

suppressPackageStartupMessages(library(data.table))


## ---------------------------------------------------------------------------
## IEDB source tables
## ---------------------------------------------------------------------------
##
## The exports carry two header rows: a group row and a name row, so a column is
## identified by the pair, e.g. `Assay :: Qualitative Measurement`. The indices
## below were read off that header. They are re-checked on every run against
## `verify`, so a schema change at IEDB fails loudly instead of silently
## producing a wrong database.
##
## data.table cannot read the single-file B-cell and MHC exports (3.2GB and
## 9.2GB; fread hits a 2GB limit and then a bus error), so the multi-file
## variants are used. Only member 00 carries the header rows.

BASE_URL <- "https://www.iedb.org/downloader.php?file_name=doc/"

SOURCES <- list(
  epitope = list(
    zip = "epitope_full_v3.zip", members = "epitope_full_v3.csv",
    cols = c(epitope_iri = 1L, object_type = 2L, name = 3L, modifications = 5L,
             source_antigen = 10L, source_organism = 16L, relation = 18L),
    verify = c("1" = "IEDB IRI", "2" = "Object Type", "3" = "Name",
               "5" = "Modifications", "16" = "Species", "18" = "Epitope Relation")
  ),
  tcell = list(
    zip = "tcell_full_v3.zip", members = "tcell_full_v3.csv",
    cols = c(ref_iri = 2L, epitope_iri = 10L, host = 44L, qualitative = 123L,
             mhc_allele = 142L, mhc_class = 146L),
    verify = c("2" = "IEDB IRI", "10" = "IEDB IRI", "44" = "Name",
               "123" = "Qualitative Measurement", "142" = "Name", "146" = "Class")
  ),
  bcell = list(
    zip = "bcell_full_v3_multi_file.zip",
    members = sprintf("bcell_full_v3_%02d.csv", 0:2),
    ## Note the name: the B-cell table says "Measure", the other two "Measurement".
    cols = c(ref_iri = 2L, epitope_iri = 10L, host = 44L, qualitative = 103L),
    verify = c("2" = "IEDB IRI", "10" = "IEDB IRI", "44" = "Name",
               "103" = "Qualitative Measure")
  ),
  mhc = list(
    zip = "mhc_ligand_full_multi_file.zip",
    members = sprintf("mhc_ligand_full_%02d.csv", 0:9),
    cols = c(ref_iri = 2L, epitope_iri = 10L, host = 44L, qualitative = 95L,
             mhc_allele = 108L, mhc_class = 112L),
    verify = c("2" = "IEDB IRI", "10" = "Epitope IRI", "44" = "Name",
               "95" = "Qualitative Measurement", "108" = "Name", "112" = "Class")
  )
)

STD_AA <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]

## Complexity thresholds. Fixed rather than exposed because they were calibrated
## against the real IEDB distribution: together they remove ~1.5% of otherwise
## eligible peptides, and the entropy cut sits at about the 1.5th percentile.
## `--include-low-complexity` turns the whole filter off.
MIN_UNIQUE_AA    <- 4L
MAX_AA_FRACTION  <- 0.50
MAX_HOMOPOLYMER  <- 4L
MIN_NORM_ENTROPY <- 0.55

## Cap on how many values a list-valued metadata column holds.
MAX_LIST_ITEMS <- 5L

## `Related Object :: Epitope Relation`. An empty value means the epitope has no
## related object, i.e. it is a straight fragment of a natural molecule. Matched
## with match() on the names, since NATURE_MAP[""] would return NA.
NATURE_MAP <- stats::setNames(
  c("NATURAL", "ANALOG", "MIMOTOPE", "NEOEPITOPE", "NEOEPITOPE", "NEOEPITOPE",
    "NEOEPITOPE"),
  c("", "analog", "mimotope", "in-frame neo-epitope", "frameshift neo-epitope",
    "fusion neo-epitope", "unspecified neo-epitope")
)


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

fmt_bytes <- function(n) {
  u <- c("B", "kB", "MB", "GB"); i <- 1L
  while (n >= 1024 && i < length(u)) { n <- n / 1024; i <- i + 1L }
  sprintf("%.1f %s", n, u[i])
}


## ---------------------------------------------------------------------------
## Command line
## ---------------------------------------------------------------------------

USAGE <- paste0(
"build_iedb_db.R ", VERSION, " -- build an AntigenSoup epitope database from IEDB.

Usage:
  Rscript scripts/build_iedb_db.R [options]

  --outdir DIR              Output directory (default: databases/iedb)
                            IEDB archives are cached in <outdir>/source and reused.
  --min-length N            Minimum peptide length (default: 8)
  --max-length N            Maximum peptide length (default: 25)
  --include-mhc-only        Also include peptides whose only positive evidence is
                            MHC binding or elution (roughly doubles the database)
  --include-low-complexity  Skip the sequence complexity filter
  --natural-only            Also drop neo-epitopes, keeping only peptides recorded
                            as fragments of a natural molecule
  --force-download          Re-download the IEDB export even if cached
  --threads N               data.table threads (default: data.table's own default)
  -h, --help                Show this message

With no options this builds the recommended database: linear 8-25mer peptides
over the 20 standard amino acids, each with at least one positive experimental
T-cell or B-cell assay, excluding analogues, mimotopes, peptides only ever seen
chemically modified, and low-complexity sequences.
")

parse_args <- function(argv) {
  o <- list(outdir = "databases/iedb", min_length = 8L, max_length = 25L,
            include_mhc_only = FALSE, include_low_complexity = FALSE,
            natural_only = FALSE, force_download = FALSE, threads = 0L)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (a %in% c("-h", "--help")) { cat(USAGE); quit(status = 0, save = "no") }
    int_arg <- function() {
      if (i == length(argv)) die(sprintf("%s requires a value", a))
      n <- suppressWarnings(as.integer(argv[i + 1L]))
      if (is.na(n)) die(sprintf("%s expects an integer, got '%s'", a, argv[i + 1L]))
      n
    }
    switch(a,
      "--outdir"                 = { if (i == length(argv)) die("--outdir requires a value")
                                     o$outdir <- argv[i + 1L]; i <- i + 2L },
      "--min-length"             = { o$min_length <- int_arg(); i <- i + 2L },
      "--max-length"             = { o$max_length <- int_arg(); i <- i + 2L },
      "--threads"                = { o$threads <- int_arg(); i <- i + 2L },
      "--include-mhc-only"       = { o$include_mhc_only <- TRUE; i <- i + 1L },
      "--include-low-complexity" = { o$include_low_complexity <- TRUE; i <- i + 1L },
      "--natural-only"           = { o$natural_only <- TRUE; i <- i + 1L },
      "--force-download"         = { o$force_download <- TRUE; i <- i + 1L },
      die(sprintf("unknown option '%s' (see --help)", a))
    )
  }
  if (o$min_length < 1L) die("--min-length must be at least 1")
  if (o$max_length < o$min_length) die("--max-length must be >= --min-length")
  o
}


## ---------------------------------------------------------------------------
## Download and read
## ---------------------------------------------------------------------------

fetch <- function(zipname, cache_dir, force) {
  dest <- file.path(cache_dir, zipname)
  if (file.exists(dest) && !force) {
    log_msg(sprintf("cached   %-32s %s", zipname, fmt_bytes(file.size(dest))))
    return(dest)
  }
  url <- paste0(BASE_URL, zipname)
  log_msg(sprintf("download %-32s ...", zipname))
  tmp <- paste0(dest, ".part")
  st <- suppressWarnings(system2("curl", c("-fsSL", "--retry", "3",
    "--connect-timeout", "30", "-o", shQuote(tmp), shQuote(url))))
  ## IEDB answers an unknown file name with HTTP 200 and an empty body, so -f
  ## does not catch it; the size check does.
  if (!identical(st, 0L) || !file.exists(tmp) || file.size(tmp) < 1000) {
    unlink(tmp)
    die(sprintf("could not download %s\n  from %s\n  curl exit %s.\n", zipname, url, st),
        "Either the network is unavailable or IEDB has renamed this export.\n",
        "You can also place the file manually in ", cache_dir)
  }
  file.rename(tmp, dest)
  log_msg(sprintf("         %-32s %s", zipname, fmt_bytes(file.size(dest))))
  dest
}

read_iedb <- function(src, zippath, label) {
  ## Read the two header rows to confirm the schema and learn the column count.
  con <- pipe(sprintf("unzip -p %s %s | head -2", shQuote(zippath),
                      shQuote(src$members[1])), "r")
  hdr <- utils::read.csv(con, header = FALSE, nrows = 2L, colClasses = "character",
                         check.names = FALSE)
  close(con)
  if (nrow(hdr) != 2L) die(sprintf("could not read the header of %s", src$zip))
  names <- as.character(hdr[2, ]); ncol <- ncol(hdr)
  idx <- as.integer(names(src$verify))
  bad <- which(idx > ncol | names[idx] != unname(src$verify))
  if (length(bad))
    die(sprintf("the IEDB '%s' export no longer matches the expected schema:\n", label),
        paste(sprintf("  column %d: expected '%s', found '%s'", idx[bad],
                      unname(src$verify)[bad],
                      ifelse(idx[bad] > ncol, "<missing>", names[idx[bad]])),
              collapse = "\n"),
        "\n\nThe column indices in SOURCES need updating before rebuilding.")

  ## fread mis-detects the column count on the headerless members, so pin it
  ## with fill = ncol.
  sel <- unname(src$cols); nm <- names(src$cols); ord <- order(sel)
  d <- rbindlist(lapply(src$members, function(m)
    fread(cmd = sprintf("unzip -p %s %s", shQuote(zippath), shQuote(m)),
          header = FALSE, sep = ",", fill = ncol,
          skip = if (m == src$members[1]) 2L else 0L,
          select = sel[ord], col.names = nm[ord],
          colClasses = "character", showProgress = FALSE)))
  setcolorder(d, nm)
  ## A mis-parsed member shows up at once as a malformed IRI.
  n_bad <- d[!grepl("^https?://www\\.iedb\\.org/epitope/[0-9]+$", epitope_iri), .N]
  if (n_bad > 0L)
    die(sprintf("%s of %s rows in the IEDB '%s' export did not parse cleanly. ",
                fmt_n(n_bad), fmt_n(nrow(d)), label),
        "Re-run with --force-download; if that does not help, the format has changed.")
  log_msg(sprintf("read     %-32s %s rows", src$zip, fmt_n(nrow(d))))
  d[, epitope_id := as.integer(sub(".*/", "", epitope_iri))][, epitope_iri := NULL]
  d[]
}


## ---------------------------------------------------------------------------
## Sequences
## ---------------------------------------------------------------------------

## IEDB writes a chemically modified epitope as "SEQUENCE + MOD(pos)", e.g.
## "AAIAAVKEEAF + METH(A10)". The residues before the " + " are the parent
## sequence, which is what a microbial protein could actually match. Residues are
## never rewritten.
parse_sequence <- function(name)
  toupper(gsub("[[:space:]]", "", sub(" \\+ .*$", "", name)))

## Per-residue counts over the standard alphabet. rowSums() equalling nchar() is
## also the alphabet check: anything else simply is not counted.
complexity <- function(seqs) {
  n <- nchar(seqs)
  cnt <- vapply(STD_AA, function(a) n - nchar(gsub(a, "", seqs, fixed = TRUE)),
                integer(length(seqs)))
  if (length(seqs) == 1L) cnt <- matrix(cnt, nrow = 1L)
  p <- cnt / n
  H <- -rowSums(p * ifelse(p > 0, log2(p), 0))
  ## Normalise by the entropy attainable at this length, not log2(20), so short
  ## and long peptides are comparable.
  homo <- rep(1L, length(seqs)); k <- 2L; idx <- seq_along(seqs)
  while (length(idx) && k <= max(n)) {
    idx <- idx[grepl(sprintf("(.)\\1{%d}", k - 1L), seqs[idx])]
    if (length(idx)) homo[idx] <- k
    k <- k + 1L
  }
  data.table(length = n,
             n_alphabet = as.integer(rowSums(cnt)),
             unique_aa = as.integer(rowSums(cnt > 0L)),
             max_aa_frac = apply(cnt, 1L, max) / n,
             norm_entropy = ifelse(log2(pmin(n, 20L)) > 0, H / log2(pmin(n, 20L)), 0),
             homopolymer = homo)
}

## Collapse a column to a bounded "; " list per sequence, ordered by frequency
## then alphabetically so the output is deterministic.
agg_list <- function(dt, col, out) {
  x <- dt[nzchar(get(col)), .N, by = c("sequence", col)]
  setnames(x, col, "value")
  if (!nrow(x)) return(setnames(data.table(sequence = character(), v = character()),
                                "v", out))
  setorderv(x, c("sequence", "N", "value"), c(1L, -1L, 1L))
  tot <- x[, .(n = .N), by = sequence]
  x[, rk := rowid(sequence)]
  r <- x[rk <= MAX_LIST_ITEMS, .(v = paste(value, collapse = "; ")), by = sequence][tot, on = "sequence"]
  r[n > MAX_LIST_ITEMS, v := paste0(v, sprintf("; ...(+%d more)", n - MAX_LIST_ITEMS))]
  setnames(r[, .(sequence, v)], "v", out)
}


## ---------------------------------------------------------------------------
## Build
## ---------------------------------------------------------------------------

main <- function(argv) {
  opt <- parse_args(argv)
  for (b in c("unzip", "curl"))
    if (!nzchar(Sys.which(b))) die(sprintf("the `%s` command is required but was not found on PATH.", b))
  if (opt$threads > 0L) setDTthreads(opt$threads)

  cache <- file.path(opt$outdir, "source")
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  built_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  ## 1. source ---------------------------------------------------------------
  log_msg("stage 1/6  acquiring IEDB export")
  zips <- vapply(SOURCES, function(s) fetch(s$zip, cache, opt$force_download), character(1))
  src <- rbindlist(lapply(names(SOURCES), function(k) {
    li <- utils::unzip(zips[[k]], list = TRUE)
    data.table(file = basename(zips[[k]]), bytes = file.size(zips[[k]]),
               md5 = unname(tools::md5sum(zips[[k]])),
               built = format(max(li$Date), "%Y-%m-%d"))
  }))
  ## The export has no version number, so its build date is the release stamp.
  release <- max(src$built)
  log_msg("IEDB release ", release)

  ## 2. epitopes -------------------------------------------------------------
  log_msg("stage 2/6  reading epitope table")
  ep <- read_iedb(SOURCES$epitope, zips[["epitope"]], "epitope")
  n_read <- nrow(ep)

  unknown <- setdiff(unique(ep$relation), names(NATURE_MAP))
  if (length(unknown))
    log_msg("WARNING: unrecognised Epitope Relation, kept as UNDETERMINED: ",
            paste(sprintf("'%s'", unknown), collapse = ", "))
  ep[, nature := NATURE_MAP[match(relation, names(NATURE_MAP))]]
  ep[is.na(nature), nature := "UNDETERMINED"]
  ep[, modified := nzchar(modifications) | grepl(" \\+ ", name)]

  ## 3. record-level QC ------------------------------------------------------
  log_msg("stage 3/6  sequence QC")
  ep[, drop := fifelse(object_type == "Non-peptidic", "non_peptide",
              fifelse(object_type != "Linear peptide", "non_linear", NA_character_))]
  ep[is.na(drop), sequence := parse_sequence(name)]
  ep[is.na(drop) & !nzchar(sequence), drop := "missing_sequence"]

  cand <- cbind(ep[is.na(drop)], complexity(ep[is.na(drop), sequence]))
  cand[, drop := fifelse(n_alphabet != length, "invalid_alphabet",
              fifelse(length < opt$min_length, "below_min_length",
              fifelse(length > opt$max_length, "above_max_length", NA_character_)))]
  dropped <- c(table(ep$drop), table(cand$drop))
  ## Carry only what the aggregation needs. Holding the full tables through
  ## stage 5 costs several GB and pushes the machine into swap.
  keep <- cand[is.na(drop), .(epitope_id, sequence, nature, modified,
                              source_organism, source_antigen)]
  rm(ep, cand); invisible(gc(FALSE))
  log_msg(sprintf("           %s of %s epitope records pass", fmt_n(nrow(keep)), fmt_n(n_read)))
  if (!nrow(keep)) die("no epitope records survived sequence QC; check --min-length/--max-length")

  ## 4. assay evidence -------------------------------------------------------
  log_msg("stage 4/6  reading assay evidence")
  assay <- rbindlist(lapply(c("tcell", "bcell", "mhc"), function(k) {
    a <- read_iedb(SOURCES[[k]], zips[[k]], k)[, class := k]
    if (!"mhc_allele" %in% names(a)) a[, `:=`(mhc_allele = "", mhc_class = "")]
    a
  }), use.names = TRUE)

  ## Every row in these tables is a curated experimental measurement; the export
  ## has no prediction-only assay records. "Qualitative Measurement" is the
  ## curated call, one of Positive, Positive-Low/Intermediate/High, Negative.
  assay[, positive := startsWith(qualitative, "Positive")]
  n_assay <- nrow(assay)

  ## 5. collapse to unique sequences -----------------------------------------
  log_msg("stage 5/6  aggregating evidence and collapsing duplicates")
  ev <- assay[keep[, .(epitope_id, sequence)], on = "epitope_id", nomatch = 0L]
  rm(assay); invisible(gc(FALSE))

  counts <- dcast(ev[qualitative != "", .N, by = .(sequence, class, positive)],
                  sequence ~ class + positive, value.var = "N", fill = 0L)
  for (col in c("tcell_TRUE", "tcell_FALSE", "bcell_TRUE", "bcell_FALSE",
                "mhc_TRUE", "mhc_FALSE"))
    if (!col %in% names(counts)) counts[, (col) := 0L]
  setnames(counts,
           c("tcell_TRUE", "tcell_FALSE", "bcell_TRUE", "bcell_FALSE", "mhc_TRUE", "mhc_FALSE"),
           c("n_positive_tcell", "n_negative_tcell", "n_positive_bcell",
             "n_negative_bcell", "n_positive_mhc", "n_negative_mhc"))

  refs <- unique(ev[, .(sequence, ref_iri)])[, .(n_references = .N), by = sequence]
  per_seq <- keep[, .(n_iedb_records = .N,
                      is_natural    = any(nature == "NATURAL"),
                      is_neoepitope = any(nature == "NEOEPITOPE"),
                      undetermined  = any(nature == "UNDETERMINED"),
                      all_modified  = all(modified)), by = sequence]

  db <- Reduce(function(a, b) merge(a, b, by = "sequence", all.x = TRUE),
               list(per_seq, counts, refs,
                    agg_list(keep, "epitope_id", "iedb_epitope_ids"),
                    agg_list(keep, "source_organism", "source_organisms"),
                    agg_list(keep, "source_antigen", "source_antigens"),
                    agg_list(ev, "host", "hosts"),
                    agg_list(ev, "mhc_class", "mhc_classes"),
                    agg_list(ev[class != "bcell"], "mhc_allele", "mhc_alleles")))
  for (col in names(db)[vapply(db, is.numeric, logical(1))]) db[is.na(get(col)), (col) := 0L]
  for (col in names(db)[vapply(db, is.character, logical(1))]) db[is.na(get(col)), (col) := ""]
  db <- cbind(db, complexity(db$sequence)[, .(length, unique_aa, max_aa_frac,
                                              norm_entropy, homopolymer)])
  log_msg(sprintf("           %s unique sequences from %s records",
                  fmt_n(nrow(db)), fmt_n(nrow(keep))))

  ## 6. categories and filtering ---------------------------------------------
  db[, n_positive_immune := n_positive_tcell + n_positive_bcell]
  db[, n_assays := n_positive_tcell + n_negative_tcell + n_positive_bcell +
                   n_negative_bcell + n_positive_mhc + n_negative_mhc]
  db[, evidence_category := fifelse(n_positive_immune > 0L, "IMMUNE_SUPPORTED",
                            fifelse(n_positive_mhc > 0L, "MHC_ONLY",
                             fifelse(n_assays > 0L, "NEGATIVE_ONLY", "NO_ASSAY")))]
  db[, evidence_strength := fifelse(evidence_category == "IMMUNE_SUPPORTED" & n_references >= 2L,
                                    "HIGH", fifelse(evidence_category == "IMMUNE_SUPPORTED",
                                                    "SUPPORTED", evidence_category))]
  ## A sequence attested anywhere as a natural fragment counts as natural, even
  ## if another study used the same string as an analogue backbone.
  db[, epitope_nature := fifelse(is_natural, "NATURAL",
                         fifelse(is_neoepitope, "NEOEPITOPE",
                          fifelse(undetermined, "UNDETERMINED", "ANALOG_OR_MIMOTOPE")))]
  db[, complexity_status := fifelse(unique_aa >= MIN_UNIQUE_AA &
                                    max_aa_frac <= MAX_AA_FRACTION &
                                    homopolymer <= MAX_HOMOPOLYMER &
                                    norm_entropy >= MIN_NORM_ENTROPY,
                                    "PASS", "LOW_COMPLEXITY")]

  ## Evaluated in order, so the reason names the first rule that rejected the
  ## peptide. Only the counts are reported; the TSV holds retained peptides only.
  allowed <- c("IMMUNE_SUPPORTED", if (opt$include_mhc_only) "MHC_ONLY")
  db[, why := NA_character_]
  db[is.na(why) & !evidence_category %in% allowed, why := tolower(evidence_category)]
  db[is.na(why) & epitope_nature == "ANALOG_OR_MIMOTOPE", why := "analog_or_mimotope"]
  if (opt$natural_only) db[is.na(why) & epitope_nature != "NATURAL", why := "not_natural"]
  db[is.na(why) & all_modified, why := "modified_only"]
  if (!opt$include_low_complexity)
    db[is.na(why) & complexity_status == "LOW_COMPLEXITY", why := "low_complexity"]

  setorder(db, sequence)            # deterministic ids, independent of input order
  inc <- db[is.na(why)]
  if (!nrow(inc)) die("no peptides passed the filters; loosen the options and retry")
  inc[, antigensoup_id := sprintf("AS_%07d", seq_len(.N))]

  ## outputs -----------------------------------------------------------------
  log_msg("stage 6/6  writing outputs")
  stem <- file.path(opt$outdir, paste0("iedb_antigensoup_", release))
  fasta <- paste0(stem, ".fasta"); tsv <- paste0(stem, ".tsv")
  report <- paste0(stem, ".report.txt")

  writeLines(paste0(">", inc$antigensoup_id, "\n", inc$sequence), fasta)
  fwrite(inc[, .(antigensoup_id, sequence, length, iedb_epitope_ids, n_iedb_records,
                 n_positive_tcell, n_negative_tcell, n_positive_bcell, n_negative_bcell,
                 n_positive_mhc, n_references, source_organisms, source_antigens,
                 hosts, mhc_classes, mhc_alleles, evidence_category, evidence_strength,
                 epitope_nature, complexity_status)],
         tsv, sep = "\t", quote = FALSE, na = "")
  write_report(report, opt, release, built_at, src, n_read, dropped, keep,
               db, inc, n_assay, fasta, tsv)

  cat(readLines(report), sep = "\n")
  log_msg("done. FASTA: ", fasta)
}


## ---------------------------------------------------------------------------
## Report
## ---------------------------------------------------------------------------

write_report <- function(path, opt, release, built_at, src, n_read, dropped, keep,
                         db, inc, n_assay, fasta, tsv) {
  L <- character(0)
  add <- function(...) L <<- c(L, paste0(...))
  kv  <- function(k, v) add(sprintf("  %-32s %s", k, v))
  tab <- function(dt, col) for (i in seq_len(nrow(dt)))
    add(sprintf("    %-20s %10s", dt[[col]][i], fmt_n(dt$N[i])))

  add("AntigenSoup IEDB database build report")
  add(strrep("=", 70)); add("")
  add("PROVENANCE")
  kv("script", paste0("build_iedb_db.R ", VERSION))
  kv("built", built_at)
  kv("IEDB release", release)
  kv("R / data.table", paste0(getRversion(), " / ", utils::packageVersion("data.table")))
  kv("source URL", BASE_URL)
  add("")
  for (i in seq_len(nrow(src)))
    add(sprintf("    %-32s %9s  %s  md5 %s", src$file[i], fmt_bytes(src$bytes[i]),
                src$built[i], src$md5[i]))
  add("")
  add("SETTINGS")
  kv("length range", sprintf("%d - %d aa", opt$min_length, opt$max_length))
  kv("accepted residues", paste(STD_AA, collapse = ""))
  kv("--include-mhc-only", opt$include_mhc_only)
  kv("--include-low-complexity", opt$include_low_complexity)
  kv("--natural-only", opt$natural_only)
  kv("complexity thresholds", sprintf("unique>=%d, max_frac<=%.2f, homopolymer<=%d, entropy>=%.2f",
     MIN_UNIQUE_AA, MAX_AA_FRACTION, MAX_HOMOPOLYMER, MIN_NORM_ENTROPY))
  add("")
  add("EPITOPE RECORDS")
  kv("read from epitope_full_v3", fmt_n(n_read))
  for (r in c("non_peptide", "non_linear", "missing_sequence", "invalid_alphabet",
              "below_min_length", "above_max_length"))
    kv(paste0("removed: ", r), fmt_n(if (r %in% names(dropped)) dropped[[r]] else 0L))
  kv("retained", fmt_n(nrow(keep)))
  add("")
  add("ASSAY EVIDENCE")
  kv("assay rows read (T/B/MHC)", fmt_n(n_assay))
  add("  The IEDB export contains curated experimental assays only, so there are")
  add("  no prediction-only records and no PREDICTION_ONLY category.")
  add("")
  add("UNIQUE SEQUENCES")
  kv("after collapsing duplicates", fmt_n(nrow(db)))
  add("  evidence category:");  tab(db[, .N, by = evidence_category][order(-N)], "evidence_category")
  add("  epitope nature:");     tab(db[, .N, by = epitope_nature][order(-N)], "epitope_nature")
  add("  complexity status:");  tab(db[, .N, by = complexity_status][order(-N)], "complexity_status")
  add("")
  add("FILTERING")
  ex <- db[!is.na(why), .N, by = why][order(-N)]
  for (i in seq_len(nrow(ex))) kv(paste0("removed: ", ex$why[i]), fmt_n(ex$N[i]))
  kv("RETAINED IN DATABASE", fmt_n(nrow(inc)))
  add("")
  add("DATABASE COMPOSITION")
  kv("positive T-cell evidence", fmt_n(inc[n_positive_tcell > 0L, .N]))
  kv("positive B-cell evidence", fmt_n(inc[n_positive_bcell > 0L, .N]))
  kv("both T-cell and B-cell", fmt_n(inc[n_positive_tcell > 0L & n_positive_bcell > 0L, .N]))
  kv("positive MHC evidence", fmt_n(inc[n_positive_mhc > 0L, .N]))
  kv("HIGH strength (>=2 references)", fmt_n(inc[evidence_strength == "HIGH", .N]))
  kv("neo-epitopes", fmt_n(inc[epitope_nature == "NEOEPITOPE", .N]))
  add("")
  add("  peptide length distribution:")
  ld <- inc[, .N, by = length][order(length)]
  for (i in seq_len(nrow(ld)))
    add(sprintf("    %3d aa %9s  %s", ld$length[i], fmt_n(ld$N[i]),
                strrep("#", max(0L, round(50 * ld$N[i] / max(ld$N))))))
  add("")
  add("  most frequent source organisms:")
  so <- head(inc[, .N, by = .(o = sub(";.*$", "", source_organisms))][order(-N)], 12)
  for (i in seq_len(nrow(so)))
    add(sprintf("    %-46s %9s", substr(ifelse(nzchar(so$o[i]), so$o[i], "<unspecified>"), 1, 46),
                fmt_n(so$N[i])))
  add("")
  add("OUTPUTS")
  kv("FASTA", fasta); kv("metadata TSV", tsv); kv("this report", path)
  add("")
  writeLines(L, path)
}

## Run only when executed as a script, so the file can be sourced for testing.
if (sys.nframe() == 0L && !interactive()) main(commandArgs(trailingOnly = TRUE))
