#!/usr/bin/env python3
import argparse
import itertools
import sys

STD_AAS = "ACDEFGHIKLMNPQRSTVWY"


def clean_pep(s: str) -> str:
    s = s.strip().upper()
    allowed = set(STD_AAS)
    return "".join(c for c in s if c in allowed)


def hamming_variants(seq: str, mismatches: int):
    """
    Yield all unique sequences within Hamming distance `mismatches` of `seq`.
    Only substitutions (no indels). Yields the original sequence too (d=0).
    """
    L = len(seq)
    aas = STD_AAS

    yield seq  # d=0: original

    for d in range(1, mismatches + 1):
        for positions in itertools.combinations(range(L), d):
            # For each combination of positions, try all AA substitutions
            substitution_choices = []
            for pos in positions:
                substitution_choices.append([aa for aa in aas if aa != seq[pos]])

            for subs in itertools.product(*substitution_choices):
                variant = list(seq)
                for pos, aa in zip(positions, subs):
                    variant[pos] = aa
                yield "".join(variant)


def parse_fasta(path: str):
    """Yield (header, seq) pairs from a FASTA file using only base Python."""
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def main():
    ap = argparse.ArgumentParser(
        description="Generate Hamming-distance variants of epitope sequences for exact matching."
    )
    ap.add_argument("epitope", help="Epitope sequence (e.g. SIINFEKL) or path to a FASTA file")
    ap.add_argument(
        "--mismatches", "-d", type=int, default=1,
        help="Maximum number of mismatches / substitutions (default: 1)"
    )
    ap.add_argument(
        "--out", "-o", default=None,
        help="Output FASTA path (default: <epitope>_variants_d<d>.fasta)"
    )
    args = ap.parse_args()

    mismatches = args.mismatches
    if mismatches < 0:
        ap.error("--mismatches must be >= 0")

    # Determine if input is a sequence or a FASTA file
    import os
    if os.path.isfile(args.epitope):
        records = [(hdr, clean_pep(seq)) for hdr, seq in parse_fasta(args.epitope)]
        records = [(hdr, seq) for hdr, seq in records if seq]
        if not records:
            print("[gen_variants] No valid sequences found in FASTA.", file=sys.stderr)
            sys.exit(1)
        # Auto-name output based on filename stem
        stem = os.path.splitext(os.path.basename(args.epitope))[0]
        default_out = f"{stem}_variants_d{mismatches}.fasta"
    else:
        seq = clean_pep(args.epitope)
        if not seq:
            ap.error(f"Could not parse a valid peptide sequence from: {args.epitope!r}")
        records = [(seq, seq)]  # use sequence itself as the "name"
        default_out = f"{seq}_variants_d{mismatches}.fasta"

    out_path = args.out if args.out else default_out

    n_records = 0
    n_variants = 0
    seen = set()  # deduplicate across epitopes

    with open(out_path, "w") as out:
        for epitope_name, seq in records:
            variants_this = 0
            for variant in hamming_variants(seq, mismatches):
                if variant in seen:
                    continue
                seen.add(variant)

                d = sum(a != b for a, b in zip(seq, variant))
                if d == 0:
                    label = "original"
                else:
                    # Encode which positions changed: e.g. A3V_K5R
                    changes = []
                    for i, (orig_aa, new_aa) in enumerate(zip(seq, variant)):
                        if orig_aa != new_aa:
                            changes.append(f"{orig_aa}{i+1}{new_aa}")
                    label = "_".join(changes)

                out.write(f">{epitope_name}|d{d}|{label}|{seq}\n{variant}\n")
                variants_this += 1

            n_variants += variants_this
            n_records += 1

    print(
        f"[gen_variants] epitopes={n_records} mismatches={mismatches} "
        f"variants={n_variants} wrote={out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
