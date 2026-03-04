#!/usr/bin/env python3
"""
bedmethyl_to_dss.py
───────────────────
Convert a modkit pileup bedMethyl file to the four-column DSS input format.

modkit pileup bedMethyl columns (0-indexed):
  0  chrom
  1  start (0-based)
  2  end
  3  modified_base_code   (e.g. "m" = 5mC, "a" = 6mA)
  4  score (unused)
  5  strand
  6  start (copy)
  7  end   (copy)
  8  color (unused)
  9  N_valid_cov          ← total coverage
  10 percent_modified     (unused; we recompute from counts)
  11 N_mod                ← methylated count
  12 N_canonical
  13 N_other_mod
  14 N_del
  15 N_fail
  16 N_diff
  17 N_nocall

DSS input columns (1-based position, no header required but we include one):
  chr   pos   N (coverage)   X (methylated count)

Only CpG-context rows with N_valid_cov > 0 are written.
Strand is collapsed: both + and − strands at the same position are summed.

Usage:
  python bedmethyl_to_dss.py input.bedMethyl output.dss.tsv
"""

import sys
import csv
from collections import defaultdict


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    in_path  = sys.argv[1]
    out_path = sys.argv[2]

    # Accumulate coverage + methylated counts per (chrom, 1-based pos).
    # Collapsing strands by summing both.
    counts: dict[tuple[str, int], list[int]] = defaultdict(lambda: [0, 0])  # [N, X]

    with open(in_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue

            chrom        = fields[0]
            start        = int(fields[1])       # 0-based
            n_valid_cov  = int(fields[9])
            n_mod        = int(fields[11])

            if n_valid_cov == 0:
                continue

            pos = start + 1                     # convert to 1-based
            counts[(chrom, pos)][0] += n_valid_cov
            counts[(chrom, pos)][1] += n_mod

    # Sort by chrom then position and write output
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["chr", "pos", "N", "X"])
        for (chrom, pos), (N, X) in sorted(counts.items(), key=lambda kv: (kv[0][0], kv[0][1])):
            writer.writerow([chrom, pos, N, X])

    print(f"[bedmethyl_to_dss] Written {len(counts):,} CpG sites → {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
