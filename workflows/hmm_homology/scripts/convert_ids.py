# workflows/hmm_homology/scripts/convert_ids.py
#
# Snakemake script: convert hmmsearch tblout protein accessions
# to Sorghum gene IDs via NCBI gene2accession.
#
# snakemake.input.tblout        -- hmmsearch --tblout output
# snakemake.input.gene2accession -- NCBI gene2accession file (plain or .gz)
# snakemake.output[0]           -- gene ID list (one per line)

import gzip
import pandas as pd
from pathlib import Path

tblout_path        = snakemake.input.tblout
gene2accession_path = snakemake.input.gene2accession
out_path           = snakemake.output[0]

# ── 1. Parse hmmsearch tblout ────────────────────────────────────────────────
# tblout format: first 3 lines are comments (#), last 10 are a footer (#).
# Column 0 is the target (protein) accession.
hit_accessions = []
with open(tblout_path, "r") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        hit_accessions.append(line.split()[0])

print(f"[convert_ids] Found {len(hit_accessions)} hmmsearch hit(s)")

if not hit_accessions:
    Path(out_path).write_text("")
    raise SystemExit(0)

# ── 2. Read gene2accession ───────────────────────────────────────────────────
gene2acc_cols = [
    "tax_id", "GeneID", "status",
    "RNA_nucleotide_accession", "RNA_nucleotide_gi",
    "protein_accession", "protein_gi",
    "genomic_nucleotide_accession", "genomic_nucleotide_gi",
    "start_position", "end_position", "orientation",
    "assembly", "mature_peptide_accession", "mature_peptide_gi", "Symbol",
]

open_fn = gzip.open if str(gene2accession_path).endswith(".gz") else open
print(f"[convert_ids] Reading gene2accession from {gene2accession_path}")

with open_fn(gene2accession_path, "rt") as f:
    gene2acc_df = pd.read_csv(
        f, sep="\t", comment="#", header=None,
        names=gene2acc_cols, low_memory=False,
    )

# Build protein_accession → GeneID lookup (strip version suffix for matching)
protein_to_gene: dict[str, int] = {}
for _, row in gene2acc_df.iterrows():
    prot = row["protein_accession"]
    gid  = row["GeneID"]
    if pd.notna(prot) and prot != "-" and pd.notna(gid) and gid != "-":
        protein_to_gene[prot] = gid
        protein_to_gene[prot.split(".")[0]] = gid  # version-stripped key

# ── 3. Map accessions → gene IDs ────────────────────────────────────────────
gene_ids  = []
unmapped  = []
for acc in hit_accessions:
    gid = protein_to_gene.get(acc) or protein_to_gene.get(acc.split(".")[0])
    if gid:
        gene_ids.append(gid)
    else:
        unmapped.append(acc)

print(f"[convert_ids] Mapped {len(gene_ids)}/{len(hit_accessions)} accessions to gene IDs")
if unmapped:
    print(f"[convert_ids] Unmapped accessions ({len(unmapped)}): {unmapped[:10]}")

# ── 4. Write output ──────────────────────────────────────────────────────────
Path(out_path).parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as out:
    for gid in sorted(set(gene_ids)):
        out.write(f"{gid}\n")

print(f"[convert_ids] Wrote {len(set(gene_ids))} unique gene IDs to {out_path}")
