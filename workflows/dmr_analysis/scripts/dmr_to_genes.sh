#!/usr/bin/env bash
#
# dmr_to_genes.sh — Intersect DMR results with a GFF3 annotation to produce
#                   a list of overlapping gene IDs.
#
# Usage:
#   dmr_to_genes.sh -d DMR_TSV -g GFF3 -o OUTPUT [-r REGION]
# /home/daffa/Work/2026/02-JSPP67/workflows/dmr_analysis/scripts/dmr_to_genes.sh -d dmr_results_15_filtered.tsv -o dmr_genes_gene_15_filtered.tsv -g /home/daffa/Work/2026/02-JSPP67/data/reference/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff
# /home/daffa/Work/2026/02-JSPP67/workflows/dmr_analysis/scripts/dmr_to_genes.sh -d dmr_results_15_filtered.tsv -r "promoter_1000" -o dmr_genes_promoter_1000_15_filtered.tsv -g /home/daffa/Work/2026/02-JSPP67/data/reference/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff
# /home/daffa/Work/2026/02-JSPP67/workflows/dmr_analysis/scripts/dmr_to_genes.sh -d dmr_results_15_filtered.tsv -r "promoter_2000" -o dmr_genes_promoter_2000_15_filtered.tsv -g /home/daffa/Work/2026/02-JSPP67/data/reference/GCF_000003195.3_Sorghum_bicolor_NCBIv3_genomic.gff
#
# Arguments:
#   -d  DMR TSV file (output of dss_dmr.R; columns: chr, start, end, ...)
#   -g  GFF3 annotation file
#   -o  Output file path (TSV of unique GeneIDs)
#   -r  Region type: "gene" (default) or "promoter_N" (e.g. "promoter_2000")
#
# Dependencies: awk, grep, sort, bedtools

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
REGION="gene"

# ── Argument parsing ──────────────────────────────────────────────────────────
usage() {
    sed -n '3,12p' "$0" | sed 's/^# \?//'
    exit 1
}

while getopts ":d:g:o:r:h" opt; do
    case $opt in
        d) DMR_TSV="$OPTARG" ;;
        g) GFF3="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        r) REGION="$OPTARG" ;;
        h) usage ;;
        :) echo "ERROR: option -$OPTARG requires an argument" >&2; usage ;;
        \?) echo "ERROR: unknown option -$OPTARG" >&2; usage ;;
    esac
done

[[ -z "${DMR_TSV:-}" ]] && { echo "ERROR: -d DMR_TSV is required" >&2; usage; }
[[ -z "${GFF3:-}"    ]] && { echo "ERROR: -g GFF3 is required"    >&2; usage; }
[[ -z "${OUTPUT:-}"  ]] && { echo "ERROR: -o OUTPUT is required"   >&2; usage; }
[[ -f "$DMR_TSV"     ]] || { echo "ERROR: DMR file not found: $DMR_TSV" >&2; exit 1; }
[[ -f "$GFF3"        ]] || { echo "ERROR: GFF3 file not found: $GFF3"   >&2; exit 1; }

# ── Derive promoter parameters from REGION ───────────────────────────────────
IS_PROMOTER=0
UPSTREAM=0
if [[ "$REGION" == promoter_* ]]; then
    IS_PROMOTER=1
    UPSTREAM="${REGION#promoter_}"
    if ! [[ "$UPSTREAM" =~ ^[0-9]+$ ]]; then
        echo "ERROR: promoter region must be 'promoter_N' with integer N (got: $REGION)" >&2
        exit 1
    fi
fi

# ── Run ───────────────────────────────────────────────────────────────────────
mkdir -p "$(dirname "$OUTPUT")"
TMP_BED="${OUTPUT}.dmr.bed"

# Convert DSS DMR TSV (chr, start, end, …) to BED3, skip header
tail -n +2 "$DMR_TSV" \
    | awk 'NF>=3 {print $1"\t"$2"\t"$3}' \
    > "$TMP_BED"

# Extract gene features; for promoter_N rewrite coords to the upstream window
grep $'\tgene\t' "$GFF3" \
    | awk -v is_prom="$IS_PROMOTER" -v up="$UPSTREAM" \
        'BEGIN{FS="\t"; OFS="\t"} {
            if (is_prom) {
                gs = $4; ge = $5; str = $7
                if (str == "+") { $4 = (gs - up > 1 ? gs - up : 1); $5 = gs - 1 }
                else             { $4 = ge + 1; $5 = ge + up }
            }
            print
        }' \
    | bedtools intersect -a stdin -b "$TMP_BED" -u \
    | awk '{
        gene_id = ""
        n = split($9, attrs, ";")
        for (i=1; i<=n; i++) {
            if (attrs[i] ~ /^Dbxref=GeneID:/) {
                sub(/^Dbxref=GeneID:/, "", attrs[i]); gene_id = attrs[i]
            }
        }
        if (gene_id != "") print gene_id
    }' \
    | sort -u \
    > "$OUTPUT"

rm -f "$TMP_BED"
echo "Done: $(wc -l < "$OUTPUT") genes written to $OUTPUT"