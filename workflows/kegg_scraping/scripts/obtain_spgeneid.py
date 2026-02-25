# workflows/kegg_scraping/scripts/obtain_spgeneid.py
#
# Extract sp:geneid pairs from a per-KO pathway gene TSV
# (output of separate_pathway_gene_files.py).
#
# Usage:
#   python obtain_spgeneid.py --input {pathway_genes_tsv} --output {spgeneid_txt}

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Extract sp:geneid pairs from a per-KO pathway gene TSV."
    )
    parser.add_argument("--input",  required=True, help="Per-KO pathway gene TSV (01-pathway_genes_from_KEGG.txt)")
    parser.add_argument("--output", required=True, help="Output txt file, one sp:geneid per line")
    args = parser.parse_args()

    input_path  = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    # TSV columns: species | gene | ncbi_geneid | ko | ncbi_proteinid
    # 'gene' (col 1) is the KEGG gene ID, which matches pep file headers (>sp:geneid)
    out = []
    with open(input_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            sp   = parts[0]          # e.g. sbi
            gid  = parts[1]          # e.g. AT1G31180 or 8055458 (KEGG gene ID)
            out.append(f"{sp}:{gid}")

    with open(output_path, 'w') as w:
        for entry in out:
            w.write(entry + '\n')

    print(f"Wrote {len(out)} sp:geneid entries to {output_path}")


if __name__ == "__main__":
    main()
