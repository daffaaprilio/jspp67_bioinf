# workflows/kegg_scraping/scripts/obtain_aa_sequences.py
#
# Extract amino acid sequences from KEGG .pep files for a list of sp:geneid pairs.
#
# Usage:
#   python obtain_aa_sequences.py --input {spgeneid_txt} --kegg-dir {kegg_dir} --output {out.faa}
#   python workflows/kegg_scraping/scripts/obtain_aa_sequences.py --input results/by_ko/map00020-K00025_kegg/map00020-K00025_kegg-01-pathway_genes_from_KEGG.txt --kegg-dir data/reference/KEGG --output results/hmm_homology/map00020-K00025_kegg/map00020-K00025_kegg-03.faa
#       

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Extract AA sequences from KEGG .pep files for a list of sp:geneid pairs."
    )
    parser.add_argument("--input",    required=True, help="sp:geneid list (one per line, 02-spgeneid.txt)")
    parser.add_argument("--kegg-dir", required=True, help="Root of local KEGG FTP mirror (contains sp/ subdirs)")
    parser.add_argument("--output",   required=True, help="Output FASTA file")
    args = parser.parse_args()

    input_path  = Path(args.input)
    kegg_dir    = Path(args.kegg_dir)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    # Read gene IDs and group by species
    gene_ids = set()
    with open(input_path) as f:
        for line in f:
            gid = line.strip()
            if gid:
                gene_ids.add(gid)   # e.g. "sbi:8055458", "ath:AT1G31180"

    species = {gid.split(':')[0] for gid in gene_ids}

    output_path.parent.mkdir(parents=True, exist_ok=True)

    found = 0
    with open(output_path, 'w') as out_f:
        for sp in sorted(species):
            pep_files = list((kegg_dir / sp).glob("*.pep"))
            if not pep_files:
                print(f"WARNING: no .pep file found for species {sp}, skipping")
                continue
            pep_file = pep_files[0]   # each species has exactly one .pep

            sp_gene_ids = {gid for gid in gene_ids if gid.startswith(sp + ':')}
            current_id  = None
            current_seq = []

            with open(pep_file) as pep_f:
                for line in pep_f:
                    if line.startswith('>'):
                        # flush previous matching record
                        if current_id and current_seq:
                            out_f.write(f">{current_id}\n")
                            out_f.write(''.join(current_seq))
                            found += 1
                        # check new header: ">sp:geneid  description"
                        header_id  = line[1:].split()[0]   # e.g. "sbi:8055458"
                        current_id = header_id if header_id in sp_gene_ids else None
                        current_seq = []
                    else:
                        if current_id:
                            current_seq.append(line)
                # flush last record
                if current_id and current_seq:
                    out_f.write(f">{current_id}\n")
                    out_f.write(''.join(current_seq))
                    found += 1

    missing = len(gene_ids) - found
    print(f"Wrote {found}/{len(gene_ids)} sequences to {output_path}" +
          (f" ({missing} IDs not found in pep files)" if missing else ""))


if __name__ == "__main__":
    main()
