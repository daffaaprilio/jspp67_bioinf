# workflows/gene_network/scripts/convert_id.py
#
# Convert hmmsearch tblout protein accessions to Sorghum gene IDs
# via NCBI gene2accession file.
#
# Usage:
#   python3 convert_id.py -m <hmm_result.tbl> -r <gene2accession> -o <output.txt>

import pandas as pd
from pathlib import Path
import argparse


def parse_hmmsearch_tblout(tblout_path):
    """
    Parse hmmsearch --tblout output file.
    Returns list of protein accessions (target names).
    
    tblout format: first 3 lines are comments (#), last 10 lines are footer.
    Column 0 is the target (protein) accession.
    """
    hit_accessions = []
    hit_evalue = []
    with open(tblout_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            hit_accessions.append(line.split()[0])
            hit_evalue.append(line.split()[4])
    return hit_accessions, hit_evalue


def load_gene2accession(gene2acc_path):
    """
    Load NCBI gene2accession file and create protein_accession -> GeneID mapping.
    Handles both versioned (XP_123.1) and unversioned (XP_123) accessions.
    """
    gene2acc_cols = [
        'tax_id', 'GeneID', 'status',
        'RNA_nucleotide_accession', 'RNA_nucleotide_gi',
        'protein_accession', 'protein_gi',
        'genomic_nucleotide_accession', 'genomic_nucleotide_gi',
        'start_position', 'end_position', 'orientation',
        'assembly', 'mature_peptide_accession', 'mature_peptide_gi', 'Symbol'
    ]
    
    print(f"Reading gene2accession file from {gene2acc_path}")
    gene2acc_df = pd.read_csv(
        gene2acc_path, 
        sep='\t', 
        comment='#',
        names=gene2acc_cols,
        low_memory=False
    )
    
    # Build protein_accession -> GeneID lookup
    protein_to_gene = {}
    for _, row in gene2acc_df.iterrows():
        protein_id = row['protein_accession']
        gene_id = row['GeneID']
        if pd.notna(protein_id) and protein_id != '-' and pd.notna(gene_id) and gene_id != '-':
            protein_to_gene[protein_id] = gene_id
            # Also store version-stripped key for matching
            protein_to_gene[protein_id.split('.')[0]] = gene_id
    
    return protein_to_gene


def convert_protein_to_gene_ids(hit_accessions, protein_to_gene):
    """
    Convert list of protein accessions to gene IDs using the mapping.
    Returns list of gene IDs and list of unmapped accessions.
    """
    gene_ids = []
    unmapped = []
    
    for protein_id in hit_accessions:
        # Try exact match first, then version-stripped
        gene_id = protein_to_gene.get(protein_id) or protein_to_gene.get(protein_id.split('.')[0])
        if gene_id:
            gene_ids.append(gene_id)
        else:
            unmapped.append(protein_id)
    
    return gene_ids, unmapped


def save_gene_ids(gene_ids, output_path):
    """
    Save in a tab separated format:
    col1: unique gene IDs. col2: 
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        for gene_id in sorted(set(gene_ids)):
            f.write(f"{gene_id}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Convert hmmsearch protein accessions to Sorghum gene IDs'
    )
    parser.add_argument(
        '-m', '--hmm-result',
        dest='hmm_result',
        required=True,
        help='Path to hmmsearch --tblout output file'
    )
    parser.add_argument(
        '-r', '--reference',
        dest='reference',
        required=True,
        help='Path to NCBI gene2accession file (filtered for Sorghum)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Path to output file for gene IDs'
    )
    
    args = parser.parse_args()
    
    # Step 1: Parse hmmsearch results
    print(f"Parsing hmmsearch results from {args.hmm_result}")
    hit_accessions = parse_hmmsearch_tblout(args.hmm_result)
    print(f"Found {len(hit_accessions)} hmmsearch hit(s)")
    
    if not hit_accessions:
        print("No hits found. Creating empty output file.")
        Path(args.output).write_text("")
        return
    
    # Step 2: Load gene2accession mapping
    protein_to_gene = load_gene2accession(args.reference)
    print(f"Loaded {len(protein_to_gene)} protein-to-gene mappings")
    
    # Step 3: Convert protein IDs to gene IDs
    print("Converting protein IDs to gene IDs")
    gene_ids, unmapped = convert_protein_to_gene_ids(hit_accessions, protein_to_gene)
    
    print(f"\nConverted {len(gene_ids)}/{len(hit_accessions)} protein IDs to gene IDs")
    if unmapped:
        print(f"Unmapped protein IDs: {len(unmapped)}")
        print(f"Examples of unmapped: {unmapped[:10]}")
    
    # Step 4: Save output
    save_gene_ids(gene_ids, args.output)
    print(f"Wrote {len(set(gene_ids))} unique gene IDs to {args.output}")


if __name__ == '__main__':
    main()

