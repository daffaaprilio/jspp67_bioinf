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
    rows = []
    with open(tblout_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            cols = line.split()
            rows.append({
                'protein_accession': cols[0],
                'evalue': cols[4],
                'description': ' '.join(cols[18:])
            })
    
    return pd.DataFrame(rows, columns=['protein_accession', 'evalue', 'description'])

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


def convert_protein_to_gene_ids(hits_df, protein_to_gene):
    """
    Convert list of protein accessions to gene IDs using the mapping.
    Returns list of gene IDs and list of unmapped accessions.
    Convert dataframe of hits (columns: protein accession, hit evalue, hit description) 
    into a more complete dataframe (namely gene IDs, protein accession, hit evalue, hit description)
    """
    # map the protein_accession column
    hits_df['gene_id'] = hits_df['protein_accession'].map(protein_to_gene)
    # reorder columns
    hits_df = hits_df[['protein_accession', 'gene_id', 'evalue', 'description']]
    mapped   = hits_df[hits_df['gene_id'].notna()].copy()
    unmapped = hits_df[hits_df['gene_id'].isna()].copy()
    
    return mapped, unmapped


def save_gene_ids(mapped, output_path):
    """
    Save the converted hits dataframe in a tab separated format
    """
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(output_path, sep='\t', header=True, index=False)

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
    hits_df = parse_hmmsearch_tblout(args.hmm_result)
    print(f"Found {len(hits_df)} hmmsearch hit(s)")
    
    if hits_df.empty:
        print("No hits found. Creating empty output file.")
        Path(args.output).write_text("")
        return
    
    # Step 2: Load gene2accession mapping
    protein_to_gene = load_gene2accession(args.reference)
    print(f"Loaded {len(protein_to_gene)} protein-to-gene mappings")
    
    # Step 3: Convert protein IDs to gene IDs
    print("Converting protein IDs to gene IDs")
    mapped, unmapped = convert_protein_to_gene_ids(hits_df, protein_to_gene)
    
    print(f"\nConverted {len(mapped)}/{len(hits_df)} protein IDs to gene IDs")
    if unmapped.empty:
        print(f"Unmapped protein IDs: {len(unmapped)}")
        print(f"Examples of unmapped: {unmapped[:10]}")
    
    # Step 4: Save output
    save_gene_ids(mapped, args.output)
    print(f"Wrote {len(set(mapped['gene_id']))} unique gene IDs to {args.output}")


if __name__ == '__main__':
    main()

