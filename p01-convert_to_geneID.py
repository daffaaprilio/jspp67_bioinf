import pandas as pd
from pathlib import Path

gene2acc_file = '/home/daffa/Work/2025/11-JSPP67/sorghum_gene2accession'

HIT_OUTPUT = '/home/daffa/Work/2025/11-JSPP67/hit_genes_for_clustering/hits.txt'
HIT_ENTREZ_OUTPUT = '/home/daffa/Work/2025/11-JSPP67/hit_genes_for_clustering/hits_EGI.txt'

# parse hit genes from HMM
hmm_hit_prpF = Path('/home/daffa/Work/2025/11-JSPP67/hmm/hmm_pipeline/prpF/prpF_homologs_in_plant-results.tbl')
with open(hmm_hit_prpF, 'r') as f:
    content = f.readlines()[3:-10]
    hit = [line.split()[0] for line in content]

# convert hit protein ID into findable gene ID in Sorghum coex data
print("Reading gene2accession file for Sorghum data")
gene2acc_df = pd.read_csv(gene2acc_file, sep='\t', comment='#', 
                 names=['tax_id', 'GeneID', 'status', 'RNA_nucleotide_accession', 'RNA_nucleotide_gi',
                        'protein_accession', 'protein_gi', 'genomic_nucleotide_accession',
                        'genomic_nucleotide_gi', 'start_position', 'end_position', 'orientation',
                        'assembly', 'mature_peptide_accession', 'mature_peptide_gi', 'Symbol'])

print("Creating mapping from protein to gene symbol")
protein_to_gene = {}
for _, row in gene2acc_df.iterrows():
    protein_id = row['protein_accession']
    gene_id = row['GeneID']
    if pd.notna(protein_id) and protein_id != '-' and pd.notna(gene_id) and gene_id != '-':
        protein_to_gene[protein_id] = gene_id

print("Converting protein ID to gene ID")
gene_ids = []
unmapped = []
for protein_id in hit:
    if protein_id in protein_to_gene:
        gene_ids.append(protein_to_gene[protein_id])
    else:
        unmapped.append(protein_id)

print(f"\nConverted {len([g for g in gene_ids if g])} out of {len(hit)} protein IDs")
if unmapped:
    print(f"Unmapped protein IDs: {len(unmapped)}")
    print("Example of unmapped protein IDs:", unmapped[:10])

# saving as tsv files
hit_output = HIT_OUTPUT
hit_entrez_output = HIT_ENTREZ_OUTPUT

def save_list(items, outfile):
    with open(outfile, 'w') as f:
        for i in items:
            f.write(f"{i}\n")

save_list(hit, hit_output)
save_list(gene_ids, hit_entrez_output)