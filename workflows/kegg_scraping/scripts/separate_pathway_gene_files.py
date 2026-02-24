# workflow/kegg_scraping/scripts/separate_pathway_gene_files.py
#
# split the output of kegg_pathway_genes.py into multiple files of the name 
# {GENEID}_KEGG/{GENEID}_KEGG-01-pathway_genes_from_KEGG.txt
#
# Usage:
#   python {this-script.py} --input {pathway_gene_file}
#

import argparse
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def parse_input_file(file):
    """
    read pathway gene file
    output df that contains sp:gene -delimiter- accessionID
    """
    with open(file, 'r') as path_gene_f:
        lines = path_gene_f.readlines()
        lines.remove('species\tgene\tko\tncbi_geneid\tncbi_proteinid\n') # remove header the hard way
        sp_gene_col = []
        acc_col =[]
        for line in lines:
            # split tabs
            sp = line.split('\t')[0]
            geneid = line.split('\t')[1]
            kegg_ortho = line.split('\t')[2]
            prot_acc = line.split('\t')[3]
            # assign to prev defined list, then eventually dataframe
            sp_gene_col.append(':'.join((sp, geneid)))
            acc_col.append(prot_acc)

        d = {"gene": sp_gene_col, "accession":acc_col}
        df = pd.DataFrame(data=d)
        
    return df

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Separates the pathway gene file (output from the kegg_pathway_genes.py script) and outputs files based on the genes"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Pathway gene file (output from the kegg_pathway_genes.py script, i.e., sbi-zma_map00660_genes.tsv)",
    )

    args = parser.parse_args()

    input = args.input
    if Path(input).exists():
        pass
    else:
        raise FileNotFoundError('Input TSV file (kegg_pathway_genes) not found')
    
    # prefix of the Output file name must be specified for the subsequent Snakemake pipeline

if __name__ == "__main__":
    main()

