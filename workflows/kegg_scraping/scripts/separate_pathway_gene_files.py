# workflow/kegg_scraping/scripts/separate_pathway_gene_files.py
#
# split the output of kegg_pathway_genes.py into multiple files of the name 
# {pathway}-{ko}_KEGG/{pathway}-{ko}_KEGG-01-pathway_genes_from_KEGG.txt
# The idea is to feed the HMM_homology pipeline various enzyme/protein sequences that
# have one same function, which is whatever the KEGG orthology function is.
#
# Usage:
#   python {this-script.py} --input {pathway_gene_file} --output {directory_for_the_output_files}
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
    outputs df that contains sp:gene -delimiter- accessionID
    """
    df = pd.read_csv(file, sep='\t')
    return df

def separate_by_orthology(df):
    """
    Group DataFrame by ko column.
    Returns dict[str, pd.DataFrame] e.g. {"K01703": sub_df, ...}
    Rows with NaN ko are dropped.
    """
    df_clean = df.dropna(subset=['ko'])
    return {ko: sub_df.reset_index(drop=True) for ko, sub_df in df_clean.groupby('ko')}

def save_outfile(sub_df, output_path):
    """
    Save a KO sub-DataFrame as a TSV file.
    Creates parent directories if they don't exist.
    Drops the pathway column so output contains: species | gene | ncbi_geneid | ko | ncbi_proteinid
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cols = [c for c in sub_df.columns if c != 'pathway']
    sub_df[cols].to_csv(output_path, sep='\t', index=False)
    

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
    parser.add_argument(
        '--output',
        required=True,
        help="Output TSV path prefix (e.g., results/kegg/)",
    )
    parser.add_argument(
        "--pathway",
        required=True,
        help="Reference pathway id (e.g. map00660)",
    )
    args = parser.parse_args()

    pathway = args.pathway
    input = args.input
    output_prefix = args.output # prefix of the Output file name must be specified for the subsequent Snakemake 

    if Path(input).exists():
        # ensure output directory exists before iterating (required for Snakemake directory() output)
        Path(output_prefix).mkdir(parents=True, exist_ok=True)
        # main executor
        input_df = parse_input_file(input)
        for ko, sub_df in separate_by_orthology(input_df).items():
            out = Path(output_prefix) / f"{pathway}-{ko}_kegg" / f"{pathway}-{ko}_kegg-01-pathway_genes_from_KEGG.txt"
            save_outfile(sub_df, out)
            print(f"Saved {pathway}-{ko} file as {out}")
    else:
        raise FileNotFoundError(f'Input TSV file {input} not found')
    
    

if __name__ == "__main__":
    main()

