# workflows/kegg_scraping/scripts/kegg_pathway_genes.py
#
# Extract all organism-specific genes belonging to a given KEGG reference
# pathway from locally stored KEGG FTP flat files.
#
# Usage:
#   python kegg_pathway_genes.py --pathway map00660 --organism sbi \
#       --kegg-dir data/reference/KEGG --output results/kegg/sbi_map00660_genes.tsv

import argparse
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def parse_kegg_list(filepath, col_names=("source", "target")):
    """
    Parse a KEGG two-column tab-separated .list file into a DataFrame.

    Namespace prefixes (e.g. 'sbi:', 'ko:', 'path:') are stripped
    automatically so that downstream joins use bare identifiers.
    """
    df = pd.read_csv(filepath, sep="\t", header=None, names=col_names)
    for col in df.columns:
        df[col] = df[col].str.split(":", n=1).str[1]
    return df


def get_pathway_genes(pathway, organism, kegg_dir):
    """
    Return a DataFrame of organism genes in *pathway* with annotations.

    Parameters
    ----------
    pathway : str
        Reference pathway id, e.g. ``map00660``.
    organism : str
        KEGG three-letter organism code, e.g. ``sbi``.
    kegg_dir : Path
        Root directory of the local KEGG FTP mirror
        (expected sub-dirs: ``<organism>/``, ``ko/``, ``map/``).

    Returns
    -------
    pd.DataFrame
        Columns: gene, ko, ncbi_geneid, ncbi_proteinid
    """
    org_dir = kegg_dir / organism
    pathway_id = pathway.replace("map", organism)  # map00660 -> sbi00660

    # 1. genes in this pathway
    print(f"[1/4] Parsing {organism}_pathway.list …")
    sbi_pathway = parse_kegg_list(
        org_dir / f"{organism}_pathway.list", ("gene", "pathway")
    )
    genes = sbi_pathway.loc[sbi_pathway["pathway"] == pathway_id, ["gene"]].copy()
    print(f"      Found {len(genes)} {organism} genes in pathway {pathway_id}")

    # 2. KO ortholog assignments
    print(f"[2/4] Parsing {organism}_ko.list …")
    sbi_ko = parse_kegg_list(org_dir / f"{organism}_ko.list", ("gene", "ko"))
    genes = genes.merge(sbi_ko, on="gene", how="left")

    # 3. NCBI Gene IDs
    print(f"[3/4] Parsing {organism}_ncbi-geneid.list …")
    sbi_geneid = parse_kegg_list(
        org_dir / f"{organism}_ncbi-geneid.list", ("gene", "ncbi_geneid")
    )
    genes = genes.merge(sbi_geneid, on="gene", how="left")

    # 4. NCBI Protein accessions
    print(f"[4/4] Parsing {organism}_ncbi-proteinid.list …")
    sbi_protid = parse_kegg_list(
        org_dir / f"{organism}_ncbi-proteinid.list", ("gene", "ncbi_proteinid")
    )
    genes = genes.merge(sbi_protid, on="gene", how="left")

    return genes


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract organism genes from a KEGG reference pathway."
    )
    parser.add_argument(
        "--pathway",
        required=True,
        help="Reference pathway id (e.g. map00660)",
    )
    parser.add_argument(
        "--organism",
        required=True,
        help="KEGG three-letter organism code (e.g. sbi)",
    )
    parser.add_argument(
        "--kegg-dir",
        default="data/reference/KEGG",
        help="Root of the local KEGG FTP mirror (default: data/reference/KEGG)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path (e.g. results/kegg/sbi_map00660_genes.tsv)",
    )
    args = parser.parse_args()

    kegg_dir = Path(args.kegg_dir)
    output = Path(args.output)

    genes = get_pathway_genes(args.pathway, args.organism, kegg_dir)

    output.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(output, sep="\t", index=False)
    print(f"\nSaved {len(genes)} genes → {output}")


if __name__ == "__main__":
    main()
