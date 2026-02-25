# workflows/kegg_scraping/scripts/kegg_pathway_genes.py
#
# Extract all species-specific genes belonging to a given KEGG reference
# pathway from locally stored KEGG FTP flat files.
#
# Usage:
#   python kegg_pathway_genes.py --pathway map00660 --species sbi \
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


def get_pathway_genes(pathway, species, kegg_dir):
    """
    Return a DataFrame of species genes in *pathway* with annotations.

    Parameters
    ----------
    pathway : str
        Reference pathway id, e.g. ``map00660``.
    species : str
        KEGG three-letter species code, e.g. ``sbi``.
    kegg_dir : Path
        Root directory of the local KEGG FTP mirror
        (expected sub-dirs: ``<species>/``, ``ko/``, ``map/``).

    Returns
    -------
    pd.DataFrame
        Columns: gene, ko, ncbi_geneid, ncbi_proteinid
    """
    org_dir = kegg_dir / species
    pathway_id = pathway.replace("map", species)  # map00660 -> sbi00660 if species is not "all". Otherwise, retain.
    
    print(f"Preparing species name: {species}")

    # 1. genes in this pathway
    print(f"[1/4] Parsing {species}_pathway.list ...")
    sp_pathway = parse_kegg_list(
        org_dir / f"{species}_pathway.list", ("gene", "pathway")
    )
    genes = sp_pathway.loc[sp_pathway["pathway"] == pathway_id, ["gene"]].copy()
    print(f"      Found {len(genes)} {species} genes in pathway {pathway_id}")

    # 2. genes in EGI format (ncbi gene id)
    print(f"[2/4] Parsing {species}_ncbi-geneid.list ...")
    sp_geneid = parse_kegg_list(
        org_dir / f"{species}_ncbi-geneid.list", ("gene", "ncbi_geneid")
    )
    genes = genes.merge(sp_geneid, on="gene", how="left")

    # 3. KO ortholog assignments
    print(f"[3/4] Parsing {species}_ko.list ...")
    sp_ko = parse_kegg_list(org_dir / f"{species}_ko.list", ("gene", "ko"))
    genes = genes.merge(sp_ko, on="gene", how="left")

    # 4. NCBI Protein accessions
    print(f"[4/4] Parsing {species}_ncbi-proteinid.list ...")
    sp_protid = parse_kegg_list(
        org_dir / f"{species}_ncbi-proteinid.list", ("gene", "ncbi_proteinid")
    )
    genes = genes.merge(sp_protid, on="gene", how="left")

    return genes

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Extract species genes from a KEGG reference pathway."
    )
    parser.add_argument(
        "--species",
        required=True,
        nargs="+",
        help=(
            "Determine homology only from provided species (KEGG three-letter "
            "species code (e.g. sbi)), rather than from all species listed in "
            'KEGG. Insert "atted-plants" if you want to insert ATTED-II plants. '
            'Insert "all" if you want to find homology from all species '
            "listed in the database (coming soon)."
        )
    )
    parser.add_argument(
        "--pathway",
        required=True,
        help="Reference pathway id (e.g. map00660)",
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
    species = args.species
    
    if "all" in species:
        raise NotImplementedError(
            'Processing for "all" species is not implemented.\n'
            'Provide one or more specific KEGG organism codes instead.'
        )
    
    if "atted-plants" in species:
        species = (
            "ath", "bdi", "bna", "brp", "cit", "cre", "ghi", "gmx",
            "mtr", "nta", "osa", "pop", "sbi", "sly", "sot", "taes",
            "vvi", "zma"
        )

    dfs = []
    for sp in species:
        df = get_pathway_genes(args.pathway, sp, kegg_dir)
        df.insert(0, "species", sp)
        dfs.append(df)
    genes = pd.concat(dfs, ignore_index=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    genes.to_csv(output, sep="\t", index=False)
    print(f"\nSaved {len(genes)} genes → {output}")


if __name__ == "__main__":
    main()