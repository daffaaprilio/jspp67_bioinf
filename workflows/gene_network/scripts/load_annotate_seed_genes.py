# workflows/gene_network/scripts/load_annotate_seed_genes.py
# 
# Load and annotate seed genes. Annotation helps identification of cluster in Cytoscape
#
# Usage: (via main snakefile) or run 
#       python3 load_annotate_seed_genes.py --graph-pickle {graph_object} --pathways map00020 map00660
# 

from pathlib import Path
import pandas as pd
import pickle
import argparse

WDIR = Path(__file__).resolve().parents[3]
HMM_RESULTS = WDIR / 'results/hmm_homology'

# ---------------------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------------------

def parse_graph_pickle_filename(graph_pickle):
    base = graph_pickle.stem
    info = base.split('-')[1]
    minZ = info.split('_')[0].replace('z', '')
    K = info.split('_')[1].replace('k', '')
    return minZ, K

def gather_top_hits(pathways, option, results_dir=HMM_RESULTS):
    """
    Scan results_dir for subdirectories matching any of the given pathways,
    load the *-07-homologous_geneID.txt file from each, and return a
    combined DataFrame of seed gene IDs.

    Parameters
    ----------
    pathways : list[str]
        Pathway prefixes to include, e.g. ['map00020', 'map00660'].
        Pass ['all'] to include every subdirectory.
    option : str
        'top_only'    : keep only the single best-evalue hit per KO
        'low_evalues' : (not yet implemented)
    results_dir : Path
        Root directory containing per-KO subdirectories.

    Returns
    -------
    pd.DataFrame with columns: gene_id, evalue, source, pathway
    """

    if option == 'low_evalues':
        raise NotImplementedError("Option 'low_evalues' will be implemented later")
    if option != 'top_only':
        raise ValueError("Options available: 'top_only' and 'low_evalues'")

    records = []

    for ko_dir in sorted(results_dir.iterdir()):
        if not ko_dir.is_dir():
            continue

        # filter by pathway prefix unless 'all' is requested
        if pathways != ['all'] and not any(ko_dir.name.startswith(p) for p in pathways):
            continue

        hit_files = list(ko_dir.glob('*-07-homologous_geneID.txt'))
        if not hit_files:
            continue

        df = pd.read_csv(hit_files[0], sep='\t')
        if df.empty:
            continue

        if option == 'top_only':
            df = df.iloc[[0]].copy()

        # derive pathway tag from directory name  (e.g. "map00020-K00025_kegg" → "map00020")
        pathway_tag = ko_dir.name.split('-')[0]
        source = ko_dir.name.replace('_kegg', '')

        df['source'] = source
        df['pathway'] = pathway_tag

        records.append(df[['gene_id', 'evalue', 'source', 'pathway']])

    if not records:
        raise FileNotFoundError(
            f"No homologous_geneID files found in {results_dir} for pathways: {pathways}"
        )

    return pd.concat(records, ignore_index=True)


def annotate_graph(G, seed_genes_df):
    """
    Add a 'pathway' node attribute to every node in G.

    Nodes present in seed_genes_df are tagged with their pathway
    (e.g. 'map00020', 'map00660'); all other nodes are tagged 'background'.

    Parameters
    ----------
    G : nx.Graph
        Co-expression graph built by build_graph.py.
    seed_genes_df : pd.DataFrame
        Output of gather_top_hits(); must contain 'gene_id' and 'pathway' columns.

    Returns
    -------
    G : nx.Graph  (modified in-place, also returned for convenience)
    """

    # default all nodes to 'background'
    nx_attrs = {node: {'pathway': 'background'} for node in G.nodes()}

    # overwrite with seed pathway labels
    # if a gene appears in multiple pathways, last write wins — keep both via comma-join
    gene_to_pathways = (
        seed_genes_df.groupby('gene_id')['pathway']
        .apply(lambda x: ','.join(sorted(set(x))))
        .to_dict()
    )

    for gene_id, pathway_label in gene_to_pathways.items():
        node = str(gene_id)   # graph nodes are strings (from coexdir_to_edgeslist)
        if node in G:
            nx_attrs[node]['pathway'] = pathway_label

    import networkx as nx
    nx.set_node_attributes(G, nx_attrs)

    n_annotated = sum(1 for n in G.nodes() if G.nodes[n]['pathway'] != 'background')
    print(f"Annotated {n_annotated} / {G.number_of_nodes()} nodes as seed genes.")

    return G


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(graph_pickle, pathways, output_filename, option):
    # load graph
    with open(graph_pickle, 'rb') as f:
        G = pickle.load(f)
    print(f"Loaded graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # gather seed genes
    seed_genes_df = gather_top_hits(pathways, option)
    print(f"Seed genes collected: {len(seed_genes_df)} rows across pathways {pathways}")

    # annotate
    G = annotate_graph(G, seed_genes_df)

    # save annotated graph
    output_filename.parent.mkdir(parents=True, exist_ok=True)
    with open(output_filename, 'wb') as f:
        pickle.dump(G, f)
    print(f"Annotated graph saved to: {output_filename}")

    # also save seed_genes_df as TSV for inspection
    tsv_out = output_filename.with_suffix('.seed_genes.tsv')
    seed_genes_df.to_csv(tsv_out, sep='\t', index=False)
    print(f"Seed genes table saved to: {tsv_out}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Load, annotate, and save seed genes on the co-expression graph')
    parser.add_argument(
        '--graph-pickle', '-g',
        type=Path,
        required=True,
        help='Path to the NetworkX graph pickle produced by build_graph.py'
    )
    parser.add_argument(
        '--pathways', '-p',
        nargs='+',
        default=['map00020', 'map00660'],
        help='Pathway prefix(es) to include as seed genes (e.g. map00020 map00660); use "all" for everything'
    )
    parser.add_argument(
        '--option', '-x',
        type=str,
        default='top_only',
        choices=['top_only', 'low_evalues'],
        help='How many hits to use per KO: "top_only" (best hit) or "low_evalues" (all significant hits)'
    )
    parser.add_argument(
        '--output', '-o',
        type=Path,
        default=None,
        help='Output path for the annotated graph pickle (.pkl)'
    )
    args = parser.parse_args()

    if args.output is None:
        minZ, K = parse_graph_pickle_filename(args.graph_pickle)
        args.output = args.graph_pickle.parent / f"sbi_G_annotated-z{minZ}_k{K}.pkl"

    main(args.graph_pickle, args.pathways, args.output, args.option)