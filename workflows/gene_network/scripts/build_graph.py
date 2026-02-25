# workflows/gene_network/scripts/build_graph.py
# 
# Build an NetworkX graph object and save it in pickle format
# to allow exploration (using py4cytoscape) in a 
# Jupyter Network interface
#
# Usage: (via main snakefile) or run 
#       python build_graph.py
# 

from pathlib import Path
import pandas as pd
import networkx as nx
import pickle
import argparse

WDIR = Path(__file__).resolve().parents[3]
sbi_coex_dir = WDIR / 'data/reference/sbi_coex'


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def coexdir_to_edgeslist(coex_dir, K, minZ):
    '''
    create a function that opens an ATTED-II style gene coexpression directory,
    read all the contents temporarily through a python dataframe,
    filter only necessary edges (important to slim down the edge-list size)
    store the final content as an edge-list
    u-v has z weight: gene u and gene v is co-expressed with a value of z (z-score)
    '''
    edges = []

    for file in coex_dir.glob('*'):
        u = str(file.stem)
        # temporary df, rewritten every iteration of for loop
        df = pd.read_csv(file, sep='\t', header=None, names=['v', 'z'])
        df['v'] = df['v'].astype(str)

        # from each file, filter top K coexpressed genes & with a certain min Z score
        df = df.nlargest(K, 'z')
        df = df[df['z'] >= minZ]    

        # store to a list
        for v, z in zip(df['v'], df['z']):
            edges.append((u, v, float(z)))
    
    return edges

def build_graph(edges):
    '''
    initiate an empty undirected graph, then load the edges list
    use the correct method (add_weighted_edges_from)
    '''
    G = nx.Graph()
    # G.add_weighted_edges_from(edges)
    for u, v, z, in edges:
        if G.has_edge(u, v):
            if z > G[u][v]['weight']:
                G[u][v]['weight'] = z
        else:
            G.add_edge(u, v, weight=z)

    return G

def save_object(pyobj, pklobj):
    with open(pklobj, 'wb') as f:
        pickle.dump(pyobj, f)

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(coex_dir, output_filename, K, minZ):
    edges = coexdir_to_edgeslist(coex_dir, K, minZ)
    G = build_graph(edges)

    save_object(G, output_filename)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create coexpression network')
    parser.add_argument(
        '--gene-no', 
        '-k', 
        type=int, 
        default=10,
        help='From each file in the coex dir, filter up to K number of genes to be included in the network'
    )
    parser.add_argument(
        '--z-score', 
        '-z', 
        type=float, 
        default=4,
        help='From each file in the coex dir, filter only genes coexpressed with min z score'
    )
    parser.add_argument(
        '--coex-dir', 
        '-c', 
        type=Path, 
        default=sbi_coex_dir,
        help='Path to ATTED-II style coexpression directory'
    )
    parser.add_argument(
        '--output', 
        '-o', 
        type=Path, 
        default=None,
        help='Output path for the graph pickle (.pkl)'
    )
    args = parser.parse_args()

    K = args.gene_no
    minZ = args.z_score
    output = args.output or WDIR / f'results/gene_network/sbi_G-z{minZ}_k{K}.pkl'

    main(args.coex_dir, output, K, minZ)