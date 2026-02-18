from pathlib import Path
import os
import pandas as pd
import networkx as nx
import pickle
from datetime import datetime
import argparse

def coexdir_to_edgeslist(coex_dir):
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
    
    return None

def _main(coex_dir, output_dirname):
    edges = coexdir_to_edgeslist(coex_dir)
    G = build_graph(edges)

    output_path = Path(output_dirname)
    output_path.mkdir(parents=True, exist_ok=True)

    save_object(G, output_path / 'sbi_G_prpF.pkl')
    save_object(edges, output_path / 'sbi_edges.pkl')

sbi_coex_dir = Path('/Users/daffa/Documents/Work/jspp67_bioinf/sbi_coex')
K = 10
minZ = float(4)

if __name__ == '__main__':
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    default_output = Path(f'TAA_cluster-1st_attempt/{timestamp}')

    parser = argparse.ArgumentParser(description='Create coexpression network')
    parser.add_argument('--coex-dir', '-c', type=Path, default=sbi_coex_dir,
                        help='Path to ATTED-II style coexpression directory')
    parser.add_argument('--output-dir', '-o', type=Path, default=default_output,
                        help=f'Output directory for graph and edges pickle (default folder name is {default_output})')
    parser.add_argument('--gene-no', '-k', type=int, default=10,
                        help='From each file in the coex dir, filter up to K number of genes to be included in the network')
    parser.add_argument('--z-score', '-z', type=float, default=4,
                        help='From each file in the coex dir, filter only genes coexpressed with min z score')
    args = parser.parse_args()

    _main(args.coex_dir, args.output_dir)
    # _main(sbi_coex_dir, output_dir)