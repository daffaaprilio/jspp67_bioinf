from pathlib import Path
import os
import pandas as pd
import networkx as nx
import pickle
from datetime import datetime

sbi_coex_dir = Path('/home/daffa/Work/2025/11-JSPP67/sbi_coex')
K = 10
minZ = float(4)

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

def main(coex_dir, output_dirname):
    edges = coexdir_to_edgeslist(coex_dir)
    G = build_graph(edges)
    
    os.makedirs(output_dirname, exist_ok=True)
    
    save_object(G, f'{output_dirname}/sbi_G_prpF.pkl')
    save_object(edges, f'{output_dirname}/sbi_edges.pkl')

if __name__ == '__main__':
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir = f'TAA_cluster-1st_attempt/{timestamp}'
    main(sbi_coex_dir, output_dir)