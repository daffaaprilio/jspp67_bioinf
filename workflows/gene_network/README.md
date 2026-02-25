# Co-expression network construction

## Steps
**Step 1.** Compile top hit gene IDs from `results/gene_network/map*-K*_kegg/map*-K*_kegg-01-homologous_geneID.txt`. Annotate these gene IDs with KEGG pathway (`map_id`) and orthology (`ko_id`) information.

**Step 2.** Load the co-expression data as a list of edges.

**Step 3.** Load these edges into graph object.

**Step 4.**

**Step last.** Validation

## Steps
**Step 1.** Build the full co-expression graph (apply filtering parameters, i.e., `K=10`, `minZ=4.0`)

**Step 2.** 
