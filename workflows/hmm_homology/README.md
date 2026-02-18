# HMM Homology workflow
In some cases, the gene that we want to find its coexpressed genes does not have annotation in the target organism. This workflow allows the identification of gene(s) homologous to the query gene and available in the target organism.

## Required package
Installed from `conda`:
```shell
# Snakemake
conda install bioconda::snakemake
# Bioinformatics tool
## BLAST
conda install bioconda::blast
```