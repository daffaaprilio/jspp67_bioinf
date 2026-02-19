# HMM Homology workflow
In some cases, the gene that we want to find its coexpressed genes does not have annotation in the target organism. This workflow allows the identification of gene(s) homologous to the query gene and available in the target organism.

## Required package
Installed in a conda environment `sbi`:
```shell
# ensure mamba is available
conda install -n base -c conda-forge mamba -y
# install required bioinformatics package
mamba install -n sbi -c conda-forge -c bioconda snakemake mafft blast hmmer
```