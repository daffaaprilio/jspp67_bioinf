# HMM Homology workflow
In some cases, the gene that we want to find its coexpressed genes does not have annotation in the target organism. This workflow allows the identification of gene(s) homologous to the query gene and available in the target organism.

## Required package
Installed in a conda environment `sbi`:
```shell
# ensure mamba is available
conda install -n base -c conda-forge mamba -y
# install required bioinformatics package
mamba install -n sbi -c conda-forge -c bioconda snakemake mafft blast hmmer entrez-direct
```

## Extended pipeline
```shell
# obtain accession sequences from results/hmm_homology/tbrB_refseq/tbrB_refseq-07-homologous_geneID.txt
cut -f1 results/hmm_homology/tbrB_refseq/tbrB_refseq-07-homologous_geneID.txt > results/hmm_homology/tbrB_refseq/tbrB_refseq-08-homologous_sequences_accid.txt
seqkit grep -f results/hmm_homology/tbrB_refseq/tbrB_refseq-08-homologous_sequences_accid.txt data/reference/GCF_000003195.3_Sorghum_bicolor_NCBIv3_protein.faa > results/hmm_homology/tbrB_refseq/tbrB_refseq-09-homologous_sequences.faa
# perform multiple sequence alignment
mafft --auto results/hmm_homology/tbrB_refseq/tbrB_refseq-09-homologous_sequences.faa > results/hmm_homology/tbrB_refseq/tbrB_refseq-10-aligned_homologous_sequences.faa
```

```shell
```