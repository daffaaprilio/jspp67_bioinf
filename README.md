# JSPP67 Bioinformatics
This repository contains scripts written to support the poster presentation in the 67th Annual Meeting of the Japan Society of Plant Physiology (JSPP67). 
To be specific, codes written in this repository will revolve around the utilization of the sorghum gene co-expression database.

Currently, this repository contains the following module(s). Other modules are coming soon.
- HMM-based homology detection
- KEGG scraping utility to obtain genes of interest
- Gene network cluster creation

## To run each module separately
Recommended for current stage: network inference method.
```shell
# cd to top working directory
cd /Users/daffa/Documents/Work/jspp67_bioinf
# download KEGG FTP data
snakemake -c 8 -s workflows/kegg_scraping/kegg_ftp.smk --configfile configs/config-kegg_ftp.yaml
cd data/reference/KEGG # remove snakemake timestamp
find . -type f -name '*.snakemake_timestamp' -exec rm {} +
find . -type f -name '*.tar.gz' -exec sh -c 'tar -xzvf "$1" -C "$(dirname "$1")"' _ {} \; # unzip all .tar.gz files first then the remaining .gz files
find . -type f -name '*.gz' -exec gunzip {} \;
# KEGG scraping (obtain neighboring pathways)
snakemake -c 8 -s workflows/kegg_scraping/Snakefile --configfile configs/config.yaml
# build network objects
snakemake -c 8 -s workflows/gene_network/Snakefile --configfile configs/config.yaml -np
```
Alternative method: find previously annotated TAA-related genes. Requires 
obtaining amino acid sequences from reviewing relevant literatures, then 
save them in the `data/input` directory.
```shell

```
