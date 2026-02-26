# HMM Homology workflow
Fortunately, we have access to databases such as KEGG, that is proven to be useful in the case of searching for the gene consituents of our pathway of intereset: the TAA biosynthetic pathway.

Although some key genes might not be available yet in the database, we can infer the TAA biosynthetic pathway cluster by looking at its adjacent pathway cluster, which in this case, the TCA cycle and the C5 metabolic pathway.

This workflow consists of code(s) to automatically obtain genes related to TCA cyle & C5-dibasic branched acid metabolic pathway from KEGG FTP.

## Data Preparation
Run this command on working directory (`config['wdir']`).
```shell
# download & prepare ko.tar.gz, KEGG species files (up to now still limited to ATTED-II 17 plants species), and map.tar.gz
snakemake -c 6 -s workflows/kegg_scraping/kegg_ftp.smk
# remove snakemake timestamp
cd data/reference/KEGG
find . -type f -name '*.snakemake_timestamp' -exec rm {} +
# unzip all .tar.gz files first then the remaining .gz files
find . -type f -name '*.tar.gz' -exec sh -c 'tar -xzvf "$1" -C "$(dirname "$1")"' _ {} \;
find . -type f -name '*.gz' -exec gunzip {} \;
```
