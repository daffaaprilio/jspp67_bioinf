#!/bin/zsh

# conv for Conversion Reference Data
conv=/home/daffa/Work/2025/11-JSPP67/ncbi_data/gene2accession
out=/home/daffa/Work/2025/11-JSPP67/sorghum_gene2accession

# filter to only include sorghum
cat $conv | awk -F'\t' '$1 == 4558' > $out
