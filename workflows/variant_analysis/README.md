# Variant analysis pipeline
## SV (Structural Variant) detection
### Genome assembly
Currently, there are 4 draft genome assemblies (SBC4, 10, 11, 23) and a reference genome (BTx623). 

The overall quality statistics of the 4 assemblies and the reference:
- QUAST report (contains metric, i.e., size, N50, #contigs, GC%). Private file (directory `/Users/daffaaprilio/Documents/Work/Sbi-TAA_genomics/assembly_eval/QUAST/sorghum_assemblies/*.html`, viewing via browser recommended).
- BUSCO completeness (`/Users/daffaaprilio/Documents/Work/Sbi-TAA_genomics/assembly_eval/BUSCO`).

![QUAST report](./images/quast_report.png)

## SNP/indel detection
Install required software/packages:
```shell
mamba install -n sbi -c conda-forge -c bioconda ncurses samtools -y
```
Steps:
```shell
# Index reference genome
samtools faidx data/reference/asm_BTx623.fna
```