# DMR analysis

## File preparation
```shell
# obtained SAM from matsu, convert to BAM
# base command
samtools view --bam out.bam in.sam
# implementation
cd data/reads/
for file in read_SBC*.sam; do samtools view -b -o "${file:r}.bam" "$file"; done
```