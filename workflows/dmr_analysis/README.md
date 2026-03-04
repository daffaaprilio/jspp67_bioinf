# DMR Analysis

Differentially Methylated Region (DMR) detection from **ONT native methylation** BAMs produced by **dorado duplex** with `4mC_5mC,6mA` modification calling (MM/ML tags).

> **This is NOT bisulfite sequencing.** Use `modkit pileup` to extract methylation from MM/ML tags — not Bismark or MethylDackel.

## Samples

| Sample | BAM location | Method |
|--------|-------------|--------|
| r0066  | `$dmr_bam_dir/r0066.bam` | dorado duplex `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` |

> Add more samples by dropping BAMs into `dmr_bam_dir` and updating `configs/config.yaml`.

## Tool chain

```
modkit pileup → bedmethyl_to_dss.py → DSS (R: DMLtest + callDMR)
```

| Step | Rule | Output |
|------|------|--------|
| Methylation pileup | `extract_methylation` | `results/dmr_analysis/{sample}/{sample}_CpG.bedMethyl` |
| Format conversion | `bedmethyl_to_dss` | `results/dmr_analysis/{sample}/{sample}_CpG.dss.tsv` |
| DMR calling | `call_dmr` | `results/dmr_analysis/{comparison}/dmr_results.tsv` |

## Requirements

```shell
mamba install -n sbi -c conda-forge -c bioconda modkit samtools -y
Rscript -e "BiocManager::install(c('DSS', 'data.table'))"
Rscript -e "install.packages('optparse')"
```

## Usage

### Standalone (pileup only — works with 1 sample)
```shell
snakemake -c 8 \
  -s workflows/dmr_analysis/Snakefile \
  --configfile configs/config.yaml
```

### Adding more samples and running DMR calling

1. Drop new BAM into `dmr_bam_dir` (naming: `{sample}.bam`)
2. Add sample name to `dmr_samples` in `configs/config.yaml`
3. Add a pairwise comparison to `dmr_comparisons`:
   ```yaml
   dmr_comparisons:
     - name: r0066_vs_r0067
       sample_a: r0066
       sample_b: r0067
   ```
4. Re-run Snakemake.

## DMR parameters (tunable in `configs/config.yaml`)

| Parameter | Default | Meaning |
|-----------|---------|----------|
| `dss_delta` | `0.1` | Min \|ΔMeth\| to call a DMR |
| `dss_p_threshold` | `0.05` | DML p-value threshold |
| `dss_min_cpg` | `3` | Min CpG sites per DMR |
| `dss_min_len` | `50` | Min DMR length (bp) |

## Outputs

```
results/dmr_analysis/
├── r0066/
│   ├── r0066_CpG.bedMethyl          ← modkit pileup output
│   └── r0066_CpG.dss.tsv            ← reformatted for DSS (chr pos N X)
└── {comparison}/
    └── dmr_results.tsv              ← DSS callDMR results
```

## Prototyping notebook

`workflows/dmr_analysis/notebook.ipynb` — step-by-step prototype for r0066.