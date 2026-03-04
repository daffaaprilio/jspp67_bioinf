#!/usr/bin/env Rscript
# dss_dmr.R
# ─────────────────────────────────────────────────────────────────────────────
# Pairwise DMR detection using the DSS Bioconductor package.
#
# Input:  Two four-column TSV files (chr, pos, N, X) produced by
#         bedmethyl_to_dss.py; no replicates required (single-sample mode).
# Output: TSV with DSS callDMR results.
#
# Required R packages:
#   Rscript -e "if (!requireNamespace('BiocManager')) install.packages('BiocManager')"
#   Rscript -e "BiocManager::install(c('DSS', 'data.table'))"
#
# Usage (called by Snakemake):
#   Rscript dss_dmr.R \
#     --sample_a A_CpG.dss.tsv \
#     --sample_b B_CpG.dss.tsv \
#     --comparison A_vs_B \
#     --delta 0.1 \
#     --p_threshold 0.05 \
#     --min_cpg 3 \
#     --min_len 50 \
#     --output dmr_results.tsv
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(optparse)
  library(DSS)
  library(data.table)
})

# ── CLI arguments ─────────────────────────────────────────────────────────────

option_list <- list(
  make_option("--sample_a",    type = "character", help = "DSS TSV for sample A"),
  make_option("--sample_b",    type = "character", help = "DSS TSV for sample B"),
  make_option("--comparison",  type = "character", help = "Comparison name (for logging)"),
  make_option("--delta",       type = "double",   default = 0.1,  help = "Min |Δmeth| to call DMR"),
  make_option("--p_threshold", type = "double",   default = 0.05, help = "DML p-value threshold"),
  make_option("--min_cpg",     type = "integer",  default = 3,    help = "Min CpGs per DMR"),
  make_option("--min_len",     type = "integer",  default = 50,   help = "Min DMR length (bp)"),
  make_option("--output",      type = "character", help = "Output TSV path")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ── Validate inputs ───────────────────────────────────────────────────────────

for (arg in c("sample_a", "sample_b", "output")) {
  if (is.null(opt[[arg]])) stop(sprintf("Missing required argument: --%s", arg))
}

message("[dss_dmr] Comparison: ", opt$comparison)
message("[dss_dmr] Sample A:   ", opt$sample_a)
message("[dss_dmr] Sample B:   ", opt$sample_b)

# ── Load data ─────────────────────────────────────────────────────────────────

read_dss_input <- function(path) {
  dt <- fread(path, header = TRUE, sep = "\t",
              colClasses = list(character = "chr", integer = c("pos", "N", "X")))
  setnames(dt, c("chr", "pos", "N", "X"))
  dt[N > 0]            # drop zero-coverage sites
}

dat_a <- read_dss_input(opt$sample_a)
dat_b <- read_dss_input(opt$sample_b)

message(sprintf("[dss_dmr] Sites loaded — A: %d, B: %d", nrow(dat_a), nrow(dat_b)))

# ── Build BSseq objects ───────────────────────────────────────────────────────

make_bs <- function(dt, sample_name) {
  makeBSseqData(
    list(as.data.frame(dt)),
    sampleNames = sample_name
  )
}

bs_a <- make_bs(dat_a, basename(opt$sample_a))
bs_b <- make_bs(dat_b, basename(opt$sample_b))

# ── DML test (no replication — use smoothing) ─────────────────────────────────

message("[dss_dmr] Running DMLtest …")
dml_result <- DMLtest(
  bs_a, bs_b,
  smoothing           = TRUE,
  smoothing.span      = 500,
  equal.disp          = FALSE
)

# ── Call DMRs ─────────────────────────────────────────────────────────────────

message(sprintf(
  "[dss_dmr] Calling DMRs (delta=%.2f, p<%.3f, minCpG=%d, minLen=%d) …",
  opt$delta, opt$p_threshold, opt$min_cpg, opt$min_len
))

dmrs <- callDMR(
  dml_result,
  delta     = opt$delta,
  p.threshold = opt$p_threshold,
  minCG     = opt$min_cpg,
  minlen    = opt$min_len,
  dis.merge = 100          # merge adjacent DMRs within 100 bp
)

message(sprintf("[dss_dmr] Found %d DMRs.", nrow(dmrs)))

# ── Write output ──────────────────────────────────────────────────────────────

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

if (nrow(dmrs) > 0) {
  fwrite(as.data.table(dmrs), opt$output, sep = "\t", quote = FALSE)
} else {
  # Write empty file with header so downstream rules don't fail
  header_cols <- c("chr", "start", "end", "length", "nCG",
                   "meanMethy1", "meanMethy2", "diff.Methy", "areaStat")
  fwrite(data.table(matrix(ncol = length(header_cols), nrow = 0,
                            dimnames = list(NULL, header_cols))),
         opt$output, sep = "\t", quote = FALSE)
  message("[dss_dmr] No DMRs found; empty file written.")
}

message("[dss_dmr] Done → ", opt$output)
