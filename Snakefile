# Root Snakefile — orchestrates all workflow modules.

configfile: "configs/config.yaml"


# ── Modules ──────────────────────────────────────────────────────────────────

module hmm_homology:
    snakefile: "workflows/hmm_homology/Snakefile"
    config: config

use rule * from hmm_homology as hmm_homology_*


module variant_analysis:
    snakefile: "workflows/variant_analysis/Snakefile"
    config: config

use rule * from variant_analysis as variant_analysis_*


module dmr_analysis:
    snakefile: "workflows/dmr_analysis/Snakefile"
    config: config

use rule * from dmr_analysis as dmr_analysis_*


# ── Default target ────────────────────────────────────────────────────────────

rule all:
    input:
        rules.hmm_homology_all.input,
        rules.variant_analysis_all.input,
        rules.dmr_analysis_all.input,
