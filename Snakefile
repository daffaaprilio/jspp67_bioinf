# Root Snakefile — orchestrates all workflow modules.

configfile: "config.yaml"


# ── Modules ──────────────────────────────────────────────────────────────────

module hmm_homology:
    snakefile: "workflows/hmm_homology/Snakefile"
    config: config["hmm_homology"]

use rule * from hmm_homology as hmm_homology_*


# ── Default target ────────────────────────────────────────────────────────────

rule all:
    input:
        rules.hmm_homology_all.input,
