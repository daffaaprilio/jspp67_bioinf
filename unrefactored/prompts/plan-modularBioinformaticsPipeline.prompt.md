# Plan: Modular Snakemake Pipeline Repository

## Overview

Create a **workflow-centric modular repository** using Snakemake's native modularity features. Your current pipeline has good foundations (Snakemake, config files, separated concerns) but suffers from hardcoded paths, disconnected stages, and code duplication. The new layout will organize workflows by biological function, share common utilities, and use a central configuration system.

## Recommended Repository Structure

```
repository_root/
├── Snakefile                    # Main entry point, orchestrates modules
├── config.yaml                  # Global configuration
├── README.md                    # Quick start, module overview
├── CHANGELOG.md                 # Track changes
├── .gitignore
├── profiles/
│   └── default/
│       └── config.yaml          # Snakemake execution settings
├── workflows/                   # Module-based organization
│   ├── hmm_homology/           # HMM search module (from current Snakefile)
│   │   ├── Snakefile
│   │   ├── config.yaml
│   │   ├── README.md
│   │   └── scripts/
│   │       └── convert_ids.py  # from p01
│   ├── coexpression/           # Network module (from p02)
│   │   ├── Snakefile
│   │   ├── config.yaml
│   │   ├── README.md
│   │   └── scripts/
│   │       └── build_network.py
│   └── variant_calling/        # Future module
│       ├── Snakefile
│       ├── config.yaml
│       └── README.md
├── common/                      # Shared resources layer
│   ├── scripts/
│   │   ├── parsers.py          # ID conversion, file parsing utilities
│   │   └── snakemake_utils.py  # Helper functions for rules
│   ├── rules/
│   │   └── id_conversion.smk   # Reusable Snakemake rules
│   ├── envs/
│   │   ├── blast.yaml          # Conda environment definitions
│   │   └── hmmer.yaml
│   └── schemas/
│       └── config.schema.yaml  # Validate configurations
├── data/                        # Data organization convention
│   ├── reference/              # ncbi_data, taxdump
│   ├── input/                  # User-provided queries
│   └── raw/                    # hmm/seed, hmm/target
├── results/                    # Pipeline outputs
│   ├── hmm_homology/
│   │   └── {timestamp}/
│   └── coexpression/
│       └── {timestamp}/
├── resources/                  # Downloaded databases (large files)
├── .test/                      # Test data and configs
│   ├── data/
│   └── test.yaml
├── test.sh                     # Run test workflow
└── docs/                       # Documentation
    ├── modules/
    │   ├── hmm_homology.md
    │   ├── coexpression.md
    │   └── variant_calling.md
    └── adding_modules.md
```

## Implementation Steps

### 1. Create module-based directory structure

Organize workflows by biological function:
- `workflows/hmm_homology/` - HMM search module (from current Snakefile)
- `workflows/coexpression/` - Network module (from p02)
- `workflows/variant_calling/` - Future module

Each module contains:
- `Snakefile` - Module-specific workflow rules
- `config.yaml` - Module-specific configuration
- `scripts/` - Module-specific Python scripts
- `README.md` - Module documentation

### 2. Create shared resources layer

Establish `common/` directory for reusable components:
- `common/scripts/parsers.py` - ID conversion, file parsing utilities
- `common/scripts/snakemake_utils.py` - Helper functions for rules
- `common/rules/id_conversion.smk` - Reusable Snakemake rules
- `common/envs/blast.yaml` - Conda environment definitions
- `common/envs/hmmer.yaml`
- `common/schemas/config.schema.yaml` - Configuration validation

### 3. Implement main entry point

At repository root:
- **Snakefile (main)**: Imports module workflows, orchestrates dependencies
- **config.yaml (global)**: Paths, organism info, module toggles
- **profiles/default/config.yaml**: Snakemake execution settings (cores, etc.)

Example root Snakefile structure:
```snakemake
module hmm_homology:
    snakefile: "workflows/hmm_homology/Snakefile"
    config: config["hmm_homology"]

use rule * from hmm_homology as hmm_*

module coexpression:
    snakefile: "workflows/coexpression/Snakefile"
    config: config["coexpression"]

use rule * from coexpression as coex_*

rule all:
    input:
        rules.hmm_search_target.output,
        "results/coexpression/network.gpickle"
```

### 4. Refactor current Snakefile HMM rules

Move HMM workflow to `workflows/hmm_homology/`:
- Extract rules (lines 14-102 of current Snakefile) to module Snakefile
- Replace hardcoded paths (lines 6-9) with `config["paths"]["target"]`, etc.
- Create module-specific config with:
  - `seed_fasta`
  - `target_proteome`
  - `prefix`
  - `evalue_threshold`
  - `coverage_threshold`

### 5. Refactor p01 and p02 scripts

**ID Conversion (p01-convert_to_geneID.py)**:
- Move ID conversion logic to `common/scripts/parsers.py` as reusable functions
- Create Snakemake rule wrapper in `workflows/hmm_homology/Snakefile` that calls conversion
- Remove hardcoded paths (p01 lines 4-7)

**Coexpression Network (p02-create_network.py)**:
- Move to `workflows/coexpression/scripts/build_network.py`
- Remove hardcoded paths (p02 line 67: `/Users/daffa/...`)
- Add Snakemake rules to `workflows/coexpression/Snakefile`
- Keep function structure (`coexdir_to_edgeslist()`, `build_graph()`)

### 6. Set up data organization convention

Establish clear data directory structure:
```
data/
  reference/         # ncbi_data, taxdump
  input/            # User-provided queries
  raw/              # hmm/seed, hmm/target
results/
  {module_name}/
    {timestamp}/    # Timestamped runs
resources/          # Downloaded databases (large files)
```

### 7. Create integration in root Snakefile

Use Snakemake's `module` directive to import and use module workflows:
- Import each module with namespace
- Use `use rule * from {module} as {prefix}_*` pattern
- Define main `rule all` that chains module outputs
- Allow dependency flow between modules

### 8. Establish configuration hierarchy

Create root config.yaml with module-specific sections:
```yaml
global:
  organism: "Sorghum_bicolor"
  ncbi_data_dir: "data/reference/ncbi_data"
  
hmm_homology:
  enabled: true
  seed: "data/input/prpF.faa"
  target: "data/reference/GCF_000003195.3_protein.faa"
  evalue: 1e-5
  coverage: 0.5
  
coexpression:
  enabled: true
  coex_data_dir: "data/raw/sbi_coex"
  min_z_score: 4
  k_neighbors: 10
  
variant_calling:
  enabled: false  # Future module
```

### 9. Add documentation structure

Create comprehensive documentation:
- **README.md**: Quick start, module overview, installation
- **docs/modules/{module_name}.md**: Detailed module documentation
- **docs/adding_modules.md**: Guide for adding new modules
- **CHANGELOG.md**: Track changes
- Each workflow subdirectory gets its own README explaining:
  - Purpose
  - Inputs/outputs
  - Configuration options
  - Example usage

### 10. Create example/test data setup

Establish testing infrastructure:
- `.test/` directory with small test datasets
- `test.yaml` config for CI testing
- Shell script to run test workflow: `bash test.sh`
- Example commands for running individual modules

## Verification Steps

1. **Run existing HMM workflow**: 
   ```bash
   snakemake --configfile config.yaml --cores 4
   ```

2. **Test module in isolation**:
   ```bash
   snakemake -s workflows/hmm_homology/Snakefile --config seed=test.faa
   ```

3. **Verify ID conversion integrates**: 
   - Check that gene IDs are automatically produced after HMM results
   - No manual script execution required

4. **Run full pipeline**: 
   - HMM search → ID conversion → coexpression network (automatic chaining)
   - Verify data flows between modules

5. **Add a new dummy module**: 
   - Follow the pattern
   - Ensure it integrates without modifying existing modules
   - Test enabling/disabling via config

## Design Decisions

### Snakemake modules over Python package
- Fits your workflow-focused use case
- No need for pip installation complexity
- Lab members can clone and run directly
- Easier to modify and extend for non-programmers

### Biological function organization
- Each `workflows/*` folder is self-contained by biological purpose
- Easier for lab members to understand than technical organization
- Clear mental model: one folder = one analysis type
- Reduces cognitive load when working on specific analysis

### Shared `common/` layer
- Avoids code duplication across modules
- Keeps modules loosely coupled
- Shared utilities (ID parsers, file readers) in one place
- Easy to version and test common code

### Timestamped results
- Preserves run history for comparison
- You're already doing this in p02 (lines 71-77)
- Prevents accidental overwrites
- Facilitates reproducibility

### Config-driven module enabling
- Can turn modules on/off without modifying code
- Flexible pipeline composition
- Easy to skip expensive steps during development
- Clear declaration of what will run

### Root Snakefile as orchestrator
- Single entry point for full pipeline
- Modules can still run standalone for development/testing
- Clear dependency chain visible in one place
- Follows Snakemake best practices

## Current Issues Being Addressed

### Hardcoded paths (FIXED)
- **Before**: Paths hardcoded in Snakefile (lines 6-9), p01 (lines 4-7), p02 (line 67)
- **After**: All paths from config, portable across systems

### Disconnected stages (FIXED)
- **Before**: Three stages must be manually run (Snakemake → p01 → p02)
- **After**: Integrated workflow, automatic dependency tracking

### Code duplication (FIXED)
- **Before**: Archive Snakefiles are 99% identical
- **After**: Single module definition, reusable for multiple genes

### Limited configuration (FIXED)
- **Before**: config.yaml only has 3 fields
- **After**: Comprehensive config with module-specific sections

### Poor reusability (FIXED)
- **Before**: Can't easily run for different organisms
- **After**: Change config, run anywhere

### No module abstraction (FIXED)
- **Before**: Monolithic pipeline
- **After**: Clear module boundaries, easy to extend

## Next Steps After Implementation

1. **Migrate existing data**: Move current data to new directory structure
2. **Test with real data**: Run full pipeline on actual Sorghum prpF example
3. **Add variant calling module**: Implement using existing ASM drafts
4. **Create CI/CD**: Set up automated testing on small datasets
5. **Document for lab**: Write user guide for lab members
6. **Add more modules**: Phylogenetics, expression analysis, etc.

## Benefits of This Approach

- **Modularity**: Add new modules without touching existing code
- **Maintainability**: Clear separation of concerns, easy to debug
- **Reproducibility**: Config-driven, documented, version-controlled
- **Flexibility**: Run full pipeline or individual modules
- **Collaboration**: Lab members can work on different modules independently
- **Scalability**: Easy to add new organisms, genes, or analysis types
- **Best practices**: Follows Snakemake conventions, uses conda environments
