# ABFE Workflow Settings Documentation

This document describes all configuration options for the Absolute Binding Free
Energy (ABFE) workflow.

## Overview

The ABFE workflow calculates absolute binding free energies using alchemical
decoupling. Unlike RBFE (Relative Binding Free Energy) which perturbs between
two ligands, ABFE completely decouples a single ligand from its environment
through three alchemical stages:

**Bound leg (protein-ligand complex):**
1. **Restrain**: Apply Boresch restraints to confine the ligand
2. **Discharge**: Remove electrostatic interactions
3. **Vanish**: Remove van der Waals interactions

**Free leg (ligand in solvent):**
1. **Discharge**: Remove electrostatic interactions
2. **Vanish**: Remove van der Waals interactions

The binding free energy is computed as:

```
ΔG_bind = ΔG_free,discharge + ΔG_free,vanish
        - ΔG_correction
        - ΔG_bound,vanish
        - ΔG_bound,discharge
        - ΔG_bound,restrain
```

Where `ΔG_correction` is the analytical correction for releasing the Boresch
restraints.

## Configuration File Structure

The configuration is defined in `config/config_abfe.yml` and validated against
the JSON schema in `config/schemas/config_abfe.schema.yml`.

### Input Files

```yaml
# Glob pattern for ligand input files (SDF format)
ligands: "workflow/inputs/ligands/*.sdf"

# Glob pattern for protein input files (parameterised)
# Expects .rst7 and .prm7 files from AmberTools
protein_files: "workflow/inputs/protein/protein.*"

# Working directory for all workflow outputs
working_directory: "abfe_runs"

# Number of threads for simulations
simulation_threads: 5
```

### Setup Settings

```yaml
setup_settings:
  # Force field for ligand parameterisation
  # Options: gaff, gaff2, openff_unconstrained-2.0.0, etc.
  ligand_forcefield: "gaff2"

  # Water model for solvation
  # Options: tip3p, tip4p, spce, etc.
  water_model: "tip3p"

  # Ion concentration in mol/L for charge neutralisation
  ion_concentration: 0.15

  # Box padding in nm (added to bounding box)
  box_length: 5.0

  # Box geometry type
  # Options: cubic, triclinic, orthorhombic
  box_type: "cubic"
```

### Restraint Search Settings

```yaml
restraint_search_settings:
  # Runtime for unrestrained equilibration before restraint search
  search_runtime: "1ns"

  # Temperature for restraint search simulation
  temperature: "300K"

  # Method for restraint selection
  # Currently only "BSS" (BioSimSpace built-in) is supported
  method: "BSS"

  # Type of restraint
  # Currently only "Boresch" is supported
  restraint_type: "Boresch"
```

### Minimisation and Equilibration Stages

The workflow supports flexible protocol definition. Each stage can be either
a minimisation or equilibration step.

**Free leg stages** (`min/eq-stages-free`):
```yaml
min/eq-stages-free:
  minimisation1:
    engine: Gromacs
    minimisation-steps: 1000
    restraint-string: heavy        # Restrain heavy atoms
    restraint-force-constant: 5.0  # kcal/(mol·Å²)

  equilibration1:
    engine: Gromacs
    timestep: 1fs
    runtime: 15ps
    temperature: 300K
    thermostat-time-constant: 0.5ps
    restraint-string: heavy
    restraint-force-constant: 5.0

  # ... additional stages
```

**Bound leg stages** (`min/eq-stages-bound`):
Similar structure, but may include backbone restraint stages for protein
stability.

### Production Settings

```yaml
production-settings:
  # Number of independent replicas for error estimation
  num_replicas: 3

  # Hydrogen mass repartitioning factor
  # Use 3 for GROMACS with 4fs timestep, 1.5 for SOMD
  hmr_factor: 3

  # Lambda schedules for each ABFE stage
  lambda_schedules:
    # Restrain stage: few windows needed (smooth transformation)
    restrain:
      - 0.0
      - 0.25
      - 0.5
      - 0.75
      - 1.0

    # Discharge stage: moderate windows for electrostatic decoupling
    discharge:
      - 0.0
      - 0.125
      - 0.25
      - 0.375
      - 0.5
      - 0.625
      - 0.75
      - 0.875
      - 1.0

    # Vanish stage: dense windows needed for LJ softcore
    vanish:
      - 0.0
      - 0.05
      - 0.1
      # ... (typically 20+ windows)
      - 0.95
      - 1.0

  # Settings for FREE leg production
  free-leg-settings:
    runtime: 2ns           # Per lambda window
    timestep: 4fs          # Requires HMR
    temperature: 300K
    pressure: 1bar
    report-interval: 250   # Steps between energy reports
    restart-interval: 25000  # Steps between checkpoints

  # Settings for BOUND leg production
  bound-leg-settings:
    runtime: 2ns
    timestep: 4fs
    temperature: 300K
    pressure: 1bar
    report-interval: 250
    restart-interval: 25000
```

### Analysis Settings

```yaml
analysis-settings:
  # Plotting backend
  # Options: native (matplotlib), cinnabar
  backend: native

  # Path to experimental results file (optional)
  # Format: CSV with columns "ligand,value,error"
  experimental-results: workflow/inputs/exp_data.csv

  # Units of experimental data
  # Options: kcal/mol, kJ/mol, Ki_uM
  experimental-units: kcal/mol
```

## Running the Workflow

### Basic Execution

```bash
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml --use-conda --cores 8
```

### With SLURM Cluster

```bash
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml \
    --executor slurm \
    --profile profiles \
    --use-conda \
    --cores 32 \
    --jobs 100 \
    --resources gpu=4
```

### Monitoring Progress and Results

```bash
# Check status of all workflow stages
snakemake status -s workflow/Snakefile --configfile config/config_abfe.yml --cores 1

# Show only binding free energy results (completed ligands)
snakemake results -s workflow/Snakefile --configfile config/config_abfe.yml --cores 1
```

The `status` command shows:
- Preparation stage completion (setup, min/eq, restraints, ABFE prep)
- Production stage completion for each ligand/replica
- Per-stage free energy values and current binding energy estimates

The `results` command shows only the final binding free energy summary.

### Cleaning Workflow Outputs

```bash
# Clean all outputs
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml clean

# Clean only analysis (keep simulations)
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml clean_analysis

# Clean production (keep equilibration)
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml clean_production
```

## Output Structure

```
working_directory/
├── setup/                      # Parameterised and solvated systems
│   ├── {ligand}_free.bss
│   └── {ligand}_bound.bss
├── preparation/                # Minimisation and equilibration stages
│   ├── minimisation/
│   │   ├── 1/
│   │   ├── 2/
│   │   └── ...
│   ├── equilibration/
│   │   ├── 1/
│   │   ├── 2/
│   │   └── ...
│   └── final/                  # Equilibrated systems ready for production
│       ├── {ligand}_free.bss
│       └── {ligand}_bound.bss
├── restraints/                 # Boresch restraint parameters
│   ├── {ligand}_restraint.json
│   └── {ligand}_correction.txt
├── abfe_prepared/              # Systems ready for production
│   ├── {ligand}_free.bss
│   └── {ligand}_bound.bss
├── production/                 # Production simulation outputs
│   └── {ligand}/
│       ├── bound_{replica}/
│       │   ├── restrain/
│       │   │   └── lambda_*/
│       │   ├── discharge/
│       │   │   └── lambda_*/
│       │   └── vanish/
│       │       └── lambda_*/
│       └── free_{replica}/
│           ├── discharge/
│           │   └── lambda_*/
│           └── vanish/
│               └── lambda_*/
├── analysis/                   # Analysis outputs
│   ├── {ligand}/
│   │   ├── bound_{replica}/
│   │   │   ├── restrain/pmf.csv
│   │   │   ├── discharge/pmf.csv
│   │   │   └── vanish/pmf.csv
│   │   └── free_{replica}/
│   │       ├── discharge/pmf.csv
│   │       └── vanish/pmf.csv
│   ├── detailed_abfe_results.csv
│   ├── final_abfe_results.csv
│   └── abfe_comparison.png
└── logs/                       # Log files
```

## Tips and Best Practices

1. **Lambda schedules**: The vanish stage typically needs more lambda windows
   (20+) due to the softcore LJ potential. The restrain stage needs fewer
   windows (4-6) since it's a smooth harmonic potential.

2. **Convergence**: Check the overlap matrices in the analysis output. Good
   overlap (>0.03 off-diagonal elements) indicates sufficient sampling between
   adjacent lambda windows.

3. **Restraint selection**: The restraint search simulation should be long
   enough to sample representative ligand fluctuations (typically 1-5 ns).

4. **HMR**: Hydrogen mass repartitioning allows 4fs timesteps. Use factor=3
   for GROMACS.

5. **Replicas**: Use at least 3 replicas to estimate sampling uncertainty.
   The final error reported is the maximum of propagated MBAR errors and
   replica standard deviation.
