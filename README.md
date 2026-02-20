# Ligand FEP Workflows

Snakemake-based ligand free-energy perturbation workflows using BioSimSpace and sire, with support for both Relative (RBFE) and Absolute (ABFE) binding free energy calculations.

Supports GROMACS and SOMD2 as production engines, with modular minimisation/equilibration staging.

## Quick Start

```bash
# 1. Create and activate environment
conda env create -f environment.yml
conda activate workflow

# 2. Run RBFE workflow
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml \
    --profile profiles/ --executor slurm --resources gpu=3

# 3. Or run ABFE workflow
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml \
    --profile profiles/ --executor slurm --resources gpu=3
```

Both workflows share a single `workflow/Snakefile` entry point. The `method` key in the config file (`rbfe` or `abfe`) determines which rules are loaded.

## Installation

### Environment

```bash
conda env create -f environment.yml
conda activate workflow
```

This installs BioSimSpace, SOMD2/sire, GROMACS (with CUDA), Snakemake, AmberTools, and the SLURM executor plugin.

### CUDA / GPU Compatibility

The conda-installed GROMACS and OpenMM both require compatible CUDA versions. If GROMACS fails with CUDA errors but OpenMM works:

1. Remove the conda GROMACS: `conda remove gromacs`
2. Build GROMACS from source with your CUDA toolkit:

   ```bash
   cmake .. -DGMX_GPU=CUDA -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX
   make -j$(nproc) && make install
   ```

### Input Files

1. **Ligand files**: SDF format with 3D coordinates in the correct binding pose, one per file in `workflow/inputs/ligands/`
2. **Protein files**: Parameterised AMBER format (`.prm7` and `.rst7`) in `workflow/inputs/protein/`

## Configuration

RBFE and ABFE each have their own config file:

- `config/config_rbfe.yml` -- RBFE settings (includes `method: rbfe`)
- `config/config_abfe.yml` -- ABFE settings (includes `method: abfe`)

Schema validation is defined in `config/schemas/`. Full documentation:

- `doc/rbfe_settings.md` -- RBFE config reference
- `doc/abfe_settings.md` -- ABFE config reference

### Production Engines

Both workflows support multiple engines:

- **GROMACS** (`engine: gromacs`): Uses BioSimSpace for setup and GROMACS for MD. Requires HMR factor in config.
- **SOMD2** (`engine: somd2`): Uses sire/OpenMM. Handles HMR internally. Hamiltonian replica exchange (`runner: repex`) is the default. For standard uncoupled windows set `runner: standard`.

The engine can be overridden on the command line: `--config engine=somd2`

## RBFE Workflow

Relative binding free energy calculations comparing ligand pairs through a perturbation network.

### Stages

| Stage | Description |
|-------|-------------|
| **network_prep** | Generate perturbation network using LOMAP |
| **setup** | Parameterise ligands (GAFF2), solvate free and bound systems |
| **min/eq** | Minimise and equilibrate end states (configurable multi-stage protocol) |
| **prepare** | Create merged (alchemical) systems, apply HMR |
| **production** | Production MD per lambda window per replica |
| **analysis** | MBAR analysis, network-wide DDG estimation |

### Network Preparation

If `{working_directory}/network/network.dat` already exists, the network preparation step is skipped. This lets you manually define or edit the perturbation network before running the workflow.

## ABFE Workflow

Absolute binding free energy calculations by decoupling individual ligands from their environment.

### Stages

| Stage | Description |
|-------|-------------|
| **setup** | Parameterise ligand, solvate free and bound systems |
| **min/eq** | Minimise and equilibrate (same configurable protocol as RBFE) |
| **restraint_search** | Find optimal Boresch restraints for bound leg |
| **prepare** | Mark ligand for decoupling, apply HMR (GROMACS only) |
| **production** | Production MD with decoupling schedule |
| **analysis** | MBAR analysis, DG_bind via thermodynamic cycle |

### Boresch Restraints

The bound leg uses Boresch restraints (1 distance, 2 angles, 3 dihedrals) to maintain ligand orientation during decoupling. The analytical correction for releasing these restraints is applied in the final analysis.

## Running on SLURM

A SLURM profile is provided at `profiles/config.yaml`. It is pre-configured for a 3-GPU workstation and will need adjusting for your cluster. The key settings to review:

| Setting | Default | Description |
|---------|---------|-------------|
| `slurm_account` | `default` | Your SLURM account/allocation name |
| `mem_mb` | `5000` | Memory per job (MB) — increase for large systems |
| `runtime` | `4320m` (72h) | Wall time per job — production runs ~52h for 21 lambda windows |
| `gpu` | `1` | GPUs per job |
| `cpus_per_gpu` | `5` | CPUs allocated per GPU — should be >= `simulation_threads` in your config |

```bash
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml \
    --profile profiles/ --executor slurm --resources gpu=3
```

The `--resources gpu=N` flag limits the total number of concurrent GPU jobs. Set this to match the number of GPUs available on your machine or allocation.

## Monitoring

```bash
# Full status
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml status --cores 1

# Results only
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml results --cores 1
```

## Min/Eq Protocol

Minimisation and equilibration stages are fully configurable. Stages execute in the order defined in the config file and can be added or removed freely. The default protocol is based on [Roe & Brooks (2020)](https://doi.org/10.1063/5.0013849).

## Useful Commands

```bash
# Dry run
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml -n

# Unlock after interrupted run
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml --unlock

# Clean analysis only
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml clean_analysis --cores 1

# Clean everything
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml clean --cores 1

# Visualise rule graph (workflow structure)
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml --rulegraph | dot -Tpng > rulegraph.png

# Visualise full DAG (every job)
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml --dag | dot -Tpng > dag.png
```

## File Structure

```
workflow/
  Snakefile                          # Unified entry point (reads method from config)
  rules/
    common.smk                       # Conditional schema validation
    ligprep.smk                      # Shared: parameterise and solvate
    mineq_stages_bound.smk           # Shared: dynamic min/eq for bound leg
    mineq_stages_free.smk            # Shared: dynamic min/eq for free leg
    rbfe/
      network_prep.smk               # RBFE: perturbation network generation
      prepare.smk                    # RBFE: merged system creation
      production.smk                 # RBFE: production MD
      analyse.smk                    # RBFE: MBAR analysis and DDG
    abfe/
      restraint_search.smk           # ABFE: Boresch restraint optimisation
      prepare.smk                    # ABFE: decoupling + HMR (GROMACS)
      production.smk                 # ABFE: GROMACS production
      production_somd2.smk           # ABFE: SOMD2 production
      analyse.smk                    # ABFE: MBAR analysis and DG_bind
  scripts/
    param_solvate.py                 # Shared: parameterisation and solvation
    minimisation.py                  # Shared: minimisation via BSS
    equilibration.py                 # Shared: equilibration via BSS
    rbfe/
      prepare_network.py             # LOMAP network generation
      rbfe_prep.py                   # Alchemical system setup
      production.py                  # Production MD
      analyse_single_leg.py          # Per-leg MBAR analysis
      final_analysis.py              # Network-wide DDG
      status.py                      # Progress monitoring
    abfe/
      restraint_search.py            # Restraint search
      abfe_prep.py                   # Decoupling + HMR
      production.py                  # GROMACS production
      production_somd2.py            # SOMD2 production
      analyse_leg.py                 # Per-leg analysis
      final_analysis.py              # Thermodynamic cycle
      status.py                      # Progress monitoring
  inputs/
    ligands/                         # SDF files
    protein/                         # AMBER prm7/rst7
config/
  config_rbfe.yml                    # RBFE configuration
  config_abfe.yml                    # ABFE configuration
  schemas/
    config_rbfe.schema.yml           # RBFE schema
    config_abfe.schema.yml           # ABFE schema
```
