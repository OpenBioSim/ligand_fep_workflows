# RBFE Workflow Settings Guide

This guide describes the configuration options for the relative binding free energy (RBFE) workflow. The config file is passed to Snakemake via `--configfile` and validated against `config/schemas/config_rbfe.schema.yml`.

A working example config is provided at `config/config_rbfe.yml`. A minimal test config is at `config/config_rbfe_test.yml`.

---

## Top-Level Settings

These fields are required at the top level of the config file.

```yaml
method: rbfe                                      # Must be "rbfe"
ligands: "workflow/inputs/ligands/*.sdf"          # Glob for ligand SDF files
protein_files: "workflow/inputs/protein/protein.*"  # Glob for protein input files
working_directory: "output/rbfe"                  # Base directory for all outputs
simulation_threads: 5                             # CPU threads per simulation job
```

- `method` must be `rbfe`. This tells the Snakefile which ruleset and schema to use.
- `ligands` and `protein_files` accept glob patterns. All matched files are used.
- `working_directory` is the root under which setup, network, production, and analysis outputs are written.
- `simulation_threads` sets the number of CPU threads available to each simulation job. Relevant for GROMACS and OpenMM.

---

## network_settings

Controls how the perturbation network is built using LOMAP scores.

```yaml
network_settings:
  lambda_windows: 17        # Lambda windows for easy perturbations (LOMAP >= threshold)
  lomap_threshold: 0.4      # Score threshold separating easy from difficult perturbations
  diff_lambda_windows: 21   # Lambda windows for difficult perturbations (LOMAP < threshold)
```

LOMAP scores the structural similarity of each perturbation pair on a scale of 0 to 1. Higher similarity means a smoother alchemical path is possible and fewer lambda windows are required.

- Perturbations with LOMAP score >= `lomap_threshold` use `lambda_windows`.
- Perturbations with LOMAP score < `lomap_threshold` use `diff_lambda_windows`.
- `diff_lambda_windows` should typically be larger than `lambda_windows` to improve sampling for difficult perturbations.

For SOMD2, these per-perturbation window counts are passed directly to the runner. For GROMACS, the lambda schedules (see below) determine the per-leg lambda points, and the number of windows is derived from the schedule length.

---

## setup_settings

Controls system preparation: parameterisation, solvation, and box setup.

```yaml
setup_settings:
  ligand_forcefield: "gaff2"    # Ligand force field ("gaff2" recommended)
  water_model: "tip3p"          # Solvent water model
  ion_concentration: 0.15       # Ionic concentration in mol/L
  box_length: 5.0               # Minimum distance from solute to box edge in nm
  box_type: "cubic"             # Box geometry: "cubic", "triclinic", or "orthorhombic"
```

- `ligand_forcefield` is passed to BioSimSpace for small molecule parameterisation. `gaff2` is the recommended default.
- `ion_concentration` adds counterions to neutralise the system plus additional salt to reach the given concentration.
- `box_length` sets the minimum padding between the solute and the box edge (in nm). Larger values increase system size and simulation cost.
- `box_type` sets the periodic box geometry. `cubic` is simplest; `triclinic` can be more efficient for elongated systems.

---

## Minimisation and Equilibration Stages

Two separate stage lists are required: one for the free (solvated ligand) leg and one for the bound (protein-ligand complex) leg. Each stage key must begin with `minimisation` or `equilibration` followed by an integer (e.g. `minimisation1`, `equilibration2`). Stages are executed in sorted key order.

Supported engines for min/eq: `gromacs`, `openmm`, `amber`.

### min/eq-stages-free

```yaml
min/eq-stages-free:
  minimisation1:
    engine: gromacs
    minimisation-steps: 5000
    restraint-string: "heavy"        # Named selection: "heavy", "backbone", "all"
    restraint-force-constant: 10.0   # kcal/(mol·Å²)

  equilibration1:
    engine: gromacs
    runtime: 100ps
    timestep: 2fs
    temperature-start: 0K
    temperature-end: 300K
    thermostat-time-constant: 1ps
    pressure: 1bar
    restraint-string: "heavy"
    restraint-force-constant: 10.0

  equilibration2:
    engine: gromacs
    runtime: 200ps
    timestep: 2fs
    temperature: 300K
    pressure: 1bar
```

### min/eq-stages-bound

```yaml
min/eq-stages-bound:
  minimisation1:
    engine: gromacs
    minimisation-steps: 5000
    restraint-string: "backbone"
    restraint-force-constant: 10.0

  equilibration1:
    engine: gromacs
    runtime: 100ps
    timestep: 2fs
    temperature-start: 0K
    temperature-end: 300K
    thermostat-time-constant: 1ps
    pressure: 1bar
    restraint-string: "backbone"
    restraint-force-constant: 10.0

  equilibration2:
    engine: gromacs
    runtime: 200ps
    timestep: 2fs
    temperature: 300K
    pressure: 1bar
```

### Stage Options

**Minimisation keys:**

| Key | Required | Description |
|-----|----------|-------------|
| `engine` | yes | `gromacs`, `openmm`, or `amber` |
| `minimisation-steps` | no | Number of energy minimisation steps |
| `restraint-string` | no | Named atom selection to restrain: `heavy`, `backbone`, `all` |
| `restraint-indices` | no | Explicit atom index list (alternative to `restraint-string`) |
| `restraint-force-constant` | no | Force constant in kcal/(mol·Å²) |

**Equilibration keys:**

| Key | Required | Description |
|-----|----------|-------------|
| `engine` | yes | `gromacs`, `openmm`, or `amber` |
| `runtime` | yes | Simulation time (e.g. `100ps`, `1ns`) |
| `timestep` | no | Integration timestep (e.g. `2fs`) |
| `temperature` | no | Target temperature (e.g. `300K`) |
| `temperature-start` | no | Starting temperature for annealing (e.g. `0K`) |
| `temperature-end` | no | Ending temperature for annealing (e.g. `300K`) |
| `thermostat-time-constant` | no | Thermostat coupling time constant (e.g. `1ps`) |
| `pressure` | no | Target pressure (e.g. `1bar`). Omit for NVT |
| `restraint-string` | no | Named atom selection to restrain |
| `restraint-indices` | no | Explicit atom index list |
| `restraint-force-constant` | no | Force constant in kcal/(mol·Å²) |

When both `temperature-start` and `temperature-end` are set, a linear temperature ramp is applied over the stage runtime. Use `temperature` alone for constant-temperature equilibration.

---

## production-settings

Controls production FEP simulations. Supported engines are `somd2` and `gromacs`. Amber is not supported for RBFE production.

```yaml
production-settings:
  num_replicas: 3       # Independent replicas per perturbation/leg
  engine: somd2         # "somd2" or "gromacs" (case-insensitive)
```

- `num_replicas` sets the number of independent repeats for each perturbation leg. A minimum of 3 is recommended for meaningful error estimates.
- The final reported error is the maximum of the propagated MBAR error and the standard deviation across replicas.

### SOMD2 Settings

Used when `engine: somd2`.

```yaml
production-settings:
  num_replicas: 3
  engine: somd2

  somd2-settings:
    runtime: 2ns
    timestep: 4fs
    temperature: 300K
    pressure: 1bar
    use-modified-dummies: true      # Apply ghostly dummy atom corrections
    cutoff_type: RF                 # "RF" (reaction field) or "PME"
    cutoff: 12A
    energy_frequency: 1ps           # Energy output frequency (keep fine for MBAR)
    frame_frequency: 500ps          # Trajectory frame output frequency
    checkpoint_frequency: 500ps     # Checkpoint write frequency
    integrator: langevin_middle
    shift_delta: 1.5A               # Soft-core shift delta parameter
    perturbable_constraint: h_bonds_not_heavy_perturbed
    runner: repex                   # "repex" (replica exchange) or "standard"
    gpus_per_job: 1                 # GPUs allocated per job
```

Key SOMD2 options:

- `use-modified-dummies`: Applies corrections to dummy (ghost) atom bonded terms using the [ghostly](https://github.com/OpenBioSim/ghostly) library. Recommended for most perturbations.
- `cutoff_type`: `RF` (reaction field) is faster and suitable for most use cases. `PME` is more rigorous but computationally more expensive.
- `energy_frequency`: Controls how often energies are written for all lambda windows. Keep this fine-grained (1ps or similar) for accurate MBAR free energy estimates.
- `runner`: `repex` runs replica exchange Monte Carlo between lambda windows and requires one OpenMM context per window. Ensure sufficient total GPU VRAM. Use `standard` for single-GPU or memory-constrained setups.
- `gpus_per_job`: Number of GPUs per job. For `repex` with many windows, increasing this can improve throughput.
- Lambda window counts are determined per perturbation by `network_settings`, not set here.
- There is no `equilibration_time` field for SOMD2 RBFE; pre-production equilibration is handled by the explicit min/eq stages.

### GROMACS Settings

Used when `engine: gromacs`.

```yaml
production-settings:
  num_replicas: 3
  engine: gromacs

  gromacs-settings:
    hmr_factor: 3                   # Hydrogen mass repartitioning factor

    lambda_schedules:
      bound:
        bonded: [0.0, 0.5, 1.0, 1.0, 1.0]
        coul:   [0.0, 0.5, 1.0, 1.0, 1.0]
        vdw:    [0.0, 0.0, 0.0, 0.5, 1.0]
      free:
        coul:   [0.0, 0.5, 1.0, 1.0, 1.0]
        vdw:    [0.0, 0.0, 0.0, 0.5, 1.0]

    free-leg-settings:
      runtime: 2ns
      timestep: 4fs
      temperature: 300K
      pressure: 1bar
      use-modified-dummies: true
      report-interval: 1ps          # Energy/trajectory reporting interval
      restart-interval: 100ps       # Checkpoint write interval

    bound-leg-settings:
      runtime: 2ns
      timestep: 4fs
      temperature: 300K
      pressure: 1bar
      use-modified-dummies: true
      report-interval: 1ps
      restart-interval: 100ps
```

Key GROMACS options:

- `hmr_factor`: Hydrogen mass repartitioning factor. Set to 3 when using a 4fs timestep. HMR is applied during preparation; the production run uses `hmr: false` internally.
- `lambda_schedules`: Defines how bonded, coulomb, and vdW interactions are scaled across lambda windows for each leg. All lists within a schedule must have the same length. The bound leg can include a `bonded` schedule; the free leg typically does not.
- `report-interval` and `restart-interval` are specified in time units and converted to steps internally based on the timestep.
- `use-modified-dummies`: Same ghostly library corrections as SOMD2. Recommended for most perturbations.
- Free and bound legs can have different runtimes if needed.

---

## analysis-settings

Controls how free energy results are computed and plotted.

```yaml
analysis-settings:
  backend: native                                    # "native" or "cinnabar"
  # experimental-results: workflow/inputs/exp_data.csv  # optional
  # experimental-units: kcal/mol                        # optional
```

- `backend`: `native` uses the built-in matplotlib-based analysis. `cinnabar` uses the cinnabar package if installed, which provides additional network-aware analysis.
- `experimental-results`: Optional path to a CSV file containing measured binding affinities. If provided, predicted and experimental values are compared in the output plots.
- `experimental-units`: Units for the experimental data. Accepted values are `kcal/mol`, `kJ/mol`, or `Ki_uM`. Required if `experimental-results` is set.

---

## Running the Workflow

```bash
# Build the perturbation network and prepare all systems (stops before production)
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml \
    --until prepare_rbfe --cores 8

# Run the full workflow (prep + production + analysis)
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml --cores 8

# Run with SLURM
snakemake -s workflow/Snakefile --configfile config/config_rbfe.yml \
    --executor slurm --profile profiles/config.yaml \
    --jobs 100 --resources gpu=4

# Check the status of all perturbations
snakemake status -s workflow/Snakefile --configfile config/config_rbfe.yml --cores 1

# Print only the final binding free energy results
snakemake results -s workflow/Snakefile --configfile config/config_rbfe.yml --cores 1

# Remove all outputs
snakemake clean -s workflow/Snakefile --configfile config/config_rbfe.yml

# Remove only analysis outputs (re-run analysis without re-running simulations)
snakemake clean_analysis -s workflow/Snakefile --configfile config/config_rbfe.yml

# Remove only production outputs
snakemake clean_production -s workflow/Snakefile --configfile config/config_rbfe.yml
```

---

## Output Directory Structure

All outputs are written under `working_directory` (e.g. `output/rbfe/`).

```
{working_directory}/
├── setup/                          # Parameterised and solvated systems
│   ├── {ligand}_free.bss
│   └── {ligand}_bound.bss
├── network/                        # Perturbation network files
│   ├── network.dat                 # Edge list with per-perturbation lambda window counts
│   ├── network.png                 # Network visualisation
│   └── lomap_scores.csv
├── rbfe_prepared/                  # Alchemically merged and equilibrated systems
│   ├── free/
│   │   └── {lig1}~{lig2}.bss
│   └── bound/
│       └── {lig1}~{lig2}.bss
├── production/                     # Production simulation outputs
│   └── {engine}/                   # e.g. somd2/ or gromacs/
│       └── {lig1}~{lig2}/
│           ├── bound_{replica}/
│           └── free_{replica}/
├── analysis/                       # Analysis outputs
│   └── {engine}/
│       └── {lig1}~{lig2}/
│           ├── bound_{replica}/pmf.csv
│           ├── free_{replica}/pmf.csv
│           ├── rbfe_results.csv
│           └── rbfe_network.png
└── logs/                           # Snakemake and SLURM job logs
```

---

## Tips and Notes

**Lambda windows**: `lambda_windows` applies to perturbations scoring at or above `lomap_threshold`. `diff_lambda_windows` applies to those below. Use a larger value for `diff_lambda_windows` to improve convergence on difficult perturbations.

**GROMACS lambda schedules**: Each list in `lambda_schedules` must have the same length. The length determines the number of lambda windows for that leg. The bound leg includes a `bonded` schedule; the free leg does not (bonded terms do not change in the solvated ligand).

**Modified dummies**: `use-modified-dummies: true` enables corrections to dummy atom bonded terms via the [ghostly](https://github.com/OpenBioSim/ghostly) library. This is recommended for most perturbations and helps avoid artifacts from ghost atoms retaining their bonded interactions.

**SOMD2 runner**: `repex` performs replica exchange Monte Carlo moves between lambda windows. It allocates one OpenMM context per window at startup, so total GPU VRAM must be sufficient to hold all windows simultaneously. For single-GPU or memory-limited jobs, use `runner: standard`.

**SOMD2 lambda windows**: For RBFE with SOMD2, the number of lambda windows per perturbation is controlled by `network_settings.lambda_windows` and `network_settings.diff_lambda_windows`, not by `somd2-settings`. There is no `num_lambda` field in the SOMD2 block.

**HMR with GROMACS**: When using `timestep: 4fs`, set `hmr_factor: 3`. HMR is applied during the preparation step; production runs then use `hmr: false` internally. Setting `hmr_factor: 1` disables HMR.

**Replicas and error estimation**: Use at least 3 replicas (`num_replicas: 3`). The final reported uncertainty is the maximum of the propagated MBAR error and the standard deviation of the mean across replicas.

**Production engine**: Only `somd2` and `gromacs` are supported for RBFE production. `amber` is available for minimisation and equilibration stages only.
