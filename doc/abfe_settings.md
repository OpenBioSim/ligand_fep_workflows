# ABFE Workflow Settings

This document describes all configuration options for the Absolute Binding Free
Energy (ABFE) workflow. The config file is `config/config_abfe.yml` and is
validated against `config/schemas/config_abfe.schema.yml`. A test configuration
is provided at `config/config_abfe_test.yml`.

## Overview

The ABFE workflow calculates absolute binding free energies using alchemical
decoupling with Boresch restraints. Each ligand is simulated independently
on two legs:

**Bound leg** (protein-ligand complex): The ligand is held in place with Boresch
restraints while its electrostatic and van der Waals interactions are decoupled
from the environment in a single 21-window schedule.

**Free leg** (ligand in solvent): The same decoupling is performed without the
protein.

The binding free energy is:

```
DG_bind = DG_free - DG_bound + correction
```

where `correction` is the analytical standard-state correction for releasing
the Boresch restraints (typically around -10 to -12 kcal/mol). It is computed
automatically during the restraint search and stored alongside the restraint
parameters.

---

## Top-Level Fields

The `method` key is required and must be set to `abfe`. Without it the workflow
defaults to RBFE and validates against the wrong schema.

```yaml
# Must be present and set to abfe
method: abfe

# Glob pattern matching ligand SDF files
ligands: "workflow/inputs/ligands/*.sdf"

# Glob pattern matching parameterised protein input files.
# Expects a .rst7 (coordinates) and .prm7 (topology) from AmberTools.
protein_files: "workflow/inputs/protein/protein.*"

# Root directory for all workflow outputs
working_directory: "output/abfe"

# Number of CPU threads passed to GROMACS during minimisation and equilibration
simulation_threads: 5
```

---

## setup_settings

Controls ligand parameterisation and solvation for both legs.

```yaml
setup_settings:
  # Force field for ligand parameterisation.
  # Common choices: gaff2, gaff, openff_unconstrained-2.0.0
  ligand_forcefield: "gaff2"

  # Water model for solvation
  # Common choices: tip3p, tip4p, spce
  water_model: "tip3p"

  # Ion concentration in mol/L for charge neutralisation
  ion_concentration: 0.15

  # Padding added to the bounding box on each side, in nm
  box_length: 5.0

  # Simulation box shape
  # Options: cubic, triclinic, orthorhombic
  box_type: "cubic"
```

---

## restraint_search_settings

Controls the unrestrained simulation used to identify Boresch restraint
parameters for the bound leg. The restraint search runs an equilibration
simulation, then uses BioSimSpace to select anchor atoms and fit harmonic
parameters. The standard-state correction for releasing the restraints is
computed analytically and saved to `restraints/{ligand}/{ligand}_correction.txt`.

```yaml
restraint_search_settings:
  # Runtime for the restraint search simulation.
  # 1ns is sufficient for a test run; use 3-5ns for production.
  search_runtime: "1ns"

  # Temperature for the restraint search simulation
  temperature: "300K"

```

---

## min/eq-stages-free and min/eq-stages-bound

These sections define the minimisation and equilibration protocol applied to the
free-leg system (ligand in solvent) and the bound-leg system (protein-ligand
complex) respectively.

Each key must begin with `minimisation` or `equilibration` followed by any
alphanumeric characters or underscores (e.g. `minimisation1`, `equilibration2`,
`minimisation_norestraint`). Stages are run in the order they appear in the
file. Only `gromacs` is supported as the engine for ABFE min/eq stages.

### Minimisation stage keys

| Key | Required | Description |
|-----|----------|-------------|
| `engine` | yes | MD engine. Must be `gromacs`. |
| `minimisation-steps` | no | Number of energy minimisation steps. |
| `restraint-string` | no | Atom selection to restrain: `heavy`, `backbone`, or omit for no restraint. |
| `restraint-indices` | no | List of atom indices to restrain (alternative to `restraint-string`). |
| `restraint-force-constant` | no | Force constant in kcal/(mol·A^2). |

### Equilibration stage keys

| Key | Required | Description |
|-----|----------|-------------|
| `engine` | yes | MD engine. Must be `gromacs`. |
| `runtime` | yes | Simulation duration, e.g. `15ps`, `1ns`. |
| `timestep` | no | Integration timestep, e.g. `1fs`, `2fs`. |
| `temperature` | no | Target temperature, e.g. `300K`. Overrides `temperature-start`/`temperature-end`. |
| `temperature-start` | no | Starting temperature for a temperature ramp. |
| `temperature-end` | no | Ending temperature for a temperature ramp. |
| `thermostat-time-constant` | no | Thermostat coupling time, e.g. `1ps`, `0.5ps`. |
| `pressure` | no | Target pressure (enables NPT), e.g. `1bar`. |
| `restraint-string` | no | Atom selection to restrain: `heavy`, `backbone`, or omit for no restraint. |
| `restraint-indices` | no | List of atom indices to restrain. |
| `restraint-force-constant` | no | Force constant in kcal/(mol·A^2). |

### Free leg example

The free leg only contains the ligand in water. The protocol ramps down
positional restraints across several minimisation and NVT/NPT stages.

```yaml
min/eq-stages-free:
  minimisation1:
    engine: gromacs
    minimisation-steps: 1000
    restraint-string: heavy       # Restrain all non-hydrogen atoms
    restraint-force-constant: 5.0

  equilibration1:
    engine: gromacs
    timestep: 1fs
    runtime: 15ps
    temperature: 300K
    thermostat-time-constant: 0.5ps
    restraint-string: heavy
    restraint-force-constant: 5.0

  minimisation2:
    engine: gromacs
    minimisation-steps: 1000
    restraint-string: heavy
    restraint-force-constant: 2.0

  minimisation3:
    engine: gromacs
    minimisation-steps: 1000
    restraint-string: heavy
    restraint-force-constant: 0.1

  minimisation4:
    engine: gromacs
    minimisation-steps: 1000      # No restraints — free minimisation

  equilibration2:
    engine: gromacs
    timestep: 1fs
    runtime: 5ps
    temperature: 300K
    pressure: 1bar                # NPT from here
    thermostat-time-constant: 1ps
    restraint-string: heavy
    restraint-force-constant: 1.0

  equilibration3:
    engine: gromacs
    timestep: 1fs
    runtime: 5ps
    temperature: 300K
    pressure: 1bar
    thermostat-time-constant: 1ps
    restraint-string: heavy
    restraint-force-constant: 0.5

  equilibration4:
    engine: gromacs
    timestep: 2fs
    runtime: 10ps
    temperature: 300K
    pressure: 1bar
    thermostat-time-constant: 1ps  # No restraints — unrestrained NPT
```

### Bound leg example

The bound leg includes the protein. Extra equilibration stages with backbone
restraints are recommended to stabilise the protein structure before the
restraints are fully released.

```yaml
min/eq-stages-bound:
  minimisation1:
    engine: gromacs
    minimisation-steps: 1000
    restraint-string: heavy
    restraint-force-constant: 5.0

  equilibration1:
    engine: gromacs
    timestep: 1fs
    runtime: 15ps
    temperature: 300K
    thermostat-time-constant: 0.5ps
    restraint-string: heavy
    restraint-force-constant: 5.0

  minimisation2:
    engine: gromacs
    minimisation-steps: 1000
    restraint-string: heavy
    restraint-force-constant: 2.0

  minimisation3:
    engine: gromacs
    minimisation-steps: 1000
    restraint-string: heavy
    restraint-force-constant: 0.1

  minimisation4:
    engine: gromacs
    minimisation-steps: 1000

  equilibration2:
    engine: gromacs
    timestep: 1fs
    runtime: 5ps
    temperature: 300K
    pressure: 1bar
    thermostat-time-constant: 1ps
    restraint-string: heavy
    restraint-force-constant: 1.0

  equilibration3:
    engine: gromacs
    timestep: 1fs
    runtime: 5ps
    temperature: 300K
    pressure: 1bar
    thermostat-time-constant: 1ps
    restraint-string: heavy
    restraint-force-constant: 0.5

  equilibration4:
    engine: gromacs
    timestep: 1fs
    runtime: 10ps
    temperature: 300K
    pressure: 1bar
    thermostat-time-constant: 1ps
    restraint-string: backbone    # Switch to backbone-only restraints
    restraint-force-constant: 0.5

  equilibration5:
    engine: gromacs
    timestep: 2fs
    runtime: 10ps
    temperature: 300K
    pressure: 1bar
    thermostat-time-constant: 1ps # No restraints — final unrestrained NPT
```

---

## production-settings

Controls the alchemical free energy production simulations. Two engines are
supported: `somd2` and `gromacs`. The engine is specified at the top of this
section and the corresponding sub-section provides engine-specific settings.

```yaml
production-settings:
  # Number of independent replicas per ligand per leg.
  # The final binding free energy error is the maximum of the propagated
  # MBAR error and the standard deviation across replicas.
  # Use at least 3.
  num_replicas: 3

  # Production MD engine. Options: somd2, gromacs (case-insensitive)
  engine: somd2
```

### GROMACS production settings

When `engine: gromacs`, settings are read from the `gromacs-settings` block.
GROMACS production runs a short per-lambda NVT + NPT equilibration internally
before the production MD.

The lambda schedule is specified as three columns (`bonded`, `coul`, `vdw`)
for 21 windows. Windows 0-10 ramp bonded and coul together while vdw stays at
zero; windows 11-20 hold bonded and coul at 1.0 and ramp vdw. The free leg
does not have bonded restraint terms, so only `coul` and `vdw` are needed.

```yaml
  gromacs-settings:
    # Hydrogen mass repartitioning factor.
    # Use 3 with a 4fs timestep to allow stable integration.
    hmr_factor: 3

    lambda_schedules:
      bound:
        bonded: [0.0, 0.010, 0.025, 0.05, 0.075, 0.10, 0.15, 0.22, 0.35, 0.75, 1.0,
                 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        coul:   [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        vdw:    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
      free:
        coul:   [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        vdw:    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    free-leg-settings:
      runtime: 2ns
      timestep: 4fs          # Requires hmr_factor: 3
      temperature: 300K
      pressure: 1bar
      report-interval: 1ps   # Time string, converted to steps internally
      restart-interval: 500ps

    bound-leg-settings:
      runtime: 2ns
      timestep: 4fs
      temperature: 300K
      pressure: 1bar
      report-interval: 1ps
      restart-interval: 500ps
```

### SOMD2 production settings

When `engine: somd2`, settings are read from the `somd2-settings` block. SOMD2
uses OpenMM internally and applies HMR automatically; there is no `hmr_factor`
to set. Per-window equilibration is handled internally via `equilibration_time`.

The lambda schedule uses a decharge then annihilate approach. `num_lambda`
specifies the number of windows per stage, so the total number of lambda
evaluations is `2 * num_lambda`.

```yaml
  somd2-settings:
    runtime: 2ns
    timestep: 4fs
    temperature: 298K

    # Long-range electrostatics treatment
    # Options: PME (particle mesh Ewald), RF (reaction field)
    cutoff_type: PME
    cutoff: 10A

    # Windows per alchemical stage (decharge + annihilate).
    # Total lambda evaluations = 2 * num_lambda.
    num_lambda: 21

    # Frequency for writing energy data
    energy_frequency: 1ps

    # Frequency for writing trajectory frames
    frame_frequency: 500ps

    # Frequency for writing checkpoint files
    checkpoint_frequency: 500ps

    # Integrator for dynamics
    integrator: langevin_middle

    # Soft-core shift parameter for LJ interactions (Angstrom)
    shift_delta: 2.25A

    # Constraint scheme for perturbing hydrogen bonds
    perturbable_constraint: h_bonds_not_heavy_perturbed

    # Per-window equilibration time before production data collection begins.
    # This is handled entirely within SOMD2 — no separate equilibration rule runs.
    equilibration_time: 20ps

    # Sampling strategy.
    # repex: replica exchange — all lambda windows run simultaneously across
    #        multiple GPUs. Requires sufficient GPU VRAM for all contexts.
    # standard: independent simulations per window, suitable for single-GPU setups.
    runner: repex

    # Number of GPUs allocated per production job.
    # For repex, increase this to accommodate the total memory required by all
    # lambda windows across both stages.
    gpus_per_job: 3
```

---

## analysis-settings

```yaml
analysis-settings:
  # Plotting backend.
  # native: uses matplotlib, always available.
  # cinnabar: uses the cinnabar package (must be installed separately).
  backend: native

  # Optional: path to a CSV file of experimental binding affinities.
  # Format: columns named "ligand", "value", "error".
  # Uncomment to enable comparison plots in the analysis output.
  # experimental-results: workflow/inputs/exp_data.csv

  # Units of the experimental data (input only; outputs are always kcal/mol).
  # Options: kcal/mol, kJ/mol, Ki_uM
  # experimental-units: kcal/mol
```

---

## Running the Workflow

```bash
# Basic execution on a single machine
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml --cores 8

# With SLURM
snakemake -s workflow/Snakefile --configfile config/config_abfe.yml \
    --executor slurm --profile profiles/config.yaml \
    --jobs 100 --resources gpu=4
```

### Monitoring progress

```bash
# Show per-ligand stage completion and estimated free energies
snakemake status -s workflow/Snakefile --configfile config/config_abfe.yml --cores 1

# Show only the final binding free energy summary table
snakemake results -s workflow/Snakefile --configfile config/config_abfe.yml --cores 1
```

### Cleaning outputs

```bash
# Remove all outputs
snakemake clean -s workflow/Snakefile --configfile config/config_abfe.yml

# Remove only analysis outputs (keeps production simulations)
snakemake clean_analysis -s workflow/Snakefile --configfile config/config_abfe.yml

# Remove production outputs (keeps equilibration)
snakemake clean_production -s workflow/Snakefile --configfile config/config_abfe.yml
```

---

## Output Directory Structure

```
{working_directory}/
├── setup/                          # Parameterised and solvated systems
│   ├── {ligand}_free.bss
│   └── {ligand}_bound.bss
├── preparation/                    # Min/eq outputs
│   ├── minimisation/
│   │   └── {N}/                    # Stage number
│   ├── equilibration/
│   │   └── {N}/
│   └── final/                      # Equilibrated systems ready for production
│       ├── {ligand}_free.bss
│       └── {ligand}_bound.bss
├── restraints/                     # Boresch restraint search outputs
│   └── {ligand}/
│       ├── restraint_search_sim/   # Restraint search simulation files
│       ├── {ligand}_restraint.json # Selected restraint parameters
│       └── {ligand}_correction.txt # Standard-state correction (kcal/mol)
├── production/                     # Production simulation outputs
│   └── {engine}/                   # somd2/ or gromacs/
│       └── {ligand}/
│           ├── bound_{replica}/    # Replica index: 0, 1, 2, ...
│           └── free_{replica}/
├── analysis/                       # Analysis outputs
│   └── {engine}/
│       └── {ligand}/
│           ├── bound_{replica}/pmf.csv
│           ├── free_{replica}/pmf.csv
│           ├── detailed_abfe_results.csv
│           ├── final_abfe_results.csv
│           └── abfe_comparison.png
└── logs/                           # Snakemake and job log files
```

---

## Tips and Best Practices

1. **Restraint search duration**: The default `1ns` in the test config is too
   short for production work. Use at least 3-5 ns to ensure the ligand samples
   representative bound-state conformations before anchor atoms are selected.

2. **SOMD2 repex and GPU memory**: The `repex` runner loads all lambda windows
   simultaneously (one OpenMM context per window). For `num_lambda: 21`, this
   means 42 contexts total across both stages. Ensure that the combined VRAM
   across the allocated GPUs is sufficient. Use `runner: standard` for
   single-GPU setups or where memory is limited.

3. **HMR for GROMACS**: Set `hmr_factor: 3` when using a 4 fs timestep.
   SOMD2 manages HMR internally and does not need this setting.

4. **SOMD2 equilibration**: SOMD2 handles per-window equilibration internally
   via `equilibration_time`. No separate equilibration Snakemake rule is needed
   or run for SOMD2 production.

5. **GROMACS production equilibration**: When using GROMACS for production, a
   short NVT heat and NPT equilibration runs inside each production job before
   the production MD. This is separate from the global `min/eq-stages-bound`
   and `min/eq-stages-free` protocol.

6. **Lambda schedule design (GROMACS)**: Windows 0-10 ramp bonded and coul
   together while vdw stays at zero. Windows 11-20 hold bonded and coul at 1.0
   and ramp vdw from 0.1 to 1.0. This ordering avoids inserting a charged
   ligand into the environment and reduces the variance in dV/dlambda at the
   vdw end states.

7. **Number of replicas**: Use at least 3 replicas. The reported uncertainty
   is the maximum of the propagated MBAR errors and the standard deviation
   across replicas, so a minimum of 3 is needed to estimate the latter.

8. **Standard-state correction sign**: The correction stored in
   `{ligand}_correction.txt` is negative (approximately -10 to -12 kcal/mol
   for typical Boresch restraints). It is added to the calculated
   `DG_free - DG_bound` as shown in the formula above.
