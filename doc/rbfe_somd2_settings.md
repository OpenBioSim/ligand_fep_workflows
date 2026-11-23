# SOMD2 RBFE Settings

`method: rbfe_somd2` runs relative binding free energies with SOMD2 only, using the protocol from the OpenBioSim benchmarks. It differs from `method: rbfe` in how the end states are built:

- Each edge's ligands are mapped and merged first, then the merged molecule is solvated once per leg. Both end states therefore share one topology, including the solvent.
- Both end states of the merged system are equilibrated. SOMD2 starts the windows nearest each end state from that end state's structure (multiconformational seeding), so each edge only runs in one direction.
- Production uses replica exchange with adaptive allocation of sampling, REST2, and Monte Carlo flips of terminal rings. GCMC is available, but off by default.

The full configuration is in `config/config_rbfe_somd2.yml`, and a fast test version in `config/config_rbfe_somd2_test.yml`.

## Stages

| Rule | Output | Description |
|---|---|---|
| `prep_network` | `network/network.dat` | LOMAP network, one direction per edge. Skipped if `network.dat` already exists. |
| `parameterise_ligand` | `parameterised/{ligand}.bss` | Ligand parameters and charges from `ligand_settings.forcefield`. |
| `prepare_edge` | `edges/{edge}/{bound,free}.bss` | Atom mapping, merging, and solvation of each leg. The mapping is saved in `edge.json`. |
| `equilibrate_lambda_0` | `edges/{edge}/{leg}_lambda_0.bss` | Equilibration of the lambda = 0 end state. |
| `equilibrate_lambda_1` | `edges/{edge}/{leg}_lambda_1.bss` | Equilibration of the lambda = 1 end state, then matching to the lambda = 0 box. |
| `production` | `production/somd2/{edge}/{leg}_{replica}/` | SOMD2, with the lambda = 1 end state as `perturbed_system`. |
| `analyse_single_target`, `collate_rbfe_analysis` | `analysis/somd2/` | MBAR per leg, then ΔΔG per edge, as for `method: rbfe`. |

All replicas of every edge finish replica N before any starts replica N + 1, so that each edge has a result as early as possible. `snakemake status` and `snakemake results` work as for `method: rbfe`.

## Settings

### Inputs

| Key | Description |
|---|---|
| `ligands` | Glob pattern for the ligand SDF files, one ligand per file. |
| `protein_files` | Glob pattern for the parameterised protein files, e.g. prm7 and rst7. Crystal waters, ions and cofactors in them are kept. |
| `working_directory` | Directory for all outputs. |
| `simulation_threads` | CPU threads per job. |

### `ligand_settings`

| Key | Default | Description |
|---|---|---|
| `forcefield` | `openff_unconstrained-2.3.0` | Ligand force field, which also assigns the partial charges. |

### `mapping_settings`

| Key | Default | Description |
|---|---|---|
| `mcs_strategy` | `full` | BioSimSpace's MCS strategy: `full`, `heavy` (heavy atoms first, then hydrogens), or `both` (the larger mapping). Values other than `full` need a BioSimSpace with the `mcs_strategy` option. |
| `allow_ring_breaking` | `false` | Allow mappings that break rings. |
| `allow_ring_size_change` | `false` | Allow mappings that change ring sizes. |

If BioSimSpace warns that relaxing ring matching gives a larger common core, the edge is mapped again with it relaxed. Ghost atoms that land on top of another atom are moved slightly, which has no effect at the end state where they are ghosts.

### `setup_settings`

| Key | Default | Description |
|---|---|---|
| `water_model` | `tip3p` | Water model. |
| `ion_concentration` | `0.15` | NaCl concentration in mol/L. Systems are also neutralised. |
| `box_type` | `dodecahedron` | `cubic`, `dodecahedron` or `octahedron`. |
| `bound_padding` | `1.0 nm` | Minimum distance from the solute to the box edge in the bound leg. |
| `free_padding` | `1.5 nm` | The same for the free leg. |

### `equilibration_settings`

| Key | Default | Description |
|---|---|---|
| `ghost_modifications` | `true` | Apply Ghostly's modifications to ghost atoms during equilibration, as SOMD2 does in production. |
| `match_box` | `true` | Rescale the lambda = 1 box to the lambda = 0 one, then run the `match_box` stages at constant volume, so both end states share a box. Needed for GCMC. |
| `bound`, `free` | | The protocol for each leg: `defaults` for every stage, the `stages`, and the `match_box` stages. |

Each stage has a `type` (`minimisation` or `dynamics`) and optionally `runtime`, `timestep`, `temperature`, `pressure`, `restraint` (`heavy`, `backbone`, or a Sire selection) and `k`, the restraint constant in E = k r². A dynamics stage that becomes unstable runs briefly without constraints on the perturbable molecule, then is retried.

### `production-settings`

| Key | Default | Description |
|---|---|---|
| `num_replicas` | `3` | Independent replicas of each leg. |
| `gpus_per_job` | `1` | GPUs per job. Replica exchange spreads the windows over them. |
| `fallback_timestep` | `2 fs` | A leg that becomes unstable is rerun from scratch with this timestep, keeping the failed run as `{leg}_{replica}.unstable_<timestep>`. `null` disables it. |
| `somd2` | | SOMD2 options, passed to SOMD2 as they are. Any SOMD2 option can be added. |
| `leg-overrides` | `free: {gcmc: false}` | SOMD2 options for one leg only. |

The default SOMD2 options:

| Option | Value | Notes |
|---|---|---|
| `runtime` | `5 ns` | Per window. |
| `timestep` | `4 fs` | With SOMD2's default hydrogen mass repartitioning. |
| `temperature` | `298.15 K` | |
| `pressure` | `1 bar` | SOMD2 runs at constant volume instead when GCMC is on. |
| `num_lambda` | `11` | |
| `replica_exchange` | `true` | |
| `adaptive_repex` | `true` | Concentrates sampling on the windows that contribute most to the uncertainty. Needs SOMD2 2026.3.0 or later. |
| `num_energy_neighbours` | `3` | Energies only for the nearest windows. |
| `terminal_flip_frequency` | `1 ps` | |
| `rest2_scale` | `4.0` | REST2 on SOMD2's automatically selected region. |
| `rest2_profile` | `exponential` | |
| `gcmc` | `false` | Set to `true` for GCMC water sampling in the bound leg. |
| `gcmc_frequency` | `5 ps` | Used when GCMC is on. |
| `coalchemical_restraint_dist` | `10 A` | Keeps the co-alchemical ion of charge-changing edges away from the ligand. |
| `checkpoint_frequency`, `frame_frequency` | `500 ps` | |

An interrupted leg restarts from its checkpoint when the workflow is run again. To extend completed legs, increase `runtime`, delete the legs' `.done` files, and run the workflow again.

### `analysis-settings`

As for `method: rbfe`: `backend`, `experimental-results`, `experimental-units` and `temperature`, which converts experimental binding constants to free energies.
