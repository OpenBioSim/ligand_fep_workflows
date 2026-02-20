# RBFE Workflow Config Schema

- [1. Property `RBFE Workflow Config Schema > ligands`](#ligands)
- [2. Property `RBFE Workflow Config Schema > protein_files`](#protein_files)
- [3. Property `RBFE Workflow Config Schema > working_directory`](#working_directory)
- [4. Property `RBFE Workflow Config Schema > simulation_threads`](#simulation_threads)
- [5. Property `RBFE Workflow Config Schema > network_settings`](#network_settings)
  - [5.1. Property `RBFE Workflow Config Schema > network_settings > lambda_windows`](#network_settings_lambda_windows)
  - [5.2. Property `RBFE Workflow Config Schema > network_settings > lomap_threshold`](#network_settings_lomap_threshold)
  - [5.3. Property `RBFE Workflow Config Schema > network_settings > diff_lambda_windows`](#network_settings_diff_lambda_windows)
- [6. Property `RBFE Workflow Config Schema > setup_settings`](#setup_settings)
  - [6.1. Property `RBFE Workflow Config Schema > setup_settings > ligand_forcefield`](#setup_settings_ligand_forcefield)
  - [6.2. Property `RBFE Workflow Config Schema > setup_settings > water_model`](#setup_settings_water_model)
  - [6.3. Property `RBFE Workflow Config Schema > setup_settings > ion_concentration`](#setup_settings_ion_concentration)
  - [6.4. Property `RBFE Workflow Config Schema > setup_settings > box_length`](#setup_settings_box_length)
  - [6.5. Property `RBFE Workflow Config Schema > setup_settings > box_type`](#setup_settings_box_type)
- [7. Property `RBFE Workflow Config Schema > min/eq-stages-free`](#min/eq-stages-free)
  - [7.1. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$`](#min/eq-stages-free_pattern1)
    - [7.1.1. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > engine`](#min/eq-stages-free_pattern1_engine)
    - [7.1.2. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > minimisation-steps`](#min/eq-stages-free_pattern1_minimisation-steps)
    - [7.1.3. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-string`](#min/eq-stages-free_pattern1_restraint-string)
    - [7.1.4. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices`](#min/eq-stages-free_pattern1_restraint-indices)
      - [7.1.4.1. RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items](#min/eq-stages-free_pattern1_restraint-indices_items)
    - [7.1.5. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-force-constant`](#min/eq-stages-free_pattern1_restraint-force-constant)
  - [7.2. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$`](#min/eq-stages-free_pattern2)
    - [7.2.1. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > engine`](#min/eq-stages-free_pattern2_engine)
    - [7.2.2. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > runtime`](#min/eq-stages-free_pattern2_runtime)
    - [7.2.3. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > timestep`](#min/eq-stages-free_pattern2_timestep)
    - [7.2.4. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > temperature-start`](#min/eq-stages-free_pattern2_temperature-start)
    - [7.2.5. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > temperature-end`](#min/eq-stages-free_pattern2_temperature-end)
    - [7.2.6. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > temperature`](#min/eq-stages-free_pattern2_temperature)
    - [7.2.7. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > thermostat-time-constant`](#min/eq-stages-free_pattern2_thermostat-time-constant)
    - [7.2.8. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > pressure`](#min/eq-stages-free_pattern2_pressure)
    - [7.2.9. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-string`](#min/eq-stages-free_pattern2_restraint-string)
    - [7.2.10. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices`](#min/eq-stages-free_pattern2_restraint-indices)
      - [7.2.10.1. RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items](#min/eq-stages-free_pattern2_restraint-indices_items)
    - [7.2.11. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-force-constant`](#min/eq-stages-free_pattern2_restraint-force-constant)
- [8. Property `RBFE Workflow Config Schema > min/eq-stages-bound`](#min/eq-stages-bound)
  - [8.1. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$`](#min/eq-stages-bound_pattern1)
    - [8.1.1. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > engine`](#min/eq-stages-bound_pattern1_engine)
    - [8.1.2. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > minimisation-steps`](#min/eq-stages-bound_pattern1_minimisation-steps)
    - [8.1.3. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-string`](#min/eq-stages-bound_pattern1_restraint-string)
    - [8.1.4. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices`](#min/eq-stages-bound_pattern1_restraint-indices)
      - [8.1.4.1. RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items](#min/eq-stages-bound_pattern1_restraint-indices_items)
    - [8.1.5. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-force-constant`](#min/eq-stages-bound_pattern1_restraint-force-constant)
  - [8.2. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$`](#min/eq-stages-bound_pattern2)
    - [8.2.1. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > engine`](#min/eq-stages-bound_pattern2_engine)
    - [8.2.2. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > runtime`](#min/eq-stages-bound_pattern2_runtime)
    - [8.2.3. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > timestep`](#min/eq-stages-bound_pattern2_timestep)
    - [8.2.4. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > temperature-start`](#min/eq-stages-bound_pattern2_temperature-start)
    - [8.2.5. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > temperature-end`](#min/eq-stages-bound_pattern2_temperature-end)
    - [8.2.6. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > temperature`](#min/eq-stages-bound_pattern2_temperature)
    - [8.2.7. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > thermostat-time-constant`](#min/eq-stages-bound_pattern2_thermostat-time-constant)
    - [8.2.8. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > pressure`](#min/eq-stages-bound_pattern2_pressure)
    - [8.2.9. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-string`](#min/eq-stages-bound_pattern2_restraint-string)
    - [8.2.10. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices`](#min/eq-stages-bound_pattern2_restraint-indices)
      - [8.2.10.1. RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items](#min/eq-stages-bound_pattern2_restraint-indices_items)
    - [8.2.11. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-force-constant`](#min/eq-stages-bound_pattern2_restraint-force-constant)
- [9. Property `RBFE Workflow Config Schema > production-settings`](#production-settings)
  - [9.1. Property `RBFE Workflow Config Schema > production-settings > num_replicas`](#production-settings_num_replicas)
  - [9.2. Property `RBFE Workflow Config Schema > production-settings > engine`](#production-settings_engine)
  - [9.3. Property `RBFE Workflow Config Schema > production-settings > hmr_factor`](#production-settings_hmr_factor)
  - [9.4. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings`](#production-settings_free-leg-settings)
    - [9.4.1. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > runtime`](#production-settings_free-leg-settings_runtime)
    - [9.4.2. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > timestep`](#production-settings_free-leg-settings_timestep)
    - [9.4.3. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > temperature`](#production-settings_free-leg-settings_temperature)
    - [9.4.4. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > pressure`](#production-settings_free-leg-settings_pressure)
    - [9.4.5. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-string`](#production-settings_free-leg-settings_restraint-string)
    - [9.4.6. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-indices`](#production-settings_free-leg-settings_restraint-indices)
      - [9.4.6.1. RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-indices > restraint-indices items](#production-settings_free-leg-settings_restraint-indices_items)
    - [9.4.7. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-force-constant`](#production-settings_free-leg-settings_restraint-force-constant)
    - [9.4.8. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > use-modified-dummies`](#production-settings_free-leg-settings_use-modified-dummies)
    - [9.4.9. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > report-interval`](#production-settings_free-leg-settings_report-interval)
    - [9.4.10. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restart-interval`](#production-settings_free-leg-settings_restart-interval)
  - [9.5. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings`](#production-settings_bound-leg-settings)
    - [9.5.1. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > runtime`](#production-settings_bound-leg-settings_runtime)
    - [9.5.2. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > timestep`](#production-settings_bound-leg-settings_timestep)
    - [9.5.3. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > temperature`](#production-settings_bound-leg-settings_temperature)
    - [9.5.4. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > pressure`](#production-settings_bound-leg-settings_pressure)
    - [9.5.5. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-string`](#production-settings_bound-leg-settings_restraint-string)
    - [9.5.6. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-indices`](#production-settings_bound-leg-settings_restraint-indices)
      - [9.5.6.1. RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-indices > restraint-indices items](#production-settings_bound-leg-settings_restraint-indices_items)
    - [9.5.7. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-force-constant`](#production-settings_bound-leg-settings_restraint-force-constant)
    - [9.5.8. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > use-modified-dummies`](#production-settings_bound-leg-settings_use-modified-dummies)
    - [9.5.9. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > report-interval`](#production-settings_bound-leg-settings_report-interval)
    - [9.5.10. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restart-interval`](#production-settings_bound-leg-settings_restart-interval)
- [10. Property `RBFE Workflow Config Schema > analysis-settings`](#analysis-settings)
  - [10.1. Property `RBFE Workflow Config Schema > analysis-settings > backend`](#analysis-settings_backend)
  - [10.2. Property `RBFE Workflow Config Schema > analysis-settings > experimental-results`](#analysis-settings_experimental-results)
  - [10.3. Property `RBFE Workflow Config Schema > analysis-settings > experimental-units`](#analysis-settings_experimental-units)

**Title:** RBFE Workflow Config Schema

|                           |                  |
| ------------------------- | ---------------- |
| **Type**                  | `object`         |
| **Required**              | No               |
| **Additional properties** | Any type allowed |

**Description:** Defines and validates the allowed stages and settings for the RBFE molecular dynamics workflow.

| Property                                       | Pattern | Type    | Deprecated | Definition | Title/Description                                                                                                                                          |
| ---------------------------------------------- | ------- | ------- | ---------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------- |
| + [ligands](#ligands )                         | No      | string  | No         | -          | Glob pattern for ligand input files.                                                                                                                       |
| + [protein_files](#protein_files )             | No      | string  | No         | -          | Glob pattern for protein input files. Assumes that this is coordinates and topology only, do not try to load multiple types of each, or multiple proteins. |
| + [working_directory](#working_directory )     | No      | string  | No         | -          | Working directory for the workflow. This will contain all outputs for all simulations, as well as analysis.                                                |
| + [simulation_threads](#simulation_threads )   | No      | integer | No         | -          | Number of threads to use for simulations.                                                                                                                  |
| + [network_settings](#network_settings )       | No      | object  | No         | -          | -                                                                                                                                                          |
| + [setup_settings](#setup_settings )       | No      | object  | No         | -          | -                                                                                                                                                          |
| + [min/eq-stages-free](#min/eq-stages-free )   | No      | object  | No         | -          | A dictionary of simulation stages. Keys must start with 'minimisation' or 'equilibration'.                                                                 |
| + [min/eq-stages-bound](#min/eq-stages-bound ) | No      | object  | No         | -          | A dictionary of simulation stages. Keys must start with 'minimisation' or 'equilibration'.                                                                 |
| + [production-settings](#production-settings ) | No      | object  | No         | -          | -                                                                                                                                                          |
| + [analysis-settings](#analysis-settings )     | No      | object  | No         | -          | -                                                                                                                                                          |
| - [](#additionalProperties )                   | No      | object  | No         | -          | -                                                                                                                                                          |

## <a name="ligands"></a>1. Property `RBFE Workflow Config Schema > ligands`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | Yes      |

**Description:** Glob pattern for ligand input files.

## <a name="protein_files"></a>2. Property `RBFE Workflow Config Schema > protein_files`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | Yes      |

**Description:** Glob pattern for protein input files. Assumes that this is coordinates and topology only, do not try to load multiple types of each, or multiple proteins.

## <a name="working_directory"></a>3. Property `RBFE Workflow Config Schema > working_directory`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | Yes      |

**Description:** Working directory for the workflow. This will contain all outputs for all simulations, as well as analysis.

## <a name="simulation_threads"></a>4. Property `RBFE Workflow Config Schema > simulation_threads`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | Yes       |

**Description:** Number of threads to use for simulations.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

## <a name="network_settings"></a>5. Property `RBFE Workflow Config Schema > network_settings`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

| Property                                                        | Pattern | Type    | Deprecated | Definition | Title/Description                                                                  |
| --------------------------------------------------------------- | ------- | ------- | ---------- | ---------- | ---------------------------------------------------------------------------------- |
| - [lambda_windows](#network_settings_lambda_windows )           | No      | integer | No         | -          | Number of windows for regular transformations (those above the LOMAP threshold).   |
| - [lomap_threshold](#network_settings_lomap_threshold )         | No      | number  | No         | -          | LOMAP threshold that defines difficult transformations.                            |
| - [diff_lambda_windows](#network_settings_diff_lambda_windows ) | No      | integer | No         | -          | Number of windows for difficult transformations (those below the LOMAP threshold). |

### <a name="network_settings_lambda_windows"></a>5.1. Property `RBFE Workflow Config Schema > network_settings > lambda_windows`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Number of windows for regular transformations (those above the LOMAP threshold).

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

### <a name="network_settings_lomap_threshold"></a>5.2. Property `RBFE Workflow Config Schema > network_settings > lomap_threshold`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** LOMAP threshold that defines difficult transformations.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |
| **Maximum**  | &le; 1 |

### <a name="network_settings_diff_lambda_windows"></a>5.3. Property `RBFE Workflow Config Schema > network_settings > diff_lambda_windows`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Number of windows for difficult transformations (those below the LOMAP threshold).

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

## <a name="setup_settings"></a>6. Property `RBFE Workflow Config Schema > setup_settings`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

| Property                                                    | Pattern | Type             | Deprecated | Definition | Title/Description                                               |
| ----------------------------------------------------------- | ------- | ---------------- | ---------- | ---------- | --------------------------------------------------------------- |
| - [ligand_forcefield](#setup_settings_ligand_forcefield ) | No      | string           | No         | -          | Force field to use for ligand parameterisation (e.g., 'gaff2'). |
| - [water_model](#setup_settings_water_model )             | No      | string           | No         | -          | Water model to use (e.g., 'tip3p').                             |
| - [ion_concentration](#setup_settings_ion_concentration ) | No      | number           | No         | -          | Ion concentration in mol/L.                                     |
| - [box_length](#setup_settings_box_length )               | No      | number           | No         | -          | Box length in nm.                                               |
| - [box_type](#setup_settings_box_type )                   | No      | enum (of string) | No         | -          | Type of box to use for solvation.                               |

### <a name="setup_settings_ligand_forcefield"></a>6.1. Property `RBFE Workflow Config Schema > setup_settings > ligand_forcefield`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Force field to use for ligand parameterisation (e.g., 'gaff2').

### <a name="setup_settings_water_model"></a>6.2. Property `RBFE Workflow Config Schema > setup_settings > water_model`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Water model to use (e.g., 'tip3p').

### <a name="setup_settings_ion_concentration"></a>6.3. Property `RBFE Workflow Config Schema > setup_settings > ion_concentration`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Ion concentration in mol/L.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

### <a name="setup_settings_box_length"></a>6.4. Property `RBFE Workflow Config Schema > setup_settings > box_length`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Box length in nm.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

### <a name="setup_settings_box_type"></a>6.5. Property `RBFE Workflow Config Schema > setup_settings > box_type`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | No                 |

**Description:** Type of box to use for solvation.

Must be one of:
* "cubic"
* "triclinic"
* "orthorhombic"

## <a name="min/eq-stages-free"></a>7. Property `RBFE Workflow Config Schema > min/eq-stages-free`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

**Description:** A dictionary of simulation stages. Keys must start with 'minimisation' or 'equilibration'.

| Property                                                        | Pattern | Type   | Deprecated | Definition | Title/Description |
| --------------------------------------------------------------- | ------- | ------ | ---------- | ---------- | ----------------- |
| - [^minimisation[a-zA-Z0-9_]*$](#min/eq-stages-free_pattern1 )  | Yes     | object | No         | -          | -                 |
| - [^equilibration[a-zA-Z0-9_]*$](#min/eq-stages-free_pattern2 ) | Yes     | object | No         | -          | -                 |

### <a name="min/eq-stages-free_pattern1"></a>7.1. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$`
> All properties whose name matches the regular expression
```^minimisation[a-zA-Z0-9_]*$``` ([Test](https://regex101.com/?regex=%5Eminimisation%5Ba-zA-Z0-9_%5D%2A%24))
must respect the following conditions

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | No          |
| **Additional properties** | Not allowed |

| Property                                                                             | Pattern | Type             | Deprecated | Definition | Title/Description                                                            |
| ------------------------------------------------------------------------------------ | ------- | ---------------- | ---------- | ---------- | ---------------------------------------------------------------------------- |
| + [engine](#min/eq-stages-free_pattern1_engine )                                     | No      | enum (of string) | No         | -          | The MD engine to use (Amber, Gromacs or OpenMM).                             |
| - [minimisation-steps](#min/eq-stages-free_pattern1_minimisation-steps )             | No      | integer          | No         | -          | Number of minimisation steps.                                                |
| - [restraint-string](#min/eq-stages-free_pattern1_restraint-string )                 | No      | string           | No         | -          | The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone'). |
| - [restraint-indices](#min/eq-stages-free_pattern1_restraint-indices )               | No      | array of integer | No         | -          | Indices of atoms to apply the restraint to.                                  |
| - [restraint-force-constant](#min/eq-stages-free_pattern1_restraint-force-constant ) | No      | number           | No         | -          | Force constant for the restraint, kcal_per_mol / angstrom^2.                 |

#### <a name="min/eq-stages-free_pattern1_engine"></a>7.1.1. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > engine`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | Yes                |

**Description:** The MD engine to use (Amber, Gromacs or OpenMM).

Must be one of:
* "OpenMM"
* "Amber"
* "Gromacs"

#### <a name="min/eq-stages-free_pattern1_minimisation-steps"></a>7.1.2. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > minimisation-steps`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Number of minimisation steps.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

#### <a name="min/eq-stages-free_pattern1_restraint-string"></a>7.1.3. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-string`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').

#### <a name="min/eq-stages-free_pattern1_restraint-indices"></a>7.1.4. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `array of integer` |
| **Required** | No                 |

**Description:** Indices of atoms to apply the restraint to.

|                      | Array restrictions |
| -------------------- | ------------------ |
| **Min items**        | N/A                |
| **Max items**        | N/A                |
| **Items unicity**    | False              |
| **Additional items** | False              |
| **Tuple validation** | See below          |

| Each item of this array must be                                                 | Description |
| ------------------------------------------------------------------------------- | ----------- |
| [restraint-indices items](#min/eq-stages-free_pattern1_restraint-indices_items) | -           |

##### <a name="min/eq-stages-free_pattern1_restraint-indices_items"></a>7.1.4.1. RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="min/eq-stages-free_pattern1_restraint-force-constant"></a>7.1.5. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^minimisation[a-zA-Z0-9_]*$ > restraint-force-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Force constant for the restraint, kcal_per_mol / angstrom^2.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

### <a name="min/eq-stages-free_pattern2"></a>7.2. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$`
> All properties whose name matches the regular expression
```^equilibration[a-zA-Z0-9_]*$``` ([Test](https://regex101.com/?regex=%5Eequilibration%5Ba-zA-Z0-9_%5D%2A%24))
must respect the following conditions

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | No          |
| **Additional properties** | Not allowed |

| Property                                                                             | Pattern | Type             | Deprecated | Definition | Title/Description                                                                          |
| ------------------------------------------------------------------------------------ | ------- | ---------------- | ---------- | ---------- | ------------------------------------------------------------------------------------------ |
| + [engine](#min/eq-stages-free_pattern2_engine )                                     | No      | enum (of string) | No         | -          | The MD engine to use (Amber, Gromacs or OpenMM).                                           |
| + [runtime](#min/eq-stages-free_pattern2_runtime )                                   | No      | string           | No         | -          | Duration of the simulation (e.g., '10ps', '100ns', '1.5ps').                               |
| - [timestep](#min/eq-stages-free_pattern2_timestep )                                 | No      | string           | No         | -          | Time step for the simulation (e.g., '2fs', '2.5fs').                                       |
| - [temperature-start](#min/eq-stages-free_pattern2_temperature-start )               | No      | string           | No         | -          | Starting temperature for equilibration (e.g., '300K', '300.5K').                           |
| - [temperature-end](#min/eq-stages-free_pattern2_temperature-end )                   | No      | string           | No         | -          | Ending temperature for equilibration (e.g., '300K', '300.5K').                             |
| - [temperature](#min/eq-stages-free_pattern2_temperature )                           | No      | string           | No         | -          | Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures. |
| - [thermostat-time-constant](#min/eq-stages-free_pattern2_thermostat-time-constant ) | No      | string           | No         | -          | Thermostat time constant for the equilibration (e.g., '1ps', '0.5ps').                     |
| - [pressure](#min/eq-stages-free_pattern2_pressure )                                 | No      | string           | No         | -          | Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.       |
| - [restraint-string](#min/eq-stages-free_pattern2_restraint-string )                 | No      | string           | No         | -          | The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').               |
| - [restraint-indices](#min/eq-stages-free_pattern2_restraint-indices )               | No      | array of integer | No         | -          | Indices of atoms to apply the restraint to.                                                |
| - [restraint-force-constant](#min/eq-stages-free_pattern2_restraint-force-constant ) | No      | number           | No         | -          | Force constant for the restraint, kcal_per_mol / angstrom^2.                               |

#### <a name="min/eq-stages-free_pattern2_engine"></a>7.2.1. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > engine`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | Yes                |

**Description:** The MD engine to use (Amber, Gromacs or OpenMM).

Must be one of:
* "OpenMM"
* "Amber"
* "Gromacs"

#### <a name="min/eq-stages-free_pattern2_runtime"></a>7.2.2. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > runtime`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | Yes      |

**Description:** Duration of the simulation (e.g., '10ps', '100ns', '1.5ps').

| Restrictions                      |                                                                                                                                |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(ps\|ns)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28ps%7Cns%29%24) |

#### <a name="min/eq-stages-free_pattern2_timestep"></a>7.2.3. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > timestep`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Time step for the simulation (e.g., '2fs', '2.5fs').

| Restrictions                      |                                                                                                                       |
| --------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(fs)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28fs%29%24) |

#### <a name="min/eq-stages-free_pattern2_temperature-start"></a>7.2.4. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > temperature-start`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Starting temperature for equilibration (e.g., '300K', '300.5K').

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="min/eq-stages-free_pattern2_temperature-end"></a>7.2.5. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > temperature-end`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Ending temperature for equilibration (e.g., '300K', '300.5K').

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="min/eq-stages-free_pattern2_temperature"></a>7.2.6. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > temperature`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures.

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="min/eq-stages-free_pattern2_thermostat-time-constant"></a>7.2.7. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > thermostat-time-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Thermostat time constant for the equilibration (e.g., '1ps', '0.5ps').

| Restrictions                      |                                                                                                                                |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(ps\|ns)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28ps%7Cns%29%24) |

#### <a name="min/eq-stages-free_pattern2_pressure"></a>7.2.8. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > pressure`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.

| Restrictions                      |                                                                                                                                    |
| --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(bar\|atm)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28bar%7Catm%29%24) |

#### <a name="min/eq-stages-free_pattern2_restraint-string"></a>7.2.9. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-string`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').

#### <a name="min/eq-stages-free_pattern2_restraint-indices"></a>7.2.10. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `array of integer` |
| **Required** | No                 |

**Description:** Indices of atoms to apply the restraint to.

|                      | Array restrictions |
| -------------------- | ------------------ |
| **Min items**        | N/A                |
| **Max items**        | N/A                |
| **Items unicity**    | False              |
| **Additional items** | False              |
| **Tuple validation** | See below          |

| Each item of this array must be                                                 | Description |
| ------------------------------------------------------------------------------- | ----------- |
| [restraint-indices items](#min/eq-stages-free_pattern2_restraint-indices_items) | -           |

##### <a name="min/eq-stages-free_pattern2_restraint-indices_items"></a>7.2.10.1. RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="min/eq-stages-free_pattern2_restraint-force-constant"></a>7.2.11. Property `RBFE Workflow Config Schema > min/eq-stages-free > ^equilibration[a-zA-Z0-9_]*$ > restraint-force-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Force constant for the restraint, kcal_per_mol / angstrom^2.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

## <a name="min/eq-stages-bound"></a>8. Property `RBFE Workflow Config Schema > min/eq-stages-bound`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

**Description:** A dictionary of simulation stages. Keys must start with 'minimisation' or 'equilibration'.

| Property                                                         | Pattern | Type   | Deprecated | Definition | Title/Description |
| ---------------------------------------------------------------- | ------- | ------ | ---------- | ---------- | ----------------- |
| - [^minimisation[a-zA-Z0-9_]*$](#min/eq-stages-bound_pattern1 )  | Yes     | object | No         | -          | -                 |
| - [^equilibration[a-zA-Z0-9_]*$](#min/eq-stages-bound_pattern2 ) | Yes     | object | No         | -          | -                 |

### <a name="min/eq-stages-bound_pattern1"></a>8.1. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$`
> All properties whose name matches the regular expression
```^minimisation[a-zA-Z0-9_]*$``` ([Test](https://regex101.com/?regex=%5Eminimisation%5Ba-zA-Z0-9_%5D%2A%24))
must respect the following conditions

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | No          |
| **Additional properties** | Not allowed |

| Property                                                                              | Pattern | Type             | Deprecated | Definition | Title/Description                                                            |
| ------------------------------------------------------------------------------------- | ------- | ---------------- | ---------- | ---------- | ---------------------------------------------------------------------------- |
| + [engine](#min/eq-stages-bound_pattern1_engine )                                     | No      | enum (of string) | No         | -          | The MD engine to use (Amber, Gromacs or OpenMM).                             |
| - [minimisation-steps](#min/eq-stages-bound_pattern1_minimisation-steps )             | No      | integer          | No         | -          | Number of minimisation steps.                                                |
| - [restraint-string](#min/eq-stages-bound_pattern1_restraint-string )                 | No      | string           | No         | -          | The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone'). |
| - [restraint-indices](#min/eq-stages-bound_pattern1_restraint-indices )               | No      | array of integer | No         | -          | Indices of atoms to apply the restraint to.                                  |
| - [restraint-force-constant](#min/eq-stages-bound_pattern1_restraint-force-constant ) | No      | number           | No         | -          | Force constant for the restraint, kcal_per_mol / angstrom^2.                 |

#### <a name="min/eq-stages-bound_pattern1_engine"></a>8.1.1. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > engine`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | Yes                |

**Description:** The MD engine to use (Amber, Gromacs or OpenMM).

Must be one of:
* "OpenMM"
* "Amber"
* "Gromacs"

#### <a name="min/eq-stages-bound_pattern1_minimisation-steps"></a>8.1.2. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > minimisation-steps`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Number of minimisation steps.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

#### <a name="min/eq-stages-bound_pattern1_restraint-string"></a>8.1.3. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-string`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').

#### <a name="min/eq-stages-bound_pattern1_restraint-indices"></a>8.1.4. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `array of integer` |
| **Required** | No                 |

**Description:** Indices of atoms to apply the restraint to.

|                      | Array restrictions |
| -------------------- | ------------------ |
| **Min items**        | N/A                |
| **Max items**        | N/A                |
| **Items unicity**    | False              |
| **Additional items** | False              |
| **Tuple validation** | See below          |

| Each item of this array must be                                                  | Description |
| -------------------------------------------------------------------------------- | ----------- |
| [restraint-indices items](#min/eq-stages-bound_pattern1_restraint-indices_items) | -           |

##### <a name="min/eq-stages-bound_pattern1_restraint-indices_items"></a>8.1.4.1. RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="min/eq-stages-bound_pattern1_restraint-force-constant"></a>8.1.5. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^minimisation[a-zA-Z0-9_]*$ > restraint-force-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Force constant for the restraint, kcal_per_mol / angstrom^2.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

### <a name="min/eq-stages-bound_pattern2"></a>8.2. Pattern Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$`
> All properties whose name matches the regular expression
```^equilibration[a-zA-Z0-9_]*$``` ([Test](https://regex101.com/?regex=%5Eequilibration%5Ba-zA-Z0-9_%5D%2A%24))
must respect the following conditions

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | No          |
| **Additional properties** | Not allowed |

| Property                                                                              | Pattern | Type             | Deprecated | Definition | Title/Description                                                                          |
| ------------------------------------------------------------------------------------- | ------- | ---------------- | ---------- | ---------- | ------------------------------------------------------------------------------------------ |
| + [engine](#min/eq-stages-bound_pattern2_engine )                                     | No      | enum (of string) | No         | -          | The MD engine to use (Amber, Gromacs or OpenMM).                                           |
| + [runtime](#min/eq-stages-bound_pattern2_runtime )                                   | No      | string           | No         | -          | Duration of the simulation (e.g., '10ps', '100ns', '1.5ps').                               |
| - [timestep](#min/eq-stages-bound_pattern2_timestep )                                 | No      | string           | No         | -          | Time step for the simulation (e.g., '2fs', '2.5fs').                                       |
| - [temperature-start](#min/eq-stages-bound_pattern2_temperature-start )               | No      | string           | No         | -          | Starting temperature for equilibration (e.g., '300K', '300.5K').                           |
| - [temperature-end](#min/eq-stages-bound_pattern2_temperature-end )                   | No      | string           | No         | -          | Ending temperature for equilibration (e.g., '300K', '300.5K').                             |
| - [temperature](#min/eq-stages-bound_pattern2_temperature )                           | No      | string           | No         | -          | Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures. |
| - [thermostat-time-constant](#min/eq-stages-bound_pattern2_thermostat-time-constant ) | No      | string           | No         | -          | Thermostat time constant for the equilibration (e.g., '1ps', '0.5ps').                     |
| - [pressure](#min/eq-stages-bound_pattern2_pressure )                                 | No      | string           | No         | -          | Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.       |
| - [restraint-string](#min/eq-stages-bound_pattern2_restraint-string )                 | No      | string           | No         | -          | The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').               |
| - [restraint-indices](#min/eq-stages-bound_pattern2_restraint-indices )               | No      | array of integer | No         | -          | Indices of atoms to apply the restraint to.                                                |
| - [restraint-force-constant](#min/eq-stages-bound_pattern2_restraint-force-constant ) | No      | number           | No         | -          | Force constant for the restraint, kcal_per_mol / angstrom^2.                               |

#### <a name="min/eq-stages-bound_pattern2_engine"></a>8.2.1. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > engine`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | Yes                |

**Description:** The MD engine to use (Amber, Gromacs or OpenMM).

Must be one of:
* "OpenMM"
* "Amber"
* "Gromacs"

#### <a name="min/eq-stages-bound_pattern2_runtime"></a>8.2.2. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > runtime`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | Yes      |

**Description:** Duration of the simulation (e.g., '10ps', '100ns', '1.5ps').

| Restrictions                      |                                                                                                                                |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(ps\|ns)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28ps%7Cns%29%24) |

#### <a name="min/eq-stages-bound_pattern2_timestep"></a>8.2.3. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > timestep`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Time step for the simulation (e.g., '2fs', '2.5fs').

| Restrictions                      |                                                                                                                       |
| --------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(fs)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28fs%29%24) |

#### <a name="min/eq-stages-bound_pattern2_temperature-start"></a>8.2.4. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > temperature-start`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Starting temperature for equilibration (e.g., '300K', '300.5K').

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="min/eq-stages-bound_pattern2_temperature-end"></a>8.2.5. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > temperature-end`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Ending temperature for equilibration (e.g., '300K', '300.5K').

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="min/eq-stages-bound_pattern2_temperature"></a>8.2.6. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > temperature`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures.

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="min/eq-stages-bound_pattern2_thermostat-time-constant"></a>8.2.7. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > thermostat-time-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Thermostat time constant for the equilibration (e.g., '1ps', '0.5ps').

| Restrictions                      |                                                                                                                                |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(ps\|ns)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28ps%7Cns%29%24) |

#### <a name="min/eq-stages-bound_pattern2_pressure"></a>8.2.8. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > pressure`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.

| Restrictions                      |                                                                                                                                    |
| --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(bar\|atm)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28bar%7Catm%29%24) |

#### <a name="min/eq-stages-bound_pattern2_restraint-string"></a>8.2.9. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-string`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').

#### <a name="min/eq-stages-bound_pattern2_restraint-indices"></a>8.2.10. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `array of integer` |
| **Required** | No                 |

**Description:** Indices of atoms to apply the restraint to.

|                      | Array restrictions |
| -------------------- | ------------------ |
| **Min items**        | N/A                |
| **Max items**        | N/A                |
| **Items unicity**    | False              |
| **Additional items** | False              |
| **Tuple validation** | See below          |

| Each item of this array must be                                                  | Description |
| -------------------------------------------------------------------------------- | ----------- |
| [restraint-indices items](#min/eq-stages-bound_pattern2_restraint-indices_items) | -           |

##### <a name="min/eq-stages-bound_pattern2_restraint-indices_items"></a>8.2.10.1. RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-indices > restraint-indices items

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="min/eq-stages-bound_pattern2_restraint-force-constant"></a>8.2.11. Property `RBFE Workflow Config Schema > min/eq-stages-bound > ^equilibration[a-zA-Z0-9_]*$ > restraint-force-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Force constant for the restraint, kcal_per_mol / angstrom^2.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

## <a name="production-settings"></a>9. Property `RBFE Workflow Config Schema > production-settings`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

| Property                                                         | Pattern | Type             | Deprecated | Definition | Title/Description                                   |
| ---------------------------------------------------------------- | ------- | ---------------- | ---------- | ---------- | --------------------------------------------------- |
| + [num_replicas](#production-settings_num_replicas )             | No      | integer          | No         | -          | Number of replicas for each perturbation.           |
| + [engine](#production-settings_engine )                         | No      | enum (of string) | No         | -          | The MD engine to use (Amber, Gromacs or OpenMM).    |
| - [hmr_factor](#production-settings_hmr_factor )                 | No      | number           | No         | -          | Hydrogen mass repartitioning factor, typically 1.5. |
| + [free-leg-settings](#production-settings_free-leg-settings )   | No      | object           | No         | -          | -                                                   |
| + [bound-leg-settings](#production-settings_bound-leg-settings ) | No      | object           | No         | -          | -                                                   |

### <a name="production-settings_num_replicas"></a>9.1. Property `RBFE Workflow Config Schema > production-settings > num_replicas`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | Yes       |

**Description:** Number of replicas for each perturbation.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

### <a name="production-settings_engine"></a>9.2. Property `RBFE Workflow Config Schema > production-settings > engine`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | Yes                |

**Description:** The MD engine to use (Amber, Gromacs or OpenMM).

Must be one of:
* "Somd"
* "Amber"
* "Gromacs"
* "Somd2"

### <a name="production-settings_hmr_factor"></a>9.3. Property `RBFE Workflow Config Schema > production-settings > hmr_factor`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Hydrogen mass repartitioning factor, typically 1.5.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

### <a name="production-settings_free-leg-settings"></a>9.4. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

| Property                                                                                       | Pattern | Type             | Deprecated | Definition | Title/Description                                                                                                                                                            |
| ---------------------------------------------------------------------------------------------- | ------- | ---------------- | ---------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| + [runtime](#production-settings_free-leg-settings_runtime )                                   | No      | string           | No         | -          | Duration of the simulation (e.g., '2ns').                                                                                                                                    |
| - [timestep](#production-settings_free-leg-settings_timestep )                                 | No      | string           | No         | -          | Time step for the simulation (e.g., '4fs').                                                                                                                                  |
| - [temperature](#production-settings_free-leg-settings_temperature )                           | No      | string           | No         | -          | Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures.                                                                                   |
| - [pressure](#production-settings_free-leg-settings_pressure )                                 | No      | string           | No         | -          | Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.                                                                                         |
| - [restraint-string](#production-settings_free-leg-settings_restraint-string )                 | No      | string           | No         | -          | The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').                                                                                                 |
| - [restraint-indices](#production-settings_free-leg-settings_restraint-indices )               | No      | array of integer | No         | -          | Indices of atoms to apply the restraint to.                                                                                                                                  |
| - [restraint-force-constant](#production-settings_free-leg-settings_restraint-force-constant ) | No      | number           | No         | -          | Force constant for the restraint, kcal_per_mol / angstrom^2.                                                                                                                 |
| - [use-modified-dummies](#production-settings_free-leg-settings_use-modified-dummies )         | No      | boolean          | No         | -          | Whether to apply modifications described here: https://doi.org/10.1021/acs.jctc.0c01328 to dummy atoms. Requires the Ghostly package: https://github.com/OpenBioSim/ghostly. |
| - [report-interval](#production-settings_free-leg-settings_report-interval )                   | No      | integer          | No         | -          | Interval for reporting simulation data, in steps (e.g., 250 for 1ps if timestep is 4fs).                                                                                     |
| - [restart-interval](#production-settings_free-leg-settings_restart-interval )                 | No      | integer          | No         | -          | Interval for restarting the simulation, in steps (e.g., 25000 for 100ps if timestep is 4fs).                                                                                 |

#### <a name="production-settings_free-leg-settings_runtime"></a>9.4.1. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > runtime`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | Yes      |

**Description:** Duration of the simulation (e.g., '2ns').

| Restrictions                      |                                                                                            |
| --------------------------------- | ------------------------------------------------------------------------------------------ |
| **Must match regular expression** | ```^[0-9]+(ps\|ns)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28ps%7Cns%29%24) |

#### <a name="production-settings_free-leg-settings_timestep"></a>9.4.2. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > timestep`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Time step for the simulation (e.g., '4fs').

| Restrictions                      |                                                                                   |
| --------------------------------- | --------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(fs)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28fs%29%24) |

#### <a name="production-settings_free-leg-settings_temperature"></a>9.4.3. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > temperature`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures.

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="production-settings_free-leg-settings_pressure"></a>9.4.4. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > pressure`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.

| Restrictions                      |                                                                                                                                    |
| --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(bar\|atm)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28bar%7Catm%29%24) |

#### <a name="production-settings_free-leg-settings_restraint-string"></a>9.4.5. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-string`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').

#### <a name="production-settings_free-leg-settings_restraint-indices"></a>9.4.6. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-indices`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `array of integer` |
| **Required** | No                 |

**Description:** Indices of atoms to apply the restraint to.

|                      | Array restrictions |
| -------------------- | ------------------ |
| **Min items**        | N/A                |
| **Max items**        | N/A                |
| **Items unicity**    | False              |
| **Additional items** | False              |
| **Tuple validation** | See below          |

| Each item of this array must be                                                           | Description |
| ----------------------------------------------------------------------------------------- | ----------- |
| [restraint-indices items](#production-settings_free-leg-settings_restraint-indices_items) | -           |

##### <a name="production-settings_free-leg-settings_restraint-indices_items"></a>9.4.6.1. RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-indices > restraint-indices items

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="production-settings_free-leg-settings_restraint-force-constant"></a>9.4.7. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restraint-force-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Force constant for the restraint, kcal_per_mol / angstrom^2.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="production-settings_free-leg-settings_use-modified-dummies"></a>9.4.8. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > use-modified-dummies`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |
| **Default**  | `false`   |

**Description:** Whether to apply modifications described here: https://doi.org/10.1021/acs.jctc.0c01328 to dummy atoms. Requires the Ghostly package: https://github.com/OpenBioSim/ghostly.

#### <a name="production-settings_free-leg-settings_report-interval"></a>9.4.9. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > report-interval`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Interval for reporting simulation data, in steps (e.g., 250 for 1ps if timestep is 4fs).

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

#### <a name="production-settings_free-leg-settings_restart-interval"></a>9.4.10. Property `RBFE Workflow Config Schema > production-settings > free-leg-settings > restart-interval`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Interval for restarting the simulation, in steps (e.g., 25000 for 100ps if timestep is 4fs).

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

### <a name="production-settings_bound-leg-settings"></a>9.5. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

| Property                                                                                        | Pattern | Type             | Deprecated | Definition | Title/Description                                                                                                                                                            |
| ----------------------------------------------------------------------------------------------- | ------- | ---------------- | ---------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| + [runtime](#production-settings_bound-leg-settings_runtime )                                   | No      | string           | No         | -          | Duration of the simulation (e.g., '2ns').                                                                                                                                    |
| - [timestep](#production-settings_bound-leg-settings_timestep )                                 | No      | string           | No         | -          | Time step for the simulation (e.g., '4fs').                                                                                                                                  |
| - [temperature](#production-settings_bound-leg-settings_temperature )                           | No      | string           | No         | -          | Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures.                                                                                   |
| - [pressure](#production-settings_bound-leg-settings_pressure )                                 | No      | string           | No         | -          | Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.                                                                                         |
| - [restraint-string](#production-settings_bound-leg-settings_restraint-string )                 | No      | string           | No         | -          | The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').                                                                                                 |
| - [restraint-indices](#production-settings_bound-leg-settings_restraint-indices )               | No      | array of integer | No         | -          | Indices of atoms to apply the restraint to.                                                                                                                                  |
| - [restraint-force-constant](#production-settings_bound-leg-settings_restraint-force-constant ) | No      | number           | No         | -          | Force constant for the restraint, kcal_per_mol / angstrom^2.                                                                                                                 |
| - [use-modified-dummies](#production-settings_bound-leg-settings_use-modified-dummies )         | No      | boolean          | No         | -          | Whether to apply modifications described here: https://doi.org/10.1021/acs.jctc.0c01328 to dummy atoms. Requires the Ghostly package: https://github.com/OpenBioSim/ghostly. |
| - [report-interval](#production-settings_bound-leg-settings_report-interval )                   | No      | integer          | No         | -          | Interval for reporting simulation data, in steps (e.g., 250 for 1ps if timestep is 4fs).                                                                                     |
| - [restart-interval](#production-settings_bound-leg-settings_restart-interval )                 | No      | integer          | No         | -          | Interval for restarting the simulation, in steps (e.g., 25000 for 100ps if timestep is 4fs).                                                                                 |

#### <a name="production-settings_bound-leg-settings_runtime"></a>9.5.1. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > runtime`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | Yes      |

**Description:** Duration of the simulation (e.g., '2ns').

| Restrictions                      |                                                                                            |
| --------------------------------- | ------------------------------------------------------------------------------------------ |
| **Must match regular expression** | ```^[0-9]+(ps\|ns)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28ps%7Cns%29%24) |

#### <a name="production-settings_bound-leg-settings_timestep"></a>9.5.2. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > timestep`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Time step for the simulation (e.g., '4fs').

| Restrictions                      |                                                                                   |
| --------------------------------- | --------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(fs)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28fs%29%24) |

#### <a name="production-settings_bound-leg-settings_temperature"></a>9.5.3. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > temperature`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Temperature for the simulation (e.g., '300K', '300.5K'). Overrides any other temperatures.

| Restrictions                      |                                                                                                                     |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(K)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28K%29%24) |

#### <a name="production-settings_bound-leg-settings_pressure"></a>9.5.4. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > pressure`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Pressure for the simulation (e.g., '1bar', '1.0bar'). Overrides any other pressures.

| Restrictions                      |                                                                                                                                    |
| --------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------- |
| **Must match regular expression** | ```^[0-9]+(\.[0-9]+)?(bar\|atm)$``` [Test](https://regex101.com/?regex=%5E%5B0-9%5D%2B%28%5C.%5B0-9%5D%2B%29%3F%28bar%7Catm%29%24) |

#### <a name="production-settings_bound-leg-settings_restraint-string"></a>9.5.5. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-string`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** The restraint applied during this stage (e.g., 'none', 'heavy', 'backbone').

#### <a name="production-settings_bound-leg-settings_restraint-indices"></a>9.5.6. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-indices`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `array of integer` |
| **Required** | No                 |

**Description:** Indices of atoms to apply the restraint to.

|                      | Array restrictions |
| -------------------- | ------------------ |
| **Min items**        | N/A                |
| **Max items**        | N/A                |
| **Items unicity**    | False              |
| **Additional items** | False              |
| **Tuple validation** | See below          |

| Each item of this array must be                                                            | Description |
| ------------------------------------------------------------------------------------------ | ----------- |
| [restraint-indices items](#production-settings_bound-leg-settings_restraint-indices_items) | -           |

##### <a name="production-settings_bound-leg-settings_restraint-indices_items"></a>9.5.6.1. RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-indices > restraint-indices items

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="production-settings_bound-leg-settings_restraint-force-constant"></a>9.5.7. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restraint-force-constant`

|              |          |
| ------------ | -------- |
| **Type**     | `number` |
| **Required** | No       |

**Description:** Force constant for the restraint, kcal_per_mol / angstrom^2.

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 0 |

#### <a name="production-settings_bound-leg-settings_use-modified-dummies"></a>9.5.8. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > use-modified-dummies`

|              |           |
| ------------ | --------- |
| **Type**     | `boolean` |
| **Required** | No        |
| **Default**  | `false`   |

**Description:** Whether to apply modifications described here: https://doi.org/10.1021/acs.jctc.0c01328 to dummy atoms. Requires the Ghostly package: https://github.com/OpenBioSim/ghostly.

#### <a name="production-settings_bound-leg-settings_report-interval"></a>9.5.9. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > report-interval`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Interval for reporting simulation data, in steps (e.g., 250 for 1ps if timestep is 4fs).

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

#### <a name="production-settings_bound-leg-settings_restart-interval"></a>9.5.10. Property `RBFE Workflow Config Schema > production-settings > bound-leg-settings > restart-interval`

|              |           |
| ------------ | --------- |
| **Type**     | `integer` |
| **Required** | No        |

**Description:** Interval for restarting the simulation, in steps (e.g., 25000 for 100ps if timestep is 4fs).

| Restrictions |        |
| ------------ | ------ |
| **Minimum**  | &ge; 1 |

## <a name="analysis-settings"></a>10. Property `RBFE Workflow Config Schema > analysis-settings`

|                           |             |
| ------------------------- | ----------- |
| **Type**                  | `object`    |
| **Required**              | Yes         |
| **Additional properties** | Not allowed |

| Property                                                           | Pattern | Type             | Deprecated | Definition | Title/Description                                                                                                 |
| ------------------------------------------------------------------ | ------- | ---------------- | ---------- | ---------- | ----------------------------------------------------------------------------------------------------------------- |
| - [backend](#analysis-settings_backend )                           | No      | enum (of string) | No         | -          | Backend to use for analysis. If cinnabar is not installed, native will be used.                                   |
| - [experimental-results](#analysis-settings_experimental-results ) | No      | string           | No         | -          | Path to the file containing experimental results.                                                                 |
| - [experimental-units](#analysis-settings_experimental-units )     | No      | enum (of string) | No         | -          | Units used in the experimental results. Note that these are input units only, outputs will always be in kcal/mol. |

### <a name="analysis-settings_backend"></a>10.1. Property `RBFE Workflow Config Schema > analysis-settings > backend`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | No                 |

**Description:** Backend to use for analysis. If cinnabar is not installed, native will be used.

Must be one of:
* "cinnabar"
* "native"

### <a name="analysis-settings_experimental-results"></a>10.2. Property `RBFE Workflow Config Schema > analysis-settings > experimental-results`

|              |          |
| ------------ | -------- |
| **Type**     | `string` |
| **Required** | No       |

**Description:** Path to the file containing experimental results.

### <a name="analysis-settings_experimental-units"></a>10.3. Property `RBFE Workflow Config Schema > analysis-settings > experimental-units`

|              |                    |
| ------------ | ------------------ |
| **Type**     | `enum (of string)` |
| **Required** | No                 |

**Description:** Units used in the experimental results. Note that these are input units only, outputs will always be in kcal/mol.

Must be one of:
* "Ki_uM"
* "kcal/mol"
* "kJ/mol"

----------------------------------------------------------------------------------------------------------------------------
Generated using [json-schema-for-humans](https://github.com/coveooss/json-schema-for-humans) on 2025-08-21 at 16:07:07 +0100
