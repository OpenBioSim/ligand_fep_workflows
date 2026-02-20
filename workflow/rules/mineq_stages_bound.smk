"""
Snakemake rules for minimisation and equilibration of bound leg in ABFE.

These rules are dynamically generated from the configuration file to allow
flexible protocol definition. Each stage (minimisation or equilibration)
becomes a separate Snakemake rule that is executed sequentially.

The bound leg requires protein-specific considerations like backbone restraints
during early equilibration stages.
"""

from pathlib import Path
import re


# Access bound leg stages from config
_bound_stages = config["min/eq-stages-bound"]
_bound_stage_names = list(_bound_stages.keys())


def _get_script_name_bound(stage_name: str) -> str:
    """
    Determine the appropriate script based on stage type.

    Args:
        stage_name: Name of the stage (e.g., 'minimisation1', 'equilibration2')

    Returns:
        Script filename (minimisation.py or equilibration.py)

    Raises:
        ValueError: If stage name doesn't start with 'minimisation' or 'equilibration'
    """
    match = re.match(r"^(minimisation|equilibration)", stage_name)
    if not match:
        raise ValueError(f"Invalid stage name: {stage_name}")
    return f"{match.group(1)}.py"


def _stage_output_bound(stage: str) -> Path:
    """
    Get the output directory for a given bound leg stage.

    Produces a nested structure: preparation/{type}/{number}/
    e.g. minimisation1 -> preparation/minimisation/1/
         equilibration3 -> preparation/equilibration/3/

    Args:
        stage: Stage name

    Returns:
        Path to the stage output directory
    """
    match = re.match(r"^(minimisation|equilibration)(\d+)$", stage)
    if not match:
        raise ValueError(f"Cannot parse stage name: {stage}")
    stage_type = match.group(1)
    stage_num = match.group(2)
    return Path(f"{config['working_directory']}/preparation/{stage_type}/{stage_num}")


# Generate rules dynamically for each stage
for _i_bound, _stage_bound in enumerate(_bound_stage_names):
    # Determine if this is the final stage (output goes to special directory)
    _is_last_stage_bound = _i_bound == len(_bound_stage_names) - 1

    _params_bound = _bound_stages[_stage_bound]
    _prefix_bound = _get_script_name_bound(_stage_bound)

    # Build stage-specific arguments
    if _stage_bound.startswith("minimisation"):
        _extra_arg_bound = f"--minimisation-steps {_params_bound['minimisation-steps']}"
    elif _stage_bound.startswith("equilibration"):
        _extra_arg_bound = f"--runtime {_params_bound['runtime']}"
        if _params_bound.get("temperature-start") is not None:
            _extra_arg_bound += f" --temperature-start {_params_bound['temperature-start']}"
        if _params_bound.get("temperature-end") is not None:
            _extra_arg_bound += f" --temperature-end {_params_bound['temperature-end']}"
        if _params_bound.get("temperature") is not None:
            _extra_arg_bound += f" --temperature {_params_bound['temperature']}"
        if _params_bound.get("pressure") is not None:
            _extra_arg_bound += f" --pressure {_params_bound['pressure']}"
        if _params_bound.get("thermostat-time-constant") is not None:
            _extra_arg_bound += (
                f" --thermostat-time-constant {_params_bound['thermostat-time-constant']}"
            )
    else:
        _extra_arg_bound = ""

    # Define input function - first stage takes from setup, others from previous stage
    if _i_bound == 0:

        def _input_func_bound(wildcards):
            return Path(
                f"{config['working_directory']}/setup/{wildcards.ligand}_bound.bss"
            )

    else:
        _prev_stage_bound = _bound_stage_names[_i_bound - 1]

        def _input_func_bound(wildcards, prev_stage=_prev_stage_bound):
            return Path(f"{_stage_output_bound(prev_stage)}/{wildcards.ligand}_bound.bss")

    # Define the rule
    rule:
        name:
            f"run_{_stage_bound}_bound"
        input:
            prev=_input_func_bound
        output:
            f"{config['working_directory']}/preparation/final/{{ligand}}_bound.bss"
            if _is_last_stage_bound
            else f"{_stage_output_bound(_stage_bound)}/{{ligand}}_bound.bss"
        threads:
            config["simulation_threads"]
        resources:
            gpu=1
        params:
            output_dir=(
                Path(f"{config['working_directory']}/preparation/final")
                if _is_last_stage_bound
                else _stage_output_bound(_stage_bound)
            ),
            engine_arg=(
                f"--engine {_params_bound.get('engine').strip().lower()}" if _params_bound.get("engine") else ""
            ),
            restraint_string_arg=(
                f"--restraint-string {_params_bound.get('restraint-string')}"
                if _params_bound.get("restraint-string")
                else ""
            ),
            restraint_indices_arg=(
                f"--restraint-indices {_params_bound.get('restraint-indices')}"
                if _params_bound.get("restraint-indices")
                else ""
            ),
            restraint_force_constant_arg=(
                f"--restraint-force-constant {_params_bound.get('restraint-force-constant')}"
                if _params_bound.get("restraint-force-constant")
                else ""
            ),
            script=f"workflow/scripts/{_prefix_bound}",
            extra_arg=_extra_arg_bound,
        shell:
            """
            echo "Running {rule} for {wildcards.ligand}"
            python {params.script} \
                --input {input.prev} \
                --output-directory {params.output_dir} \
                --outfile-name {wildcards.ligand}_bound \
                {params.engine_arg} \
                {params.extra_arg} \
                {params.restraint_string_arg} \
                {params.restraint_indices_arg} \
                {params.restraint_force_constant_arg}
            """
