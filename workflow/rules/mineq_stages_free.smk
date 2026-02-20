"""
Snakemake rules for minimisation and equilibration of free leg in ABFE.

These rules are dynamically generated from the configuration file.
The free leg contains only the ligand in solvent, so backbone restraints
are not applicable.
"""

from pathlib import Path
import re


# Access free leg stages from config
_free_stages = config["min/eq-stages-free"]
_free_stage_names = list(_free_stages.keys())


def _get_script_name_free(stage_name: str) -> str:
    """
    Determine the appropriate script based on stage type.

    Args:
        stage_name: Name of the stage

    Returns:
        Script filename
    """
    match = re.match(r"^(minimisation|equilibration)", stage_name)
    if not match:
        raise ValueError(f"Invalid stage name: {stage_name}")
    return f"{match.group(1)}.py"


def _stage_output_free(stage: str) -> Path:
    """Get the output directory for a given free leg stage.

    Produces a nested structure: preparation/{type}/{number}/
    e.g. minimisation1 -> preparation/minimisation/1/
         equilibration3 -> preparation/equilibration/3/
    """
    match = re.match(r"^(minimisation|equilibration)(\d+)$", stage)
    if not match:
        raise ValueError(f"Cannot parse stage name: {stage}")
    stage_type = match.group(1)
    stage_num = match.group(2)
    return Path(f"{config['working_directory']}/preparation/{stage_type}/{stage_num}")


# Generate rules dynamically for each stage
for _i_free, _stage_free in enumerate(_free_stage_names):
    _is_last_stage_free = _i_free == len(_free_stage_names) - 1

    _params_free = _free_stages[_stage_free]
    _prefix_free = _get_script_name_free(_stage_free)

    # Build stage-specific arguments
    if _stage_free.startswith("minimisation"):
        _extra_arg_free = f"--minimisation-steps {_params_free['minimisation-steps']}"
    elif _stage_free.startswith("equilibration"):
        _extra_arg_free = f"--runtime {_params_free['runtime']}"
        if _params_free.get("temperature-start") is not None:
            _extra_arg_free += f" --temperature-start {_params_free['temperature-start']}"
        if _params_free.get("temperature-end") is not None:
            _extra_arg_free += f" --temperature-end {_params_free['temperature-end']}"
        if _params_free.get("temperature") is not None:
            _extra_arg_free += f" --temperature {_params_free['temperature']}"
        if _params_free.get("pressure") is not None:
            _extra_arg_free += f" --pressure {_params_free['pressure']}"
        if _params_free.get("thermostat-time-constant") is not None:
            _extra_arg_free += (
                f" --thermostat-time-constant {_params_free['thermostat-time-constant']}"
            )
    else:
        _extra_arg_free = ""

    # Define input function
    if _i_free == 0:

        def _input_func_free(wildcards):
            return Path(
                f"{config['working_directory']}/setup/{wildcards.ligand}_free.bss"
            )

    else:
        _prev_stage_free = _free_stage_names[_i_free - 1]

        def _input_func_free(wildcards, prev_stage=_prev_stage_free):
            return Path(f"{_stage_output_free(prev_stage)}/{wildcards.ligand}_free.bss")

    # Define the rule
    rule:
        name:
            f"run_{_stage_free}_free"
        input:
            prev=_input_func_free
        output:
            f"{config['working_directory']}/preparation/final/{{ligand}}_free.bss"
            if _is_last_stage_free
            else f"{_stage_output_free(_stage_free)}/{{ligand}}_free.bss"
        threads:
            config["simulation_threads"]
        resources:
            gpu=1
        params:
            output_dir=(
                Path(f"{config['working_directory']}/preparation/final")
                if _is_last_stage_free
                else _stage_output_free(_stage_free)
            ),
            engine_arg=(
                f"--engine {_params_free.get('engine').strip().lower()}" if _params_free.get("engine") else ""
            ),
            restraint_string_arg=(
                f"--restraint-string {_params_free.get('restraint-string')}"
                if _params_free.get("restraint-string")
                else ""
            ),
            restraint_indices_arg=(
                f"--restraint-indices {_params_free.get('restraint-indices')}"
                if _params_free.get("restraint-indices")
                else ""
            ),
            restraint_force_constant_arg=(
                f"--restraint-force-constant {_params_free.get('restraint-force-constant')}"
                if _params_free.get("restraint-force-constant")
                else ""
            ),
            script=f"workflow/scripts/{_prefix_free}",
            extra_arg=_extra_arg_free,
        shell:
            """
            echo "Running {rule} for {wildcards.ligand}"
            python {params.script} \
                --input {input.prev} \
                --output-directory {params.output_dir} \
                --outfile-name {wildcards.ligand}_free \
                {params.engine_arg} \
                {params.extra_arg} \
                {params.restraint_string_arg} \
                {params.restraint_indices_arg} \
                {params.restraint_force_constant_arg}
            """
