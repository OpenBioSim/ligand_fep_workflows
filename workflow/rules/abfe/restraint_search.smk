"""
Snakemake rule for Boresch restraint search in ABFE workflow.

This rule runs a short unrestrained MD simulation of the protein-ligand complex,
then analyses the trajectory to find optimal Boresch restraints. Boresch restraints
define 6 degrees of freedom (1 distance, 2 angles, 3 dihedrals) using 3 anchor
points each in the protein and ligand.

The restraint correction term is also computed and saved for use in the
final free energy calculation.
"""

from pathlib import Path


# Get restraint search settings
restraint_settings = config.get("restraint_search_settings", {})


rule restraint_search:
    """
    Find optimal Boresch restraints for ABFE bound leg.

    This rule:
        1. Runs a short unrestrained equilibration of the complex
        2. Analyses the trajectory to find stable anchor points
        3. Calculates force constants from native fluctuations
        4. Computes the analytical correction for restraint release
        5. Saves restraint parameters for production simulations

    The output restraint.json file contains all information needed
    to apply Boresch restraints during the bound leg calculation.
    """
    input:
        system=Path(f"{config['working_directory']}/preparation/final")
        / "{ligand}_bound.bss",
    output:
        restraint=Path(f"{config['working_directory']}/restraints")
        / "{ligand}_restraint.json",
        correction=Path(f"{config['working_directory']}/restraints")
        / "{ligand}_correction.txt",
    log:
        Path(f"{config['working_directory']}/logs") / "{ligand}_restraint_search.log",
    threads:
        config["simulation_threads"]
    resources:
        gpu=1
    params:
        script=Path("workflow/scripts/abfe/restraint_search.py"),
        output_directory=Path(f"{config['working_directory']}/restraints"),
        search_runtime=restraint_settings.get("search_runtime", "1ns"),
        temperature=restraint_settings.get("temperature", "300K"),
    shell:
        """
        echo "Running restraint search for {wildcards.ligand}"
        python {params.script} \
            --input {input.system} \
            --output-directory {params.output_directory} \
            --ligand-name {wildcards.ligand} \
            --runtime {params.search_runtime} \
            --temperature {params.temperature} \
            2>&1 | tee {log}
        """
