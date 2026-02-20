"""
Snakemake rule for ABFE system preparation.

This rule prepares systems for alchemical decoupling by:
    1. Marking the ligand for decoupling using BSS.Align.decouple()
    2. Applying hydrogen mass repartitioning (HMR) for 4fs timesteps
    3. Saving prepared systems for both bound and free legs

Unlike RBFE which creates merged molecules between two ligands, ABFE
marks a single ligand for complete decoupling from the environment.
"""

from pathlib import Path


rule abfe_prep_bound:
    """
    Prepare bound leg system for ABFE decoupling.

    The ligand is marked for decoupling and the restraint parameters
    are incorporated into the system for the restrain stage.
    """
    input:
        system=Path(f"{config['working_directory']}/preparation/final")
        / "{ligand}_bound.bss",
        restraint=Path(f"{config['working_directory']}/restraints")
        / "{ligand}_restraint.json",
    output:
        prepared=Path(f"{config['working_directory']}/abfe_prepared")
        / "{ligand}_bound.bss",
    log:
        Path(f"{config['working_directory']}/logs") / "{ligand}_abfe_prep_bound.log",
    threads:
        config["simulation_threads"]
    resources:
        mem_mb=5000
    params:
        script=Path("workflow/scripts/abfe/abfe_prep.py"),
        output_directory=Path(f"{config['working_directory']}/abfe_prepared"),
        hmr_factor=config["production-settings"].get("gromacs-settings", {}).get("hmr_factor", 3),
    shell:
        """
        echo "Preparing bound leg for {wildcards.ligand}"
        python {params.script} \
            --input {input.system} \
            --output-directory {params.output_directory} \
            --ligand-name {wildcards.ligand} \
            --leg bound \
            --restraint-file {input.restraint} \
            --hmr-factor {params.hmr_factor} \
            2>&1 | tee {log}
        """


rule abfe_prep_free:
    """
    Prepare free leg system for ABFE decoupling.

    Only the ligand decoupling is applied; no restraints are needed
    for the free leg since the ligand is not interacting with a receptor.
    """
    input:
        system=Path(f"{config['working_directory']}/preparation/final")
        / "{ligand}_free.bss",
    output:
        prepared=Path(f"{config['working_directory']}/abfe_prepared")
        / "{ligand}_free.bss",
    log:
        Path(f"{config['working_directory']}/logs") / "{ligand}_abfe_prep_free.log",
    threads:
        config["simulation_threads"]
    resources:
        mem_mb=5000
    params:
        script=Path("workflow/scripts/abfe/abfe_prep.py"),
        output_directory=Path(f"{config['working_directory']}/abfe_prepared"),
        hmr_factor=config["production-settings"].get("gromacs-settings", {}).get("hmr_factor", 3),
    shell:
        """
        echo "Preparing free leg for {wildcards.ligand}"
        python {params.script} \
            --input {input.system} \
            --output-directory {params.output_directory} \
            --ligand-name {wildcards.ligand} \
            --leg free \
            --hmr-factor {params.hmr_factor} \
            2>&1 | tee {log}
        """
