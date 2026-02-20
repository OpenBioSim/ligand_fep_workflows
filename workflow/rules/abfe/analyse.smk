"""
Snakemake rules for ABFE analysis.

The analysis workflow consists of:
    1. Per-leg analysis: Extract PMF from each leg's lambda windows
    2. Final analysis: Combine legs and apply correction

The binding free energy is computed as:
    DG_bind = DG_free - DG_correction - DG_bound
"""

from pathlib import Path

_engine = config.get("engine", config["production-settings"].get("engine", "gromacs")).strip().lower()


def get_all_leg_pmfs(wildcards) -> list[str]:
    """
    Get all PMF files needed for final analysis.

    Returns paths to all leg PMF files for all ligands and replicas.
    """
    num_replicas = config["production-settings"]["num_replicas"]
    working_dir = config["working_directory"]

    pmf_files = []
    for ligand in LIGANDS:
        for replica in range(num_replicas):
            # Bound leg
            pmf_files.append(
                f"{working_dir}/analysis/{_engine}/{ligand}/bound_{replica}/pmf.csv"
            )
            # Free leg
            pmf_files.append(
                f"{working_dir}/analysis/{_engine}/{ligand}/free_{replica}/pmf.csv"
            )
    return pmf_files


# Per-leg analysis rules
# ======================

rule analyse_leg:
    """
    Analyse a single ABFE leg to extract PMF.

    Uses MBAR to estimate the free energy change across lambda windows
    and generates overlap matrix for convergence assessment.
    """
    priority: 10
    input:
        done=Path(
            f"{config['working_directory']}/production/{_engine}/{{ligand}}/{{leg}}_{{replica}}/.done"
        ),
    output:
        pmf=Path(
            f"{config['working_directory']}/analysis/{_engine}/{{ligand}}/{{leg}}_{{replica}}/pmf.csv"
        ),
        overlap=Path(
            f"{config['working_directory']}/analysis/{_engine}/{{ligand}}/{{leg}}_{{replica}}/overlap.npy"
        ),
    log:
        Path(f"{config['working_directory']}/logs")
        / "{ligand}_analysis_{leg}_{replica}.log",
    threads:
        config["simulation_threads"]
    resources:
        gpu=0
    params:
        script=Path("workflow/scripts/abfe/analyse_leg.py"),
        input_directory=lambda wc: Path(
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/{wc.leg}_{wc.replica}"
        ),
        output_directory=lambda wc: Path(
            f"{config['working_directory']}/analysis/{_engine}/{wc.ligand}/{wc.leg}_{wc.replica}"
        ),
        temperature=lambda wc: (
            config["production-settings"].get("somd2-settings", {}).get("temperature", "298K")
            if _engine == "somd2"
            else config["production-settings"].get("gromacs-settings", {}).get("bound-leg-settings", {}).get("temperature", "300K")
        ),
    shell:
        """
        echo "Analysing {wildcards.ligand} {wildcards.leg} replica {wildcards.replica}"
        python {params.script} \
            --input-directory {params.input_directory} \
            --output-directory {params.output_directory} \
            --temperature {params.temperature} \
            --plot-overlap-matrix \
            --plot-pmf \
            2>&1 | tee {log}
        """


# Final analysis rule
# ===================

rule collate_abfe_analysis:
    """
    Collate all ABFE results and compute binding free energies.

    This rule:
        1. Combines leg PMFs for each ligand
        2. Applies restraint correction
        3. Calculates DG_bind for each ligand
        4. Averages over replicas with error estimation
        5. Optionally compares to experimental data
    """
    priority: 20
    input:
        pmf_files=get_all_leg_pmfs,
        corrections=expand(
            f"{config['working_directory']}/restraints/{{ligand}}_correction.txt",
            ligand=LIGANDS,
        ),
    output:
        results=Path(f"{config['working_directory']}/analysis/{_engine}/final_abfe_results.csv"),
    log:
        Path(f"{config['working_directory']}/logs/final_analysis.log"),
    threads:
        config["simulation_threads"]
    resources:
        gpu=0
    params:
        script=Path("workflow/scripts/abfe/final_analysis.py"),
        output_directory=Path(f"{config['working_directory']}/analysis/{_engine}"),
        num_replicas=config["production-settings"]["num_replicas"],
        exp_results=(
            config["analysis-settings"].get("experimental-results")
            if config["analysis-settings"].get("experimental-results")
            else ""
        ),
        exp_units=config["analysis-settings"].get("experimental-units", "kcal/mol"),
        backend=config["analysis-settings"].get("backend", "native"),
    run:
        import json

        # Build ligands list
        ligands_str = json.dumps(LIGANDS)

        # Build command
        cmd = f"""
            python {params.script} \
                --analysis-directory {params.output_directory} \
                --restraints-directory {config['working_directory']}/restraints \
                --output-directory {params.output_directory} \
                --ligands '{ligands_str}' \
                --num-replicas {params.num_replicas} \
                --plotting-backend {params.backend}
            """

        if params.exp_results:
            cmd += f" --experimental-results {params.exp_results}"
            cmd += f" --experimental-units {params.exp_units}"

        shell(f"{cmd} 2>&1 | tee {log}")
