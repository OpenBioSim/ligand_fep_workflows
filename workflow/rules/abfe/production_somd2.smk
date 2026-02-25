"""
Snakemake rules for ABFE production simulations using SOMD2.

SOMD2 (sire/OpenMM) handles minimisation, equilibration, and production
in a single step via its Runner/RepexRunner API. It takes equilibrated
systems from the preparation stage and applies sire-native decoupling,
so the abfe_prep step (BSS decoupling + HMR) is not needed.

The lambda schedule uses a decharge -> annihilate approach:
    Stage 1 (decharge): ramp charges off, ramp restraints on
    Stage 2 (annihilate): ramp vdW off, charges off, restraints on

Each leg runs for multiple replicas to estimate the free energy change
and associated uncertainty.
"""

from pathlib import Path

_engine = config.get("engine", config["production-settings"].get("engine", "gromacs")).strip().lower()


# Bound leg production (SOMD2)
# ============================

rule somd2_production_bound:
    """
    SOMD2 production for bound leg.

    Loads the equilibrated system, applies sire-native decoupling,
    constructs the lambda schedule, and runs SOMD2 with Boresch restraints.

    Replica N waits for replica N-1's production to complete,
    so all ligands get a first result before any starts a second replica.
    """
    priority: 2
    input:
        system=Path(f"{config['working_directory']}/preparation/final")
        / "{ligand}_bound.bss",
        restraint=Path(f"{config['working_directory']}/restraints")
        / "{ligand}_restraint.json",
        prev_replica=lambda wc: [] if int(wc.replica) == 0 else [
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/bound_{int(wc.replica) - 1}/.done",
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/free_{int(wc.replica) - 1}/.done",
        ],
    output:
        done=Path(
            f"{config['working_directory']}/production/{_engine}/{{ligand}}/bound_{{replica}}/.done"
        ),
    threads:
        config["simulation_threads"]
    resources:
        gpu=config["production-settings"].get("somd2-settings", {}).get("gpus_per_job", 1),
        tasks_per_gpu=0
    log:
        Path(f"{config['working_directory']}/logs")
        / "{ligand}_somd2_production_bound_{replica}.log",
    params:
        script=Path("workflow/scripts/abfe/production_somd2.py"),
        runtime=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "runtime", "2ns"
        ),
        timestep=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "timestep", "4fs"
        ),
        temperature=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "temperature", "298K"
        ),
        cutoff_type=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "cutoff_type", "PME"
        ),
        cutoff=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "cutoff", "10A"
        ),
        num_lambda=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "num_lambda", 21
        ),
        energy_frequency=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "energy_frequency", "1ps"
        ),
        frame_frequency=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "frame_frequency", "500ps"
        ),
        checkpoint_frequency=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "checkpoint_frequency", "500ps"
        ),
        integrator=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "integrator", "langevin_middle"
        ),
        shift_delta=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "shift_delta", "2.25A"
        ),
        perturbable_constraint=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "perturbable_constraint", "h_bonds_not_heavy_perturbed"
        ),
        equilibration_time=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "equilibration_time", "20ps"
        ),
        runner=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "runner", "repex"
        ),
        perturbation_type=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "perturbation_type", "annihilate"
        ),
        output_directory=lambda wc: Path(
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/bound_{wc.replica}"
        ),
    shell:
        """
        echo "Running SOMD2 bound leg for {wildcards.ligand} replica {wildcards.replica}"
        python {params.script} \
            --input {input.system} \
            --output-directory {params.output_directory} \
            --leg bound \
            --restraint-file {input.restraint} \
            --runtime {params.runtime} \
            --timestep {params.timestep} \
            --temperature {params.temperature} \
            --cutoff-type {params.cutoff_type} \
            --cutoff {params.cutoff} \
            --num-lambda {params.num_lambda} \
            --energy-frequency {params.energy_frequency} \
            --frame-frequency {params.frame_frequency} \
            --checkpoint-frequency {params.checkpoint_frequency} \
            --integrator {params.integrator} \
            --shift-delta {params.shift_delta} \
            --perturbable-constraint {params.perturbable_constraint} \
            --equilibration-time {params.equilibration_time} \
            --runner {params.runner} \
            --perturbation-type {params.perturbation_type} \
            2>&1 | tee {log}
        touch {output.done}
        """


# Free leg production (SOMD2)
# ============================

rule somd2_production_free:
    """
    SOMD2 production for free leg.

    Loads the equilibrated system, applies sire-native decoupling,
    constructs the lambda schedule, and runs SOMD2 without restraints.

    Replica N waits for replica N-1's production to complete,
    so all ligands get a first result before any starts a second replica.
    """
    priority: 2
    input:
        system=Path(f"{config['working_directory']}/preparation/final")
        / "{ligand}_free.bss",
        prev_replica=lambda wc: [] if int(wc.replica) == 0 else [
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/bound_{int(wc.replica) - 1}/.done",
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/free_{int(wc.replica) - 1}/.done",
        ],
    output:
        done=Path(
            f"{config['working_directory']}/production/{_engine}/{{ligand}}/free_{{replica}}/.done"
        ),
    threads:
        config["simulation_threads"]
    resources:
        gpu=config["production-settings"].get("somd2-settings", {}).get("gpus_per_job", 1),
        tasks_per_gpu=0
    log:
        Path(f"{config['working_directory']}/logs")
        / "{ligand}_somd2_production_free_{replica}.log",
    params:
        script=Path("workflow/scripts/abfe/production_somd2.py"),
        runtime=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "runtime", "2ns"
        ),
        timestep=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "timestep", "4fs"
        ),
        temperature=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "temperature", "298K"
        ),
        cutoff_type=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "cutoff_type", "PME"
        ),
        cutoff=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "cutoff", "10A"
        ),
        num_lambda=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "num_lambda", 21
        ),
        energy_frequency=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "energy_frequency", "1ps"
        ),
        frame_frequency=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "frame_frequency", "500ps"
        ),
        checkpoint_frequency=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "checkpoint_frequency", "500ps"
        ),
        integrator=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "integrator", "langevin_middle"
        ),
        shift_delta=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "shift_delta", "2.25A"
        ),
        perturbable_constraint=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "perturbable_constraint", "h_bonds_not_heavy_perturbed"
        ),
        equilibration_time=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "equilibration_time", "20ps"
        ),
        runner=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "runner", "repex"
        ),
        perturbation_type=lambda wc: config["production-settings"].get("somd2-settings", {}).get(
            "perturbation_type", "annihilate"
        ),
        output_directory=lambda wc: Path(
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/free_{wc.replica}"
        ),
    shell:
        """
        echo "Running SOMD2 free leg for {wildcards.ligand} replica {wildcards.replica}"
        python {params.script} \
            --input {input.system} \
            --output-directory {params.output_directory} \
            --leg free \
            --runtime {params.runtime} \
            --timestep {params.timestep} \
            --temperature {params.temperature} \
            --cutoff-type {params.cutoff_type} \
            --cutoff {params.cutoff} \
            --num-lambda {params.num_lambda} \
            --energy-frequency {params.energy_frequency} \
            --frame-frequency {params.frame_frequency} \
            --checkpoint-frequency {params.checkpoint_frequency} \
            --integrator {params.integrator} \
            --shift-delta {params.shift_delta} \
            --perturbable-constraint {params.perturbable_constraint} \
            --equilibration-time {params.equilibration_time} \
            --runner {params.runner} \
            --perturbation-type {params.perturbation_type} \
            2>&1 | tee {log}
        touch {output.done}
        """
