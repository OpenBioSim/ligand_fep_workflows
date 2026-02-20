from pathlib import Path


def create_python_script_call(wc, input, leg):
    output_directory = Path(f"{config['working_directory']}/production/{_engine}/{wc.ligand1}~{wc.ligand2}/{leg}_{wc.replica_number}")
    args = [f"--engine {_engine}"]

    if _engine == "somd2":
        cfg = config["production-settings"]["somd2-settings"]
        args.append(f"--runtime {cfg['runtime']}")
        if cfg.get("timestep"):
            args.append(f"--timestep {cfg['timestep']}")
        if cfg.get("temperature"):
            args.append(f"--temperature {cfg['temperature']}")
        if cfg.get("pressure"):
            args.append(f"--pressure {cfg['pressure']}")
        if cfg.get("use-modified-dummies", False):
            args.append("--use-modified-dummies")
        if cfg.get("cutoff_type"):
            args.append(f"--cutoff-type {cfg['cutoff_type']}")
        if cfg.get("cutoff"):
            args.append(f"--cutoff {cfg['cutoff']}")
        if cfg.get("energy_frequency"):
            args.append(f"--energy-frequency {cfg['energy_frequency']}")
        if cfg.get("frame_frequency"):
            args.append(f"--frame-frequency {cfg['frame_frequency']}")
        if cfg.get("checkpoint_frequency"):
            args.append(f"--checkpoint-frequency {cfg['checkpoint_frequency']}")
        if cfg.get("integrator"):
            args.append(f"--integrator {cfg['integrator']}")
        if cfg.get("shift_delta"):
            args.append(f"--shift-delta {cfg['shift_delta']}")
        if cfg.get("perturbable_constraint"):
            args.append(f"--perturbable-constraint {cfg['perturbable_constraint']}")
        runner = cfg.get("runner", "repex")
        args.append(f"--runner {runner}")
    else:
        cfg = config["production-settings"]["gromacs-settings"][f"{leg}-leg-settings"]
        args.append(f"--runtime {cfg['runtime']}")
        if cfg.get("timestep"):
            args.append(f"--timestep {cfg['timestep']}")
        if cfg.get("temperature"):
            args.append(f"--temperature {cfg['temperature']}")
        if cfg.get("pressure"):
            args.append(f"--pressure {cfg['pressure']}")
        if cfg.get("restraint-string") is not None:
            args.append(f"--restraint-string {cfg['restraint-string']}")
        if cfg.get("restraint-indices") is not None:
            args.append("--restraint-indices " + " ".join(map(str, cfg["restraint-indices"])))
        if cfg.get("restraint-force-constant") is not None:
            args.append(f"--restraint-force-constant {cfg['restraint-force-constant']}")
        if cfg.get("use-modified-dummies", False):
            args.append("--use-modified-dummies")
        if cfg.get("report-interval"):
            args.append(f"--report-interval {cfg['report-interval']}")
        if cfg.get("restart-interval"):
            args.append(f"--restart-interval {cfg['restart-interval']}")

    args.append(f"--network-location {config['working_directory']}/network")

    return f"""
    python workflow/scripts/rbfe/production.py --input {input.file} --output-directory {output_directory} {" ".join(args)}
    """


def _run_gromacs_stages(output_directory):
    """Run GROMACS minimisation, heating, equilibration, and production stages."""
    outdir_path_min = Path(output_directory) / "minimisation"
    lambda_values = [d.name.split('_')[1] for d in outdir_path_min.glob("lambda_*") if d.is_dir()]
    if not lambda_values:
        raise FileNotFoundError(f"No minimisation directories found in {outdir_path_min}")
    lambda_values.sort(key=float)

    print("Minimising")
    for lambda_value in lambda_values:
        d = f"{output_directory}/minimisation/lambda_{lambda_value}"
        shell(f"gmx grompp -f {d}/gromacs.mdp -c {d}/gromacs_ref.gro -p {d}/gromacs.top -o {d}/gromacs.tpr 2>&1 | tee {d}/grompp.log")
        shell(f"gmx mdrun -ntmpi 1 -deffnm {d}/gromacs 2>&1 | tee {d}/mdrun.log")

    print("Heating")
    for lambda_value in lambda_values:
        d = f"{output_directory}/heat/lambda_{lambda_value}"
        prev_gro = f"{output_directory}/minimisation/lambda_{lambda_value}/gromacs.gro"
        shell(f"gmx grompp -f {d}/gromacs.mdp -c {prev_gro} -p {d}/gromacs.top -o {d}/gromacs.tpr 2>&1 | tee {d}/grompp.log")
        shell(f"gmx mdrun -ntmpi 1 -deffnm {d}/gromacs 2>&1 | tee {d}/mdrun.log")

    print("Equilibrating")
    for lambda_value in lambda_values:
        d = f"{output_directory}/eq/lambda_{lambda_value}"
        prev_gro = f"{output_directory}/heat/lambda_{lambda_value}/gromacs.gro"
        shell(f"gmx grompp -f {d}/gromacs.mdp -c {prev_gro} -p {d}/gromacs.top -o {d}/gromacs.tpr 2>&1 | tee {d}/grompp.log")
        shell(f"gmx mdrun -ntmpi 1 -deffnm {d}/gromacs 2>&1 | tee {d}/mdrun.log")

    print("Running production")
    for lambda_value in lambda_values:
        d = f"{output_directory}/lambda_{lambda_value}"
        prev_gro = f"{output_directory}/eq/lambda_{lambda_value}/gromacs.gro"
        shell(f"gmx grompp -f {d}/gromacs.mdp -c {prev_gro} -p {d}/gromacs.top -o {d}/gromacs.tpr 2>&1 | tee {d}/grompp.log")
        shell(f"gmx mdrun -ntmpi 1 -deffnm {d}/gromacs 2>&1 | tee {d}/mdrun.log")

    # Remove intermediate directories as they confuse the analysis
    shell(f"rm -rf {output_directory}/minimisation {output_directory}/heat {output_directory}/eq")


rule production_bound:
    priority: 2
    input:
        file = Path(f"{config['working_directory']}/rbfe_prepared/bound/{{ligand1}}~{{ligand2}}.bss"),
        prev_replica = lambda wc: [] if int(wc.replica_number) == 0 else [
            f"{config['working_directory']}/production/{_engine}/{wc.ligand1}~{wc.ligand2}/bound_{int(wc.replica_number) - 1}/.done"
        ],
    output:
        done = Path(f"{config['working_directory']}/production/{_engine}/{{ligand1}}~{{ligand2}}/bound_{{replica_number}}/.done")
    threads: config["simulation_threads"]
    resources:
        gpu=config["production-settings"].get("somd2-settings", {}).get("gpus_per_job", 1) if _engine == "somd2" else 1
    log:
        Path(f"{config['working_directory']}/logs/{{ligand1}}~{{ligand2}}_production_bound_{{replica_number}}.log")
    run:
        python_cmd = create_python_script_call(wildcards, input, "bound")
        shell(python_cmd)
        if _engine == "gromacs":
            output_directory = str(Path(f"{config['working_directory']}/production/{_engine}/{wildcards.ligand1}~{wildcards.ligand2}/bound_{wildcards.replica_number}"))
            _run_gromacs_stages(output_directory)
        shell(f"touch {output.done}")


rule production_free:
    priority: 1
    input:
        file = Path(f"{config['working_directory']}/rbfe_prepared/free/{{ligand1}}~{{ligand2}}.bss"),
        prev_replica = lambda wc: [] if int(wc.replica_number) == 0 else [
            f"{config['working_directory']}/production/{_engine}/{wc.ligand1}~{wc.ligand2}/free_{int(wc.replica_number) - 1}/.done"
        ],
    output:
        done = Path(f"{config['working_directory']}/production/{_engine}/{{ligand1}}~{{ligand2}}/free_{{replica_number}}/.done")
    threads: config["simulation_threads"]
    resources:
        gpu=config["production-settings"].get("somd2-settings", {}).get("gpus_per_job", 1) if _engine == "somd2" else 1
    log:
        Path(f"{config['working_directory']}/logs/{{ligand1}}~{{ligand2}}_production_free_{{replica_number}}.log")
    run:
        python_cmd = create_python_script_call(wildcards, input, "free")
        shell(python_cmd)
        if _engine == "gromacs":
            output_directory = str(Path(f"{config['working_directory']}/production/{_engine}/{wildcards.ligand1}~{wildcards.ligand2}/free_{wildcards.replica_number}"))
            _run_gromacs_stages(output_directory)
        shell(f"touch {output.done}")
