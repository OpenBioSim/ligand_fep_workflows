from pathlib import Path
import pandas as pd

_gromacs_runner = config["production-settings"].get("gromacs-settings", {}).get("runner", "standard").strip().lower()
_repex_frequency = config["production-settings"].get("gromacs-settings", {}).get("repex-frequency", 1000)
_nex = config["production-settings"].get("gromacs-settings", {}).get("nex", 1000000)
_oversubscribe = config["production-settings"].get("gromacs-settings", {}).get("oversubscribe", True)
_gromacs_gpus_per_job = config["production-settings"].get("gromacs-settings", {}).get("gpus_per_job", 1)
_gromacs_oversubscription_factor = config["production-settings"].get("gromacs-settings", {}).get("oversubscription_factor", 1)
_gromacs_mps_pct = config["production-settings"].get("gromacs-settings", {}).get("mps_active_thread_percentage")
_integrator = config["production-settings"].get("gromacs-settings", {}).get("integrator", "sd").strip().lower()


_rbfe_network_cache = None


def _gromacs_num_lambda_for_edge(ligand1: str, ligand2: str) -> int:
    """Number of lambda windows (= MPI ranks for repex) for a specific RBFE edge.

    Unlike ABFE (one fixed window count from a static config schedule), RBFE's
    network.dat gives each edge its own num_lambda -- read directly rather than
    discovering it dynamically from directory globbing inside _run_gromacs_stages,
    since Snakemake resolves `resources:` before the rule body runs.
    """
    global _rbfe_network_cache
    if _rbfe_network_cache is None:
        network_file = Path(f"{config['working_directory']}/network/network.dat")
        _rbfe_network_cache = pd.read_csv(
            str(network_file),
            sep=r"\s+",
            names=["ligand1", "ligand2", "num_lambda", "lambda_windows"],
        )
    row = _rbfe_network_cache[
        (_rbfe_network_cache["ligand1"] == ligand1) & (_rbfe_network_cache["ligand2"] == ligand2)
    ]
    if row.empty:
        raise WorkflowError(f"No network.dat entry found for edge {ligand1}~{ligand2}.")
    return int(row["num_lambda"].iloc[0])


def _get_rbfe_pairs():
    """Read all ligand pairs from the network file (requires network prep to have run)."""
    network_file = Path(f"{config['working_directory']}/network/network.dat")
    data = pd.read_csv(
        str(network_file),
        sep=r"\s+",
        names=["ligand1", "ligand2", "num_lambda", "lambda_windows"],
    )
    return [f"{row['ligand1']}~{row['ligand2']}" for _, row in data.iterrows()]


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
        if cfg.get("num_energy_neighbours") is not None:
            args.append(f"--num-energy-neighbours {cfg['num_energy_neighbours']}")
        if cfg.get("use_dispersion_correction", False):
            args.append("--use-dispersion-correction")
        if cfg.get("oversubscription_factor", 1) > 1:
            args.append(f"--oversubscription-factor {cfg['oversubscription_factor']}")
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
        if config["production-settings"].get("restart", False):
            args.append("--restart")
        args.append(f"--network-location {config['working_directory']}/network")
    elif _engine == "amber":
        cfg = config["production-settings"].get("amber-settings", {})
        amber_leg_cfg = cfg.get(f"{leg}-leg-settings", cfg)
        args.append(f"--runtime {amber_leg_cfg['runtime']}")
        if amber_leg_cfg.get("timestep"):
            args.append(f"--timestep {amber_leg_cfg['timestep']}")
        if amber_leg_cfg.get("temperature"):
            args.append(f"--temperature {amber_leg_cfg['temperature']}")
        if amber_leg_cfg.get("pressure"):
            args.append(f"--pressure {amber_leg_cfg['pressure']}")
        if amber_leg_cfg.get("use-modified-dummies", False):
            args.append("--use-modified-dummies")
        if amber_leg_cfg.get("report-interval"):
            args.append(f"--report-interval {amber_leg_cfg['report-interval']}")
        if amber_leg_cfg.get("restart-interval"):
            args.append(f"--restart-interval {amber_leg_cfg['restart-interval']}")
        amber_runner = cfg.get("runner", "standard")
        args.append(f"--runner {amber_runner}")
        if amber_runner == "repex":
            args.append(f"--repex-frequency {cfg.get('repex-frequency', 1000)}")
            if cfg.get("exe"):
                args.append(f"--amber-exe {cfg['exe']}")
        args.append(f"--network-location {config['working_directory']}/network")
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
        args.append(f"--runner {_gromacs_runner}")
        args.append(f"--repex-frequency {_repex_frequency}")
        if _oversubscribe:
            args.append("--oversubscribe")
        args.append(f"--network-location {config['working_directory']}/network")

    return f"""
    python workflow/scripts/rbfe/production.py --input {input.file} --output-directory {output_directory} {" ".join(args)}
    """


def _run_gromacs_stages(output_directory, repex=False, repex_frequency=1000):
    """Run GROMACS minimisation, heating, equilibration, and production stages.

    When repex=True the production step runs all lambda windows together via
    ``gmx mdrun -multidir -replex`` instead of independent per-window runs.
    Min/heat/eq always run per-window regardless of the runner mode.
    """
    # Restart path: skip min/heat/eq and continue production from checkpoint.
    # production.py (BSS setup_only) has already regenerated gromacs.mdp with
    # the new nsteps before this function is called.
    if config["production-settings"].get("restart", False):
        if repex:
            raise NotImplementedError("Restart is not yet supported for GROMACS repex.")
        prod_path = Path(output_directory)
        lambda_values = sorted(
            [d.name.split("_")[1] for d in prod_path.glob("lambda_*") if d.is_dir()],
            key=float,
        )
        if lambda_values:
            print("Restarting GROMACS production from checkpoint")
            for lv in lambda_values:
                d = prod_path / f"lambda_{lv}"
                cpt_arg = f"-t {d}/gromacs.cpt" if (d / "gromacs.cpt").exists() else ""
                shell(f"gmx grompp -f {d}/gromacs.mdp -c {d}/gromacs.gro {cpt_arg} -p {d}/gromacs.top -o {d}/gromacs.tpr 2>&1 | tee {d}/grompp.log")
                shell(f"gmx mdrun -ntmpi 1 -deffnm {d}/gromacs 2>&1 | tee {d}/mdrun.log")
            return

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

    if repex:
        # HREX production: all lambda windows share one gmx mdrun -multidir invocation.
        # BSS has already written the shared topology to output_directory/gromacs.top
        # and per-lambda MDPs to output_directory/lambda_*/gromacs.mdp.
        # Re-run grompp for each lambda using equilibrated coordinates, then launch
        # all windows together with -replex.
        print("Running HREX production (grompp per lambda, then multidir mdrun)")
        shared_top = Path(f"{output_directory}/gromacs.top")
        for lambda_value in lambda_values:
            lam_dir = f"{output_directory}/lambda_{lambda_value}"
            eq_gro = f"{output_directory}/eq/lambda_{lambda_value}/gromacs.gro"
            top_file = str(shared_top) if shared_top.exists() else f"{lam_dir}/gromacs.top"
            shell(
                f"gmx grompp -f {lam_dir}/gromacs.mdp -c {eq_gro} -p {top_file} "
                f"-o {lam_dir}/gromacs.tpr -maxwarn 1 2>&1 | tee {lam_dir}/grompp.log"
            )
        n_replicas = len(lambda_values)
        multidir = " ".join(f"lambda_{lv}" for lv in lambda_values)
        shell(
            f"cd {output_directory} && "
            # mpi=srun (see resources: in the production_bound/production_free rules)
            # makes the slurm-jobstep executor skip its own nested srun wrapper and run
            # this shell command directly -- but the plugin unconditionally strips
            # SLURM_* env vars first, including ones mpirun's own PRRTE runtime needs
            # to detect the allocation (SLURM_NODELIST, SLURM_TASKS_PER_NODE). Re-export
            # them here -- values are trivial and static for this single-node cluster.
            f"export SLURM_NODELIST=$(hostname) SLURM_TASKS_PER_NODE={n_replicas} SLURM_JOB_NUM_NODES=1 SLURM_NNODES=1 && "
            f"mpirun {'--oversubscribe ' if _oversubscribe else ''}"
            f"-mca opal_cuda_support 1 -x OMP_NUM_THREADS=1 "
            f"{'-x CUDA_MPS_ACTIVE_THREAD_PERCENTAGE=' + str(_gromacs_mps_pct) + ' ' if _gromacs_mps_pct else ''}"
            f"-np {n_replicas} gmx_mpi mdrun -deffnm gromacs "
            # -nbfe (separate free-energy nonbonded GPU kernel) requires GROMACS 2026+.
            f"-nb gpu -nbfe gpu "
            f"{'-pme gpu ' if _gromacs_oversubscription_factor > 1 or _integrator == 'md' else ''}-bonded gpu {'-update gpu ' if _integrator == 'md' else ''}-nstlist 100 -cpt -1 "
            f"-c gromacs_out.gro -multidir {multidir} -replex {repex_frequency} -nex {_nex} "
            f"2>&1 | tee mdrun.log"
        )
    else:
        print("Running production")
        for lambda_value in lambda_values:
            d = f"{output_directory}/lambda_{lambda_value}"
            prev_gro = f"{output_directory}/eq/lambda_{lambda_value}/gromacs.gro"
            shell(f"gmx grompp -f {d}/gromacs.mdp -c {prev_gro} -p {d}/gromacs.top -o {d}/gromacs.tpr 2>&1 | tee {d}/grompp.log")
            shell(f"gmx mdrun -ntmpi 1 -deffnm {d}/gromacs 2>&1 | tee {d}/mdrun.log")

    # Remove intermediate directories as they confuse the analysis
    shell(f"rm -rf {output_directory}/minimisation {output_directory}/heat {output_directory}/eq")


# Replica barrier
# ================
#
# Ensures ALL edges complete both legs of replica N before any edge starts
# replica N+1. This guarantees at least one result per edge before additional
# replicas are run — prioritising results over throughput.

rule replica_barrier:
    input:
        bound=lambda wc: expand(
            f"{config['working_directory']}/production/{_engine}/{{pair}}/bound_{wc.replica}/.done",
            pair=_get_rbfe_pairs(),
        ),
        free=lambda wc: expand(
            f"{config['working_directory']}/production/{_engine}/{{pair}}/free_{wc.replica}/.done",
            pair=_get_rbfe_pairs(),
        ),
    output:
        touch(
            Path(f"{config['working_directory']}/production/{_engine}/.replica_{{replica}}_barrier")
        ),
    priority: 3


rule production_bound:
    priority: 2
    input:
        file = Path(f"{config['working_directory']}/rbfe_prepared/bound/{{ligand1}}~{{ligand2}}.bss"),
        prev_replica = lambda wc: [] if int(wc.replica_number) == 0 else [
            f"{config['working_directory']}/production/{_engine}/.replica_{int(wc.replica_number) - 1}_barrier",
        ],
    output:
        done = Path(f"{config['working_directory']}/production/{_engine}/{{ligand1}}~{{ligand2}}/bound_{{replica_number}}/.done")
    threads: config["simulation_threads"]
    resources:
        gpu=config["production-settings"].get("somd2-settings", {}).get("gpus_per_job", 1) if _engine == "somd2" else (_gromacs_gpus_per_job if _gromacs_runner == "repex" else 1),
        mpi=lambda wc: "srun" if (_engine == "gromacs" and _gromacs_runner == "repex") else None,
        tasks_per_node=lambda wc: _gromacs_num_lambda_for_edge(wc.ligand1, wc.ligand2) if (_engine == "gromacs" and _gromacs_runner == "repex") else None,
        cpus_per_task=lambda wc: 1 if (_engine == "gromacs" and _gromacs_runner == "repex") else None
    log:
        Path(f"{config['working_directory']}/logs/{{ligand1}}~{{ligand2}}_production_bound_{{replica_number}}.log")
    run:
        python_cmd = create_python_script_call(wildcards, input, "bound")
        shell(python_cmd)
        if _engine == "gromacs":
            output_directory = str(Path(f"{config['working_directory']}/production/{_engine}/{wildcards.ligand1}~{wildcards.ligand2}/bound_{wildcards.replica_number}"))
            _run_gromacs_stages(output_directory, repex=(_gromacs_runner == "repex"), repex_frequency=_repex_frequency)
        shell(f"touch {output.done}")


rule production_free:
    priority: 2
    input:
        file = Path(f"{config['working_directory']}/rbfe_prepared/free/{{ligand1}}~{{ligand2}}.bss"),
        prev_replica = lambda wc: [] if int(wc.replica_number) == 0 else [
            f"{config['working_directory']}/production/{_engine}/.replica_{int(wc.replica_number) - 1}_barrier",
        ],
    output:
        done = Path(f"{config['working_directory']}/production/{_engine}/{{ligand1}}~{{ligand2}}/free_{{replica_number}}/.done")
    threads: config["simulation_threads"]
    resources:
        gpu=config["production-settings"].get("somd2-settings", {}).get("gpus_per_job", 1) if _engine == "somd2" else (_gromacs_gpus_per_job if _gromacs_runner == "repex" else 1),
        mpi=lambda wc: "srun" if (_engine == "gromacs" and _gromacs_runner == "repex") else None,
        tasks_per_node=lambda wc: _gromacs_num_lambda_for_edge(wc.ligand1, wc.ligand2) if (_engine == "gromacs" and _gromacs_runner == "repex") else None,
        cpus_per_task=lambda wc: 1 if (_engine == "gromacs" and _gromacs_runner == "repex") else None
    log:
        Path(f"{config['working_directory']}/logs/{{ligand1}}~{{ligand2}}_production_free_{{replica_number}}.log")
    run:
        python_cmd = create_python_script_call(wildcards, input, "free")
        shell(python_cmd)
        if _engine == "gromacs":
            output_directory = str(Path(f"{config['working_directory']}/production/{_engine}/{wildcards.ligand1}~{wildcards.ligand2}/free_{wildcards.replica_number}"))
            _run_gromacs_stages(output_directory, repex=(_gromacs_runner == "repex"), repex_frequency=_repex_frequency)
        shell(f"touch {output.done}")
