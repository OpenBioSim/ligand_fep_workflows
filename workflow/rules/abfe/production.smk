"""
Snakemake rules for ABFE equilibration and production simulations.

ABFE uses a unified single-leg protocol where bonded/coul/vdw
lambdas are all controlled in a single DataFrame schedule. This uses
`perturbation_type="full"` which applies bonded soft-core by default.

Equilibration: BioSimSpace setup -> GROMACS minimisation -> NVT heating -> NPT eq
Production: GROMACS NPT production MD (using equilibrated coordinates)

Bound leg: 21 windows (bonded+coul ramp, then vdw ramp)
Free leg: 21 windows (coul ramp, then vdw ramp)

Each leg runs at multiple lambda values and for multiple replicas
to estimate the free energy change and associated uncertainty.
"""

from pathlib import Path

_engine = config.get("engine", config["production-settings"].get("engine", "gromacs")).strip().lower()


# Bound leg equilibration
# =======================

rule equilibrate_bound:
    """
    Equilibration for bound leg using unified protocol.

    Sets up BioSimSpace simulations and runs GROMACS minimisation,
    NVT heating, and NPT equilibration for each lambda window.
    Boresch restraints are applied throughout.

    Replica N waits for replica N-1's production to complete,
    so all ligands get a first result before any starts a second replica.
    """
    input:
        system=Path(f"{config['working_directory']}/abfe_prepared")
        / "{ligand}_bound.bss",
        restraint=Path(f"{config['working_directory']}/restraints")
        / "{ligand}_restraint.json",
        prev_replica=lambda wc: [] if int(wc.replica) == 0 else [
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/bound_{int(wc.replica) - 1}/.done"
        ],
    output:
        done=Path(
            f"{config['working_directory']}/equilibration/{{ligand}}/bound_{{replica}}/.done"
        ),
    threads:
        config["simulation_threads"]
    resources:
        gpu=1
    log:
        Path(f"{config['working_directory']}/logs")
        / "{ligand}_equilibrate_bound_{replica}.log",
    params:
        script=Path("workflow/scripts/abfe/production.py"),
        lambda_schedule=lambda wc: config["production-settings"]["gromacs-settings"]["lambda_schedules"]["bound"],
        runtime=lambda wc: config["production-settings"]["gromacs-settings"]["bound-leg-settings"].get(
            "runtime", "2ns"
        ),
        timestep=lambda wc: config["production-settings"]["gromacs-settings"]["bound-leg-settings"].get(
            "timestep", "4fs"
        ),
        temperature=lambda wc: config["production-settings"]["gromacs-settings"]["bound-leg-settings"].get(
            "temperature", "300K"
        ),
        pressure=lambda wc: config["production-settings"]["gromacs-settings"]["bound-leg-settings"].get(
            "pressure", "1bar"
        ),
        report_interval=lambda wc: config["production-settings"][
            "gromacs-settings"
        ]["bound-leg-settings"].get("report-interval", "1ps"),
        restart_interval=lambda wc: config["production-settings"][
            "gromacs-settings"
        ]["bound-leg-settings"].get("restart-interval", "500ps"),
    run:
        import json
        from pathlib import Path

        output_directory = str(Path(str(output.done)).resolve().parent)

        # Write lambda schedule to a temp file (avoids brace escaping issues with shell())
        schedule_file = Path(output_directory) / "_lambda_schedule.json"
        Path(output_directory).mkdir(parents=True, exist_ok=True)
        with open(schedule_file, "w") as f:
            json.dump(params.lambda_schedule, f)

        # BioSimSpace setup phase (creates min/heat/eq/production directories)
        shell(
            f"echo 'Setting up bound leg for {wildcards.ligand} replica {wildcards.replica}' && "
            f"python {params.script} "
            f"--input {input.system} "
            f"--output-directory {output_directory} "
            f"--ligand-name {wildcards.ligand} "
            f"--leg bound "
            f"--lambda-schedule-file {schedule_file} "
            f"--restraint-file {input.restraint} "
            f"--runtime {params.runtime} "
            f"--timestep {params.timestep} "
            f"--temperature {params.temperature} "
            f"--pressure {params.pressure} "
            f"--report-interval {params.report_interval} "
            f"--restart-interval {params.restart_interval} "
            f"2>&1 | tee {log}"
        )

        # Clean up schedule file
        schedule_file.unlink(missing_ok=True)

        # Discover lambda values from minimisation directories
        outdir_path_min = Path(output_directory) / "minimisation"
        lambda_values = [d.name.split('_')[1] for d in outdir_path_min.glob("lambda_*") if d.is_dir()]
        if not lambda_values:
            raise FileNotFoundError(f"No minimisation directories found in {outdir_path_min}")
        lambda_values.sort(key=float)

        # Run GROMACS equilibration stages
        # All gmx commands cd into the target directory first so that
        # auxiliary files (mdout.mdp, step*.pdb crash dumps) stay contained.
        print("Minimising")
        for lambda_value in lambda_values:
            d = f"{output_directory}/minimisation/lambda_{lambda_value}"
            shell(f"cd {d} && gmx grompp -f gromacs.mdp -c gromacs_ref.gro -p gromacs.top -o gromacs.tpr -maxwarn 1 2>&1 | tee grompp.log")
            shell(f"cd {d} && gmx mdrun -ntmpi 1 -deffnm gromacs 2>&1 | tee mdrun.log")

        print("Heating")
        for lambda_value in lambda_values:
            d = f"{output_directory}/heat/lambda_{lambda_value}"
            prev_gro = f"{output_directory}/minimisation/lambda_{lambda_value}/gromacs.gro"
            shell(f"cd {d} && gmx grompp -f gromacs.mdp -c {prev_gro} -p gromacs.top -o gromacs.tpr -maxwarn 1 2>&1 | tee grompp.log")
            shell(f"cd {d} && gmx mdrun -ntmpi 1 -deffnm gromacs 2>&1 | tee mdrun.log")

        print("Equilibrating")
        for lambda_value in lambda_values:
            d = f"{output_directory}/eq/lambda_{lambda_value}"
            prev_gro = f"{output_directory}/heat/lambda_{lambda_value}/gromacs.gro"
            shell(f"cd {d} && gmx grompp -f gromacs.mdp -c {prev_gro} -p gromacs.top -o gromacs.tpr -maxwarn 1 2>&1 | tee grompp.log")
            shell(f"cd {d} && gmx mdrun -ntmpi 1 -deffnm gromacs 2>&1 | tee mdrun.log")

        # Clean up minimisation and heating directories (keep eq/ and lambda_*/ for production)
        shell(f"rm -rf {output_directory}/minimisation {output_directory}/heat")

        # Mark as complete
        shell(f"touch {output.done}")


# Free leg equilibration
# ======================

rule equilibrate_free:
    """
    Equilibration for free leg using unified protocol.

    Sets up BioSimSpace simulations and runs GROMACS minimisation,
    NVT heating, and NPT equilibration for each lambda window.
    No restraints needed for free leg.

    Replica N waits for replica N-1's production to complete,
    so all ligands get a first result before any starts a second replica.
    """
    input:
        system=Path(f"{config['working_directory']}/abfe_prepared")
        / "{ligand}_free.bss",
        prev_replica=lambda wc: [] if int(wc.replica) == 0 else [
            f"{config['working_directory']}/production/{_engine}/{wc.ligand}/free_{int(wc.replica) - 1}/.done"
        ],
    output:
        done=Path(
            f"{config['working_directory']}/equilibration/{{ligand}}/free_{{replica}}/.done"
        ),
    threads:
        config["simulation_threads"]
    resources:
        gpu=1
    log:
        Path(f"{config['working_directory']}/logs")
        / "{ligand}_equilibrate_free_{replica}.log",
    params:
        script=Path("workflow/scripts/abfe/production.py"),
        lambda_schedule=lambda wc: config["production-settings"]["gromacs-settings"]["lambda_schedules"]["free"],
        runtime=lambda wc: config["production-settings"]["gromacs-settings"]["free-leg-settings"].get(
            "runtime", "2ns"
        ),
        timestep=lambda wc: config["production-settings"]["gromacs-settings"]["free-leg-settings"].get(
            "timestep", "4fs"
        ),
        temperature=lambda wc: config["production-settings"]["gromacs-settings"]["free-leg-settings"].get(
            "temperature", "300K"
        ),
        pressure=lambda wc: config["production-settings"]["gromacs-settings"]["free-leg-settings"].get(
            "pressure", "1bar"
        ),
        report_interval=lambda wc: config["production-settings"][
            "gromacs-settings"
        ]["free-leg-settings"].get("report-interval", "1ps"),
        restart_interval=lambda wc: config["production-settings"][
            "gromacs-settings"
        ]["free-leg-settings"].get("restart-interval", "500ps"),
    run:
        import json
        from pathlib import Path

        output_directory = str(Path(str(output.done)).resolve().parent)

        # Write lambda schedule to a temp file (avoids brace escaping issues with shell())
        schedule_file = Path(output_directory) / "_lambda_schedule.json"
        Path(output_directory).mkdir(parents=True, exist_ok=True)
        with open(schedule_file, "w") as f:
            json.dump(params.lambda_schedule, f)

        # BioSimSpace setup phase (creates min/heat/eq/production directories)
        shell(
            f"echo 'Setting up free leg for {wildcards.ligand} replica {wildcards.replica}' && "
            f"python {params.script} "
            f"--input {input.system} "
            f"--output-directory {output_directory} "
            f"--ligand-name {wildcards.ligand} "
            f"--leg free "
            f"--lambda-schedule-file {schedule_file} "
            f"--runtime {params.runtime} "
            f"--timestep {params.timestep} "
            f"--temperature {params.temperature} "
            f"--pressure {params.pressure} "
            f"--report-interval {params.report_interval} "
            f"--restart-interval {params.restart_interval} "
            f"2>&1 | tee {log}"
        )

        # Clean up schedule file
        schedule_file.unlink(missing_ok=True)

        # Discover lambda values from minimisation directories
        outdir_path_min = Path(output_directory) / "minimisation"
        lambda_values = [d.name.split('_')[1] for d in outdir_path_min.glob("lambda_*") if d.is_dir()]
        if not lambda_values:
            raise FileNotFoundError(f"No minimisation directories found in {outdir_path_min}")
        lambda_values.sort(key=float)

        # Run GROMACS equilibration stages
        # All gmx commands cd into the target directory first so that
        # auxiliary files (mdout.mdp, step*.pdb crash dumps) stay contained.
        print("Minimising")
        for lambda_value in lambda_values:
            d = f"{output_directory}/minimisation/lambda_{lambda_value}"
            shell(f"cd {d} && gmx grompp -f gromacs.mdp -c gromacs_ref.gro -p gromacs.top -o gromacs.tpr -maxwarn 1 2>&1 | tee grompp.log")
            shell(f"cd {d} && gmx mdrun -ntmpi 1 -deffnm gromacs 2>&1 | tee mdrun.log")

        print("Heating")
        for lambda_value in lambda_values:
            d = f"{output_directory}/heat/lambda_{lambda_value}"
            prev_gro = f"{output_directory}/minimisation/lambda_{lambda_value}/gromacs.gro"
            shell(f"cd {d} && gmx grompp -f gromacs.mdp -c {prev_gro} -p gromacs.top -o gromacs.tpr -maxwarn 1 2>&1 | tee grompp.log")
            shell(f"cd {d} && gmx mdrun -ntmpi 1 -deffnm gromacs 2>&1 | tee mdrun.log")

        print("Equilibrating")
        for lambda_value in lambda_values:
            d = f"{output_directory}/eq/lambda_{lambda_value}"
            prev_gro = f"{output_directory}/heat/lambda_{lambda_value}/gromacs.gro"
            shell(f"cd {d} && gmx grompp -f gromacs.mdp -c {prev_gro} -p gromacs.top -o gromacs.tpr -maxwarn 1 2>&1 | tee grompp.log")
            shell(f"cd {d} && gmx mdrun -ntmpi 1 -deffnm gromacs 2>&1 | tee mdrun.log")

        # Clean up minimisation and heating directories (keep eq/ and lambda_*/ for production)
        shell(f"rm -rf {output_directory}/minimisation {output_directory}/heat")

        # Mark as complete
        shell(f"touch {output.done}")


# Bound leg production
# ====================

rule production_bound:
    """
    Production for bound leg using equilibrated coordinates.

    Runs GROMACS production MD for each lambda window using the
    equilibrated system. Input files (MDP, topology) come from the
    equilibration directory; starting coordinates from NPT equilibration.
    """
    priority: 2
    input:
        eq_done=Path(
            f"{config['working_directory']}/equilibration/{{ligand}}/bound_{{replica}}/.done"
        ),
    output:
        done=Path(
            f"{config['working_directory']}/production/{_engine}/{{ligand}}/bound_{{replica}}/.done"
        ),
    threads:
        config["simulation_threads"]
    resources:
        gpu=1
    log:
        Path(f"{config['working_directory']}/logs")
        / "{ligand}_production_bound_{replica}.log",
    run:
        from pathlib import Path

        eq_dir = Path(str(input.eq_done)).resolve().parent
        prod_dir = Path(str(output.done)).resolve().parent
        prod_dir.mkdir(parents=True, exist_ok=True)

        # Discover lambda values from equilibration setup directories
        lambda_values = [
            d.name.split('_')[1]
            for d in eq_dir.glob("lambda_*")
            if d.is_dir() and d.name != "_lambda_schedule.json"
        ]
        if not lambda_values:
            raise FileNotFoundError(f"No lambda directories found in {eq_dir}")
        lambda_values.sort(key=float)

        print(f"Running production for {len(lambda_values)} lambda windows")
        for lambda_value in lambda_values:
            lam_prod = prod_dir / f"lambda_{lambda_value}"
            lam_prod.mkdir(exist_ok=True)

            eq_setup_dir = eq_dir / f"lambda_{lambda_value}"
            eq_gro = eq_dir / "eq" / f"lambda_{lambda_value}" / "gromacs.gro"

            # Run grompp from eq setup directory (so topology includes resolve correctly)
            # Output tpr goes to production directory
            shell(
                f"cd {eq_setup_dir} && "
                f"gmx grompp -f gromacs.mdp -c {eq_gro} -p gromacs.top "
                f"-o {lam_prod}/gromacs.tpr -maxwarn 1 "
                f"> {lam_prod}/grompp.log 2>&1"
            )
            # Run mdrun from production directory (crash files stay contained)
            shell(f"cd {lam_prod} && gmx mdrun -ntmpi 1 -deffnm gromacs > mdrun.log 2>&1")

        # Mark as complete
        shell(f"touch {output.done}")


# Free leg production
# ===================

rule production_free:
    """
    Production for free leg using equilibrated coordinates.

    Runs GROMACS production MD for each lambda window using the
    equilibrated system. Input files (MDP, topology) come from the
    equilibration directory; starting coordinates from NPT equilibration.
    """
    priority: 1
    input:
        eq_done=Path(
            f"{config['working_directory']}/equilibration/{{ligand}}/free_{{replica}}/.done"
        ),
    output:
        done=Path(
            f"{config['working_directory']}/production/{_engine}/{{ligand}}/free_{{replica}}/.done"
        ),
    threads:
        config["simulation_threads"]
    resources:
        gpu=1
    log:
        Path(f"{config['working_directory']}/logs")
        / "{ligand}_production_free_{replica}.log",
    run:
        from pathlib import Path

        eq_dir = Path(str(input.eq_done)).resolve().parent
        prod_dir = Path(str(output.done)).resolve().parent
        prod_dir.mkdir(parents=True, exist_ok=True)

        # Discover lambda values from equilibration setup directories
        lambda_values = [
            d.name.split('_')[1]
            for d in eq_dir.glob("lambda_*")
            if d.is_dir() and d.name != "_lambda_schedule.json"
        ]
        if not lambda_values:
            raise FileNotFoundError(f"No lambda directories found in {eq_dir}")
        lambda_values.sort(key=float)

        print(f"Running production for {len(lambda_values)} lambda windows")
        for lambda_value in lambda_values:
            lam_prod = prod_dir / f"lambda_{lambda_value}"
            lam_prod.mkdir(exist_ok=True)

            eq_setup_dir = eq_dir / f"lambda_{lambda_value}"
            eq_gro = eq_dir / "eq" / f"lambda_{lambda_value}" / "gromacs.gro"

            # Run grompp from eq setup directory (so topology includes resolve correctly)
            # Output tpr goes to production directory
            shell(
                f"cd {eq_setup_dir} && "
                f"gmx grompp -f gromacs.mdp -c {eq_gro} -p gromacs.top "
                f"-o {lam_prod}/gromacs.tpr -maxwarn 1 "
                f"> {lam_prod}/grompp.log 2>&1"
            )
            # Run mdrun from production directory (crash files stay contained)
            shell(f"cd {lam_prod} && gmx mdrun -ntmpi 1 -deffnm gromacs > mdrun.log 2>&1")

        # Mark as complete
        shell(f"touch {output.done}")
