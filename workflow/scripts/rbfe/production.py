import argparse
import math
from pathlib import Path

import BioSimSpace as BSS
import pandas as pd
import sire as sr


def main():
    parser = argparse.ArgumentParser(description="Production script")
    parser.add_argument(
        "--input", type=str, help="Full path to the input file", required=True
    )
    parser.add_argument(
        "--output-directory",
        type=str,
        help="Full path to the directory in which the output is to be saved",
        required=True,
    )
    parser.add_argument(
        "--engine",
        type=str,
        help="The engine to use for production",
        choices=["somd", "gromacs", "amber", "somd2"],
    )
    parser.add_argument(
        "--runtime",
        type=str,
        help="Number of total runtime per window.",
        default="1ns",
    )
    parser.add_argument(
        "--timestep",
        type=str,
        help="Timestep to use for production, only use 4fs if HMR has been applied.",
        default="2fs",
    )
    parser.add_argument(
        "--temperature",
        type=str,
        help="The temperature for the production, e.g., '300K'.",
        default="300K",
    )
    parser.add_argument(
        "--pressure",
        type=str,
        help="The pressure for the production, e.g., '1bar'. If not specified, no pressure coupling is applied.",
        default="1atm",
    )
    parser.add_argument(
        "--restraint-string",
        type=str,
        choices=["backbone", "heavy", "all"],
        help="The atoms to restrain during production, apply this argument OR `restraint-indices`, not both.\
             If neither are specified the no restraints will be applied.",
        default=None,
    )
    parser.add_argument(
        "--restraint-indices",
        nargs="*",
        type=int,
        help="A list of indices of atoms to restrain during production. Apply this argument OR `restraint-string`, not both.\
             If neither are specified no restraints will be applied.",
        default=None,
    )
    parser.add_argument(
        "--restraint-force-constant",
        type=float,
        help="Force constant for the restraint, in kcal_per_mol / angstrom^2. Default is 5.0.",
        default=5.0,
    )
    parser.add_argument(
        "--use-modified-dummies",
        action="store_true",
        help="Whether to apply modifications to dummy atoms",
        default=False,
    )
    parser.add_argument(
        "--report-interval",
        type=str,
        help="Interval at which to save frames and energies (e.g., '1ps'). Converted to steps using the timestep.",
        default="1ps",
    )
    parser.add_argument(
        "--restart-interval",
        type=str,
        help="Interval at which to save restart files (e.g., '500ps'). Converted to steps using the timestep.",
        default="500ps",
    )
    parser.add_argument(
        "--network-location",
        type=str,
        help="Full path to the network file",
        default="network/network.dat",
    )
    parser.add_argument(
        "--cutoff-type",
        type=str,
        default="rf",
        help="Electrostatics cutoff type for SOMD2 (RF or PME).",
    )
    parser.add_argument(
        "--cutoff",
        type=str,
        default="12A",
        help="Cutoff distance for SOMD2 (e.g., '12A').",
    )
    parser.add_argument(
        "--energy-frequency",
        type=str,
        default="1ps",
        help="Energy sampling interval for SOMD2 (e.g., '1ps').",
    )
    parser.add_argument(
        "--frame-frequency",
        type=str,
        default=None,
        help="Trajectory frame saving interval for SOMD2. Defaults to the full runtime (single frame).",
    )
    parser.add_argument(
        "--checkpoint-frequency",
        type=str,
        default=None,
        help="Checkpoint saving interval for SOMD2. Defaults to the full runtime (single checkpoint).",
    )
    parser.add_argument(
        "--integrator",
        type=str,
        default="langevin_middle",
        help="Integrator type for SOMD2.",
    )
    parser.add_argument(
        "--shift-delta",
        type=str,
        default="1.5A",
        help="Soft-core shift delta parameter for SOMD2.",
    )
    parser.add_argument(
        "--perturbable-constraint",
        type=str,
        default="h_bonds_not_heavy_perturbed",
        help="Constraint type for perturbable molecules in SOMD2.",
    )
    parser.add_argument(
        "--runner",
        type=str,
        choices=["repex", "standard"],
        default="repex",
        help="Runner type for SOMD2: 'repex' (replica exchange, default) or 'standard' (independent windows).",
    )

    args = parser.parse_args()

    if args.restraint_string is not None and args.restraint_indices is not None:
        raise ValueError(
            "You cannot specify both 'restraint-string' and 'restraint-indices'. Please choose one."
        )

    # first we need to convert the inputs into the appropriate BSS objects
    try:
        runtime = BSS.Types.Time(args.runtime)
    except ValueError:
        raise ValueError(
            f"Invalid runtime '{args.runtime}'. Please specify a valid time duration (e.g., '10ps', '1ns')."
        )
    try:
        timestep = BSS.Types.Time(args.timestep)
    except ValueError:
        raise ValueError(
            f"Invalid timestep '{args.timestep}'. Please specify a valid time (e.g., '2fs', '1fs')."
        )
    try:
        temperature = BSS.Types.Temperature(args.temperature)
    except ValueError:
        raise ValueError(
            f"Invalid temperature '{args.temperature}'. Please specify a valid temperature (e.g., '300K')."
        )
    if args.pressure is not None:
        try:
            pressure = BSS.Types.Pressure(args.pressure)
        except ValueError:
            raise ValueError(
                f"Invalid pressure '{args.pressure}'. Please specify a valid pressure (e.g., '1bar')."
            )
    else:
        pressure = None

    # Convert time-based intervals to integer steps
    report_interval = math.ceil(sr.u(args.report_interval) / sr.u(args.timestep))
    restart_interval = math.ceil(sr.u(args.restart_interval) / sr.u(args.timestep))

    if args.restraint_indices:
        restraint = args.restraint_indices
        force_constant = args.restraint_force_constant
    elif args.restraint_string:
        restraint = args.restraint_string
        force_constant = args.restraint_force_constant
    else:
        restraint = None
        force_constant = 0.0

    # first we need to read the network file
    network_dir = Path(f"{args.network_location}/network.dat")
    network_df = pd.read_csv(
        network_dir,
        sep=r"\s+",
        header=None,
        names=["lig1", "lig2", "num_lam", "lam_array"],
    )
    # now figure out the lig1 and lig2 values from our input filename
    input_pert = Path(args.input).stem
    lig1, lig2 = input_pert.split("~")
    # now find the relevant line in `network`
    settings = network_df[
        (network_df["lig1"] == lig1) & (network_df["lig2"] == lig2)
    ].iloc[0]
    # get lambda array (this should just be range(0,1,num_lam), but in theory this way we can support  wider range of schedules)
    lam_arr = settings["lam_array"]
    # convert from string to list of floats
    lam_arr = [float(x) for x in lam_arr.strip("[]").split(",")]

    # Now split into three - one for amber and somd, another for gromacs, and a third for somd2.
    engine = args.engine.strip().lower()
    working_dir = str(Path(args.output_directory))
    if engine not in ["somd2", "gromacs"]:

        sire_system = sr.stream.load(args.input)
        # now setup dommy mods if true
        if args.use_modified_dummies:

            from loguru import logger

            # create an output folder for the logs
            log_dir = Path(working_dir) / "ghostly_logs"
            log_dir.mkdir(parents=True, exist_ok=True)

            logger.remove()
            logger.add(log_dir / f"{lig1}~{lig2}_ghostly.log", level="DEBUG")

            from ghostly import modify

            sire_system, _ = modify(sire_system)

        # now to BSS
        system = BSS._SireWrappers.System(sire_system._system)
        if system is None:
            raise ValueError(
                f"Failed to load system from {args.input}. Please check the file."
            )
        protocol = BSS.Protocol.FreeEnergyProduction(
            lam_vals=lam_arr,
            timestep=timestep,
            runtime=runtime,
            temperature=temperature,
            pressure=pressure,
            report_interval=report_interval,
            restart_interval=restart_interval,
            restraint=restraint,
            force_constant=force_constant,
        )
        if engine == "amber":
            process_production = BSS.FreeEnergy.Relative(
                system,
                protocol,
                engine="amber",
                work_dir=working_dir,
                is_gpu=True,
            )
        elif engine == "somd":
            process_production = BSS.FreeEnergy.Relative(
                system,
                protocol,
                engine="somd",
                work_dir=working_dir,
                extra_options={"minimise": "True", "gpu": "0"},
                no_dummy_modifications=args.use_modified_dummies,  # This arg needs to be True if we already applied ghostly modifications
            )
        process_production.run()
        process_production.wait()
    elif engine == "gromacs":
        sire_system = sr.stream.load(args.input)
        # now setup dommy mods if true
        if args.use_modified_dummies:

            from loguru import logger

            # create an output folder for the logs
            log_dir = Path(working_dir) / "ghostly_logs"
            log_dir.mkdir(parents=True, exist_ok=True)

            logger.remove()
            logger.add(log_dir / f"{lig1}~{lig2}_ghostly.log", level="DEBUG")

            from ghostly import modify

            sire_system, _ = modify(sire_system)

        # now to BSS
        system = BSS._SireWrappers.System(sire_system._system)
        if system is None:
            raise ValueError(
                f"Failed to load system from {args.input}. Please check the file."
            )
        # Gromacs needs a pre-production minimisation and equilibration step
        # We can't do this in python, so we use BSS to setup the directories,
        # then our snakemake rule will handle the actual running of Gromacs
        min_protocol = BSS.Protocol.FreeEnergyMinimisation(lam_vals=lam_arr)
        BSS.FreeEnergy.Relative(
            system,
            min_protocol,
            engine="gromacs",
            work_dir=working_dir + "/minimisation",
            extra_args={
                "-ntmpi": "1",
            },
            ignore_warnings=True,
            setup_only=True,
        )
        nvt_protocol = BSS.Protocol.FreeEnergyEquilibration(
            lam_vals=lam_arr,
            pressure=None,
            temperature_start=BSS.Types.Temperature("100K"),
            temperature_end=BSS.Types.Temperature("300K"),
            timestep=BSS.Types.Time("2fs"),
            runtime=BSS.Types.Time("10ps"),
            report_interval=2500,  # as little as possible
            restart_interval=2500,
        )
        BSS.FreeEnergy.Relative(
            system,
            nvt_protocol,
            engine="gromacs",
            work_dir=working_dir + "/heat",
            extra_args={
                "-ntmpi": "1",
            },
            ignore_warnings=True,
            setup_only=True,
        )
        npt_protocol = BSS.Protocol.FreeEnergyEquilibration(
            lam_vals=lam_arr,
            pressure=1 * BSS.Units.Pressure.atm,
            temperature=BSS.Types.Temperature("300K"),
            timestep=BSS.Types.Time("2fs"),
            runtime=BSS.Types.Time("20ps"),
            report_interval=2500,
            restart_interval=2500,
        )
        BSS.FreeEnergy.Relative(
            system,
            npt_protocol,
            engine="gromacs",
            work_dir=working_dir + "/eq",
            extra_args={
                "-ntmpi": "1",
            },
            ignore_warnings=True,
            setup_only=True,
        )
        protocol_production = BSS.Protocol.FreeEnergyProduction(
            lam_vals=lam_arr,
            timestep=timestep,
            runtime=runtime,
            temperature=temperature,
            pressure=pressure,
            report_interval=report_interval,
            restart_interval=restart_interval,
            restraint=restraint,
            force_constant=force_constant,
        )
        process_production = BSS.FreeEnergy.Relative(
            system,
            protocol_production,
            engine="gromacs",
            work_dir=working_dir,
            extra_args={
                "-ntmpi": "1",
            },
            ignore_warnings=True,
            setup_only=True,
        )
    else:
        from somd2.config import Config
        from somd2.runner import RepexRunner, Runner

        system = sr.stream.load(args.input)

        # first step is to create restraints
        if restraint is not None:
            # the only universally supported restaints are positional, so these are the only ones we will worry about
            # Use the in-built BSS convenience function to find the indices if a string is given
            if type(restraint) is str:
                system_bss = BSS._SireWrappers.System(system._system)
                restraint_indices = system_bss.getRestraintAtoms(
                    restraint, is_absolute=True
                )
                # don't need to keep the BSS system
                del system_bss
            elif type(restraint) is list:
                # if restraint is a list, we assume it is a list of indices
                restraint_indices = restraint
            else:
                # if no restraint is specified, we do not apply any restraints
                restraint_indices = None
            # now create the sire restraints
            if restraint_indices is not None:
                # convert the force constant to the appropriate unit
                force_constant = sr.u(f"{force_constant} kcal_per_mol / angstrom^2")
                # create positional restraints
                restraints = sr.restraints.positional(
                    system, restraint_indices, k=force_constant
                )
        else:
            restraints = None

        runner_type = args.runner.strip().lower()
        frame_frequency = args.frame_frequency if args.frame_frequency else args.runtime
        checkpoint_frequency = (
            args.checkpoint_frequency if args.checkpoint_frequency else args.runtime
        )

        # somd2 has in-built type conversion, so we can pass arguments as strings
        config = Config(
            lambda_values=lam_arr,
            output_directory=working_dir,
            runtime=args.runtime,
            temperature=args.temperature,
            pressure=args.pressure,
            cutoff_type=args.cutoff_type,
            cutoff=args.cutoff,
            timestep=args.timestep,
            restraints=restraints,
            ghost_modifications=args.use_modified_dummies,
            energy_frequency=args.energy_frequency,
            frame_frequency=frame_frequency,
            checkpoint_frequency=checkpoint_frequency,
            integrator=args.integrator,
            shift_delta=args.shift_delta,
            perturbable_constraint=args.perturbable_constraint,
            hmr=False,  # HMR is applied in the pre-production step
            overwrite=True,  # Allow reruns without manual cleanup
            restart=True,  # Allow restarting from checkpoints
            replica_exchange=(runner_type == "repex"),
        )

        if runner_type == "repex":
            # RepexRunner uses threads, safe to call run() directly
            runner = RepexRunner(system, config)
            runner.run()
            return None
        else:
            # Standard Runner uses the spawn start method, so runner.run()
            # must be called inside the if __name__ == "__main__" guard.
            runner = Runner(system, config)
            return runner

    return None


if __name__ == "__main__":
    runner = main()
    if runner is not None:
        runner.run()

        # Verify that SOMD2 produced output — runner.run() can silently fail
        output_path = Path(runner._config.output_directory)
        energy_files = list(output_path.glob("energy_traj_*.parquet"))
        if not energy_files:
            raise RuntimeError(
                f"SOMD2 produced no energy trajectory files in {output_path}. "
                f"Check {output_path / 'log.txt'} for errors."
            )
