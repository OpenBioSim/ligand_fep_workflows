import argparse
from pathlib import Path

import BioSimSpace as BSS
from _common import runProcess

parser = argparse.ArgumentParser(description="Minimisation script")
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
    "--outfile-name",
    type=str,
    help="Name of the output file (without extension).",
    required=True,
)
parser.add_argument(
    "--engine",
    type=str,
    help="The engine to use for minimisation",
    choices=["openmm", "gromacs", "amber"],
)
parser.add_argument(
    "--runtime",
    type=str,
    help="Number of minimisation steps to perform",
    default="10ps",
)
parser.add_argument(
    "--timestep",
    type=str,
    help="The timestep for the equilibration, e.g., '2fs'.",
    default="2fs",
)
parser.add_argument(
    "--temperature-start",
    type=str,
    help="The starting temperature for the equilibration, e.g., '300K'.",
    default="300K",
)
parser.add_argument(
    "--temperature-end",
    type=str,
    help="The ending temperature for the equilibration, e.g., '300K'.",
    default="300K",
)
parser.add_argument(
    "--temperature",
    type=str,
    help="The temperature for the equilibration, e.g., '300K'. If specified, overrides `temperature-start` and `temperature-end`.",
    default=None,
)
parser.add_argument(
    "--thermostat-time-constant",
    type=str,
    help="The thermostat time constant for the equilibration, e.g., '1ps'.",
    default="1ps",
)
parser.add_argument(
    "--pressure",
    type=str,
    help="The pressure for the equilibration, e.g., '1bar'.",
    default=None,
)
parser.add_argument(
    "--restraint-string",
    type=str,
    choices=["backbone", "heavy", "all"],
    help="The atoms to restrain during minimisation, apply this argument OR `restraint-indices`, not both.\
         If neither are specified the no restraints will be applied.",
    default=None,
)
parser.add_argument(
    "--restraint-indices",
    nargs="*",
    type=int,
    help="A list of indices of atoms to restrain during minimisation. Apply this argument OR `restraint-string`, not both.\
         If neither are specified no restraints will be applied.",
    default=None,
)
parser.add_argument(
    "--restraint-force-constant",
    type=float,
    help="Force constant for the restraint, in kcal_per_mol / angstrom^2. Default is 5.0.",
    default=5.0,
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
    temperature_start = BSS.Types.Temperature(args.temperature_start)
except ValueError:
    raise ValueError(
        f"Invalid temperature start '{args.temperature_start}'. Please specify a valid temperature (e.g., '300K')."
    )
try:
    temperature_end = BSS.Types.Temperature(args.temperature_end)
except ValueError:
    raise ValueError(
        f"Invalid temperature end '{args.temperature_end}'. Please specify a valid temperature (e.g., '300K')."
    )
if args.temperature is not None:
    try:
        temperature = BSS.Types.Temperature(args.temperature)
    except ValueError:
        raise ValueError(
            f"Invalid temperature '{args.temperature}'. Please specify a valid temperature (e.g., '300K')."
        )
else:
    temperature = None
if args.pressure is not None:
    try:
        pressure = BSS.Types.Pressure(args.pressure)
    except ValueError:
        raise ValueError(
            f"Invalid pressure '{args.pressure}'. Please specify a valid pressure (e.g., '1bar')."
        )
else:
    pressure = None
if args.thermostat_time_constant is not None:
    try:
        thermostat_time_constant = BSS.Types.Time(args.thermostat_time_constant)
    except ValueError:
        raise ValueError(
            f"Invalid thermostat time constant '{args.thermostat_time_constant}'. Please specify a valid time (e.g., '1ps')."
        )
else:
    thermostat_time_constant = None


if args.restraint_indices:
    restraint = args.restraint_indices
    force_constant = args.restraint_force_constant
elif args.restraint_string:
    restraint = args.restraint_string
    force_constant = args.restraint_force_constant
else:
    restraint = None
    force_constant = 0.0  # Can't be `None`, needs to be float or unit or string


system = BSS.Stream.load(args.input)

# make the output directory if it doesn't exist
output_directory = Path(args.output_directory)
output_directory.mkdir(parents=True, exist_ok=True)

protocol = BSS.Protocol.Equilibration(
    timestep=timestep,
    runtime=runtime,
    temperature_start=temperature_start,
    temperature_end=temperature_end,
    temperature=temperature,
    thermostat_time_constant=thermostat_time_constant,
    pressure=pressure,
    restraint=restraint,
    force_constant=force_constant,
)
system = runProcess(
    system=system,
    protocol=protocol,
    engine=args.engine,
)
# Save the system to the output directory
output_file = output_directory / f"{args.outfile_name}"
BSS.Stream.save(system, str(output_file))
