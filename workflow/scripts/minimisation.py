import argparse
from pathlib import Path

import BioSimSpace as BSS
from _common import runProcess
from BioSimSpace import _Exceptions

parser = argparse.ArgumentParser(description="Minimisation script")
parser.add_argument(
    "--input", type=str, help="Full path to the input file", required=True
)
parser.add_argument(
    "--outfile-name",
    type=str,
    help="Name of the output file (without extension).",
    required=True,
)
parser.add_argument(
    "--output-directory",
    type=str,
    help="Full path to the output directory.",
    required=True,
)
parser.add_argument(
    "--engine",
    type=str,
    help="The engine to use for minimisation",
    choices=["openmm", "gromacs", "amber"],
    default="gromacs",
)
parser.add_argument(
    "--minimisation-steps",
    type=int,
    help="Number of minimisation steps to perform",
    default=10000,
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

if args.restraint_indices:
    restraint = args.restraint_indices
    force_constant = args.restraint_force_constant
elif args.restraint_string:
    restraint = args.restraint_string
    force_constant = args.restraint_force_constant
else:
    restraint = None
    force_constant = 0.0  # Can't be `None`, needs to be float or unit or string

try:
    system = BSS.Stream.load(args.input)
except Exception as e:
    raise RuntimeError(f"Failed to load input file: {e}")

# make the outdir if it doesn't exist
output_dir = Path(args.output_directory)
output_dir.mkdir(parents=True, exist_ok=True)


protocol = BSS.Protocol.Minimisation(
    steps=args.minimisation_steps,
    restraint=restraint,
    force_constant=force_constant,
)
system = runProcess(
    system=system,
    protocol=protocol,
    engine=args.engine,
)
# Save the system to the output directory
output_file = output_dir / args.outfile_name
BSS.Stream.save(system, str(output_file))
