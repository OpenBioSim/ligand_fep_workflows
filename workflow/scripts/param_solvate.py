import argparse
import sys
from pathlib import Path

import BioSimSpace as BSS

parser = argparse.ArgumentParser(
    description="Parameterise and solvate both free and bound legs."
)
parser.add_argument(
    "--ligand-file",
    type=str,
    help="Path to the ligand file (e.g. mol2, sdf, pdb)",
    required=True,
)
parser.add_argument(
    "--protein-files",
    type=str,
    nargs="+",
    help="Path to parameterised protein files (rst7 and prm7)",
    required=True,
)
parser.add_argument(
    "--ligand-forcefield",
    type=str,
    help="The name of the force field to use for parameterisation.",
    choices=BSS.Parameters.forceFields(),
    default="gaff2",
)
parser.add_argument(
    "--water-model",
    type=str,
    help="The name of the water model to use for solvation.",
    choices=BSS.Solvent.waterModels(),
    default="tip3p",
)
parser.add_argument(
    "--box-length",
    type=float,
    help="Length of the box edges. Applied in x, y and z directions. In nanometres.",
    default=5.0,
)
parser.add_argument(
    "--box-type",
    type=str,
    help="Box type to use for both bound and free legs",
    choices=BSS.Box.boxTypes(),
    default="cubic",
)
parser.add_argument(
    "--ion-conc",
    type=float,
    help="Concentration of ions to add to the system (in mol/L)",
    default=0.15,
)
parser.add_argument(
    "--output-directory",
    type=str,
    help="Directory to save the output files",
    default="output",
)
args = parser.parse_args()

# first grab the name of the ligand
ligand_name = Path(args.ligand_file).stem
outdir = Path(args.output_directory)

# Check if the outputs we need already exist
if (outdir / (ligand_name + "_free.bss")).exists() and (
    outdir / (ligand_name + "_bound.bss")
).exists():
    print(
        f"Output files for {ligand_name} already exist in {outdir}. Skipping parameterisation and solvation."
    )
    sys.exit(0)

# Create the output directory if it doesn't exist
outdir.mkdir(parents=True, exist_ok=True)

# load the unparameterised ligand
lig_o = BSS.IO.readMolecules(args.ligand_file)
# parameterise the ligand
lig_p = BSS.Parameters.parameterise(
    lig_o[0],
    args.ligand_forcefield,
).getMolecule()

# read the protein files
protein = BSS.IO.readMolecules(args.protein_files)
# combine to form the bound system
system_bound = lig_p + protein

box_len = args.box_length * BSS.Types.Length("1nm")

# find the bounding box for the ligand only
box_min, box_max = lig_p.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + box_len for x in box_size]

boxtype_dict = {
    "cubic": BSS.Box.cubic,
    "rhombicDodecahedronHexagon": BSS.Box.rhombicDodecahedronHexagon,
    "rhombicDodecahedronSquare": BSS.Box.rhombicDodecahedronSquare,
    "truncatedOctahedron": BSS.Box.truncatedOctahedron,
}
boxtype_func = boxtype_dict[args.box_type]
box, angles = boxtype_func(max(box_sizes))
# print(f"box of free leg: {box}, angles: {angles}")
lig_p_solvated = BSS.Solvent.solvate(
    args.water_model,
    molecule=lig_p,
    box=box,
    angles=angles,
    ion_conc=args.ion_conc,
)
BSS.Stream.save(lig_p_solvated, str(outdir / f"{ligand_name}_free"))

# again for protein-ligand system
box_min, box_max = system_bound.getAxisAlignedBoundingBox()
box_size = [y - x for x, y in zip(box_min, box_max)]
box_sizes = [x + box_len for x in box_size]

box, angles = boxtype_func(max(box_sizes))
system_solvated = BSS.Solvent.solvate(
    args.water_model,
    molecule=system_bound,
    box=box,
    angles=angles,
    ion_conc=args.ion_conc,
)
BSS.Stream.save(system_solvated, str(outdir / f"{ligand_name}_bound"))
