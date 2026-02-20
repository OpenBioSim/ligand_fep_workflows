# Creates merged molecules for a specified list of perturbations, both free and bound legs
from pathlib import Path

import BioSimSpace as BSS
import configargparse
from BioSimSpace import _Exceptions


# TODO: Add a check of production settings and warn if 4fs timestep is used without HMR.
def prep_fep_free(
    ligand1_name: str,
    ligand2_name: str,
    ligand1_files: list | Path | str,
    ligand2_files: list | Path | str,
    output_location: str | Path,
    md_engine: str = "gromacs",
    hmr_factor: float = 1.5,
) -> None:
    md_engine = md_engine.strip().lower()

    if md_engine not in ["gromacs", "somd", "somd2", "amber"]:
        raise ValueError("md_engine must be one of: gromacs, somd, somd2, amber")

    # first lets make the output directory
    outpath = Path(output_location)
    outpath_full = outpath / "free"
    outpath_full.mkdir(parents=True, exist_ok=True)
    # print(f"Working directory: {outpath_full}")

    # load ligand files to BSS, if they're lists assume that they are some kind of prm and rst files
    if isinstance(ligand1_files, list) and isinstance(ligand2_files, list):
        ligand_1_sys = BSS.IO.readMolecules(ligand1_files)
        ligand_2_sys = BSS.IO.readMolecules(ligand2_files)
    # if they're not lists, assume that they are stream files
    else:
        ligand_1_sys = BSS.Stream.load(ligand1_files)
        ligand_2_sys = BSS.Stream.load(ligand2_files)

    if hmr_factor > 1.0:
        ligand_1_sys.repartitionHydrogenMass(factor=hmr_factor, water="no")
        ligand_2_sys.repartitionHydrogenMass(factor=hmr_factor, water="no")

    # Guessing that the first molecule in the system is the ligand.
    ligand_1 = ligand_1_sys[0]
    ligand_2 = ligand_2_sys[0]

    # Generate the mapping, AMBER needs special treatment.
    if md_engine != "amber":
        mapping = BSS.Align.matchAtoms(
            ligand_1,
            ligand_2,
            complete_rings_only=True,
        )
    else:
        mapping = BSS.Align.matchAtoms(
            ligand_1,
            ligand_2,
            complete_rings_only=True,
            prune_perturbed_constraints=True,
            prune_crossing_constraints=True,
        )
    inv_mapping = {v: k for k, v in mapping.items()}
    ligand_2_a = BSS.Align.rmsdAlign(ligand_2, ligand_1, inv_mapping)

    # Generate merged molecule.
    # print("Merging..")
    merged_ligs = BSS.Align.merge(ligand_1, ligand_2_a, mapping)

    ligand_1_sys.removeMolecules(ligand_1)
    ligand_1_sys.addMolecules(merged_ligs)
    system_free = ligand_1_sys
    # print("Saving merged system..")
    BSS.Stream.save(
        system_free,
        str(outpath_full / f"{ligand1_name}~{ligand2_name}"),
    )
    return None


def prep_fep_bound(
    ligand1_name: str,
    ligand2_name: str,
    ligand1_files: list,
    ligand2_files: list,
    output_location: str | Path,
    md_engine: str = "gromacs",
    hmr_factor: float = 1.5,
) -> None:
    md_engine = md_engine.strip().lower()

    if md_engine not in ["gromacs", "somd", "somd2", "amber"]:
        raise ValueError("md_engine must be one of: gromacs, somd, somd2, amber")

    # first lets make the output directory
    outpath = Path(output_location)
    outpath_full = outpath / "bound"
    outpath_full.mkdir(parents=True, exist_ok=True)
    # print(f"Working directory: {outpath_full}")

    # load ligand files to BSS, if they're lists assume that they are some kind of prm and rst files
    if isinstance(ligand1_files, list) and isinstance(ligand2_files, list):
        system_1 = BSS.IO.readMolecules(ligand1_files)
        system_2 = BSS.IO.readMolecules(ligand2_files)
    # if they're not lists, assume that they are stream files
    else:
        system_1 = BSS.Stream.load(ligand1_files)
        system_2 = BSS.Stream.load(ligand2_files)

    if hmr_factor > 1.0:
        system_1.repartitionHydrogenMass(factor=hmr_factor, water="no")
        system_2.repartitionHydrogenMass(factor=hmr_factor, water="no")

    # Extract ligands and protein. Do this based on nAtoms and nResidues, as sometimes
    # the order of molecules is switched, so we can't use index alone.
    # bugfix in BSS makes the below redundant but keeping this in to be 100% sure we're getting the correct structures.
    system_ligand_1 = None
    protein = None
    n_residues = [mol.nResidues() for mol in system_1]
    n_atoms = [mol.nAtoms() for mol in system_1]
    for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):
        if n_resi == 1 and n_at > 5:
            system_ligand_1 = system_1.getMolecule(i)
            print(f"Found ligand 1: {system_ligand_1}")
        elif n_resi > 1:
            protein = system_1.getMolecule(i)
            print(f"Found protein: {protein}")
        else:
            pass

    # loop over molecules in system to extract the ligand
    system_ligand_2 = None

    n_residues = [mol.nResidues() for mol in system_2]
    n_atoms = [mol.nAtoms() for mol in system_2]
    for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):
        # grab the system's ligand and the protein. ignore the waters.
        if n_resi == 1 and n_at > 5:
            system_ligand_2 = system_2.getMolecule(i)
            print(f"Found ligand 2: {system_ligand_2}")
        else:
            pass

    if system_ligand_1 and system_ligand_2 and protein:
        # print("Using molecules ligand_1, ligand_2, protein:")
        # print(system_ligand_1, system_ligand_2, protein)
        pass
    else:
        unable_to_find = []
        if not system_ligand_1:
            unable_to_find.append("ligand 1")
        if not system_ligand_2:
            unable_to_find.append("ligand 2")
        if not protein:
            unable_to_find.append("protein")
        raise _Exceptions.AlignmentError(
            f"Could not extract ligands or protein from input systems. Unable to find {' and '.join(unable_to_find)}."
        )

    # Generate the mapping, AMBER needs special treatment.
    if md_engine != "amber":
        mapping = BSS.Align.matchAtoms(
            system_ligand_1,
            system_ligand_2,
            complete_rings_only=True,
        )
    else:
        mapping = BSS.Align.matchAtoms(
            system_ligand_1,
            system_ligand_2,
            complete_rings_only=True,
            prune_perturbed_constraints=True,
            prune_crossing_constraints=True,
        )
    inv_mapping = {v: k for k, v in mapping.items()}

    # print("Aligning..")
    system_ligand_2_a = BSS.Align.rmsdAlign(
        system_ligand_2, system_ligand_1, inv_mapping
    )

    # Generate merged molecule.
    # print("Merging..")
    system_merged_ligs = BSS.Align.merge(system_ligand_1, system_ligand_2_a, mapping)

    system_1.removeMolecules(system_ligand_1)
    system_1.addMolecules(system_merged_ligs)
    system_bound = system_1
    # print("Saving merged system..")
    BSS.Stream.save(
        system_bound,
        str(outpath_full / f"{ligand1_name}~{ligand2_name}"),
    )
    return None


if __name__ == "__main__":
    parser = configargparse.ArgParser(
        description="Creates merged molecules for a specified list of perturbations, both free and bound legs."
    )
    parser.add_argument(
        "--ligand1_name",
        type=str,
        required=True,
        help="Name of the first ligand.",
    )
    parser.add_argument(
        "--ligand2_name",
        type=str,
        required=True,
        help="Name of the second ligand.",
    )
    parser.add_argument(
        "--ligand1_files",
        type=str,
        nargs="+",
        required=True,
        help="File(s) for the first ligand.",
    )
    parser.add_argument(
        "--ligand2_files",
        type=str,
        nargs="+",
        required=True,
        help="File(s) for the second ligand.",
    )
    parser.add_argument(
        "--output_location",
        type=str,
        required=True,
        help="Output location for merged files.",
    )
    parser.add_argument(
        "--leg",
        type=str,
        default="free",
        choices=["free", "bound"],
        help="Specify whether to prepare free or bound leg.",
    )
    parser.add_argument(
        "--md_engine",
        type=str,
        default="gromacs",
        choices=["gromacs", "somd", "somd2", "amber"],
        help="Molecular dynamics engine to use.",
    )
    parser.add_argument(
        "--hmr-factor",
        type=float,
        default=1.5,
        help="Hydrogen mass repartitioning factor.",
    )

    args = parser.parse_args()

    if args.leg == "free":
        prep_fep_free(
            ligand1_name=args.ligand1_name,
            ligand2_name=args.ligand2_name,
            ligand1_files=args.ligand1_files,
            ligand2_files=args.ligand2_files,
            output_location=args.output_location,
            md_engine=args.md_engine,
            hmr_factor=args.hmr_factor,
        )
    elif args.leg == "bound":
        prep_fep_bound(
            ligand1_name=args.ligand1_name,
            ligand2_name=args.ligand2_name,
            ligand1_files=args.ligand1_files,
            ligand2_files=args.ligand2_files,
            output_location=args.output_location,
            md_engine=args.md_engine,
            hmr_factor=args.hmr_factor,
        )
    else:
        raise ValueError("Invalid leg type specified. Use 'free' or 'bound'.")
