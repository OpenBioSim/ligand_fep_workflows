"""Parameterise a ligand from an SDF file, keeping its atom order and coordinates.
The force field assigns the partial charges."""

import argparse
import tempfile
from pathlib import Path

import BioSimSpace as BSS

LIGAND_RESNAME = "LIG"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligand-file", required=True)
    parser.add_argument("--output", required=True, help="Output stream file, e.g. ligand.bss")
    parser.add_argument("--forcefield", default="openff_unconstrained-2.3.0")
    args = parser.parse_args()

    ligand = BSS.IO.readMolecules(args.ligand_file)[0]
    with tempfile.TemporaryDirectory() as tmp:
        mol = BSS.Parameters.parameterise(ligand, args.forcefield, work_dir=tmp).getMolecule()

    elements = [a.element().num_protons() for a in mol._sire_object.atoms()]
    if elements != [a.element().num_protons() for a in ligand._sire_object.atoms()]:
        raise ValueError("Parameterisation changed the atom order")

    cursor = mol._sire_object.cursor()
    for residue in cursor.residues():
        residue.name = LIGAND_RESNAME
    mol._sire_object = cursor.molecule().commit()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    BSS.Stream.save(mol.toSystem(), str(output.with_suffix("")))
    print(f"Saved {output}: net charge {round(mol.charge().value())}")


if __name__ == "__main__":
    main()
