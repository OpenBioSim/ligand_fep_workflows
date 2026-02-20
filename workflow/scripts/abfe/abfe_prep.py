#!/usr/bin/env python
"""
ABFE system preparation script.

This script prepares equilibrated systems for alchemical free energy
calculations by marking the ligand for decoupling. Unlike RBFE which
creates merged molecules between two ligands, ABFE marks a single ligand
for complete decoupling from the environment.

The preparation includes:
    1. Loading the equilibrated system
    2. Identifying the ligand molecule
    3. Marking it for alchemical decoupling using BSS.Align.decouple()
    4. Applying hydrogen mass repartitioning (HMR) for longer timesteps
    5. Saving the prepared system

For bound leg systems, restraint parameters are also loaded for later
application during the restrain stage of production.

Usage:
    python abfe_prep.py --input system.bss --output-directory prepared/
                        --ligand-name ejm42 --leg bound
                        --restraint-file restraints/ejm42_restraint.json

Author: ABFE Workflow
"""

import argparse
import json
import sys
from pathlib import Path

import BioSimSpace.Sandpit.Exscientia as BSS


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Prepare systems for ABFE alchemical decoupling.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the equilibrated system (.bss file).",
    )
    parser.add_argument(
        "--output-directory",
        type=str,
        required=True,
        help="Directory to save the prepared system.",
    )
    parser.add_argument(
        "--ligand-name",
        type=str,
        required=True,
        help="Name of the ligand for output file naming.",
    )
    parser.add_argument(
        "--leg",
        type=str,
        choices=["bound", "free"],
        required=True,
        help="Which leg of the calculation (bound or free).",
    )
    parser.add_argument(
        "--restraint-file",
        type=str,
        default=None,
        help="Path to restraint parameters JSON (required for bound leg).",
    )
    parser.add_argument(
        "--hmr-factor",
        type=float,
        default=3.0,
        help="Hydrogen mass repartitioning factor (1.0 = no HMR).",
    )
    return parser.parse_args()


def find_ligand(system: BSS._SireWrappers.System) -> BSS._SireWrappers.Molecule:
    """
    Find the ligand molecule in the system.

    The ligand is identified as a molecule with exactly 1 residue and
    more than 5 atoms (to exclude single ions).

    Args:
        system: BioSimSpace system

    Returns:
        The ligand molecule

    Raises:
        ValueError: If no ligand can be identified
    """
    for mol in system:
        n_residues = mol.nResidues()
        n_atoms = mol.nAtoms()
        # Ligand has 1 residue but many atoms
        if n_residues == 1 and n_atoms > 5:
            print(f"Found ligand: {n_atoms} atoms, {n_residues} residue(s)")
            return mol

    raise ValueError(
        "Could not identify ligand in system. "
        "Expected molecule with 1 residue and >5 atoms."
    )


def decouple_ligand(
    system: BSS._SireWrappers.System,
    ligand: BSS._SireWrappers.Molecule,
) -> BSS._SireWrappers.System:
    """
    Mark the ligand for alchemical decoupling.

    This uses BSS.Align.decouple() to mark the ligand for complete
    decoupling from its environment. The decoupled ligand has its
    intramolecular interactions preserved while intermolecular
    interactions can be turned off via lambda scaling.

    Args:
        system: Complete system
        ligand: Ligand molecule to decouple

    Returns:
        System with decoupled ligand
    """
    print("Marking ligand for alchemical decoupling...")

    # Decouple the ligand
    decoupled_ligand = BSS.Align.decouple(ligand)

    # Update the system with the decoupled ligand
    system.updateMolecules(decoupled_ligand)

    return system


def apply_hmr(
    system: BSS._SireWrappers.System,
    factor: float,
) -> BSS._SireWrappers.System:
    """
    Apply hydrogen mass repartitioning to the system.

    HMR transfers mass from heavy atoms to bonded hydrogens, allowing
    larger integration timesteps (up to 4fs with factor=3).

    Args:
        system: System to modify
        factor: HMR factor (mass multiplier for hydrogens)

    Returns:
        System with HMR applied
    """
    if factor > 1.0:
        print(f"Applying hydrogen mass repartitioning (factor={factor})...")
        system.repartitionHydrogenMass(factor=factor, water="no")
    return system


def load_restraint_parameters(restraint_file: str) -> dict:
    """
    Load restraint parameters from JSON file.

    Args:
        restraint_file: Path to restraint JSON file

    Returns:
        Dictionary of restraint parameters
    """
    with open(restraint_file, "r") as f:
        params = json.load(f)
    print(f"Loaded restraint parameters from {restraint_file}")
    return params


def main():
    """Main entry point for ABFE preparation."""
    args = parse_args()

    # Validate arguments
    if args.leg == "bound" and args.restraint_file is None:
        print("Warning: No restraint file provided for bound leg.")
        print("Restraints will need to be set up during production.")

    # Create output directory
    output_dir = Path(args.output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the equilibrated system
    print(f"Loading system from {args.input}...")
    system = BSS.Stream.load(args.input)
    print(f"Loaded system with {system.nMolecules()} molecules")

    # Find the ligand
    ligand = find_ligand(system)

    # Mark ligand for decoupling
    system = decouple_ligand(system, ligand)

    # Apply HMR
    system = apply_hmr(system, args.hmr_factor)

    # For bound leg, load and store restraint parameters
    if args.leg == "bound" and args.restraint_file:
        restraint_params = load_restraint_parameters(args.restraint_file)
        # The restraint parameters are stored in the JSON file and will be
        # read during production to set up the Boresch restraints

    # Save the prepared system
    output_file = output_dir / f"{args.ligand_name}_{args.leg}"
    BSS.Stream.save(system, str(output_file))
    print(f"Saved prepared system to {output_file}.bss")

    print(f"\nABFE preparation complete for {args.ligand_name} ({args.leg} leg)")


if __name__ == "__main__":
    main()
