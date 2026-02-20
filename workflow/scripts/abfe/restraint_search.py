#!/usr/bin/env python
"""
Restraint search script for ABFE workflow.

This script performs a short unrestrained MD simulation of the protein-ligand
complex and analyses the trajectory to find optimal Boresch restraints. The
restraints are used in the bound leg to maintain the ligand orientation as it
is decoupled from the protein.

Boresch restraints define 6 external degrees of freedom using:
    - 1 distance (r): between anchor atoms in protein and ligand
    - 2 angles (theta_A, theta_B): at the protein and ligand anchors
    - 3 dihedrals (phi_A, phi_B, phi_C): along the anchor chain

The analytical correction for releasing these restraints is computed and
saved for use in the final free energy calculation.

Usage:
    python restraint_search.py --input system.bss --output-directory restraints/
                               --ligand-name ejm42 --runtime 1ns

Author: ABFE Workflow
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import BioSimSpace.Sandpit.Exscientia as BSS


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

    Args:
        system: Complete system
        ligand: Ligand molecule to decouple

    Returns:
        System with decoupled ligand
    """
    print("Marking ligand for alchemical decoupling...")
    decoupled_ligand = BSS.Align.decouple(ligand)
    system.updateMolecules(decoupled_ligand)
    return system


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Find optimal Boresch restraints for ABFE calculations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the equilibrated protein-ligand system (.bss file).",
    )
    parser.add_argument(
        "--output-directory",
        type=str,
        required=True,
        help="Directory to save restraint parameters and correction.",
    )
    parser.add_argument(
        "--ligand-name",
        type=str,
        required=True,
        help="Name of the ligand for output file naming.",
    )
    parser.add_argument(
        "--runtime",
        type=str,
        default="1ns",
        help="Duration of the restraint search simulation.",
    )
    parser.add_argument(
        "--temperature",
        type=str,
        default="300K",
        help="Temperature for the simulation.",
    )
    return parser.parse_args()


def run_restraint_search_simulation(
    system: BSS._SireWrappers.System,
    runtime: BSS.Types.Time,
    temperature: BSS.Types.Temperature,
    work_dir: Path,
) -> tuple[BSS._SireWrappers.System, BSS.Trajectory.Trajectory]:
    """
    Run a short unrestrained MD simulation for restraint search.

    Args:
        system: Equilibrated protein-ligand system
        runtime: Duration of simulation
        temperature: Simulation temperature
        work_dir: Working directory for simulation files

    Returns:
        Tuple of (final system, trajectory)
    """
    print(f"Running restraint search simulation for {runtime}...")

    protocol = BSS.Protocol.Equilibration(
        runtime=runtime,
        temperature=temperature,
        pressure=BSS.Types.Pressure(1, "atm"),
        restraint=None,
    )

    sim_dir = str(work_dir / "restraint_search_sim")
    process = BSS.Process.Gromacs(
        system, protocol, work_dir=sim_dir, ignore_warnings=True
    )
    process.start()
    process.wait()

    if process.isError():
        raise RuntimeError(
            f"Restraint search simulation failed: {process.stdout()}\n{process.stderr()}"
        )

    return process.getSystem(), process.getTrajectory()


def analyse_restraints(
    system: BSS._SireWrappers.System,
    trajectory: BSS.Trajectory.Trajectory,
    temperature: BSS.Types.Temperature,
    work_dir: Path,
) -> BSS.FreeEnergy.Restraint:
    """
    Analyse trajectory to find optimal Boresch restraints.

    Args:
        system: Protein-ligand system
        trajectory: MD trajectory for analysis
        temperature: Temperature for correction calculation
        work_dir: Working directory for analysis files

    Returns:
        BioSimSpace Restraint object containing anchor points and force constants
    """
    print("Analysing trajectory for optimal restraint parameters...")

    return BSS.FreeEnergy.RestraintSearch.analyse(
        work_dir=str(work_dir / "restraint_search_sim"),
        system=system,
        traj=trajectory,
        temperature=temperature,
        method="BSS",
        restraint_type="Boresch",
    )


def save_restraint_parameters(
    restraint: BSS.FreeEnergy.Restraint,
    output_dir: Path,
    ligand_name: str,
) -> dict[str, Any]:
    """
    Save restraint parameters to JSON file.

    Args:
        restraint: BioSimSpace Restraint object
        output_dir: Directory for output files
        ligand_name: Name of ligand for file naming

    Returns:
        Dictionary of restraint parameters
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract restraint dictionary (private attribute)
    restraint_dict = restraint._restraint_dict

    # Helper to get atom index from BSS Atom object
    def get_atom_index(atom):
        if atom is None:
            return None
        return atom.index()

    # Convert to serialisable format
    params = {
        "anchor_points": {
            "r1": get_atom_index(restraint_dict.get("anchor_points", {}).get("r1")),
            "r2": get_atom_index(restraint_dict.get("anchor_points", {}).get("r2")),
            "r3": get_atom_index(restraint_dict.get("anchor_points", {}).get("r3")),
            "l1": get_atom_index(restraint_dict.get("anchor_points", {}).get("l1")),
            "l2": get_atom_index(restraint_dict.get("anchor_points", {}).get("l2")),
            "l3": get_atom_index(restraint_dict.get("anchor_points", {}).get("l3")),
        },
        "equilibrium_values": {
            "r0": str(restraint_dict.get("equilibrium_values", {}).get("r0")),
            "thetaA0": str(restraint_dict.get("equilibrium_values", {}).get("thetaA0")),
            "thetaB0": str(restraint_dict.get("equilibrium_values", {}).get("thetaB0")),
            "phiA0": str(restraint_dict.get("equilibrium_values", {}).get("phiA0")),
            "phiB0": str(restraint_dict.get("equilibrium_values", {}).get("phiB0")),
            "phiC0": str(restraint_dict.get("equilibrium_values", {}).get("phiC0")),
        },
        "force_constants": {
            "kr": str(restraint_dict.get("force_constants", {}).get("kr")),
            "kthetaA": str(restraint_dict.get("force_constants", {}).get("kthetaA")),
            "kthetaB": str(restraint_dict.get("force_constants", {}).get("kthetaB")),
            "kphiA": str(restraint_dict.get("force_constants", {}).get("kphiA")),
            "kphiB": str(restraint_dict.get("force_constants", {}).get("kphiB")),
            "kphiC": str(restraint_dict.get("force_constants", {}).get("kphiC")),
        },
    }

    # Save to JSON
    output_file = output_dir / f"{ligand_name}_restraint.json"
    with open(output_file, "w") as f:
        json.dump(params, f, indent=2, default=str)

    print(f"Saved restraint parameters to {output_file}")
    return params


def save_correction(
    restraint: BSS.FreeEnergy.Restraint,
    output_dir: Path,
    ligand_name: str,
) -> float:
    """
    Calculate and save the analytical correction for restraint release.

    The correction term accounts for the free energy of releasing the
    Boresch restraints in the reference state (non-interacting ligand).

    Args:
        restraint: BioSimSpace Restraint object
        output_dir: Directory for output files
        ligand_name: Name of ligand

    Returns:
        Correction value in kcal/mol
    """
    correction = restraint.getCorrection()
    correction_value = correction.value()  # in kcal/mol

    output_file = output_dir / f"{ligand_name}_correction.txt"
    with open(output_file, "w") as f:
        f.write(f"{correction_value}\n")

    print(f"Restraint correction: {correction_value:.4f} kcal/mol")
    print(f"Saved correction to {output_file}")

    return correction_value


def main():
    """Main entry point for restraint search."""
    args = parse_args()

    # Parse units
    try:
        runtime = BSS.Types.Time(args.runtime)
    except ValueError:
        print(f"Error: Invalid runtime '{args.runtime}'")
        sys.exit(1)

    try:
        temperature = BSS.Types.Temperature(args.temperature)
    except ValueError:
        print(f"Error: Invalid temperature '{args.temperature}'")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the equilibrated system
    print(f"Loading system from {args.input}...")
    system = BSS.Stream.load(args.input)

    # Run restraint search simulation
    work_dir = output_dir / args.ligand_name
    system_final, trajectory = run_restraint_search_simulation(
        system=system,
        runtime=runtime,
        temperature=temperature,
        work_dir=work_dir,
    )

    # Decouple the ligand (required for RestraintSearch.analyse)
    ligand = find_ligand(system_final)
    system_final = decouple_ligand(system_final, ligand)

    # Analyse trajectory for optimal restraints
    restraint = analyse_restraints(
        system=system_final,
        trajectory=trajectory,
        temperature=temperature,
        work_dir=work_dir,
    )

    # Save restraint parameters
    save_restraint_parameters(
        restraint=restraint,
        output_dir=output_dir,
        ligand_name=args.ligand_name,
    )

    # Calculate and save correction
    save_correction(
        restraint=restraint,
        output_dir=output_dir,
        ligand_name=args.ligand_name,
    )

    print(f"\nRestraint search complete for {args.ligand_name}")


if __name__ == "__main__":
    main()
