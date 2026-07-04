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
import sire as sr


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
    parser.add_argument(
        "--method",
        type=str,
        choices=["bss", "sire"],
        default="bss",
        help="Restraint search backend: 'bss' (BioSimSpace/Aldeghi-style search, "
        "runs via GROMACS, shared by both engines) or 'sire' (sire-native "
        "boresch_search(), mirroring SOMD2's own built-in restraint "
        "auto-generation; SOMD2-only).",
    )
    parser.add_argument(
        "--protocol",
        type=str,
        choices=["rxrx", "aldeghi"],
        default="rxrx",
        help="sire.restraints.boresch_search() protocol (only used with --method sire).",
    )
    parser.add_argument(
        "--timestep",
        type=str,
        default="4fs",
        help="Integration timestep for the search trajectory (only used with --method sire).",
    )
    parser.add_argument(
        "--cutoff-type",
        type=str,
        default="PME",
        help="Electrostatics cutoff type for the search trajectory (only used with --method sire).",
    )
    parser.add_argument(
        "--cutoff",
        type=str,
        default="10A",
        help="Cutoff distance for the search trajectory (only used with --method sire).",
    )
    parser.add_argument(
        "--perturbable-constraint",
        type=str,
        default="h_bonds_not_heavy_perturbed",
        help="Constraint type for the perturbable ligand during the search trajectory "
        "(only used with --method sire).",
    )
    parser.add_argument(
        "--frame-frequency",
        type=str,
        default="10ps",
        help="Trajectory frame sampling interval for the search (only used with --method sire).",
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

    process = BSS.Process.Gromacs(
        system, protocol, work_dir=str(work_dir), ignore_warnings=True
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
        work_dir=str(work_dir),
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

    return _write_correction_file(correction_value, output_dir, ligand_name)


def _write_correction_file(
    correction_value: float,
    output_dir: Path,
    ligand_name: str,
) -> float:
    """
    Write a correction value (in kcal/mol) to the shared {ligand}_correction.txt
    file consumed by the engine-agnostic analysis pipeline. Shared by both the
    BSS and sire-native restraint search methods.
    """
    output_file = output_dir / f"{ligand_name}_correction.txt"
    with open(output_file, "w") as f:
        f.write(f"{correction_value}\n")

    print(f"Restraint correction: {correction_value:.4f} kcal/mol")
    print(f"Saved correction to {output_file}")

    return correction_value


def find_ligand_sire(system: "sr.system.System"):
    """
    Find the ligand molecule in a sire system.

    Uses the same heuristic as find_ligand()/production_somd2.py: the ligand
    is a molecule with exactly 1 residue and more than 5 atoms.

    Args:
        system: Sire system

    Returns:
        The ligand molecule (a sire Molecule view)

    Raises:
        ValueError: If no ligand can be identified
    """
    for i in range(system.num_molecules()):
        mol = system[i]
        if mol.num_residues() == 1 and mol.num_atoms() > 5:
            print(f"Found ligand: molecule {i}, {mol.num_atoms()} atoms")
            return mol

    raise ValueError(
        "Could not identify ligand in system. "
        "Expected molecule with 1 residue and >5 atoms."
    )


def run_native_restraint_search(
    system: "sr.system.System",
    runtime: str,
    frequency: str,
    temperature: str,
    protocol: str,
    timestep: str,
    cutoff_type: str,
    cutoff: str,
    perturbable_constraint: str,
):
    """
    Run a short unrestrained lambda=0 trajectory and use sire's native
    boresch_search() to find optimal Boresch restraints, mirroring SOMD2's
    own built-in restraint auto-generation (_generate_boresch_restraint in
    somd2.runner._base), but performed here (once per ligand, ahead of
    production) rather than per-replica, so the resulting restraint is
    shared across all bound-leg replicas and the standard-state correction
    can be written to the shared {ligand}_correction.txt file.

    Args:
        system: Sire system with the ligand already sire-natively decoupled
        runtime: Duration of the search trajectory
        frequency: Frame-saving interval for the search trajectory
        temperature: Simulation/correction temperature
        protocol: sire.restraints.boresch_search() protocol ("rxrx" or "aldeghi")
        timestep: Integration timestep
        cutoff_type: Electrostatics cutoff type
        cutoff: Cutoff distance
        perturbable_constraint: Constraint type for the perturbable ligand

    Returns:
        Tuple of (sire.mm.BoreschRestraints, correction value in kcal/mol,
        starting_structure). The starting structure is the least-strained
        trajectory frame returned by boresch_search, linked to the reference
        (lambda=0) end state and with the search trajectory frames dropped, to
        be used to seed production instead of the pre-search input.
    """
    from sire.restraints import boresch_search

    dynamics_kwargs = {
        "timestep": timestep,
        "temperature": temperature,
        "cutoff_type": cutoff_type,
        "cutoff": cutoff,
        "constraint": "h_bonds",
        "perturbable_constraint": perturbable_constraint,
        "platform": "CUDA",
        "lambda_value": 0.0,
    }

    print("Minimising before restraint search trajectory...")
    dynamics = system.dynamics(**dynamics_kwargs)
    dynamics.minimise()

    print(f"Running restraint search trajectory for {runtime}...")
    dynamics.run(
        runtime,
        energy_frequency=0,
        frame_frequency=frequency,
        save_velocities=False,
    )
    search_system = dynamics.commit()

    print(f"Analysing trajectory with boresch_search(protocol={protocol!r})...")
    restraints, correction, starting_structure = boresch_search(
        search_system, protocol=protocol, temperature=temperature
    )

    correction_value = float(correction.to(sr.units.kcal_per_mol))
    print(f"Standard state correction: {correction_value:.4f} kcal mol-1")

    # boresch_search fits the restraint equilibrium values to trajectory
    # averages, so the search input structure is generally not consistent with
    # them; seed production from the least-strained frame instead to avoid a
    # large restraint force at t=0. Link the frame (from the perturbable search
    # system) back to the reference (lambda=0) end state, matching how SOMD2
    # itself handles the boresch_search starting structure, and keep only that
    # single snapshot.
    starting_structure = sr.morph.link_to_reference(starting_structure)
    starting_structure.delete_all_frames()

    return restraints, correction_value, starting_structure


def save_native_restraint(
    restraints: "sr.mm.BoreschRestraints",
    output_dir: Path,
    ligand_name: str,
) -> Path:
    """
    Save a sire-native Boresch restraint via sire's own stream serialisation
    (no JSON massaging required, unlike the BSS-derived restraint).

    Args:
        restraints: sire.mm.BoreschRestraints object
        output_dir: Directory for output files
        ligand_name: Name of ligand for file naming

    Returns:
        Path to the saved restraint file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{ligand_name}_restraint.s3"
    sr.stream.save(restraints, str(output_file))
    print(f"Saved restraint to {output_file}")
    return output_file


def save_native_starting_structure(
    starting_structure: "sr.system.System",
    output_dir: Path,
    ligand_name: str,
) -> Path:
    """
    Save the least-strained starting structure returned by boresch_search
    (see run_native_restraint_search), so production can seed from a structure
    consistent with the restraint rather than the pre-search input. Saved as a
    sibling of the restraint file (``{ligand}_restraint_structure.s3``), which
    production_somd2.py discovers automatically.

    Args:
        starting_structure: sire System (single frame, reference-linked)
        output_dir: Directory for output files
        ligand_name: Name of ligand for file naming

    Returns:
        Path to the saved starting structure file
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{ligand_name}_restraint_structure.s3"
    sr.stream.save(starting_structure, str(output_file))
    print(f"Saved restraint starting structure to {output_file}")
    return output_file


def main():
    """Main entry point for restraint search."""
    args = parse_args()

    # Create output directory
    output_dir = Path(args.output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.method == "sire":
        _main_sire(args, output_dir)
    else:
        _main_bss(args, output_dir)

    print(f"\nRestraint search complete for {args.ligand_name}")


def _main_bss(args: argparse.Namespace, output_dir: Path):
    """Legacy restraint search: BioSimSpace/Aldeghi-style search via GROMACS."""
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


def _main_sire(args: argparse.Namespace, output_dir: Path):
    """
    Native restraint search: sire.restraints.boresch_search(), mirroring
    SOMD2's own built-in restraint auto-generation. SOMD2-only.
    """
    # Load the equilibrated system and convert to sire.
    print(f"Loading system from {args.input}...")
    bss_system = BSS.Stream.load(args.input)
    system = sr.system.System(bss_system._sire_object)

    # Find the ligand and apply sire-native decoupling, matching
    # production_somd2.py so the search trajectory sees the same
    # perturbable topology used for production.
    ligand = find_ligand_sire(system)
    ligand = sr.morph.decouple(ligand, as_new_molecule=False)
    system.update(ligand)

    # Normalise e.g. "300K" -> "300 K" so sire's unit parser accepts it,
    # matching production_somd2.py's handling of the same CLI convention.
    temp_value = float("".join(c for c in args.temperature if c.isdigit() or c == "."))
    temperature = f"{temp_value} K"

    restraints, correction_value, starting_structure = run_native_restraint_search(
        system,
        runtime=args.runtime,
        frequency=args.frame_frequency,
        temperature=temperature,
        protocol=args.protocol,
        timestep=args.timestep,
        cutoff_type=args.cutoff_type,
        cutoff=args.cutoff,
        perturbable_constraint=args.perturbable_constraint,
    )

    save_native_restraint(
        restraints=restraints,
        output_dir=output_dir,
        ligand_name=args.ligand_name,
    )

    save_native_starting_structure(
        starting_structure=starting_structure,
        output_dir=output_dir,
        ligand_name=args.ligand_name,
    )

    _write_correction_file(
        correction_value=correction_value,
        output_dir=output_dir,
        ligand_name=args.ligand_name,
    )


if __name__ == "__main__":
    main()
