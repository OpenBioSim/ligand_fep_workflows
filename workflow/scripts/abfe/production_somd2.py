#!/usr/bin/env python
"""
ABFE production simulation script using SOMD2.

This script runs alchemical free energy simulations for ABFE calculations
using SOMD2 (sire/OpenMM). It handles:
    1. Loading the equilibrated system from BSS format
    2. Converting to sire and applying sire-native decoupling
    3. Setting up the lambda schedule (decharge → annihilate)
    4. Configuring Boresch restraints (bound leg only)
    5. Running SOMD2 (minimisation + production)

SOMD2 handles HMR internally with an OpenMM-appropriate factor,
so no external HMR application is needed.

Usage:
    python production_somd2.py --input preparation/final/ejm42_bound.bss \
                               --output-directory production/ejm42/bound_0 \
                               --leg bound \
                               --restraint-file restraints/ejm42_restraint.json \
                               --runtime 2ns --num-lambda 21
"""

import argparse
import json
import sys
from pathlib import Path

import BioSimSpace.Sandpit.Exscientia as BSS
import sire as sr
from somd2.config import Config
from somd2.runner import RepexRunner, Runner


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run ABFE alchemical simulations with SOMD2.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the equilibrated system (.bss file from preparation/final/).",
    )
    parser.add_argument(
        "--output-directory",
        type=str,
        required=True,
        help="Directory for simulation outputs.",
    )
    parser.add_argument(
        "--leg",
        type=str,
        choices=["bound", "free"],
        required=True,
        help="Which leg of the calculation.",
    )
    parser.add_argument(
        "--restraint-file",
        type=str,
        default=None,
        help="Path to restraint parameters JSON (required for bound leg).",
    )
    parser.add_argument(
        "--runtime",
        type=str,
        default="2ns",
        help="Production runtime per lambda window.",
    )
    parser.add_argument(
        "--timestep",
        type=str,
        default="4fs",
        help="Integration timestep.",
    )
    parser.add_argument(
        "--temperature",
        type=str,
        default="298K",
        help="Simulation temperature.",
    )
    parser.add_argument(
        "--cutoff-type",
        type=str,
        default="PME",
        help="Electrostatics cutoff type (PME or RF).",
    )
    parser.add_argument(
        "--cutoff",
        type=str,
        default="10A",
        help="Cutoff distance.",
    )
    parser.add_argument(
        "--num-lambda",
        type=int,
        default=21,
        help="Number of lambda windows per stage.",
    )
    parser.add_argument(
        "--energy-frequency",
        type=str,
        default="1ps",
        help="Energy reporting interval.",
    )
    parser.add_argument(
        "--frame-frequency",
        type=str,
        default="500ps",
        help="Trajectory frame saving interval.",
    )
    parser.add_argument(
        "--checkpoint-frequency",
        type=str,
        default="500ps",
        help="Checkpoint saving interval.",
    )
    parser.add_argument(
        "--integrator",
        type=str,
        default="langevin_middle",
        help="Integrator type.",
    )
    parser.add_argument(
        "--shift-delta",
        type=str,
        default="2.25A",
        help="Soft-core shift delta parameter.",
    )
    parser.add_argument(
        "--perturbable-constraint",
        type=str,
        default="h_bonds_not_heavy_perturbed",
        help="Constraint type for perturbable molecules.",
    )
    parser.add_argument(
        "--equilibration-time",
        type=str,
        default="20ps",
        help="Per-lambda equilibration time before production.",
    )
    parser.add_argument(
        "--runner",
        type=str,
        choices=["repex", "standard"],
        default="repex",
        help="Runner type: 'repex' (replica exchange, default) or 'standard'.",
    )
    return parser.parse_args()


def load_boresch_restraints(system, restraint_file: str, temperature: float):
    """
    Load Boresch restraint parameters from JSON and construct sire restraints.

    The JSON stores molecule-local atom indices (from BSS atom.index()).
    We map these to system-level indices for the sire BoreschRestraint.

    Args:
        system: Sire system with all atoms
        restraint_file: Path to restraint JSON file
        temperature: Temperature in Kelvin (for standard state correction)

    Returns:
        Sire BoreschRestraints object
    """
    with open(restraint_file, "r") as f:
        params = json.load(f)

    # Find the protein and ligand molecules to map molecule-local
    # indices from the JSON to system-level indices.
    protein = None
    ligand = None
    for i in range(system.num_molecules()):
        mol = system[i]
        if mol.num_residues() == 1 and mol.num_atoms() > 5:
            ligand = mol
        elif mol.num_residues() > 1 and protein is None:
            protein = mol

    if protein is None or ligand is None:
        raise ValueError("Could not identify protein and ligand in system")

    # Map molecule-local indices to system-level indices
    ap = params["anchor_points"]
    all_atoms = system.atoms()

    receptor_idx = [
        all_atoms.find(protein.atoms()[ap["r1"]]),
        all_atoms.find(protein.atoms()[ap["r2"]]),
        all_atoms.find(protein.atoms()[ap["r3"]]),
    ]
    ligand_idx = [
        all_atoms.find(ligand.atoms()[ap["l1"]]),
        all_atoms.find(ligand.atoms()[ap["l2"]]),
        all_atoms.find(ligand.atoms()[ap["l3"]]),
    ]

    print(f"  Receptor system indices: {receptor_idx}")
    print(f"  Ligand system indices:   {ligand_idx}")

    # Parse equilibrium values
    ev = params["equilibrium_values"]
    r0 = sr.u(ev["r0"])
    theta0 = [
        sr.u(ev["thetaA0"]),
        sr.u(ev["thetaB0"]),
    ]
    phi0 = [
        sr.u(ev["phiA0"]),
        sr.u(ev["phiB0"]),
        sr.u(ev["phiC0"]),
    ]

    # Parse force constants
    fc = params["force_constants"]
    kr = sr.u(fc["kr"])
    ktheta = [
        sr.u(fc["kthetaA"]),
        sr.u(fc["kthetaB"]),
    ]
    kphi = [
        sr.u(fc["kphiA"]),
        sr.u(fc["kphiB"]),
        sr.u(fc["kphiC"]),
    ]

    # Construct the Boresch restraint using system-level indices
    b = sr.mm.BoreschRestraint(
        receptor=receptor_idx,
        ligand=ligand_idx,
        r0=r0,
        theta0=theta0,
        phi0=phi0,
        kr=kr,
        ktheta=ktheta,
        kphi=kphi,
    )

    correction = sr.restraints.get_standard_state_correction(b, temperature=temperature)
    print(f"Standard state correction: {correction}")

    return sr.mm.BoreschRestraints(b)


def build_lambda_schedule(decoupled_lig):
    """
    Build the ABFE lambda schedule from the decoupled molecule.

    Uses the decharge → annihilate approach:
    - Decharge: ramp charges off while ramping restraints on
    - Annihilate: ramp vdW off with charges off and restraints on

    Args:
        decoupled_lig: Sire molecule with decoupling applied

    Returns:
        Lambda schedule object
    """
    l = decoupled_lig.property("schedule")

    # Remove the default decouple stage
    l.remove_stage("decouple")

    # Stage 1: Decharge (remove electrostatics, ramp up restraints)
    l.add_stage("decharge", equation=l.initial())
    l.set_equation(
        stage="decharge",
        lever="charge",
        equation=l.lam() * l.final() + l.initial() * (1 - l.lam()),
    )
    l.set_equation(
        stage="decharge",
        force="restraint",
        equation=l.lam() * l.final(),
    )

    # Stage 2: Annihilate (remove vdW, charges stay off, restraints stay on)
    l.add_stage(
        "annihilate",
        equation=(-l.lam() + 1) * l.initial() + l.lam() * l.final(),
    )
    l.set_equation(
        stage="annihilate",
        lever="charge",
        equation=l.final(),
    )
    l.set_equation(
        stage="annihilate",
        force="restraint",
        equation=l.final(),
    )

    return l


def main():
    """Main entry point for SOMD2 ABFE production."""
    args = parse_args()

    # Validate bound leg has restraints
    if args.leg == "bound" and args.restraint_file is None:
        print("Error: Bound leg requires --restraint-file.")
        sys.exit(1)

    output_dir = Path(args.output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Load the equilibrated system from BSS format
    bss_system = BSS.Stream.load(args.input)
    print(f"Loaded system with {bss_system.nMolecules()} molecules")

    # Step 2: Convert to Sire format.
    system = sr.system.System(bss_system._sire_object)

    # Step 4: Find the ligand and apply sire-native decoupling
    print("Applying sire-native decoupling...")
    lig = None
    for i in range(system.num_molecules()):
        mol = system[i]
        if mol.num_residues() == 1 and mol.num_atoms() > 5:
            lig = mol
            print(f"  Ligand is molecule {i} ({mol.num_atoms()} atoms)")
            break
    if lig is None:
        raise ValueError("Could not identify ligand in sire system")
    lig = sr.morph.decouple(lig, as_new_molecule=False)
    system.update(lig)

    # Step 5: Build lambda schedule
    print("Building lambda schedule (decharge → annihilate)...")
    lam_schedule = build_lambda_schedule(lig)

    # Step 6: Set up Boresch restraints (bound leg only)
    # Parse temperature value for restraint correction
    temp_str = args.temperature
    temp_value = float("".join(c for c in temp_str if c.isdigit() or c == "."))
    boresch_restraints = None
    if args.leg == "bound" and args.restraint_file:
        print(f"Loading Boresch restraints from {args.restraint_file}...")
        boresch_restraints = load_boresch_restraints(
            system, args.restraint_file, temp_value
        )

    # Step 7: Configure SOMD2
    print("Configuring SOMD2...")
    runner_type = args.runner.strip().lower()
    config = Config(
        runtime=args.runtime,
        timestep=args.timestep,
        temperature=f"{temp_value} K",
        cutoff_type=args.cutoff_type,
        cutoff=args.cutoff,
        num_lambda=args.num_lambda,
        energy_frequency=args.energy_frequency,
        frame_frequency=args.frame_frequency,
        checkpoint_frequency=args.checkpoint_frequency,
        restraints=boresch_restraints,
        lambda_schedule=lam_schedule,
        perturbable_constraint=args.perturbable_constraint,
        output_directory=str(output_dir),
        minimise=True,
        equilibration_time=args.equilibration_time,
        integrator=args.integrator,
        shift_delta=args.shift_delta,
        platform="CUDA",
        replica_exchange=(runner_type == "repex"),
    )

    # Step 8: Run SOMD2
    print(f"Running SOMD2 ({args.leg} leg, runner={args.runner})...")
    if runner_type == "repex":
        # RepexRunner uses threads, safe to call run() directly
        runner = RepexRunner(config=config, system=system)
        runner.run()
    else:
        # Standard Runner uses the spawn start method, so runner.run()
        # must be called inside the if __name__ == "__main__" guard.
        runner = Runner(config=config, system=system)
        return runner

    # Step 9: Clean up temporary AMBER files
    for suffix in [".prm7", ".rst7"]:
        temp_file = Path(str(amber_base) + suffix)
        temp_file.unlink(missing_ok=True)

    print(f"\nSOMD2 production complete ({args.leg} leg)")


if __name__ == "__main__":
    runner = main()
    if runner is not None:
        runner.run()
