#!/usr/bin/env python
"""
ABFE production simulation script using unified single-leg protocol.

This script runs alchemical free energy simulations for ABFE calculations
using a unified approach where bonded/coul/vdw lambdas are all controlled
in a single DataFrame schedule. This uses `perturbation_type="full"` which
applies bonded soft-core by default.

For GROMACS, the script:
    1. Sets up simulation directories using BioSimSpace
    2. Runs minimisation, heating, and equilibration for each lambda
    3. Runs production for each lambda window

Usage:
    python production.py --input system.bss --output-directory production/
                         --ligand-name ejm42 --leg bound
                         --lambda-schedule '{"bonded": [...], "coul": [...], "vdw": [...]}'

Author: ABFE Workflow
"""

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Optional

import BioSimSpace.Sandpit.Exscientia as BSS
import pandas as pd
import sire as sr


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Run ABFE alchemical production simulations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the prepared system (.bss file).",
    )
    parser.add_argument(
        "--output-directory",
        type=str,
        required=True,
        help="Directory for simulation outputs.",
    )
    parser.add_argument(
        "--ligand-name",
        type=str,
        required=True,
        help="Name of the ligand.",
    )
    parser.add_argument(
        "--leg",
        type=str,
        choices=["bound", "free"],
        required=True,
        help="Which leg of the calculation.",
    )
    parser.add_argument(
        "--lambda-schedule-file",
        type=str,
        required=True,
        help="Path to JSON file containing the lambda schedule.",
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
        default="300K",
        help="Simulation temperature.",
    )
    parser.add_argument(
        "--pressure",
        type=str,
        default="1bar",
        help="Simulation pressure.",
    )
    parser.add_argument(
        "--report-interval",
        type=str,
        default="1ps",
        help="Energy reporting interval (e.g., '1ps'). Converted to steps using the timestep.",
    )
    parser.add_argument(
        "--restart-interval",
        type=str,
        default="500ps",
        help="Checkpoint interval (e.g., '500ps'). Converted to steps using the timestep.",
    )
    return parser.parse_args()


def create_restraint_from_json(
    system: BSS._SireWrappers.System,
    restraint_file: str,
    temperature: BSS.Types.Temperature,
) -> BSS.FreeEnergy.Restraint:
    """
    Create a BSS Restraint object from saved JSON parameters.

    Args:
        system: The molecular system (must contain the ligand and protein)
        restraint_file: Path to restraint JSON file
        temperature: Simulation temperature

    Returns:
        BSS.FreeEnergy.Restraint object
    """
    with open(restraint_file, "r") as f:
        params = json.load(f)

    # Find the protein and ligand molecules
    # Indices in the JSON are molecule-local, not system-global
    protein = None
    ligand = None
    for mol in system:
        if mol.isDecoupled():
            ligand = mol
        elif mol.nResidues() > 1:  # Protein has many residues
            protein = mol

    if protein is None or ligand is None:
        raise ValueError("Could not identify protein and decoupled ligand in system")

    protein_atoms = protein.getAtoms()
    ligand_atoms = ligand.getAtoms()

    # Build anchor points with actual atom objects using molecule-local indices
    anchor_points = {
        "r1": protein_atoms[params["anchor_points"]["r1"]],
        "r2": protein_atoms[params["anchor_points"]["r2"]],
        "r3": protein_atoms[params["anchor_points"]["r3"]],
        "l1": ligand_atoms[params["anchor_points"]["l1"]],
        "l2": ligand_atoms[params["anchor_points"]["l2"]],
        "l3": ligand_atoms[params["anchor_points"]["l3"]],
    }

    def parse_unit(string):
        return BSS.Types._GeneralUnit(sr.u(string))

    # Parse equilibrium values
    equilibrium_values = {
        "r0": parse_unit(params["equilibrium_values"]["r0"]),
        "thetaA0": parse_unit(params["equilibrium_values"]["thetaA0"]),
        "thetaB0": parse_unit(params["equilibrium_values"]["thetaB0"]),
        "phiA0": parse_unit(params["equilibrium_values"]["phiA0"]),
        "phiB0": parse_unit(params["equilibrium_values"]["phiB0"]),
        "phiC0": parse_unit(params["equilibrium_values"]["phiC0"]),
    }

    # Parse force constants
    force_constants = {
        "kr": parse_unit(params["force_constants"]["kr"]),
        "kthetaA": parse_unit(params["force_constants"]["kthetaA"]),
        "kthetaB": parse_unit(params["force_constants"]["kthetaB"]),
        "kphiA": parse_unit(params["force_constants"]["kphiA"]),
        "kphiB": parse_unit(params["force_constants"]["kphiB"]),
        "kphiC": parse_unit(params["force_constants"]["kphiC"]),
    }

    # Build the restraint dictionary
    restraint_dict = {
        "anchor_points": anchor_points,
        "equilibrium_values": equilibrium_values,
        "force_constants": force_constants,
    }

    # Create and return the Restraint object
    return BSS.FreeEnergy.Restraint(
        system, restraint_dict, temperature, restraint_type="Boresch"
    )


def setup_gromacs_abfe(
    system: BSS._SireWrappers.System,
    lam_vals_df: pd.DataFrame,
    output_dir: Path,
    leg: str,
    restraint: Optional[BSS.FreeEnergy.Restraint],
    runtime: BSS.Types.Time,
    timestep: BSS.Types.Time,
    temperature: BSS.Types.Temperature,
    pressure: BSS.Types.Pressure,
    report_interval: int,
    restart_interval: int,
) -> None:
    """
    Set up GROMACS ABFE simulations using BioSimSpace unified protocol.

    This creates the directory structure and input files for:
        1. Minimisation at each lambda
        2. NVT heating at each lambda
        3. NPT equilibration at each lambda
        4. Production at each lambda

    Uses `perturbation_type="full"` for unified bonded/coul/vdw transformations.

    Args:
        system: Prepared system with decoupled ligand
        lam_vals_df: DataFrame with columns for different lambda terms
        output_dir: Base output directory
        leg: Calculation leg (bound, free)
        restraint: BSS Boresch restraint object (bound leg only)
        runtime: Production runtime per window
        timestep: Integration timestep
        temperature: Simulation temperature
        pressure: Simulation pressure
        report_interval: Steps between energy reports
        restart_interval: Steps between checkpoints
    """
    print(f"Setting up GROMACS ABFE for {leg} leg (unified protocol)...")
    print(f"Lambda schedule ({len(lam_vals_df)} windows):\n{lam_vals_df}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create initial lambda state (first row of DataFrame)
    lam = pd.Series({col: lam_vals_df[col].iloc[0] for col in lam_vals_df.columns})
    print(f"Initial lambda state: {lam.to_dict()}")

    if restraint is not None:
        print("Using Boresch restraints for bound leg")

    # Setup minimisation
    print("Setting up minimisation...")
    min_protocol = BSS.Protocol.FreeEnergyMinimisation(
        lam=lam,
        lam_vals=lam_vals_df,
        perturbation_type="full",
    )
    BSS.FreeEnergy.AlchemicalFreeEnergy(
        system,
        min_protocol,
        engine="gromacs",
        work_dir=str(output_dir / "minimisation"),
        restraint=restraint,
        setup_only=True,
    )

    # Setup NVT heating
    print("Setting up NVT heating...")
    nvt_protocol = BSS.Protocol.FreeEnergyEquilibration(
        lam=lam,
        lam_vals=lam_vals_df,
        perturbation_type="full",
        pressure=None,  # NVT
        temperature_start=BSS.Types.Temperature("100K"),
        temperature_end=temperature,
        timestep=BSS.Types.Time("2fs"),
        runtime=BSS.Types.Time("10ps"),
        report_interval=2500,
        restart_interval=2500,
    )
    BSS.FreeEnergy.AlchemicalFreeEnergy(
        system,
        nvt_protocol,
        engine="gromacs",
        work_dir=str(output_dir / "heat"),
        restraint=restraint,
        setup_only=True,
    )

    # Setup NPT equilibration
    print("Setting up NPT equilibration...")
    npt_protocol = BSS.Protocol.FreeEnergyEquilibration(
        lam=lam,
        lam_vals=lam_vals_df,
        perturbation_type="full",
        pressure=pressure,
        temperature=temperature,
        timestep=BSS.Types.Time("2fs"),
        runtime=BSS.Types.Time("20ps"),
        report_interval=2500,
        restart_interval=2500,
    )
    BSS.FreeEnergy.AlchemicalFreeEnergy(
        system,
        npt_protocol,
        engine="gromacs",
        work_dir=str(output_dir / "eq"),
        restraint=restraint,
        setup_only=True,
        ignore_warnings=True,
    )

    # Setup production
    print("Setting up production...")
    prod_protocol = BSS.Protocol.FreeEnergy(
        lam=lam,
        lam_vals=lam_vals_df,
        perturbation_type="full",
        runtime=runtime,
        timestep=timestep,
        temperature=temperature,
        pressure=pressure,
        report_interval=report_interval,
        restart_interval=restart_interval,
    )
    BSS.FreeEnergy.AlchemicalFreeEnergy(
        system,
        prod_protocol,
        engine="gromacs",
        work_dir=str(output_dir),
        restraint=restraint,
        setup_only=True,
        ignore_warnings=True,
    )

    print("Setup complete.")


def main():
    """Main entry point for ABFE production."""
    args = parse_args()

    # Validate bound leg has restraints
    if args.leg == "bound" and args.restraint_file is None:
        print("Warning: No restraint file provided for bound leg.")

    # Parse lambda schedule from JSON file
    try:
        with open(args.lambda_schedule_file, "r") as f:
            lam_schedule_dict = json.load(f)
        lam_vals_df = pd.DataFrame(lam_schedule_dict)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"Error reading lambda schedule file: {e}")
        sys.exit(1)

    # Parse units
    try:
        runtime = BSS.Types.Time(args.runtime)
        timestep = BSS.Types.Time(args.timestep)
        temperature = BSS.Types.Temperature(args.temperature)
        pressure = BSS.Types.Pressure(args.pressure)
    except ValueError as e:
        print(f"Error parsing units: {e}")
        sys.exit(1)

    # Convert time-based intervals to integer steps
    report_interval = math.ceil(sr.u(args.report_interval) / sr.u(args.timestep))
    restart_interval = math.ceil(sr.u(args.restart_interval) / sr.u(args.timestep))

    # Load system
    print(f"Loading system from {args.input}...")
    system = BSS.Stream.load(args.input)

    # Create restraint object if provided (bound leg only)
    restraint = None
    if args.restraint_file and args.leg == "bound":
        print(f"Loading Boresch restraints from {args.restraint_file}...")
        restraint = create_restraint_from_json(system, args.restraint_file, temperature)
        print(f"Restraint correction: {restraint.getCorrection().value():.4f} kcal/mol")

    # Create output directory
    output_dir = Path(args.output_directory)

    # Setup GROMACS simulations (actual running is handled by Snakemake rule)
    setup_gromacs_abfe(
        system=system,
        lam_vals_df=lam_vals_df,
        output_dir=output_dir,
        leg=args.leg,
        restraint=restraint,
        runtime=runtime,
        timestep=timestep,
        temperature=temperature,
        pressure=pressure,
        report_interval=report_interval,
        restart_interval=restart_interval,
    )

    print(f"\nABFE setup complete: {args.ligand_name} {args.leg}")


if __name__ == "__main__":
    main()
