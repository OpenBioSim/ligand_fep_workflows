#!/usr/bin/env python
"""
ABFE final analysis script.

This script collates results from ABFE legs and computes the absolute
binding free energy for each ligand. The thermodynamic cycle is:

    DG_bind = DG_free - DG_bound - DG_correction

Where DG_correction is the analytical Boresch restraint correction.

The script:
    1. Reads PMF files from bound and free legs
    2. Extracts the final DG value (at lambda=1)
    3. Combines legs according to the thermodynamic cycle
    4. Averages over replicas with error estimation
    5. Optionally compares to experimental data

Usage:
    python final_analysis.py --analysis-directory analysis/
                            --restraints-directory restraints/
                            --output-directory analysis/
                            --ligands '["ejm42", "ejm43"]'
                            --num-replicas 3

Author: ABFE Workflow
"""

import argparse
import json
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate absolute binding free energies from ABFE simulations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--analysis-directory",
        type=str,
        required=True,
        help="Directory containing leg analysis results.",
    )
    parser.add_argument(
        "--restraints-directory",
        type=str,
        required=True,
        help="Directory containing restraint correction files.",
    )
    parser.add_argument(
        "--output-directory",
        type=str,
        required=True,
        help="Directory to save final results.",
    )
    parser.add_argument(
        "--ligands",
        type=str,
        required=True,
        help="JSON list of ligand names.",
    )
    parser.add_argument(
        "--num-replicas",
        type=int,
        default=3,
        help="Number of replicas per ligand.",
    )
    parser.add_argument(
        "--experimental-results",
        type=str,
        default=None,
        help="Path to experimental results CSV.",
    )
    parser.add_argument(
        "--experimental-units",
        type=str,
        choices=["kcal/mol", "kJ/mol", "Ki_uM"],
        default="kcal/mol",
        help="Units of experimental data.",
    )
    parser.add_argument(
        "--plotting-backend",
        type=str,
        choices=["native", "cinnabar"],
        default="native",
        help="Plotting backend to use.",
    )
    parser.add_argument(
        "--temperature",
        type=float,
        default=300.0,
        help="Temperature in Kelvin for unit conversions.",
    )
    return parser.parse_args()


def read_leg_dg(pmf_file: Path) -> tuple[float, float]:
    """
    Read the final DG value from a leg PMF file.

    Args:
        pmf_file: Path to pmf.csv file

    Returns:
        Tuple of (dg, error) in kcal/mol
    """
    try:
        df = pd.read_csv(pmf_file)
        # Get value at lambda=1
        final_row = df[df["lambda"] == 1.0]
        if len(final_row) == 0:
            # Try getting the last row
            final_row = df.iloc[-1:]

        dg = float(final_row["free_energy (kcal/mol)"].iloc[0])
        error = float(final_row["error (kcal/mol)"].iloc[0])
        return dg, error
    except Exception as e:
        print(f"Warning: Could not read {pmf_file}: {e}")
        return np.nan, np.nan


def read_correction(correction_file: Path) -> float:
    """
    Read the Boresch restraint correction.

    Args:
        correction_file: Path to correction.txt file

    Returns:
        Correction value in kcal/mol
    """
    try:
        with open(correction_file, "r") as f:
            return float(f.read().strip())
    except Exception as e:
        print(f"Warning: Could not read correction from {correction_file}: {e}")
        return 0.0


def calculate_binding_free_energy(
    bound_dg: float,
    free_dg: float,
    correction: float,
) -> float:
    """
    Calculate the absolute binding free energy.

    The thermodynamic cycle is:
        DG_bind = DG_free - DG_bound + DG_correction

    Args:
        bound_dg: DG for bound leg (full transformation)
        free_dg: DG for free leg (full transformation)
        correction: Analytical restraint correction

    Returns:
        Absolute binding free energy in kcal/mol
    """
    dg_bind = free_dg - bound_dg + correction
    return dg_bind


def propagate_error(*errors: float) -> float:
    """
    Propagate errors assuming independence (sum in quadrature).

    Args:
        errors: Individual error values

    Returns:
        Combined error
    """
    return np.sqrt(sum(e**2 for e in errors if not np.isnan(e)))


def convert_experimental_units(
    value: float,
    error: float,
    units: str,
    temperature: float = 300.0,
) -> tuple[float, float]:
    """
    Convert experimental values to kcal/mol.

    Args:
        value: Experimental value
        error: Experimental error
        units: Units of input values
        temperature: Temperature in Kelvin

    Returns:
        Tuple of (value, error) in kcal/mol
    """
    kB = 0.0019872041  # kcal/(mol*K)

    if units == "kcal/mol":
        return value, error
    elif units == "kJ/mol":
        return value * 0.239006, error * 0.239006
    elif units == "Ki_uM":
        # Convert Ki (in uM) to DG using: DG = RT * ln(Ki / C0)
        # where C0 = 1M = 10^6 uM
        dg = kB * temperature * np.log(value / 1e6)
        # Error propagation
        dg_upper = kB * temperature * np.log((value + error) / 1e6)
        dg_lower = kB * temperature * np.log((value - error) / 1e6)
        dg_error = abs(dg_upper - dg_lower) / 2.0
        return dg, dg_error
    else:
        raise ValueError(f"Unknown units: {units}")


def plot_results(
    df: pd.DataFrame,
    df_exp: Optional[pd.DataFrame],
    output_dir: Path,
    temperature: float = 300.0,
    exp_units: str = "kcal/mol",
) -> None:
    """
    Plot computed vs experimental binding free energies.

    Args:
        df: DataFrame with computed results
        df_exp: DataFrame with experimental results (or None)
        output_dir: Directory to save plots
        temperature: Temperature in Kelvin
        exp_units: Units of experimental data
    """
    if df_exp is None:
        print("No experimental data provided, skipping comparison plot.")
        return

    # Merge computed and experimental
    merged = df.merge(df_exp, on="ligand", how="inner")

    if len(merged) == 0:
        print("No matching ligands between computed and experimental data.")
        return

    # Convert experimental units
    exp_dg = []
    exp_err = []
    for _, row in merged.iterrows():
        dg, err = convert_experimental_units(
            row["exp_value"], row["exp_error"], exp_units, temperature
        )
        exp_dg.append(dg)
        exp_err.append(err)

    merged["exp_dg"] = exp_dg
    merged["exp_dg_error"] = exp_err

    # Create plot
    plt.figure(figsize=(8, 8))

    # Plot data points
    plt.errorbar(
        merged["exp_dg"],
        merged["DG_bind"],
        xerr=merged["exp_dg_error"],
        yerr=merged["error"],
        fmt="o",
        color="steelblue",
        capsize=3,
        markersize=8,
        label="Computed",
    )

    # Add ligand labels
    for _, row in merged.iterrows():
        plt.annotate(
            row["ligand"],
            (row["exp_dg"], row["DG_bind"]),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=8,
        )

    # Add reference lines
    lims = [
        min(merged["exp_dg"].min(), merged["DG_bind"].min()) - 1,
        max(merged["exp_dg"].max(), merged["DG_bind"].max()) + 1,
    ]
    plt.plot(lims, lims, "k--", alpha=0.5, label="y=x")

    # Add +/- 1 kcal/mol bands
    plt.fill_between(
        lims,
        [l - 1 for l in lims],
        [l + 1 for l in lims],
        alpha=0.2,
        color="gray",
        label=r"$\pm$ 1 kcal/mol",
    )

    plt.xlim(lims)
    plt.ylim(lims)
    plt.xlabel(r"Experimental $\Delta G_{bind}$ (kcal/mol)", fontsize=12)
    plt.ylabel(r"Computed $\Delta G_{bind}$ (kcal/mol)", fontsize=12)
    plt.title("ABFE Results")
    plt.legend(loc="upper left")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.savefig(output_dir / "abfe_comparison.png", dpi=300)
    plt.close()

    # Calculate statistics
    diff = merged["DG_bind"] - merged["exp_dg"]
    rmse = np.sqrt(np.mean(diff**2))
    mae = np.mean(np.abs(diff))
    r = np.corrcoef(merged["exp_dg"], merged["DG_bind"])[0, 1]

    print(f"\nStatistics:")
    print(f"  RMSE: {rmse:.2f} kcal/mol")
    print(f"  MAE:  {mae:.2f} kcal/mol")
    print(f"  R:    {r:.3f}")


def main():
    """Main entry point for final ABFE analysis."""
    args = parse_args()

    # Parse ligand list
    ligands = json.loads(args.ligands)
    print(f"Analysing {len(ligands)} ligands: {ligands}")

    analysis_dir = Path(args.analysis_directory)
    restraints_dir = Path(args.restraints_directory)
    output_dir = Path(args.output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Collect results for each ligand and replica
    results = []

    for ligand in ligands:
        print(f"\nProcessing {ligand}...")

        # Read restraint correction
        correction = read_correction(restraints_dir / f"{ligand}_correction.txt")
        print(f"  Restraint correction: {correction:.4f} kcal/mol")

        for replica in range(args.num_replicas):
            # Read bound leg
            bound_dg, bound_err = read_leg_dg(
                analysis_dir / ligand / f"bound_{replica}" / "pmf.csv"
            )

            # Read free leg
            free_dg, free_err = read_leg_dg(
                analysis_dir / ligand / f"free_{replica}" / "pmf.csv"
            )

            # Calculate binding free energy
            dg_bind = calculate_binding_free_energy(
                bound_dg,
                free_dg,
                correction,
            )

            # Propagate error
            error = propagate_error(bound_err, free_err)

            results.append(
                {
                    "ligand": ligand,
                    "replica": replica,
                    "bound_dg": bound_dg,
                    "free_dg": free_dg,
                    "correction": correction,
                    "DG_bind": dg_bind,
                    "error": error,
                }
            )

            print(
                f"  Replica {replica}: DG_bind = {dg_bind:.4f} +/- {error:.4f} kcal/mol"
            )

    # Convert to DataFrame
    df_all = pd.DataFrame(results)

    # Save detailed results
    df_all.to_csv(output_dir / "detailed_abfe_results.csv", index=False)

    # Average over replicas
    df_summary = (
        df_all.groupby("ligand")
        .agg(
            {
                "DG_bind": "mean",
                "error": lambda x: np.sqrt(np.sum(x**2)) / len(x),  # SEM
                "bound_dg": "mean",
                "free_dg": "mean",
                "correction": "first",
            }
        )
        .reset_index()
    )

    # Add standard deviation over replicas
    df_summary["DG_bind_std"] = df_all.groupby("ligand")["DG_bind"].std().values

    # Final error is max of propagated error and replica std
    df_summary["error"] = np.maximum(df_summary["error"], df_summary["DG_bind_std"])

    # Save summary results
    df_summary.to_csv(output_dir / "final_abfe_results.csv", index=False)

    print("\n=== Summary Results ===")
    print(df_summary[["ligand", "DG_bind", "error"]].to_string(index=False))

    # Load experimental data if provided
    df_exp = None
    if args.experimental_results:
        try:
            df_exp = pd.read_csv(args.experimental_results)
            # Standardize column names
            df_exp = df_exp.rename(
                columns={
                    df_exp.columns[0]: "ligand",
                    df_exp.columns[1]: "exp_value",
                    df_exp.columns[2]: "exp_error",
                }
            )
        except Exception as e:
            print(f"Warning: Could not load experimental data: {e}")

    # Generate plots
    plot_results(
        df_summary,
        df_exp,
        output_dir,
        args.temperature,
        args.experimental_units,
    )

    print(f"\nFinal results saved to {output_dir / 'final_abfe_results.csv'}")


if __name__ == "__main__":
    main()
