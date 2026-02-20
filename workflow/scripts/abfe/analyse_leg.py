#!/usr/bin/env python
"""
ABFE leg analysis script.

This script analyses a single ABFE leg (bound or free) to extract the
potential of mean force (PMF) across lambda windows using the MBAR estimator.

The output includes:
    - pmf.csv: Free energy at each lambda value
    - overlap.npy: Overlap matrix for convergence assessment
    - Optional: overlap_matrix.png and pmf.png plots

Usage:
    python analyse_leg.py --input-directory production/ligand/bound_0/
                          --output-directory analysis/ligand/bound_0/

Author: ABFE Workflow
"""

import argparse
from pathlib import Path

import BioSimSpace as BSS
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyse a single ABFE leg.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-directory",
        type=str,
        required=True,
        help="Directory containing the simulation data (lambda_* subdirs).",
    )
    parser.add_argument(
        "--output-directory",
        type=str,
        required=True,
        help="Directory to save analysis results.",
    )
    parser.add_argument(
        "--plot-overlap-matrix",
        action="store_true",
        help="Generate overlap matrix plot.",
    )
    parser.add_argument(
        "--plot-pmf",
        action="store_true",
        help="Generate PMF plot.",
    )
    parser.add_argument(
        "--estimator",
        type=str,
        choices=["MBAR", "TI"],
        default="MBAR",
        help="Free energy estimator to use.",
    )
    parser.add_argument(
        "--temperature",
        type=str,
        default="300K",
        help="Simulation temperature.",
    )
    return parser.parse_args()


def plot_pmf(df: pd.DataFrame, output_dir: Path) -> None:
    """
    Plot the potential of mean force.

    Args:
        df: DataFrame with lambda, free_energy, and error columns
        output_dir: Directory to save the plot
    """
    plt.figure(figsize=(8, 6))
    plt.errorbar(
        df["lambda"],
        df["free_energy (kcal/mol)"],
        yerr=df["error (kcal/mol)"],
        marker="o",
        capsize=3,
        color="steelblue",
        linewidth=1.5,
        markersize=6,
    )
    plt.xlabel(r"$\lambda$", fontsize=12)
    plt.ylabel("PMF (kcal/mol)", fontsize=12)
    plt.title("Potential of Mean Force")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / "pmf.png", dpi=300)
    plt.close()


def plot_overlap_matrix(
    overlap: np.ndarray,
    output_dir: Path,
    color_bar_cutoffs: list[float] = [0.03, 0.1, 0.3],
) -> None:
    """
    Plot the overlap matrix for convergence assessment.

    Good overlap (>0.03 off-diagonal) indicates sufficient sampling
    between adjacent lambda windows for reliable MBAR estimates.

    Args:
        overlap: 2D overlap matrix
        output_dir: Directory to save the plot
        color_bar_cutoffs: Cutoffs for colour coding
    """
    num_rows = len(overlap)

    # Define colour scheme
    color_bounds = [0] + color_bar_cutoffs + [1]
    box_colors = ["#FFD3E0", "#88CCEE", "#78C592", "#117733"]
    cmap = colors.ListedColormap(box_colors)
    norm = colors.BoundaryNorm(color_bounds, cmap.N)

    # Create figure
    fig_size = max(8, num_rows / 2)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=150)

    # Plot heatmap
    im = ax.imshow(overlap, cmap=cmap, norm=norm)

    # Add cell separators
    for i in range(num_rows - 1):
        ax.axhline(i + 0.5, color="white", linewidth=0.5)
        ax.axvline(i + 0.5, color="white", linewidth=0.5)

    # Add text labels
    for i in range(num_rows):
        for j in range(num_rows):
            val = overlap[i, j]
            # Choose text colour based on background
            text_color = "white" if val > 0.3 else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                fontsize=8,
                color=text_color,
            )

    # Colorbar
    cbar = fig.colorbar(
        im,
        ax=ax,
        cmap=cmap,
        norm=norm,
        boundaries=color_bounds,
        ticks=color_bounds,
        shrink=0.7,
    )
    cbar.outline.set_visible(False)

    # Labels
    ax.set_xlabel(r"$\lambda$ Index")
    ax.set_ylabel(r"$\lambda$ Index")
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()
    ax.set_xticks(range(num_rows))
    ax.set_yticks(range(num_rows))

    # Remove borders
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / "overlap_matrix.png", dpi=300)
    plt.close()


def main():
    """Main entry point for leg analysis."""
    args = parse_args()

    input_dir = Path(args.input_directory)
    output_dir = Path(args.output_directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Analysing leg: {input_dir}")

    # Parse temperature
    temperature = BSS.Types.Temperature(args.temperature)

    try:
        pmf, overlap = BSS.FreeEnergy.Relative.analyse(
            str(input_dir),
            temperature=temperature,
            estimator=args.estimator,
        )

        # Convert PMF to DataFrame
        lambda_values = [x[0] for x in pmf]
        free_energies = [x[1].value() for x in pmf]
        errors = [x[2].value() for x in pmf]

        df = pd.DataFrame(
            {
                "lambda": lambda_values,
                "free_energy (kcal/mol)": free_energies,
                "error (kcal/mol)": errors,
            }
        )

        # Get final DG value
        final_dg = free_energies[-1]
        final_error = errors[-1]
        print(f"Leg DG: {final_dg:.4f} +/- {final_error:.4f} kcal/mol")

    except Exception as e:
        print(f"Warning: Analysis failed: {e}")
        print("Writing empty results...")

        df = pd.DataFrame(
            {
                "lambda": [],
                "free_energy (kcal/mol)": [],
                "error (kcal/mol)": [],
            }
        )
        overlap = np.array([])

    # Save PMF
    df.to_csv(output_dir / "pmf.csv", index=False)
    print(f"Saved PMF to {output_dir / 'pmf.csv'}")

    # Save overlap matrix
    np.save(output_dir / "overlap.npy", overlap)

    # Generate plots if requested
    if args.plot_overlap_matrix and len(overlap) > 0:
        try:
            plot_overlap_matrix(overlap, output_dir)
            print(f"Saved overlap matrix plot")
        except Exception as e:
            print(f"Warning: Could not plot overlap matrix: {e}")

    if args.plot_pmf and len(df) > 0:
        try:
            plot_pmf(df, output_dir)
            print(f"Saved PMF plot")
        except Exception as e:
            print(f"Warning: Could not plot PMF: {e}")

    print("Leg analysis complete")


if __name__ == "__main__":
    main()
