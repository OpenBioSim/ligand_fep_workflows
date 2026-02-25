#!/usr/bin/env python
"""
RBFE workflow status and monitoring script.

This script displays the current status of the RBFE workflow, including:
- Which preparation stages are complete for each ligand
- Which FEP prep (merged) systems are ready for each edge
- Which production legs are complete for each edge/replica
- Current relative free energy estimates (DDG) for completed analyses

Usage:
    python status.py --working-directory output/rbfe/ --network-file output/rbfe/network/network.dat
    python status.py --working-directory output/rbfe/ --network-file output/rbfe/network/network.dat --results-only
"""

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Terminal width for consistent formatting
TERM_WIDTH = int(os.environ.get("TERM_WIDTH", os.environ.get("COLUMNS", 70)))

# Column widths — shared across all sections for visual consistency
COL_LABEL = 22  # Edge names ("ejm42~ejm46 R0") can be long
COL_STATUS = 10  # Simple [OK]/[  ] columns
COL_WIDE = 22  # Status + speed annotation columns


class Colors:
    """ANSI color codes for terminal output."""

    def __init__(self, enabled: bool = True):
        self.enabled = enabled

    @property
    def GREEN(self) -> str:
        return "\033[92m" if self.enabled else ""

    @property
    def RED(self) -> str:
        return "\033[91m" if self.enabled else ""

    @property
    def YELLOW(self) -> str:
        return "\033[93m" if self.enabled else ""

    @property
    def BLUE(self) -> str:
        return "\033[94m" if self.enabled else ""

    @property
    def CYAN(self) -> str:
        return "\033[96m" if self.enabled else ""

    @property
    def BOLD(self) -> str:
        return "\033[1m" if self.enabled else ""

    @property
    def DIM(self) -> str:
        return "\033[2m" if self.enabled else ""

    @property
    def RESET(self) -> str:
        return "\033[0m" if self.enabled else ""


def supports_color() -> bool:
    """Check if the terminal supports color output."""
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    if not hasattr(sys.stdout, "isatty") or not sys.stdout.isatty():
        return False
    if os.environ.get("TERM", "") == "dumb":
        return False
    return True


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Display RBFE workflow status and results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--working-directory",
        type=str,
        required=True,
        help="Working directory containing workflow outputs.",
    )
    parser.add_argument(
        "--network-file",
        type=str,
        required=True,
        help="Path to network.dat file defining perturbation edges.",
    )
    parser.add_argument(
        "--num-replicas",
        type=int,
        default=3,
        help="Number of replicas per edge.",
    )
    parser.add_argument(
        "--engine",
        type=str,
        default="gromacs",
        help="Engine name for production/analysis directory namespacing.",
    )
    parser.add_argument(
        "--results-only",
        action="store_true",
        help="Only show DDG results, skip status.",
    )
    parser.add_argument(
        "--no-color",
        action="store_true",
        help="Disable colored output.",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def _visible_len(s: str) -> int:
    """Return the visible length of a string, ignoring ANSI escape codes."""
    return len(re.sub(r"\033\[[0-9;]*m", "", s))


def pad(s: str, width: int) -> str:
    """Pad a string to a visible width, accounting for ANSI codes."""
    return s + " " * max(0, width - _visible_len(s))


def status_ok(colors: Colors) -> str:
    return f"{colors.GREEN}[OK]{colors.RESET}"


def status_missing(colors: Colors) -> str:
    return f"{colors.RED}[  ]{colors.RESET}"


def check_file_exists(path: Path, colors: Colors) -> str:
    return status_ok(colors) if path.exists() else status_missing(colors)


def print_header(title: str, colors: Colors) -> None:
    print()
    print(f"{colors.BOLD}{colors.CYAN}{'=' * TERM_WIDTH}{colors.RESET}")
    print(f"{colors.BOLD}{colors.CYAN}{title:^{TERM_WIDTH}}{colors.RESET}")
    print(f"{colors.BOLD}{colors.CYAN}{'=' * TERM_WIDTH}{colors.RESET}")


def print_table_header(columns: list[tuple[str, int]], colors: Colors) -> None:
    """Print a column header row and separator, using pad() for alignment."""
    parts = [pad(name, width) for name, width in columns[:-1]]
    parts.append(columns[-1][0])
    print(f"\n{colors.BOLD}{''.join(parts)}{colors.RESET}")
    print(f"{colors.DIM}{'-' * TERM_WIDTH}{colors.RESET}")


def row_label(edge: str, replica: int, num_replicas: int) -> str:
    """Build the first-column label: 'edge' or 'edge R0'."""
    if num_replicas == 1:
        return edge
    return f"{edge} R{replica}"


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


def read_network(network_file: str) -> pd.DataFrame:
    """Read the network.dat file to get perturbation edges."""
    return pd.read_csv(
        network_file,
        sep=r"\s+",
        names=["ligand1", "ligand2", "num_lambda", "lambda_windows"],
    )


def get_leg_dg(pmf_file: Path) -> Optional[float]:
    """Read the final DG value from a PMF file."""
    if not pmf_file.exists():
        return None
    try:
        df = pd.read_csv(pmf_file)
        if len(df) == 0:
            return None
        final_row = df[df["lambda"] == 1.0]
        if len(final_row) == 0:
            final_row = df.iloc[-1:]
        return float(final_row["free_energy (kcal/mol)"].iloc[0])
    except Exception:
        return None


def _check_production_complete(prod_dir: Path) -> bool:
    """Check if a production directory has completed output."""
    return (prod_dir / ".done").exists()


def _get_gromacs_performance(prod_dir: Path) -> Optional[float]:
    """Extract average per-window ns/day from a GROMACS production directory."""
    ns_per_day_values = []
    for lambda_dir in sorted(prod_dir.glob("lambda_*")):
        log_file = lambda_dir / "gromacs.log"
        if not log_file.exists():
            continue
        try:
            with open(log_file, "rb") as f:
                f.seek(0, 2)
                size = f.tell()
                f.seek(max(0, size - 1024))
                tail = f.read().decode("utf-8", errors="replace")
            for line in tail.splitlines():
                if line.strip().startswith("Performance:"):
                    parts = line.split()
                    ns_per_day_values.append(float(parts[1]))
                    break
        except (OSError, ValueError, IndexError):
            continue
    if ns_per_day_values:
        return sum(ns_per_day_values) / len(ns_per_day_values)
    return None


def _get_somd2_performance(prod_dir: Path) -> Optional[float]:
    """Extract per-window ns/day from a SOMD2 production directory.

    Searches backwards through the log so that re-runs with overwrite=True
    (which append new output) don't obscure the timing from the original run.
    Returns the last occurrence of 'Overall performance:'.
    """
    log_file = prod_dir / "log.txt"
    if not log_file.exists():
        return None
    chunk_size = 8192
    target = b"Overall performance:"
    try:
        with open(log_file, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            offset = size
            carry = b""
            while offset > 0:
                read_size = min(chunk_size, offset)
                offset -= read_size
                f.seek(offset)
                chunk = f.read(read_size) + carry
                pos = chunk.rfind(target)
                if pos != -1:
                    line_end = chunk.find(b"\n", pos)
                    line = chunk[pos : line_end if line_end != -1 else None]
                    parts = (
                        line.decode("utf-8", errors="replace")
                        .split("Overall performance:")[1]
                        .strip()
                        .split()
                    )
                    return float(parts[0])
                carry = chunk[: len(target) - 1]
    except (OSError, ValueError, IndexError):
        pass
    return None


def _get_performance(prod_dir: Path, engine: str) -> Optional[float]:
    """Get performance for the given engine."""
    if engine == "somd2":
        return _get_somd2_performance(prod_dir)
    return _get_gromacs_performance(prod_dir)


def _format_status_with_speed(
    is_complete: bool, perf: Optional[float], colors: Colors
) -> str:
    """Format a status indicator with optional speed annotation."""
    if not is_complete:
        return status_missing(colors)
    if perf is not None:
        return f"{status_ok(colors)} {colors.DIM}{perf:.0f} ns/d{colors.RESET}"
    return status_ok(colors)


def _format_dg(dg: Optional[float], colors: Colors) -> str:
    if dg is None:
        return status_missing(colors)
    return f"{colors.GREEN}{dg:+.2f}{colors.RESET}"


# ---------------------------------------------------------------------------
# Section printers
# ---------------------------------------------------------------------------


def print_preparation_status(
    working_dir: Path, ligands: list[str], colors: Colors
) -> None:
    print_header("PREPARATION STATUS", colors)

    cols = [
        ("Ligand", 16),
        ("Setup", COL_STATUS),
        ("Free", COL_STATUS),
        ("Bound", COL_STATUS),
    ]
    print_table_header(cols, colors)

    for ligand in ligands:
        setup_free = working_dir / "setup" / f"{ligand}_free.bss"
        setup_bound = working_dir / "setup" / f"{ligand}_bound.bss"
        setup_status = (
            status_ok(colors)
            if setup_free.exists() and setup_bound.exists()
            else status_missing(colors)
        )

        prep_free = check_file_exists(
            working_dir / "preparation" / "final" / f"{ligand}_free.bss", colors
        )
        prep_bound = check_file_exists(
            working_dir / "preparation" / "final" / f"{ligand}_bound.bss", colors
        )

        print(
            f"{pad(ligand, 16)}"
            f"{pad(setup_status, COL_STATUS)}"
            f"{pad(prep_free, COL_STATUS)}"
            f"{prep_bound}"
        )


def print_rbfe_prep_status(
    working_dir: Path, edges: pd.DataFrame, colors: Colors
) -> None:
    print_header("RBFE PREP STATUS", colors)

    cols = [
        ("Edge", COL_LABEL),
        ("Free", COL_STATUS),
        ("Bound", COL_STATUS),
    ]
    print_table_header(cols, colors)

    for _, row in edges.iterrows():
        pair = f"{row['ligand1']}~{row['ligand2']}"
        free_status = check_file_exists(
            working_dir / "rbfe_prepared" / "free" / f"{pair}.bss", colors
        )
        bound_status = check_file_exists(
            working_dir / "rbfe_prepared" / "bound" / f"{pair}.bss", colors
        )

        print(f"{pad(pair, COL_LABEL)}{pad(free_status, COL_STATUS)}{bound_status}")


def print_production_status(
    working_dir: Path,
    edges: pd.DataFrame,
    num_replicas: int,
    engine: str,
    colors: Colors,
) -> None:
    print_header("PRODUCTION STATUS", colors)

    cols = [
        ("Edge", COL_LABEL),
        ("Bound", COL_WIDE),
        ("Free", COL_WIDE),
    ]
    print_table_header(cols, colors)

    for i, (_, row) in enumerate(edges.iterrows()):
        pair = f"{row['ligand1']}~{row['ligand2']}"

        if i > 0 and num_replicas > 1:
            print()

        for replica in range(num_replicas):
            label = row_label(pair, replica, num_replicas)

            bound_dir = working_dir / "production" / engine / pair / f"bound_{replica}"
            bound_complete = _check_production_complete(bound_dir)
            bound_perf = _get_performance(bound_dir, engine) if bound_complete else None
            bound_str = _format_status_with_speed(bound_complete, bound_perf, colors)

            free_dir = working_dir / "production" / engine / pair / f"free_{replica}"
            free_complete = _check_production_complete(free_dir)
            free_perf = _get_performance(free_dir, engine) if free_complete else None
            free_str = _format_status_with_speed(free_complete, free_perf, colors)

            print(f"{pad(label, COL_LABEL)}{pad(bound_str, COL_WIDE)}{free_str}")


def print_analysis_status(
    working_dir: Path,
    edges: pd.DataFrame,
    num_replicas: int,
    engine: str,
    colors: Colors,
) -> None:
    print_header("ANALYSIS STATUS", colors)

    analysis_dir = working_dir / "analysis" / engine

    col_analysis = 16  # Narrower columns for 4-column analysis layout
    cols = [
        ("Edge", COL_LABEL),
        ("Bound", col_analysis),
        ("Free", col_analysis),
        ("DDG", col_analysis),
    ]
    print_table_header(cols, colors)

    for i, (_, row) in enumerate(edges.iterrows()):
        pair = f"{row['ligand1']}~{row['ligand2']}"

        if i > 0 and num_replicas > 1:
            print()

        replica_ddgs = []

        for replica in range(num_replicas):
            label = row_label(pair, replica, num_replicas)

            free_pmf = analysis_dir / pair / f"free_{replica}" / "pmf.csv"
            bound_pmf = analysis_dir / pair / f"bound_{replica}" / "pmf.csv"

            free_dg = get_leg_dg(free_pmf)
            bound_dg = get_leg_dg(bound_pmf)

            if free_dg is not None and bound_dg is not None:
                ddg = bound_dg - free_dg
                replica_ddgs.append(ddg)
                ddg_str = f"{colors.GREEN}{colors.BOLD}{ddg:+.2f}{colors.RESET}"
            else:
                ddg_str = status_missing(colors)

            print(
                f"{pad(label, COL_LABEL)}"
                f"{pad(_format_dg(bound_dg, colors), col_analysis)}"
                f"{pad(_format_dg(free_dg, colors), col_analysis)}"
                f"{ddg_str}"
            )

        if len(replica_ddgs) > 1:
            mean_ddg = np.mean(replica_ddgs)
            std_ddg = np.std(replica_ddgs)
            print(f"{colors.DIM}{'-' * TERM_WIDTH}{colors.RESET}")
            avg_str = f"{colors.CYAN}{colors.BOLD}{mean_ddg:+.2f} +/- {std_ddg:.2f}{colors.RESET}"
            print(
                f"{pad('', COL_LABEL)}"
                f"{pad('', col_analysis)}"
                f"{pad('', col_analysis)}"
                f"{avg_str}"
            )


def print_results_summary(
    working_dir: Path,
    edges: pd.DataFrame,
    num_replicas: int,
    engine: str,
    colors: Colors,
) -> None:
    print_header("RESULTS SUMMARY", colors)

    analysis_dir = working_dir / "analysis" / engine

    # Check if final results file exists
    final_results = analysis_dir / "final_simulation_results.csv"
    if final_results.exists():
        cols = [
            ("Edge", COL_LABEL),
            ("DDG (kcal/mol)", 20),
            ("Error", COL_STATUS),
        ]
        print_table_header(cols, colors)
        df = pd.read_csv(final_results)
        for _, row in df.iterrows():
            edge = f"{row['ligand1']}~{row['ligand2']}"
            dg_str = f"{colors.GREEN}{row['DDG']:+.3f}{colors.RESET}"
            err_str = f"{colors.YELLOW}+/-{row['error']:.3f}{colors.RESET}"
            print(f"{pad(edge, COL_LABEL)}{pad(dg_str, 20)}{err_str}")
        return

    # Otherwise calculate from available PMFs
    cols = [
        ("Edge", COL_LABEL),
        ("DDG (kcal/mol)", 20),
        ("Replicas", COL_STATUS),
    ]
    print(f"\n{colors.YELLOW}Estimates from completed legs:{colors.RESET}")
    print_table_header(cols, colors)

    for _, row in edges.iterrows():
        pair = f"{row['ligand1']}~{row['ligand2']}"
        replica_ddgs = []

        for replica in range(num_replicas):
            free_dg = get_leg_dg(analysis_dir / pair / f"free_{replica}" / "pmf.csv")
            bound_dg = get_leg_dg(analysis_dir / pair / f"bound_{replica}" / "pmf.csv")

            if free_dg is not None and bound_dg is not None:
                replica_ddgs.append(bound_dg - free_dg)

        complete_color = (
            colors.GREEN if len(replica_ddgs) == num_replicas else colors.YELLOW
        )
        count_str = f"{complete_color}{len(replica_ddgs)}/{num_replicas}{colors.RESET}"

        if replica_ddgs:
            mean_ddg = np.mean(replica_ddgs)
            if len(replica_ddgs) > 1:
                std_ddg = np.std(replica_ddgs)
                dg_str = (
                    f"{complete_color}{mean_ddg:+.3f} +/- {std_ddg:.3f}{colors.RESET}"
                )
            else:
                dg_str = f"{complete_color}{mean_ddg:+.3f}{colors.RESET}"
        else:
            dg_str = f"{colors.RED}(incomplete){colors.RESET}"

        print(f"{pad(pair, COL_LABEL)}{pad(dg_str, 20)}{count_str}")

    print(f"{colors.DIM}{'-' * TERM_WIDTH}{colors.RESET}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    use_color = supports_color() and not args.no_color
    colors = Colors(enabled=use_color)

    working_dir = Path(args.working_directory)
    network_file = Path(args.network_file)
    engine = args.engine

    if not working_dir.exists():
        print(
            f"{colors.RED}Working directory does not exist: {working_dir}{colors.RESET}"
        )
        print(
            f"{colors.YELLOW}No workflow outputs found. Run the workflow first.{colors.RESET}"
        )
        return

    if not network_file.exists():
        print(f"{colors.RED}Network file does not exist: {network_file}{colors.RESET}")
        print(f"{colors.YELLOW}Run network preparation first.{colors.RESET}")
        return

    # Read network to discover edges
    edges = read_network(args.network_file)

    # Extract unique ligands from edges
    ligands = sorted(set(edges["ligand1"].tolist() + edges["ligand2"].tolist()))

    print_header("RBFE WORKFLOW STATUS", colors)
    print(f"\n  {colors.BOLD}Working Directory:{colors.RESET} {working_dir}")
    print(f"  {colors.BOLD}Engine:{colors.RESET}            {engine}")
    print(f"  {colors.BOLD}Network:{colors.RESET}           {network_file}")
    print(f"  {colors.BOLD}Ligands:{colors.RESET}           {', '.join(ligands)}")
    print(f"  {colors.BOLD}Edges:{colors.RESET}             {len(edges)}")
    print(f"  {colors.BOLD}Replicas:{colors.RESET}          {args.num_replicas}")

    if args.results_only:
        print_results_summary(working_dir, edges, args.num_replicas, engine, colors)
    else:
        print_preparation_status(working_dir, ligands, colors)
        print_rbfe_prep_status(working_dir, edges, colors)
        print_production_status(working_dir, edges, args.num_replicas, engine, colors)
        print_analysis_status(working_dir, edges, args.num_replicas, engine, colors)
        print_results_summary(working_dir, edges, args.num_replicas, engine, colors)

    print()


if __name__ == "__main__":
    main()
