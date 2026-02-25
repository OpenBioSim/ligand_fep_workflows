#!/usr/bin/env python
"""
ABFE workflow status and monitoring script.

This script displays the current status of the ABFE workflow, including:
- Which preparation stages are complete
- Which production legs (bound/free) are complete for each ligand/replica
- Current binding free energy estimates for completed analyses

Usage:
    python status.py --working-directory abfe_runs/ --ligands '["ejm42", "ejm43"]'
    python status.py --working-directory abfe_runs/ --ligands '["ejm42"]' --results-only
    python status.py --working-directory abfe_runs/ --ligands '["ejm42"]' --no-color

Author: ABFE Workflow
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Terminal width for consistent formatting
TERM_WIDTH = int(os.environ.get("TERM_WIDTH", os.environ.get("COLUMNS", 90)))

# Column widths — shared across all sections for visual consistency
COL_LABEL = 18  # Ligand name or "ligand Rn"
COL_STATUS = 14  # Simple [OK]/[  ] or short value columns
COL_WIDE = 20  # Status + speed annotation columns


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
        description="Display ABFE workflow status and results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--working-directory",
        type=str,
        required=True,
        help="Working directory containing workflow outputs.",
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
        "--results-only",
        action="store_true",
        help="Only show binding free energy results, skip status.",
    )
    parser.add_argument(
        "--engine",
        type=str,
        default="gromacs",
        help="Engine name for production/analysis directory namespacing.",
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
    """Print a double-line section header."""
    print()
    print(f"{colors.BOLD}{colors.CYAN}{'=' * TERM_WIDTH}{colors.RESET}")
    print(f"{colors.BOLD}{colors.CYAN}{title:^{TERM_WIDTH}}{colors.RESET}")
    print(f"{colors.BOLD}{colors.CYAN}{'=' * TERM_WIDTH}{colors.RESET}")


def print_table_header(columns: list[tuple[str, int]], colors: Colors) -> None:
    """Print a column header row and separator, using pad() for alignment."""
    parts = [pad(name, width) for name, width in columns[:-1]]
    parts.append(columns[-1][0])  # last column: no padding
    print(f"\n{colors.BOLD}{''.join(parts)}{colors.RESET}")
    print(f"{colors.DIM}{'-' * TERM_WIDTH}{colors.RESET}")


def row_label(ligand: str, replica: int, num_replicas: int) -> str:
    """Build the first-column label: 'ligand' or 'ligand R0'."""
    if num_replicas == 1:
        return ligand
    return f"{ligand} R{replica}"


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------


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


def calculate_binding_energy(
    bound_dg: Optional[float],
    free_dg: Optional[float],
    correction: Optional[float],
) -> Optional[float]:
    if bound_dg is None or free_dg is None or correction is None:
        return None
    return free_dg - bound_dg + correction


def get_correction(restraints_dir: Path, ligand: str) -> Optional[float]:
    correction_file = restraints_dir / f"{ligand}_correction.txt"
    if correction_file.exists():
        try:
            with open(correction_file, "r") as f:
                return float(f.read().strip())
        except Exception:
            pass
    return None


def _check_production_complete(prod_dir: Path, engine: str) -> bool:
    """Check if a production directory has completed output."""
    return (prod_dir / ".done").exists()


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
                # Search for target in this chunk
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
                # Carry the start of this chunk to avoid missing a match split across chunks
                carry = chunk[: len(target) - 1]
    except (OSError, ValueError, IndexError):
        pass
    return None


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


def _get_performance(prod_dir: Path, engine: str) -> Optional[float]:
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
    working_dir: Path, ligands: list[str], engine: str, colors: Colors
) -> None:
    print_header("PREPARATION STATUS", colors)

    is_somd2 = engine == "somd2"

    if is_somd2:
        cols = [
            ("Ligand", COL_LABEL),
            ("Setup", COL_STATUS),
            ("Free", COL_STATUS),
            ("Bound", COL_STATUS),
            ("Restraint", COL_STATUS),
        ]
    else:
        cols = [
            ("Ligand", COL_LABEL),
            ("Setup", COL_STATUS),
            ("Free", COL_STATUS),
            ("Bound", COL_STATUS),
            ("Restraint", COL_STATUS),
            ("ABFE-Prep", COL_STATUS),
        ]
    print_table_header(cols, colors)

    for ligand in ligands:
        setup_free = check_file_exists(
            working_dir / "setup" / f"{ligand}_free.bss", colors
        )
        setup_bound = check_file_exists(
            working_dir / "setup" / f"{ligand}_bound.bss", colors
        )
        setup_status = (
            status_ok(colors)
            if "[OK]" in setup_free and "[OK]" in setup_bound
            else status_missing(colors)
        )

        prep_free = check_file_exists(
            working_dir / "preparation" / "final" / f"{ligand}_free.bss", colors
        )
        prep_bound = check_file_exists(
            working_dir / "preparation" / "final" / f"{ligand}_bound.bss", colors
        )
        restraint = check_file_exists(
            working_dir / "restraints" / f"{ligand}_restraint.json", colors
        )

        cells = [
            pad(ligand, COL_LABEL),
            pad(setup_status, COL_STATUS),
            pad(prep_free, COL_STATUS),
            pad(prep_bound, COL_STATUS),
        ]
        if is_somd2:
            cells.append(restraint)
        else:
            cells.append(pad(restraint, COL_STATUS))
            abfe_free = check_file_exists(
                working_dir / "abfe_prepared" / f"{ligand}_free.bss", colors
            )
            abfe_bound = check_file_exists(
                working_dir / "abfe_prepared" / f"{ligand}_bound.bss", colors
            )
            abfe_prep = (
                status_ok(colors)
                if "[OK]" in abfe_free and "[OK]" in abfe_bound
                else status_missing(colors)
            )
            cells.append(abfe_prep)
        print("".join(cells))


def print_production_status(
    working_dir: Path,
    ligands: list[str],
    num_replicas: int,
    engine: str,
    colors: Colors,
) -> None:
    is_somd2 = engine == "somd2"

    if is_somd2:
        print_header("PRODUCTION STATUS", colors)
        cols = [
            ("Ligand", COL_LABEL),
            ("Bound", COL_WIDE),
            ("Free", COL_WIDE),
        ]
    else:
        print_header("EQUILIBRATION & PRODUCTION STATUS", colors)
        cols = [
            ("Ligand", COL_LABEL),
            ("Eq-Bound", COL_STATUS),
            ("Eq-Free", COL_STATUS),
            ("Bound", COL_WIDE),
            ("Free", COL_WIDE),
        ]
    print_table_header(cols, colors)

    for i, ligand in enumerate(ligands):
        if i > 0 and num_replicas > 1:
            print()  # blank line between ligand groups

        for replica in range(num_replicas):
            label = row_label(ligand, replica, num_replicas)

            # Production status (engine-namespaced)
            bound_dir = (
                working_dir / "production" / engine / ligand / f"bound_{replica}"
            )
            bound_complete = _check_production_complete(bound_dir, engine)
            bound_perf = _get_performance(bound_dir, engine) if bound_complete else None
            bound_str = _format_status_with_speed(bound_complete, bound_perf, colors)

            free_dir = working_dir / "production" / engine / ligand / f"free_{replica}"
            free_complete = _check_production_complete(free_dir, engine)
            free_perf = _get_performance(free_dir, engine) if free_complete else None
            free_str = _format_status_with_speed(free_complete, free_perf, colors)

            cells = [pad(label, COL_LABEL)]

            if not is_somd2:
                # Equilibration status (GROMACS only)
                eq_bound_done = (
                    working_dir
                    / "equilibration"
                    / ligand
                    / f"bound_{replica}"
                    / ".done"
                )
                eq_bound = (
                    status_ok(colors)
                    if eq_bound_done.exists()
                    else status_missing(colors)
                )
                eq_free_done = (
                    working_dir / "equilibration" / ligand / f"free_{replica}" / ".done"
                )
                eq_free = (
                    status_ok(colors)
                    if eq_free_done.exists()
                    else status_missing(colors)
                )
                cells.append(pad(eq_bound, COL_STATUS))
                cells.append(pad(eq_free, COL_STATUS))

            cells.append(pad(bound_str, COL_WIDE))
            cells.append(free_str)
            print("".join(cells))


def print_analysis_status(
    working_dir: Path,
    ligands: list[str],
    num_replicas: int,
    engine: str,
    colors: Colors,
) -> None:
    print_header("ANALYSIS STATUS", colors)

    analysis_dir = working_dir / "analysis" / engine
    restraints_dir = working_dir / "restraints"

    cols = [
        ("Ligand", COL_LABEL),
        ("Bound", COL_STATUS),
        ("Free", COL_STATUS),
        ("Correction", COL_STATUS),
        ("DG_bind", COL_STATUS),
    ]
    print_table_header(cols, colors)

    for i, ligand in enumerate(ligands):
        if i > 0 and num_replicas > 1:
            print()

        correction = get_correction(restraints_dir, ligand)
        correction_str = (
            status_missing(colors) if correction is None else f"{correction:+.2f}"
        )
        replica_dgs = []

        for replica in range(num_replicas):
            label = row_label(ligand, replica, num_replicas)

            bound_pmf = analysis_dir / ligand / f"bound_{replica}" / "pmf.csv"
            free_pmf = analysis_dir / ligand / f"free_{replica}" / "pmf.csv"

            bound_dg = get_leg_dg(bound_pmf)
            free_dg = get_leg_dg(free_pmf)
            dg_bind = calculate_binding_energy(bound_dg, free_dg, correction)

            if dg_bind is not None:
                replica_dgs.append(dg_bind)
                dg_bind_str = f"{colors.GREEN}{colors.BOLD}{dg_bind:+.2f}{colors.RESET}"
            else:
                dg_bind_str = status_missing(colors)

            print(
                f"{pad(label, COL_LABEL)}"
                f"{pad(_format_dg(bound_dg, colors), COL_STATUS)}"
                f"{pad(_format_dg(free_dg, colors), COL_STATUS)}"
                f"{pad(correction_str, COL_STATUS)}"
                f"{dg_bind_str}"
            )

        # Average row for multi-replica
        if len(replica_dgs) > 1:
            mean_dg = np.mean(replica_dgs)
            std_dg = np.std(replica_dgs)
            print(f"{colors.DIM}{'-' * TERM_WIDTH}{colors.RESET}")
            avg_str = f"{colors.CYAN}{colors.BOLD}{mean_dg:+.2f} +/- {std_dg:.2f}{colors.RESET}"
            print(
                f"{pad('', COL_LABEL)}"
                f"{pad('', COL_STATUS)}"
                f"{pad('', COL_STATUS)}"
                f"{avg_str}"
            )


def print_results_summary(
    working_dir: Path,
    ligands: list[str],
    num_replicas: int,
    engine: str,
    colors: Colors,
) -> None:
    print_header("RESULTS SUMMARY", colors)

    analysis_dir = working_dir / "analysis" / engine
    restraints_dir = working_dir / "restraints"

    # Check if final results file exists
    final_results = analysis_dir / "final_abfe_results.csv"
    if final_results.exists():
        cols = [
            ("Ligand", COL_LABEL),
            ("DG_bind (kcal/mol)", 20),
            ("Error", COL_STATUS),
        ]
        print_table_header(cols, colors)
        df = pd.read_csv(final_results)
        for _, row in df.iterrows():
            dg_str = f"{colors.GREEN}{row['DG_bind']:+.3f}{colors.RESET}"
            err_str = f"{colors.YELLOW}+/-{row['error']:.3f}{colors.RESET}"
            print(f"{pad(row['ligand'], COL_LABEL)}{pad(dg_str, 20)}{err_str}")
        return

    # Otherwise calculate from available PMFs
    cols = [
        ("Ligand", COL_LABEL),
        ("DG_bind (kcal/mol)", 20),
        ("Replicas", COL_STATUS),
    ]
    print(f"\n{colors.YELLOW}Estimates from completed legs:{colors.RESET}")
    print_table_header(cols, colors)

    for ligand in ligands:
        correction = get_correction(restraints_dir, ligand)
        replica_dgs = []

        for replica in range(num_replicas):
            bound_dg = get_leg_dg(
                analysis_dir / ligand / f"bound_{replica}" / "pmf.csv"
            )
            free_dg = get_leg_dg(analysis_dir / ligand / f"free_{replica}" / "pmf.csv")
            dg_bind = calculate_binding_energy(bound_dg, free_dg, correction)
            if dg_bind is not None:
                replica_dgs.append(dg_bind)

        complete_color = (
            colors.GREEN if len(replica_dgs) == num_replicas else colors.YELLOW
        )
        count_str = f"{complete_color}{len(replica_dgs)}/{num_replicas}{colors.RESET}"

        if replica_dgs:
            mean_dg = np.mean(replica_dgs)
            if len(replica_dgs) > 1:
                std_dg = np.std(replica_dgs)
                dg_str = (
                    f"{complete_color}{mean_dg:+.3f} +/- {std_dg:.3f}{colors.RESET}"
                )
            else:
                dg_str = f"{complete_color}{mean_dg:+.3f}{colors.RESET}"
        else:
            dg_str = f"{colors.RED}(incomplete){colors.RESET}"

        print(f"{pad(ligand, COL_LABEL)}{pad(dg_str, 20)}{count_str}")

    print(f"{colors.DIM}{'-' * TERM_WIDTH}{colors.RESET}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    use_color = supports_color() and not args.no_color
    colors = Colors(enabled=use_color)

    ligands = json.loads(args.ligands)
    working_dir = Path(args.working_directory)
    engine = args.engine.strip().lower()

    if not working_dir.exists():
        print(
            f"{colors.RED}Working directory does not exist: {working_dir}{colors.RESET}"
        )
        print(
            f"{colors.YELLOW}No workflow outputs found. Run the workflow first.{colors.RESET}"
        )
        return

    print_header("ABFE WORKFLOW STATUS", colors)
    print(f"\n  {colors.BOLD}Working Directory:{colors.RESET} {working_dir}")
    print(f"  {colors.BOLD}Engine:{colors.RESET}            {engine}")
    print(f"  {colors.BOLD}Ligands:{colors.RESET}           {', '.join(ligands)}")
    print(f"  {colors.BOLD}Replicas:{colors.RESET}          {args.num_replicas}")

    if args.results_only:
        print_results_summary(working_dir, ligands, args.num_replicas, engine, colors)
    else:
        print_preparation_status(working_dir, ligands, engine, colors)
        print_production_status(working_dir, ligands, args.num_replicas, engine, colors)
        print_analysis_status(working_dir, ligands, args.num_replicas, engine, colors)
        print_results_summary(working_dir, ligands, args.num_replicas, engine, colors)

    print()


if __name__ == "__main__":
    main()
