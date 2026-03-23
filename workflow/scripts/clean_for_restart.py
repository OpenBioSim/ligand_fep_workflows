#!/usr/bin/env python
"""
Clean production .done files and analysis outputs to allow restart or extension.

This script should be run BEFORE updating the runtime in the config and
re-running Snakemake. It removes:
  - production/.done files (so Snakemake reruns production)
  - replica barrier markers
  - analysis outputs (pmf.csv, overlap.npy, plots, summary CSVs)

Equilibration .done files are intentionally preserved so that ABFE GROMACS
min/heat/eq stages are not repeated.

Usage:
    python workflow/scripts/clean_for_restart.py --config config/config_rbfe.yml
    python workflow/scripts/clean_for_restart.py --config config/config_abfe.yml --dry-run
"""

import argparse
from pathlib import Path

import yaml


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove production .done files and analysis outputs for restart.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to workflow config YAML (config_rbfe.yml or config_abfe.yml).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print files that would be removed without deleting anything.",
    )
    return parser.parse_args()


def remove(path: Path, dry_run: bool) -> None:
    if dry_run:
        print(f"  would remove: {path}")
    else:
        path.unlink()
        print(f"  removed: {path}")


def main() -> None:
    args = parse_args()

    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    working_dir = Path(cfg["working_directory"])
    engine = cfg["production-settings"].get("engine", "gromacs").strip().lower()
    prod_dir = working_dir / "production" / engine
    analysis_dir = working_dir / "analysis" / engine

    print(f"Working directory : {working_dir}")
    print(f"Engine            : {engine}")
    if args.dry_run:
        print("Dry run — no files will be deleted.\n")
    else:
        print()

    removed = 0

    # Production .done files
    for path in sorted(prod_dir.glob("**/.done")):
        remove(path, args.dry_run)
        removed += 1

    # Replica barrier markers
    for path in sorted(prod_dir.glob(".replica_*_barrier")):
        remove(path, args.dry_run)
        removed += 1

    # Per-leg analysis outputs
    for pattern in ("**/pmf.csv", "**/overlap.npy", "**/*.png"):
        for path in sorted(analysis_dir.glob(pattern)):
            remove(path, args.dry_run)
            removed += 1

    # Final / detailed summary CSVs in the engine-level directory
    for pattern in ("final_*.csv", "detailed_*.csv"):
        for path in sorted(analysis_dir.glob(pattern)):
            remove(path, args.dry_run)
            removed += 1

    verb = "Would remove" if args.dry_run else "Removed"
    print(f"\n{verb} {removed} file(s).")

    if not args.dry_run and removed > 0:
        print(
            "\nNext steps:\n"
            "  1. Update runtime (and any other settings) in your config.\n"
            "  2. Set  restart: true  under production-settings.\n"
            "  3. Re-run Snakemake — production and analysis will be rerun;\n"
            "     equilibration stages are unchanged.\n"
            "  4. After completion, set  restart: false  to prevent accidental restarts."
        )


if __name__ == "__main__":
    main()
