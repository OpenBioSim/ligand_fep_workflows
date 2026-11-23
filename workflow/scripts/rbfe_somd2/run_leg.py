"""Run one SOMD2 leg from equilibrated end states, seeding the windows nearest each
end state from its own equilibrated structure. A leg that was interrupted restarts
from its checkpoint. A leg that becomes unstable is rerun from scratch with the
fallback timestep, keeping the failed run alongside."""

import argparse
import re
import shutil
import subprocess
from pathlib import Path

import yaml

# How Sire reports an unstable simulation.
UNSTABLE = re.compile(r"NaN|has become unstable")


def run_somd2(system, perturbed_system, config_file, work):
    # SOMD2 runs in the leg's directory, so the paths must be absolute.
    command = [
        "somd2",
        str(Path(system).resolve()),
        "--perturbed-system",
        str(Path(perturbed_system).resolve()),
        "--config",
        str(Path(config_file).resolve()),
        "--output-directory",
        str(Path(work).resolve()),
    ]
    # Replica exchange keeps its state in repex_state.pkl, other runs in checkpoint files.
    restart = (work / "repex_state.pkl").exists() or any(work.glob("checkpoint_*.s3"))
    command += ["--restart"] if restart else ["--overwrite"]
    print(" ".join(command), flush=True)
    # Crash reports are written to the working directory.
    work.mkdir(parents=True, exist_ok=True)
    subprocess.run(command, check=True, cwd=work)


def is_unstable(work):
    log = work / "log.txt"
    return log.exists() and bool(UNSTABLE.search(log.read_text()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", required=True, help="Equilibrated lambda = 0 system")
    parser.add_argument("--perturbed-system", required=True, help="Equilibrated lambda = 1 system")
    parser.add_argument("--config", required=True, help="SOMD2 config file for the leg")
    parser.add_argument("--output-directory", required=True)
    parser.add_argument("--fallback-timestep", default=None, help="e.g. '2 fs'. Omit to disable.")
    args = parser.parse_args()

    work = Path(args.output_directory)
    marker = work.with_name(f"{work.name}.fallback_timestep")
    config_file = work.with_name(f"{work.name}.somd2.yaml")
    config = yaml.safe_load(Path(args.config).read_text())
    if marker.exists():
        config["timestep"] = marker.read_text().strip()
    config_file.write_text(yaml.safe_dump(config))

    try:
        run_somd2(args.system, args.perturbed_system, config_file, work)
    except subprocess.CalledProcessError:
        fallback = args.fallback_timestep
        if fallback is None or str(config.get("timestep")) == fallback or not is_unstable(work):
            raise
        timestep = str(config.get("timestep")).replace(" ", "")
        print(f"Unstable with a {timestep} timestep. Rerunning with {fallback}.", flush=True)
        failed = work.with_name(f"{work.name}.unstable_{timestep}")
        shutil.rmtree(failed, ignore_errors=True)
        work.rename(failed)
        marker.write_text(f"{fallback}\n")
        config["timestep"] = fallback
        config_file.write_text(yaml.safe_dump(config))
        run_somd2(args.system, args.perturbed_system, config_file, work)

    if not any(work.glob("energy_traj_*.parquet")):
        raise RuntimeError(f"SOMD2 produced no energy trajectories in {work}. See {work / 'log.txt'}.")


if __name__ == "__main__":
    main()
