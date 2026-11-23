"""Equilibrate one end state of a merged system with Sire, using a staged protocol
of restrained minimisation and dynamics. The lambda = 0 and lambda = 1 outputs are
the input and 'perturbed_system' for SOMD2, which seeds the windows nearest each end
state from them."""

import argparse
from pathlib import Path

import numpy as np
import sire as sr
import yaml
from sire.legacy.IO import setCoordinates

# The SOMD2 defaults.
CONSTRAINT = "h-bonds"
PERTURBABLE_CONSTRAINT = "h-bonds-not-heavy-perturbed"

# Constraints can fail for perturbable molecules with many ghost atoms. If a
# dynamics stage blows up, run briefly without them and try again.
INSTABILITY_ERRORS = ("NaN", "kinetic energy has exceeded")
FALLBACK_RUNTIME = "10 ps"
FALLBACK_TIMESTEP = "1 fs"


def restraint_atoms(mols, kind):
    if kind == "heavy":
        solutes = [m for m in mols.molecules("not water") if m.num_atoms() > 1]
        selection = "not element H"
    elif kind == "backbone":
        solutes = [m for m in mols.molecules("protein")]
        selection = "atomname N, CA, C, O"
    else:
        return mols[kind]

    if not solutes:
        raise ValueError(f"No atoms found for '{kind}' restraints")
    atoms = solutes[0].atoms(selection)
    for mol in solutes[1:]:
        atoms += mol.atoms(selection)
    return atoms


def run_dynamics(mols, stage, options):
    d = mols.dynamics(
        timestep=stage["timestep"],
        temperature=stage["temperature"],
        pressure=stage.get("pressure"),
        **options,
    )
    d.run(stage["runtime"], save_frequency=0, auto_fix_minimise=False)
    return d.commit()


def run_stage(mols, stage, defaults, lambda_value, platform):
    stage = {**defaults, **stage}
    kind = stage.get("type")

    restraints = None
    if stage.get("restraint") is not None:
        restraints = sr.restraints.positional(
            mols, atoms=restraint_atoms(mols, stage["restraint"]), k=stage["k"]
        )

    options = {
        "lambda_value": lambda_value,
        "restraints": restraints,
        "constraint": stage.get("constraint", CONSTRAINT),
        "perturbable_constraint": stage.get("perturbable_constraint", PERTURBABLE_CONSTRAINT),
        "cutoff": stage.get("cutoff", "9 A"),
        "cutoff_type": "PME",
    }
    if platform is not None:
        options["platform"] = platform

    if kind == "minimisation":
        options["constraint"] = stage.get("constraint", "none")
        options["perturbable_constraint"] = stage.get("perturbable_constraint", "none")
        m = mols.minimisation(**options)
        m.run(max_iterations=stage.get("max_iterations", 10000))
        return m.commit()

    if kind == "dynamics":
        try:
            return run_dynamics(mols, stage, options)
        except Exception as e:
            if not any(text in str(e) for text in INSTABILITY_ERRORS):
                raise
            print(
                f"Dynamics failed ({str(e).splitlines()[0]}). Running "
                f"{stage.get('fallback_runtime', FALLBACK_RUNTIME)} without constraints "
                "on perturbable molecules, then retrying."
            )
        unconstrained = {
            **stage,
            "runtime": stage.get("fallback_runtime", FALLBACK_RUNTIME),
            "timestep": stage.get("fallback_timestep", FALLBACK_TIMESTEP),
        }
        mols = run_dynamics(mols, unconstrained, {**options, "perturbable_constraint": "none"})
        return run_dynamics(mols, stage, options)

    raise ValueError(f"Unknown stage type '{kind}'")


def with_coordinates(system, coords, state, space):
    """Return a copy of 'system' with 'coords' as the coordinates of end state
    'state', keeping those of the other end state."""
    result = sr.system.System(setCoordinates(system._system, coords.tolist()))
    result.set_space(space)

    # setCoordinates writes to 'coordinates0' of perturbable molecules.
    for mol in result.molecules("property is_perturbable"):
        number = f"molnum {mol.number().value()}"
        cursor = mol.cursor()
        cursor[f"coordinates{state}"] = mol.property("coordinates0")
        cursor[f"coordinates{1 - state}"] = system[number].molecule().property(
            f"coordinates{1 - state}"
        )
        result.update(cursor.commit())
    return result


def rescale(system, space):
    """Scale the molecule centres of the lambda = 1 end state of 'system' so that it
    fits in 'space', which has the same shape but a different volume."""
    scale = (space.volume().value() / system.space().volume().value()) ** (1 / 3)
    print(f"Rescaling the lambda = 1 box by {scale:.5f} to match lambda = 0")

    coords = sr.io.get_coords_array(sr.morph.link_to_perturbed(system))
    counts = np.array([mol.num_atoms() for mol in system.molecules()])
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
    centres = np.add.reduceat(coords, starts) / counts[:, None]
    coords += (scale - 1) * np.repeat(centres, counts, axis=0)
    return with_coordinates(system, coords, 1, space)


def reduce_box(system):
    """Put a triclinic box in the reduced form OpenMM requires."""
    space = system.space()
    if hasattr(space, "reduce"):
        space.reduce(0)
        system.set_space(space)


def modify_ghosts(system):
    """Return a copy of 'system' with Ghostly's modifications to the bonded terms of
    ghost atoms, as SOMD2 applies them. Only the copy is simulated, so the saved
    system is unmodified and SOMD2 can apply its own."""
    from ghostly import modify

    try:
        return modify(system)[0]
    except Exception:
        # Angle optimisation can occasionally fail, so try once more, as SOMD2 does.
        return modify(system)[0]


def run_stages(mols, stages, defaults, state, platform, name="stage"):
    for i, stage in enumerate(stages):
        print(f"lambda = {state}: {name} {i + 1}/{len(stages)} {stage}", flush=True)
        mols = run_stage(mols, stage, defaults, float(state), platform)
    return mols


def equilibrate(system, protocol, state, platform, box=None, ghost_modifications=True):
    """Equilibrate end state 'state' of 'system'. If 'box' is given, the result is
    rescaled to it and re-equilibrated at constant volume, so that both end states
    share the same box, as needed for GCMC."""
    link = sr.morph.link_to_reference if state == 0 else sr.morph.link_to_perturbed
    defaults = protocol.get("defaults", {})
    prepare = modify_ghosts if ghost_modifications else (lambda s: s)

    mols = run_stages(link(prepare(system)), protocol["stages"], defaults, state, platform)
    result = with_coordinates(system, sr.io.get_coords_array(mols), state, mols.space())

    if box is not None:
        stages = [{**stage, "pressure": None} for stage in protocol.get("match_box", [])]
        mols = run_stages(
            link(prepare(rescale(result, box))), stages, defaults, state, platform, "match_box stage"
        )
        result = with_coordinates(system, sr.io.get_coords_array(mols), state, mols.space())

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Merged system, e.g. bound.bss")
    parser.add_argument("--protocol", required=True, help="YAML file with the leg's protocol")
    parser.add_argument("--end-state", type=int, choices=[0, 1], required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--match-box",
        default=None,
        help="Equilibrated lambda = 0 system whose box lambda = 1 is rescaled to.",
    )
    parser.add_argument(
        "--ghost-modifications",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Apply Ghostly's modifications to ghost atom bonded terms, as SOMD2 does.",
    )
    parser.add_argument("--platform", default=None, help="e.g. cuda, opencl, cpu")
    args = parser.parse_args()

    protocol = yaml.safe_load(Path(args.protocol).read_text())
    system = sr.stream.load(args.input)
    reduce_box(system)

    box = sr.stream.load(args.match_box).space() if args.match_box else None
    result = equilibrate(
        system, protocol, args.end_state, args.platform, box, args.ghost_modifications
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    sr.stream.save(result, str(output))
    print(f"Saved {output}")


if __name__ == "__main__":
    main()
