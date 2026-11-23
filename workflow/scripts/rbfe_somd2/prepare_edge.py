"""Map and merge two parameterised ligands, then solvate the merged molecule for the
bound and free legs. Both end states share the merged system's topology, so they can
be equilibrated separately and used to seed SOMD2's lambda windows."""

import argparse
import glob
import json
import tempfile
import time
import warnings
from pathlib import Path

import BioSimSpace as BSS
import numpy as np
import sire as sr

# BioSimSpace warns with this text when relaxing ring matching gives a larger MCS.
RELAXED_RING_WARNING = "ringMatchesRingOnly"

# Atoms closer than this give a NaN energy. Ghosts that close are moved by
# NUDGE_DISTANCE, which has no effect at the end state where they are ghosts.
COINCIDENT_DISTANCE = 0.01
NUDGE_DISTANCE = 0.05


def map_atoms(lig_a, lig_b, mcs_strategy):
    """Map with BioSimSpace, re-mapping with ring matching relaxed if BioSimSpace
    warns that this gives a larger common core."""
    start = time.time()
    # Only pass the strategy when it isn't the default, since older BioSimSpace
    # versions don't accept it.
    options = {"mcs_strategy": mcs_strategy} if mcs_strategy != "full" else {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mapping = BSS.Align.matchAtoms(lig_a, lig_b, **options)
    mcs_kwargs = {}
    if any(RELAXED_RING_WARNING in str(w.message) for w in caught):
        mcs_kwargs = {"ringMatchesRingOnly": False}
        mapping = BSS.Align.matchAtoms(lig_a, lig_b, mcs_kwargs=mcs_kwargs, **options)
    print(f"{len(mapping)} mapped atoms in {time.time() - start:.1f} s")
    return mapping, {"mcs_strategy": mcs_strategy, "mcs_kwargs": mcs_kwargs}


def nudge_coincident_ghosts(merged):
    """Move ghost atoms that sit on another atom of the merged molecule, in the
    coordinates of each end state."""
    mol = merged._sire_object
    rng = np.random.default_rng(0)
    cursor = mol.cursor()
    for state in (0, 1):
        coords = mol.property(f"coordinates{state}").to_vector()
        x = np.array([[c.x().value(), c.y().value(), c.z().value()] for c in coords])
        ghosts = [i for i, atom in enumerate(mol.atoms()) if atom.property(f"ambertype{state}") == "du"]
        moved = set()
        for i in ghosts:
            while True:
                d = np.linalg.norm(x - x[i], axis=1)
                d[i] = np.inf
                if d.min() >= COINCIDENT_DISTANCE:
                    break
                v = rng.normal(size=3)
                x[i] += NUDGE_DISTANCE * v / np.linalg.norm(v)
                moved.add(i)
        for i in moved:
            cursor.atoms()[i][f"coordinates{state}"] = sr.maths.Vector(*x[i])
        if moved:
            print(f"Moved {len(moved)} coincident ghost atoms at lambda = {state}")
    merged._sire_object = cursor.commit()


def box_width(system, padding):
    """Box width from padding, using the same rule as openmm.app.Modeller."""
    mols = sr.system.System(system._sire_object)
    mols = sr.morph.link_to_reference(mols)
    coords = sr.io.get_coords_array(mols)
    centre = 0.5 * (coords.min(axis=0) + coords.max(axis=0))
    radius = np.linalg.norm(coords - centre, axis=1).max()
    padding = BSS.Types.Length(padding).angstroms().value()
    return max(2 * radius + padding, 2 * padding)


def solvate(system, water_model, padding, box_type, ion_conc):
    width = BSS.Types.Length(box_width(system, padding), "angstrom")
    if box_type == "cubic":
        box, angles = BSS.Box.cubic(width)
    elif box_type == "dodecahedron":
        box, angles = BSS.Box.rhombicDodecahedronSquare(width)
    elif box_type == "octahedron":
        box, angles = BSS.Box.truncatedOctahedron(width)
    else:
        raise ValueError(f"Unsupported box type '{box_type}'")
    print(f"Solvating in a {box_type} box with image distance {width}")
    with tempfile.TemporaryDirectory() as tmp:
        system = BSS.Solvent.solvate(
            water_model,
            molecule=system,
            box=box,
            angles=angles,
            ion_conc=ion_conc,
            is_neutral=True,
            work_dir=tmp,
        )
    # Rounding can leave the box marginally outside the reduced form OpenMM requires.
    system.reduceBoxVectors()
    return system


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ligand1", required=True, help="Parameterised ligand stream file")
    parser.add_argument("--ligand2", required=True, help="Parameterised ligand stream file")
    parser.add_argument("--protein-files", required=True, help="Glob pattern for the protein files")
    parser.add_argument("--output-directory", required=True)
    parser.add_argument("--water-model", default="tip3p")
    parser.add_argument("--ion-concentration", type=float, default=0.15, help="mol/L")
    parser.add_argument("--box-type", choices=["cubic", "dodecahedron", "octahedron"], default="dodecahedron")
    parser.add_argument("--bound-padding", default="1.0 nm")
    parser.add_argument("--free-padding", default="1.5 nm")
    parser.add_argument("--mcs-strategy", choices=["full", "heavy", "both"], default="full")
    parser.add_argument("--allow-ring-breaking", action="store_true")
    parser.add_argument("--allow-ring-size-change", action="store_true")
    args = parser.parse_args()

    lig_a = BSS.Stream.load(args.ligand1)[0]
    lig_b = BSS.Stream.load(args.ligand2)[0]

    mapping, details = map_atoms(lig_a, lig_b, args.mcs_strategy)
    inverse = {b: a for a, b in mapping.items()}
    lig_b = BSS.Align.rmsdAlign(lig_b, lig_a, inverse)
    merged = BSS.Align.merge(
        lig_a,
        lig_b,
        mapping,
        allow_ring_breaking=args.allow_ring_breaking,
        allow_ring_size_change=args.allow_ring_size_change,
    )
    nudge_coincident_ghosts(merged)

    protein_files = sorted(glob.glob(args.protein_files))
    if not protein_files:
        raise ValueError(f"No protein files match '{args.protein_files}'")

    out = Path(args.output_directory)
    out.mkdir(parents=True, exist_ok=True)
    for leg, padding in (("bound", args.bound_padding), ("free", args.free_padding)):
        if leg == "bound":
            system = BSS.IO.readMolecules(protein_files)
            system.addMolecules(merged)
        else:
            system = merged.toSystem()
        system = solvate(system, args.water_model, padding, args.box_type, args.ion_concentration)
        BSS.Stream.save(system, str(out / leg))

    (out / "edge.json").write_text(
        json.dumps(
            {
                "ligand1": Path(args.ligand1).stem,
                "ligand2": Path(args.ligand2).stem,
                "charge1": round(lig_a.charge().value()),
                "charge2": round(lig_b.charge().value()),
                "mapping": {str(a): b for a, b in mapping.items()},
                "mapping_details": details,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
