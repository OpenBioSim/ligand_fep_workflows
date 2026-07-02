"""
Snakemake rule for Boresch restraint search in ABFE workflow.

This rule runs a short unrestrained MD simulation of the protein-ligand complex,
then analyses the trajectory to find optimal Boresch restraints. Boresch restraints
define 6 degrees of freedom (1 distance, 2 angles, 3 dihedrals) using 3 anchor
points each in the protein and ligand.

The restraint correction term is also computed and saved for use in the
final free energy calculation.
"""

from pathlib import Path


# Get restraint search settings
restraint_settings = config.get("restraint_search_settings", {})

_engine = config.get("engine", config["production-settings"].get("engine", "gromacs")).strip().lower()
_somd2_settings = config["production-settings"].get("somd2-settings", {})
_restraint_style = _somd2_settings.get("restraint_style", "legacy")
_native_restraint = _engine == "somd2" and _restraint_style == "native"
_restraint_ext = "s3" if _native_restraint else "json"


rule restraint_search:
    """
    Find optimal Boresch restraints for ABFE bound leg.

    With restraint_style="legacy" (default), this rule runs a short
    unrestrained equilibration via GROMACS and analyses the trajectory with
    BioSimSpace's Aldeghi-style search (BSS.FreeEnergy.RestraintSearch),
    shared by both the GROMACS and SOMD2 production paths. The output
    restraint.json file contains all information needed to apply Boresch
    restraints during the bound leg calculation.

    With restraint_style="native" (SOMD2-only), it instead runs a short
    unrestrained trajectory via sire/OpenMM and analyses it with sire's
    native sire.restraints.boresch_search(), mirroring SOMD2's own built-in
    restraint auto-generation. The output restraint.s3 file is a native
    sire-serialised sire.mm.BoreschRestraints object.

    In both cases, the analytical correction for restraint release is
    computed and saved to the shared {ligand}_correction.txt file consumed
    by the engine-agnostic analysis pipeline.
    """
    input:
        system=Path(f"{config['working_directory']}/preparation/final")
        / "{ligand}_bound.bss",
    output:
        restraint=Path(f"{config['working_directory']}/restraints")
        / f"{{ligand}}_restraint.{_restraint_ext}",
        correction=Path(f"{config['working_directory']}/restraints")
        / "{ligand}_correction.txt",
    log:
        Path(f"{config['working_directory']}/logs") / "{ligand}_restraint_search.log",
    threads:
        config["simulation_threads"]
    resources:
        gpu=1
    params:
        script=Path("workflow/scripts/abfe/restraint_search.py"),
        output_directory=Path(f"{config['working_directory']}/restraints"),
        method="sire" if _native_restraint else "bss",
        search_runtime=restraint_settings.get("search_runtime", "1ns"),
        temperature=restraint_settings.get("temperature", "300K"),
        protocol=restraint_settings.get("protocol", "rxrx"),
        timestep=_somd2_settings.get("timestep", "4fs"),
        cutoff_type=_somd2_settings.get("cutoff_type", "PME"),
        cutoff=_somd2_settings.get("cutoff", "10A"),
        perturbable_constraint=_somd2_settings.get(
            "perturbable_constraint", "h_bonds_not_heavy_perturbed"
        ),
        frame_frequency=restraint_settings.get("frame_frequency", "10ps"),
    shell:
        """
        echo "Running restraint search for {wildcards.ligand} (method={params.method})"
        python {params.script} \
            --input {input.system} \
            --output-directory {params.output_directory} \
            --ligand-name {wildcards.ligand} \
            --runtime {params.search_runtime} \
            --temperature {params.temperature} \
            --method {params.method} \
            --protocol {params.protocol} \
            --timestep {params.timestep} \
            --cutoff-type {params.cutoff_type} \
            --cutoff {params.cutoff} \
            --perturbable-constraint {params.perturbable_constraint} \
            --frame-frequency {params.frame_frequency} \
            2>&1 | tee {log}
        """
