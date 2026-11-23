"""Rules for the SOMD2 RBFE protocol: each edge's ligands are merged and solvated
once per leg, both end states of the merged system are equilibrated, and SOMD2 runs
each leg in one direction with its windows seeded from the nearest end state."""

from pathlib import Path

import pandas as pd
import yaml

_wd = Path(config["working_directory"])
_ligand_paths = {Path(f).stem: f for f in LIGAND_FILES}
_ligand_settings = config.get("ligand_settings", {})
_setup = config.get("setup_settings", {})
_mapping = config.get("mapping_settings", {})
_equilibration = config["equilibration_settings"]
_production = config["production-settings"]


wildcard_constraints:
    ligand="[^~/]+",
    ligand1="[^~/]+",
    ligand2="[^~/]+",
    leg="bound|free",
    replica_number=r"\d+",


def _edge_dir(wc):
    return _wd / "edges" / f"{wc.ligand1}~{wc.ligand2}"


def _get_rbfe_pairs():
    network = checkpoints.prep_network.get().output[0]
    data = pd.read_csv(str(network), sep=r"\s+", names=["ligand1", "ligand2", "num_lambda", "lambda_windows"])
    return [f"{row['ligand1']}~{row['ligand2']}" for _, row in data.iterrows()]


def _somd2_config(leg):
    """The SOMD2 settings for a leg, with any per-leg overrides applied."""
    settings = dict(_production["somd2"])
    settings.update((_production.get("leg-overrides") or {}).get(leg, {}))
    return settings


checkpoint prep_network:
    input:
        ligands=LIGAND_FILES,
    output:
        _wd / "network" / "network.dat",
    log:
        _wd / "logs" / "network_prep.log",
    params:
        folder=Path(config["ligands"]).parent,
        num_lambda=_production["somd2"].get("num_lambda", 11),
        output_dir=_wd / "network",
    resources:
        gpu=0,
    shell:
        """
        python workflow/scripts/rbfe/prepare_network.py \
            --ligand-folder {params.folder} \
            --lambda-windows {params.num_lambda} \
            --diff-lambda-windows {params.num_lambda} \
            --one-direction \
            --output-directory {params.output_dir} \
            2>&1 | tee {log}
        """


rule parameterise_ligand:
    input:
        lambda wc: _ligand_paths[wc.ligand],
    output:
        _wd / "parameterised" / "{ligand}.bss",
    log:
        _wd / "logs" / "{ligand}_parameterise.log",
    params:
        forcefield=_ligand_settings.get("forcefield", "openff_unconstrained-2.3.0"),
    resources:
        gpu=0,
    shell:
        """
        python workflow/scripts/rbfe_somd2/parameterise_ligand.py \
            --ligand-file {input} \
            --output {output} \
            --forcefield {params.forcefield} \
            2>&1 | tee {log}
        """


rule prepare_edge:
    input:
        ligand1=_wd / "parameterised" / "{ligand1}.bss",
        ligand2=_wd / "parameterised" / "{ligand2}.bss",
    output:
        bound=_wd / "edges" / "{ligand1}~{ligand2}" / "bound.bss",
        free=_wd / "edges" / "{ligand1}~{ligand2}" / "free.bss",
        info=_wd / "edges" / "{ligand1}~{ligand2}" / "edge.json",
    log:
        _wd / "logs" / "{ligand1}~{ligand2}_prepare.log",
    params:
        output_dir=lambda wc: _edge_dir(wc),
        protein_files=config["protein_files"],
        water_model=_setup.get("water_model", "tip3p"),
        ion_concentration=_setup.get("ion_concentration", 0.15),
        box_type=_setup.get("box_type", "dodecahedron"),
        bound_padding=_setup.get("bound_padding", "1.0 nm"),
        free_padding=_setup.get("free_padding", "1.5 nm"),
        mcs_strategy=_mapping.get("mcs_strategy", "full"),
        ring_flags=" ".join(
            flag
            for flag, key in (
                ("--allow-ring-breaking", "allow_ring_breaking"),
                ("--allow-ring-size-change", "allow_ring_size_change"),
            )
            if _mapping.get(key, False)
        ),
    resources:
        gpu=0,
    shell:
        """
        python workflow/scripts/rbfe_somd2/prepare_edge.py \
            --ligand1 {input.ligand1} \
            --ligand2 {input.ligand2} \
            --protein-files '{params.protein_files}' \
            --output-directory {params.output_dir} \
            --water-model {params.water_model} \
            --ion-concentration {params.ion_concentration} \
            --box-type {params.box_type} \
            --bound-padding '{params.bound_padding}' \
            --free-padding '{params.free_padding}' \
            --mcs-strategy {params.mcs_strategy} \
            {params.ring_flags} \
            2>&1 | tee {log}
        """


def _equilibrate(wc, input, output, log, state):
    protocol = Path(output[0]).with_name(f"{wc.leg}_protocol.yaml")
    protocol.write_text(yaml.safe_dump(_equilibration[wc.leg]))
    args = [
        f"--input {input.system}",
        f"--protocol {protocol}",
        f"--end-state {state}",
        f"--output {output[0]}",
    ]
    if state == 1 and _equilibration.get("match_box", True):
        args.append(f"--match-box {input.reference}")
    if not _equilibration.get("ghost_modifications", True):
        args.append("--no-ghost-modifications")
    shell(f"python workflow/scripts/rbfe_somd2/equilibrate.py {' '.join(args)} 2>&1 | tee {log}")


rule equilibrate_lambda_0:
    input:
        system=_wd / "edges" / "{ligand1}~{ligand2}" / "{leg}.bss",
    output:
        _wd / "edges" / "{ligand1}~{ligand2}" / "{leg}_lambda_0.bss",
    log:
        _wd / "logs" / "{ligand1}~{ligand2}_equilibrate_{leg}_lambda_0.log",
    threads: config["simulation_threads"]
    resources:
        gpu=1,
    run:
        _equilibrate(wildcards, input, output, log, 0)


# lambda = 1 follows lambda = 0, whose box it is rescaled to so that both end
# states share one box, as needed for GCMC.
rule equilibrate_lambda_1:
    input:
        system=_wd / "edges" / "{ligand1}~{ligand2}" / "{leg}.bss",
        reference=_wd / "edges" / "{ligand1}~{ligand2}" / "{leg}_lambda_0.bss",
    output:
        _wd / "edges" / "{ligand1}~{ligand2}" / "{leg}_lambda_1.bss",
    log:
        _wd / "logs" / "{ligand1}~{ligand2}_equilibrate_{leg}_lambda_1.log",
    threads: config["simulation_threads"]
    resources:
        gpu=1,
    run:
        _equilibrate(wildcards, input, output, log, 1)


# Every edge finishes replica N of both legs before any edge starts replica N+1,
# so that each edge has a result as early as possible.
rule replica_barrier:
    input:
        lambda wc: expand(
            _wd / "production" / _engine / "{pair}" / "{leg}_{replica}" / ".done",
            pair=_get_rbfe_pairs(),
            leg=["bound", "free"],
            replica=wc.replica,
        ),
    output:
        touch(_wd / "production" / _engine / ".replica_{replica}_barrier"),
    resources:
        gpu=0,
    priority: 3


rule production:
    priority: 2
    input:
        system=_wd / "edges" / "{ligand1}~{ligand2}" / "{leg}_lambda_0.bss",
        perturbed_system=_wd / "edges" / "{ligand1}~{ligand2}" / "{leg}_lambda_1.bss",
        prev_replica=lambda wc: []
        if int(wc.replica_number) == 0
        else [_wd / "production" / _engine / f".replica_{int(wc.replica_number) - 1}_barrier"],
    output:
        done=_wd / "production" / _engine / "{ligand1}~{ligand2}" / "{leg}_{replica_number}" / ".done",
    log:
        _wd / "logs" / "{ligand1}~{ligand2}_production_{leg}_{replica_number}.log",
    threads: config["simulation_threads"]
    resources:
        gpu=_production.get("gpus_per_job", 1),
    run:
        work = Path(output.done).parent
        settings = work.with_name(f"{work.name}.input.yaml")
        settings.write_text(yaml.safe_dump(_somd2_config(wildcards.leg)))
        fallback = _production.get("fallback_timestep")
        fallback_arg = f"--fallback-timestep '{fallback}'" if fallback else ""
        shell(
            f"python workflow/scripts/rbfe_somd2/run_leg.py "
            f"--system {input.system} --perturbed-system {input.perturbed_system} "
            f"--config {settings} --output-directory {work} {fallback_arg} 2>&1 | tee {log}"
        )
        shell(f"touch {output.done}")
