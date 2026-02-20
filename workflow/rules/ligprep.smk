from pathlib import Path

protein_files = config["protein_files"]
setup_settings = config["setup_settings"]
_method = config.get("method", "rbfe").strip().lower()


rule setup:
    input:
        ligand=Path("workflow/inputs") / "ligands" / "{ligand}.sdf",
        # RBFE needs network prep to complete first; ABFE has no network dependency
        network=(
            [Path(f"{config['working_directory']}/network") / "network.dat"]
            if _method == "rbfe"
            else []
        ),
    output:
        setup_free_out=Path(f"{config['working_directory']}/setup")
        / "{ligand}_free.bss",
        setup_bound_out=Path(f"{config['working_directory']}/setup")
        / "{ligand}_bound.bss",
    log:
        setup_log=Path(f"{config['working_directory']}/logs")
        / "{ligand}_setup.log",
    params:
        setup_script=Path("workflow/scripts/param_solvate.py"),
        output_directory=Path(config["working_directory"]) / "setup",
        ligand_forcefield=setup_settings["ligand_forcefield"],
        water_model=setup_settings["water_model"],
        ion_concentration=setup_settings["ion_concentration"],
        box_length=setup_settings["box_length"],
        box_type=setup_settings["box_type"],
    shell:
        """
        echo "Running setup for {wildcards.ligand}"
        python {params.setup_script} \
            --ligand-file {input.ligand} \
            --protein-files {protein_files} \
            --output-directory {params.output_directory} \
            --ligand-forcefield {params.ligand_forcefield} \
            --water-model {params.water_model} \
            --ion-conc {params.ion_concentration} \
            --box-length {params.box_length} \
            --box-type {params.box_type} \
            2>&1 | tee {log.setup_log}
        """
