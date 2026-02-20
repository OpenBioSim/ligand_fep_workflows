from pathlib import Path

rbfe_prep_settings = config["production-settings"]

rule rbfe_prep:
    input:
        file1 = Path(f"{config['working_directory']}/preparation/final") / "{ligand1}_{leg}.bss",
        file2 = Path(f"{config['working_directory']}/preparation/final") / "{ligand2}_{leg}.bss"
    output:
        done = Path(f"{config['working_directory']}/rbfe_prepared/{{leg}}/{{ligand1}}~{{ligand2}}.bss")
    params:
        out_dir = Path(config["working_directory"]) / "rbfe_prepared",
        hmr_factor = rbfe_prep_settings["gromacs-settings"]["hmr_factor"],
    shell:
        """
        echo "Preparing FEP input for {wildcards.ligand1} and {wildcards.ligand2} in leg {wildcards.leg}"  # Debug info
        python workflow/scripts/rbfe/rbfe_prep.py \
        --ligand1_name {wildcards.ligand1} \
        --ligand2_name {wildcards.ligand2} \
        --ligand1_files {input.file1} \
        --ligand2_files {input.file2} \
        --output_location {params.out_dir} \
        --leg {wildcards.leg} \
        --hmr-factor {params.hmr_factor}
        """
