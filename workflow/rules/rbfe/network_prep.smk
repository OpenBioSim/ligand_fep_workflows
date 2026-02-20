# Side note: checkpoints are quite painful to use,
# they;re hard to debug and can lead to unexpected behavior,
# unfortunately, they are necessary for this workflow as we don't know what the
# file structure will look like until the network is prepared.

from pathlib import Path

ligand_folder = Path(config["ligands"]).parent

network_settings = config["network_settings"]

checkpoint prep_network:
    input:
        folder = ligand_folder,
    output:
        network_out = Path(f"{config['working_directory']}/network") / "network.dat"
    log:
        network_log = Path(f"{config['working_directory']}/logs") / "network_prep.log",
    params:
        network_script = "workflow/scripts/rbfe/prepare_network.py",
        lambda_windows = network_settings["lambda_windows"],
        lomap_threshold = network_settings["lomap_threshold"],
        diff_lambda_windows = network_settings["diff_lambda_windows"],
        output_dir = Path(f"{config['working_directory']}/network")
    shell:
        """
        echo "Preparing network from {input.folder}"  # Debug info
        python {params.network_script} \
            --ligand-folder {input.folder} \
            --lambda-windows {params.lambda_windows} \
            --lomap-threshold {params.lomap_threshold} \
            --diff-lambda-windows {params.diff_lambda_windows} \
            --output-directory {params.output_dir}
        """