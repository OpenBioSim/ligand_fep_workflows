from pathlib import Path


def create_final_analysis_call(wc, input):
    config_analysis = config["analysis-settings"]
    config_sim = config["production-settings"]
    script = Path("workflow/scripts/rbfe/final_analysis.py")
    output_directory = Path(f"{config['working_directory']}/analysis/{_engine}")
    args = []
    args.append(f"--num-replicas {config_sim.get('num_replicas', 1)}")
    if config_analysis.get("backend") is not None:
        args.append(f"--plotting-backend {config_analysis['backend']}")
    if config_analysis.get("experimental-results") is not None:
        args.append(f"--experimental-results {config_analysis['experimental-results']}")
    if config_analysis.get("experimental-units") is not None:
        args.append(f"--experimental-units {config_analysis['experimental-units']}")
    if config_analysis.get("temperature") is not None:
        args.append(f"--temperature {config_analysis['temperature']}")

    extra_args = " ".join(args)
    return f"""
    python {script} --input-files {input.files} --output-directory {output_directory} {extra_args}
    """

rule analyse_single_target:
    priority: 10
    input:
        done = Path(f"{config['working_directory']}/production/{_engine}/{{ligand1}}~{{ligand2}}/{{leg}}_{{replica_number}}/.done")
    output:
        done = Path(f"{config['working_directory']}/analysis/{_engine}/{{ligand1}}~{{ligand2}}/{{leg}}_{{replica_number}}/pmf.csv")
    threads: config["simulation_threads"]
    resources:
        gpu=0
    params:
        input_directory = lambda wildcards: Path(f"{config['working_directory']}/production/{_engine}/{wildcards.ligand1}~{wildcards.ligand2}/{wildcards.leg}_{wildcards.replica_number}"),
        output_directory = lambda wildcards: Path(f"{config['working_directory']}/analysis/{_engine}/{wildcards.ligand1}~{wildcards.ligand2}/{wildcards.leg}_{wildcards.replica_number}"),
        script = Path("workflow/scripts/rbfe/analyse_single_leg.py")
    shell:
        "python {params.script} --input-directory {params.input_directory} --output-directory {params.output_directory} --plot-overlap-matrix --plot-pmf"

rule collate_rbfe_analysis:
    priority: 20
    resources:
        gpu=0
    input:
        files = read_ligand_pairs
    output:
        collated = Path(f"{config['working_directory']}/analysis/{_engine}/final_simulation_results.csv")
    threads: config["simulation_threads"]
    run:
        python_cmd = create_final_analysis_call(wildcards, input)
        shell(python_cmd)
