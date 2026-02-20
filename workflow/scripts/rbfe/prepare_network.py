import argparse
import csv
import glob
from pathlib import Path

import BioSimSpace as BSS
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

parser = argparse.ArgumentParser(
    description="Network preparation script. This sets up the FEP network, as well as defining some settings for production runs."
)
parser.add_argument(
    "--ligand-folder",
    type=str,
    help="Path to the folder containing the ligand files (assumed to be sdf format)",
    required=True,
)
parser.add_argument(
    "--lambda-windows",
    type=int,
    help="Number of lambda windows to use for perturbations with LOMAP score above the threshold",
    default=11,
    choices=range(3, 21, 1),
    metavar="[3 - 20]",
)
parser.add_argument(
    "--lomap-threshold",
    type=float,
    help="LOMAP score threshold for perturbation selection",
    default=0.5,
    choices=list(map(lambda x: x / 100, range(10, 100, 1))),
    metavar="[0.1 - 0.99]",
)
parser.add_argument(
    "--diff-lambda-windows",
    type=int,
    help="Number of lambda windows to use for difficult perturbations - those with LOMAP score below the threshold",
    default=17,
    choices=range(3, 21, 1),
    metavar="[3 - 20]",
)
parser.add_argument(
    "--output-directory",
    type=str,
    help="Path to the output directory for the network information",
    default="network",
)
args = parser.parse_args()

# first create the output directory for the network information
output_dir = Path(args.output_directory)
output_dir.mkdir(exist_ok=True, parents=True)

ligand_folder = sorted(glob.glob(f"{Path(args.ligand_folder).resolve()}/*.sdf"))

ligands = []
ligand_names = []

for file in ligand_folder:
    print(file)
    ligand = BSS.IO.readMolecules(file)[0]
    ligands.append(ligand)
    ligand_name = Path(file).stem
    ligand_names.append(ligand_name)

transformations, lomap_scores = BSS.Align.generateNetwork(
    ligands,
    plot_network=True,
    names=ligand_names,
    work_dir=str(output_dir / "visualise_network"),
)

# generate information for use downstream
pert_network_dict = {}
transformations_named = [
    (ligand_names[transf[0]], ligand_names[transf[1]]) for transf in transformations
]
for score, transf in sorted(zip(lomap_scores, transformations_named)):
    pert_network_dict[transf] = score
    print(transf, score)

graph = nx.Graph()

# Loop over the nligands and add as nodes to the graph.
for lig in ligand_names:
    graph.add_node(lig, label=lig, labelloc="t")

# Loop over the edges in the dictionary and add to the graph.
for edge in pert_network_dict:
    graph.add_edge(edge[0], edge[1], label=(pert_network_dict[edge]))

# Plot the networkX graph.
pos = nx.kamada_kawai_layout(graph)
plt.figure(figsize=(14, 14), dpi=150)
nx.draw(
    graph,
    pos,
    edge_color="black",
    width=1,
    linewidths=1,
    node_size=1500,
    node_color="skyblue",
    font_size=12,
    labels={node: node for node in graph.nodes()},
)

nx.draw_networkx_edge_labels(
    graph, pos, edge_labels=pert_network_dict, font_color="purple", font_size=10
)

plt.savefig(output_dir / "network.png", bbox_inches="tight", dpi=300)

# now write the network to a dat file
with open(output_dir / "ligands.dat", "w") as ligands_file:
    writer = csv.writer(ligands_file)
    for lig in ligand_names:
        writer.writerow([lig])

# write perts file. Base the lambda schedule on the file generated in the previous cell.
np.set_printoptions(formatter={"float": "{: .4f}".format})

with open(output_dir / "network.dat", "w") as network_file:
    writer = csv.writer(network_file, delimiter=" ")

    for pert, lomap_score in pert_network_dict.items():
        # based on the provided (at top of notebook) lambda allocations and LOMAP threshold, decide allocation.
        if (lomap_score is None) or (lomap_score < args.lomap_threshold):
            num_lambda = args.diff_lambda_windows
        else:
            num_lambda = args.lambda_windows

        # given the number of allocated lambda windows, generate an array for parsing downstream.
        lam_array_np = np.around(np.linspace(0, 1, int(num_lambda)), decimals=5)

        # make the array into a format readable by bash.
        lam_array = (
            str(lam_array_np)
            .replace("[ ", "")
            .replace("]", "")
            .replace("  ", ",")
            .replace("\n", "")
        )

        # write out both directions for this perturbation.
        writer.writerow(
            [
                pert[0],
                pert[1],
                len(lam_array_np),
                lam_array,
            ]
        )
        writer.writerow(
            [
                pert[1],
                pert[0],
                len(lam_array_np),
                lam_array,
            ]
        )
