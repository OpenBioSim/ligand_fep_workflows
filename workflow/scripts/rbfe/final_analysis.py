import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def get_info_from_directory(path: str | Path) -> tuple:
    """
    Extracts name1, name2, leg, and replica number from paths like:
    analysis/ejm42~ejm46/bound_0/pmf.csv
    """
    path = Path(path)

    # Get the two parent directories: "ejm42~ejm46" and "bound_0"
    pair_dir = path.parents[1].name  # e.g., 'ejm42~ejm46'
    leg_dir = path.parents[0].name  # e.g., 'bound_0'

    # Extract name1 and name2
    m_pair = re.match(r"([^~]+)~([^~]+)", pair_dir)
    if not m_pair:
        raise ValueError(f"Could not parse name pair from '{pair_dir}'")
    name1, name2 = m_pair.groups()

    # Extract leg and replica number
    m_leg = re.match(r"([a-zA-Z0-9]+)_(\d+)", leg_dir)
    if not m_leg:
        raise ValueError(f"Could not parse leg/replica from '{leg_dir}'")
    leg, replica = m_leg.groups()

    return name1, name2, leg, int(replica)


def _preprocess_experimental_data(
    experimental_values: tuple, exp_units: str, temperature: float = 300
) -> tuple:
    """
    Preprocess the experimental data by normalizing the values based on the temperature.
    """
    if exp_units == "kcal/mol":
        exp_kcal = float(experimental_values[0])
        err_kcal = float(experimental_values[1])

    elif exp_units == "Ki_uM":
        # Normalize the experimental values based on the temperature
        kB = 0.0019872041  # kcal/(mol*K)
        exp_val = float(experimental_values[0])
        exp_kcal = kB * temperature * np.log(exp_val / (10**9))
        error = float(experimental_values[1])
        exp_upper = exp_val + error
        exp_lower = exp_val - error
        exp_upper_kcal = kB * temperature * np.log(exp_upper / (10**9))
        exp_lower_kcal = kB * temperature * np.log(exp_lower / (10**9))
        err_kcal = abs(exp_upper_kcal - exp_lower_kcal) / 2.0

    elif exp_units == "kJ/mol":
        exp_kcal = float(experimental_values[0]) * 0.239006
        err_kcal = float(experimental_values[1]) * 0.239006

    else:
        raise ValueError(f"Unsupported experimental units: {exp_units}")

    return exp_kcal, err_kcal


def plot_results(
    df_simulation: pd.DataFrame,
    df_exp: pd.DataFrame,
    output_directory: Path,
    temperature: float = 300,
    exp_units: str = "kcal/mol",
) -> None:
    """
    Plots the results of the analysis.
    For plotting_backend = native
    """

    # find unique ligand1-ligand2 pairs in case any are repeated
    unique_pairs_df_sim = df_simulation.drop_duplicates(subset=["ligand1", "ligand2"])

    # experimental data is given in the form of dgs.
    # convert the experimental data to a dict, with the key being the value of `ligand` and the value being a tuple of `(`value`,`error`)`
    exp_data_dict = {
        row["ligand"]: _preprocess_experimental_data(
            (row["value"], row["error"]), exp_units=exp_units, temperature=temperature
        )
        for _, row in df_exp.iterrows()
    }
    unique_pairs_df_sim["DDG_exp"] = unique_pairs_df_sim.apply(
        lambda row: exp_data_dict[row["ligand2"]][0] - exp_data_dict[row["ligand1"]][0],
        axis=1,
    )
    # now the experimental error
    unique_pairs_df_sim["error_exp"] = unique_pairs_df_sim.apply(
        lambda row: np.sqrt(
            np.power(float(exp_data_dict[row["ligand2"]][1]), 2)
            + np.power(float(exp_data_dict[row["ligand1"]][1]), 2)
        ),
        axis=1,
    )
    # now plot DDG_exp on x and DDG_sim on y with yerr=error
    plt.figure(figsize=(10, 6))
    plt.errorbar(
        x=unique_pairs_df_sim["DDG_exp"],
        y=unique_pairs_df_sim["DDG"],
        yerr=unique_pairs_df_sim["error"],
        xerr=unique_pairs_df_sim["error_exp"],
        fmt="o",
        ls="none",
        lw=0.5,
        capsize=2,
        color="black",
        zorder=5,
    )
    # plot 1/2 kcal bounds:
    plt.fill_between(
        x=[-100, 100],
        y2=[-100.25, 99.75],
        y1=[-99.75, 100.25],
        lw=0,
        zorder=-10,
        alpha=0.3,
        color="grey",
    )
    # upper bound:
    plt.fill_between(
        x=[-100, 100],
        y2=[-99.5, 100.5],
        y1=[-99.75, 100.25],
        lw=0,
        zorder=-10,
        color="grey",
        alpha=0.2,
    )
    # lower bound:
    plt.fill_between(
        x=[-100, 100],
        y2=[-100.25, 99.75],
        y1=[-100.5, 99.5],
        lw=0,
        zorder=-10,
        color="grey",
        alpha=0.2,
    )

    # set xlim and ylim according to the experimental and simulation data
    max_sim_ddg = unique_pairs_df_sim["DDG"].max()
    max_exp_ddg = unique_pairs_df_sim["DDG_exp"].max()
    min_sim_ddg = unique_pairs_df_sim["DDG"].min()
    min_exp_ddg = unique_pairs_df_sim["DDG_exp"].min()
    max_overall = max(max_sim_ddg, max_exp_ddg)
    min_overall = min(min_sim_ddg, min_exp_ddg)
    plt.xlim(min_overall * 1.3, max_overall * 1.3)
    plt.ylim(min_overall * 1.3, max_overall * 1.3)
    plt.ylabel(r"Computed $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
    plt.xlabel(r"Experimental $\Delta\Delta$G$_{bind}$ / kcal$\cdot$mol$^{-1}$")
    plt.tight_layout()
    plt.savefig(output_directory / "DDGs.png")


def plot_results_cinnabar(
    df_simulation: pd.DataFrame,
    df_exp: pd.DataFrame,
    output_directory: Path,
    cinnabar_backend: str,
    temperature: float = 300,
    exp_units: str = "kcal/mol",
) -> None:
    import csv

    # remove all entries of df_exp whose `ligand` is not present in `ligand1` or `ligand2` of df_simulation
    df_exp = df_exp[
        df_exp["ligand"].isin(df_simulation["ligand1"])
        | df_exp["ligand"].isin(df_simulation["ligand2"])
    ]
    exp_data_dict = {
        row["ligand"]: _preprocess_experimental_data(
            (row["value"], row["error"]), exp_units=exp_units, temperature=temperature
        )
        for _, row in df_exp.iterrows()
    }
    # for plotting_backend = cinnabar
    # first we need to format our data and write it to an analysis file
    with open(output_directory / "cinnabar_data.csv", "w") as cinnabar_data_file:
        writer = csv.writer(cinnabar_data_file, delimiter=",")

        # first, write the experimental data
        writer.writerow(["# Experimental block"])
        writer.writerow(["# Ligand", "expt_DDG", "expt_dDDG"])
        for lig in exp_data_dict.keys():
            writer.writerow([lig, exp_data_dict[lig][0], exp_data_dict[lig][1]])

        writer.writerow([" "])
        writer.writerow(["# Calculated block"])
        writer.writerow(
            [
                "# Ligand1",
                "Ligand2",
                "calc_DDG",
                "calc_dDDG(MBAR)",
                "calc_dDDG(additional)",
            ]
        )

        # find unique ligand1-ligand2 pairs in case any are repeated
        unique_pairs_df_sim = df_simulation.drop_duplicates(
            subset=["ligand1", "ligand2"]
        )
        for _, row in unique_pairs_df_sim.iterrows():
            writer.writerow(
                [row["ligand1"], row["ligand2"], row["DDG"], row["error"], "0.0"]
            )

    if cinnabar_backend == "old":
        network = wrangle.FEMap(output_directory / "cinnabar_data.csv")
        plotting.plot_DDGs(
            network.graph,
            title="DDGs",
            filename=str(output_directory / "DDGs.png"),
            figsize=6,
        )
        plotting.plot_DGs(
            network.graph,
            title="DGs",
            filename=str(output_directory / "DGs.png"),
            figsize=6,
        )

    elif cinnabar_backend == "new":
        network = FEMap.from_csv(output_directory / "cinnabar_data.csv")
        try:
            plotting.plot_DDGs(
                network.to_legacy_graph(),
                title="DDGs",
                filename=str(output_directory / "DDGs.png"),
                figsize=6,
            )
        except Exception as e:
            print(
                f"Could not calculate DDG values for this network, cinnabar failed with the following error: {e}"
            )

        try:
            plotting.plot_DGs(
                network.to_legacy_graph(),
                title="DGs",
                filename=str(output_directory / "DGs.png"),
                figsize=6,
            )
        except Exception as e:
            print(
                f"Could not calculate DG values for this network, cinnabar failed with the following error: {e}"
            )


parser = argparse.ArgumentParser(
    description="Final analysis of free energy perturbation results."
)
parser.add_argument(
    "--input-files",
    type=Path,
    required=True,
    help="Path to the input CSV file.",
    nargs="+",
)
parser.add_argument(
    "--num-replicas",
    type=int,
    default=3,
    help="Total number of replicas for each system.",
)
parser.add_argument(
    "--experimental-results",
    type=Path,
    help="Path to the experimental results file.",
    default=None,
)
parser.add_argument(
    "--plotting-backend",
    type=str,
    choices=["native", "cinnabar"],
    default="cinnabar",
)
parser.add_argument(
    "--temperature",
    type=str,
    default="300K",
    help="Temperature in Kelvin at which production simulations were run.",
)
parser.add_argument(
    "--experimental-units",
    type=str,
    choices=["kcal/mol", "kJ/mol", "Ki_uM"],
    default="kcal/mol",
    help="Units used in the experimental results.",
)
parser.add_argument(
    "--output-directory",
    type=Path,
    required=True,
    help="Directory to save the output plots.",
)
log = "Logging for final analysis"
args = parser.parse_args()

# make the output dir if it doesn't exist
args.output_directory.mkdir(parents=True, exist_ok=True)

# convert temperature to a float
temp = float(re.match(r"[\d.]+", args.temperature).group(0))
# We need to get the ligands, leg and replica number from the input filenames
info = {"ligand1": [], "ligand2": [], "leg": [], "replica": [], "dg": []}
for input_file in args.input_files:
    df = pd.read_csv(input_file)
    # all we need from this is the value of the pmf at lambda=1

    try:
        val = df[df["lambda"] == 1.0]
        val = val["free_energy (kcal/mol)"].astype(float)
        pmf_value = float(val.iloc[0])
    except:
        log += f"Error reading pmf value from {input_file}. Possible that results are missing or simulation crashed.\n"
        continue
    name1, name2, leg, replica = get_info_from_directory(input_file)
    info["dg"].append(pmf_value)
    info["ligand1"].append(name1)
    info["ligand2"].append(name2)
    info["leg"].append(leg)
    info["replica"].append(replica)

# convert info to dataframe for easier searching and manipulation
df_info = pd.DataFrame(info)

if df_info.empty:
    print(
        "ERROR: No valid PMF data found in any input files. Cannot compute DDG values."
    )
    print(
        "This may indicate that the analysis engine does not support the simulation output format."
    )
    sys.exit(1)

# Boolean mask for rows where ligand1 > ligand2 (needs swapping)
mask = df_info["ligand1"] > df_info["ligand2"]

# Swap where needed
df_info.loc[mask, ["ligand1", "ligand2"]] = df_info.loc[
    mask, ["ligand2", "ligand1"]
].values

# Flip dg sign where swapped
df_info.loc[mask, "dg"] *= -1
# completely arbitrary, we just want the number
df_info.loc[mask, "replica"] += args.num_replicas

# Pivot so 'leg' becomes columns ('bound' and 'free')
pivoted = df_info.pivot_table(
    index=["ligand1", "ligand2", "replica"], columns="leg", values="dg"
).reset_index()

# Check that both legs are present
for required_leg in ["bound", "free"]:
    if required_leg not in pivoted.columns:
        print(f"ERROR: No valid '{required_leg}' leg data found after pivoting.")
        print(
            f"Available legs: {[c for c in pivoted.columns if c not in ['ligand1', 'ligand2', 'replica']]}"
        )
        print("Check that the analysis script can parse the simulation output format.")
        sys.exit(1)

# Calculate bound - free for each replica
pivoted["DDG"] = pivoted["bound"] - pivoted["free"]

# Now average over replicas for each ligand1/ligand2 pair
result = pivoted.groupby(["ligand1", "ligand2"])["DDG"].mean().reset_index()
# Now the error
result["error"] = (
    pivoted.groupby(["ligand1", "ligand2"])["DDG"].std().reset_index()["DDG"]
)

# save "result" to the output directory as a csv file
pd.DataFrame(result).to_csv(
    args.output_directory / "final_simulation_results.csv", index=False
)
if args.experimental_results is not None:
    # now read the experimental results, and do the same sorting by ligand1 and ligand2 and multiplying by -1 if needed
    df_exp = pd.read_csv(args.experimental_results)
    plotting_backend = args.plotting_backend
    if plotting_backend == "cinnabar":
        # Importing cinnabar - there are two possible API versions, so try both
        try:
            from cinnabar import plotting, wrangle

            plotting_backend = "cinnabar"
            cinnabar_API = "old"
        except ImportError:
            from cinnabar import FEMap, plotting

            plotting_backend = "cinnabar"
            cinnabar_API = "new"
        except:
            print("No cinnabar found, falling back to native plotting method")
            plotting_backend = "native"

    if plotting_backend == "native":
        plot_results(
            result, df_exp, args.output_directory, exp_units=args.experimental_units
        )
    elif plotting_backend == "cinnabar":
        plot_results_cinnabar(
            result,
            df_exp,
            output_directory=args.output_directory,
            cinnabar_backend=cinnabar_API,
            exp_units=args.experimental_units,
        )
