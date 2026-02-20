import argparse
from pathlib import Path

import BioSimSpace as BSS
import matplotlib.colors as _colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_pmfs(results_dataframe: pd.DataFrame, output_dir: Path) -> None:

    results_dataframe.plot(
        x="lambda",
        y="free_energy (kcal/mol)",
        yerr="error (kcal/mol)",
        kind="line",
        marker="o",
        ax=plt.gca(),
    )
    plt.xlabel("lambda")
    plt.ylabel("PMF (kcal/mol)")
    plt.legend(fancybox=False, edgecolor="black")
    plt.tight_layout()
    plt.savefig(output_dir / "pmf.png", dpi=300)
    plt.clf()


def plotOverlapMatrix(
    overlap: dict,
    continuous_cbar: bool = False,
    color_bar_cutoffs: list = [0.03, 0.1, 0.3],
    save_directory: Path = None,
    name: str = "overlap_matrix",
) -> None:
    """
    Plot the overlap matrix from a free-energy perturbation analysis.

    Parameters
    ----------

    overlap : List of List of float, or 2D numpy array of float
        The overlap matrix.

    continuous_cbar : bool, optional, default=False
        If True, use a continuous colour bar. Otherwise, use a discrete
        set of values defined by the 'color_bar_cutoffs' argument to
        assign a colour to each element in the matrix.

    color_bar_cutoffs : List of float, optional, default=[0.03, 0.1, 0.3]
        The cutoffs to use when assigning a colour to each element in the
        matrix. This is used for both the continuous and discrete color bars.
        Can not contain more than 3 elements.
    """

    # Validate the input
    if not isinstance(overlap, (list, tuple, np.ndarray)):
        raise TypeError(
            "The 'overlap' matrix must be a list of list types, or a numpy array!"
        )

    # Try converting to a NumPy array.
    try:
        overlap = np.array(overlap)
    except:
        raise TypeError(
            "'overlap' must be of type 'np.matrix',  'np.ndarray', or a list of lists."
        )

    # Store the number of rows.
    num_rows = len(overlap)

    # Check the data in each row.
    for row in overlap:
        if not isinstance(row, (list, tuple, np.ndarray)):
            raise TypeError("The 'overlap' matrix must be a list of list types!")
        if len(row) != num_rows:
            raise ValueError("The 'overlap' matrix must be square!")
        if not all(isinstance(x, float) for x in row):
            raise TypeError("The 'overlap' matrix must contain 'float' types!")

    # Check the colour bar options
    if not isinstance(continuous_cbar, bool):
        raise TypeError("The 'continuous_cbar' option must be a boolean!")
    if not isinstance(color_bar_cutoffs, (list, tuple, np.ndarray)):
        raise TypeError(
            "The 'color_bar_cutoffs' option must be a list of floats "
            " or a numpy array when 'continuous_cbar' is False!"
        )
    if not all(isinstance(x, float) for x in color_bar_cutoffs):
        raise TypeError("The 'color_bar_cutoffs' option must be a list of floats!")
    if len(color_bar_cutoffs) > 3:
        raise ValueError(
            "The 'color_bar_cutoffs' option must contain no more than 3 elements!"
        )

    # Add 0 and 1 to the colour bar cutoffs.
    if color_bar_cutoffs is not None:
        color_bounds = [0] + color_bar_cutoffs + [1]

    # Tuple of colours and associated font colours.
    # The last and first colours are for the top and bottom of the scale
    # for the continuous colour bar, but are ignored for the discrete bar.
    all_colors = (
        ("#FBE8EB", "black"),  # Lighter pink
        ("#FFD3E0", "black"),
        ("#88CCEE", "black"),
        ("#78C592", "black"),
        ("#117733", "white"),
        ("#004D00", "white"),
    )  # Darker green

    # Set the colour map.
    if continuous_cbar:
        # Create a color map using the extended palette and positions
        box_colors = [all_colors[i][0] for i in range(len(color_bounds) + 1)]
        cmap = _colors.LinearSegmentedColormap.from_list(
            "CustomMap", list(zip(color_bounds, box_colors))
        )

        # Normalise the same way each time so that plots are always comparable.
        norm = _colors.Normalize(vmin=0, vmax=1)
    else:
        # Throw away the first and last colours.
        box_colors = [colors[0] for colors in all_colors[1:-1]]
        cmap = _colors.ListedColormap(
            [box_colors[i] for i in range(len(color_bounds) - 1)]
        )
        norm = _colors.BoundaryNorm(color_bounds, cmap.N)

    # Create the figure and axis. Use a default size for fewer than 16 windows,
    # otherwise scale the figure size to the number of windows.
    if num_rows < 16:
        fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    else:
        fig, ax = plt.subplots(figsize=(num_rows / 2, num_rows / 2), dpi=300)

    # Create the heatmap. Separate the cells with white lines.
    im = ax.imshow(overlap, cmap=cmap, norm=norm)
    for i in range(num_rows - 1):
        for j in range(num_rows - 1):
            # Make sure these are on the edges of the cells.
            ax.axhline(i + 0.5, color="white", linewidth=0.5)
            ax.axvline(j + 0.5, color="white", linewidth=0.5)

    # Label each cell with the overlap value.
    for i in range(num_rows):
        for j in range(num_rows):
            # Get the text colour based on the overlap value.
            overlap_val = overlap[i][j]
            # Get the index of first color bound greater than the overlap value.
            for idx, bound in enumerate(color_bounds):
                if bound > overlap_val:
                    break
            text_color = all_colors[1:-1][idx - 1][1]
            ax.text(
                j,
                i,
                "{:.2f}".format(overlap[i][j]),
                ha="center",
                va="center",
                fontsize=10,
                color=text_color,
            )

    # Create a colorbar. Reduce the height of the colorbar to match the figure and remove the border.
    if continuous_cbar:
        cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap, norm=norm, shrink=0.7)
    else:
        cbar = ax.figure.colorbar(
            im,
            ax=ax,
            cmap=cmap,
            norm=norm,
            boundaries=color_bounds,
            ticks=color_bounds,
            shrink=0.7,
        )
    cbar.outline.set_visible(False)

    # Set the axis labels.
    # Set the x axis at the top of the plot.
    plt.xlabel(r"$\lambda$ Index")
    ax.xaxis.set_label_position("top")
    plt.ylabel(r"$\lambda$ Index")

    ticks = [x for x in range(0, num_rows)]

    # Set ticks every lambda window.
    plt.xticks(ticks)
    ax.xaxis.tick_top()
    plt.yticks(ticks)

    # Remove the borders.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Create a tight layout to trim whitespace.
    fig.tight_layout()
    if not save_directory:
        return plt.show()
    else:
        plt.savefig(save_directory / f"{name}.png", dpi=300)
        plt.clf()
        return None


parser = argparse.ArgumentParser(description="Analysis for a single leg.")
parser.add_argument(
    "--input-directory",
    help="Input directory containing the simulation data.",
    required=True,
)
parser.add_argument(
    "--output-directory",
    help="Output directory to save the analysis results.",
    required=True,
)
parser.add_argument(
    "--plot-overlap-matrix", action="store_true", help="Save the overlap matrix."
)
parser.add_argument("--plot-pmf", action="store_true", help="Save the PMF plot.")
args = parser.parse_args()

input_dir = Path(args.input_directory)
output_dir = Path(args.output_directory)
# fist make sure the output directory exists
output_dir.mkdir(parents=True, exist_ok=True)
try:
    pmf, overlap = BSS.FreeEnergy.Relative.analyse(str(input_dir))
    # pmf is a list of 3-element tuples (lambda, free energy, error)
    # now convert to a dataframe and write
    lambda_values = [x[0] for x in pmf]
    free_energies = [x[1].value() for x in pmf]
    errors = [x[2].value() for x in pmf]

    df = pd.DataFrame(
        {
            "lambda": lambda_values,
            "free_energy (kcal/mol)": free_energies,
            "error (kcal/mol)": errors,
        }
    )


except Exception as e:
    print(f"Could not analyse the input directory {input_dir}: {e}")
    # now we want to write blank dataframes and overlaps for later logging
    df = pd.DataFrame(
        {
            "lambda": [],
            "free_energy (kcal/mol)": [],
            "error (kcal/mol)": [],
        }
    )
    output_file = output_dir / "pmf.csv"

    overlap = np.array([])
    overlap_file = output_dir / "overlap.npy"

output_file = output_dir / "pmf.csv"
df.to_csv(output_file, index=False)
if args.plot_overlap_matrix:
    try:
        plotOverlapMatrix(overlap, save_directory=output_dir, name="overlap_matrix")
    except Exception as e:
        print(f"Could not plot overlap matrix: {e}")
if args.plot_pmf:
    try:
        plot_pmfs(df, output_dir)
    except Exception as e:
        print(f"Could not plot PMF: {e}")

# overlap matrix is a numpy array, so we can save it directly
overlap_file = output_dir / "overlap.npy"
np.save(overlap_file, overlap)

print(f"Analysis completed. Results saved to {output_dir}.")
