import argparse
import os
import re

import numpy as np
import pandas as pd
import plottingstuff as ps
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
fig_env = str(os.getenv("FPE_FIGURE_ENV"))
data_env = str(os.getenv("THESIS_DATA_PATH"))

ode_dir = os.path.join(fig_env, "ode")
fpe_dir = os.path.join(fig_env, "five_var")

# file_dir = ode_dir

## Parser
parser = argparse.ArgumentParser(
    prog="Plotter",
    description="A python3 script that plots the data generated elsewhere",
)

parser.add_argument(
    "source",
    help="Which source the data should be from",
    choices=[
        "markov",
        "ode",
        "phase",
        "ssa",
    ],
    type=str,
)

parser.add_argument(
    "-m",
    "--model",
    dest="model",
    help="State which model to run",
    choices=[
        "jesper",
        "2S",
        "2L",
        "5_2",
        "5_3",
        "ode_2",
        "ode_2_2",
        "ode_3",
        "ode_3_2",
        "ode_3_2_alt",
        "ode_3_3",
        "ode_5_2",
        "ode_5_3",
    ],
    type=str,
)

# parser.add_argument(
#     "-st",
#     "--steps",
#     dest="steps",
#     help="Number of Steps",
#     type=int,
#     default=int(1e5),
# )
parser.add_argument(
    "-n",
    "--size",
    dest="size",
    help="Size of System",
    type=int,
    default=int(100),
)


parser.add_argument(
    "-b",
    "--bins",
    dest="bin_count",
    help="Number of bins",
    type=int,
    default=int(100),
)
parser.add_argument(
    "-f",
    "--filter",
    dest="filter",
    help="Filter",
)
parser.add_argument(
    "-i",
    "--index",
    dest="index",
    help="Which index",
    type=int,
)
parser.add_argument(
    "-cp",
    "--compare-plots",
    dest="compare_plots",
    help="Compare Plots",
    action="store_true",
)
# type=int,
# const=1,

parser.add_argument(
    "-o",
    "--options",
    dest="opt",
    help="List options",
    action="store_true",
)
parser.add_argument(
    "-sh",
    "--show",
    dest="show",
    help="Show Plot",
    action="store_true",
)

parser.add_argument(
    "-ns",
    "--no-save",
    dest="save",
    help="Don't save plot",
    action="store_false",
)

parser.add_argument(
    "-nb",
    "--not-both",
    dest="both",
    help="Don't plot both x and y",
    action="store_false",
)

# parser.add_argument(
#     "-ns",
#     "--no-save",
#     dest="save",
#     help="Don't save plot",
# )


parser.add_argument(
    "-ic",
    "--initial-conds",
    nargs="*",
    dest="initial_conds",
    help="Initial Conditions",
    type=int,
)

parser.add_argument(
    "-is",
    "--include-starts",
    dest="include_starts",
    help="Include starts in legend",
    action="store_true",
)
# default=[99, 1, 100],

# parser.add_argument(
#     "-ks",
#     "--parameters",
#     dest="k",
#     help="Test",
#     type=float,
#     default=1,
# )


args = parser.parse_args()

# print(args.k1)

## Compiling the defaults and the choice of parameters


data_source: str = args.source
model_count: int = args.index

is_ode = False
if args.model is not None:
    model: str = args.model
    check_ode = re.compile(r"ode").search(model)
    if check_ode is not None:
        is_ode = True


ratio_pattern: re.Pattern = re.compile(r"R(?P<ratio>b.n)R")
capture_patterns: re.Pattern = re.compile(
    r"(?P<datasource>^[^M]*)"  # phase, ssa, etc.
    r"M(?P<model>[^P]*)"  # Model specifier
    r"P(?P<metadata>[^SR]*)"  # Parameters and
    r"[^S]*"  # a stupid way to ignore the ratio
    r"S(?P<steps>[^T]*)"  # Step count
    r"I(?P<initcond>[^CT]*)C"  # Initial conditions
    r"T(?P<t>\d+)"
    r"\.(?P<filetype>\w*)$"
)


file_names = os.listdir(data_env)

## So there are two methods, search, and match
## It seems that search is what works for what I want,
## match is more for if you are decomposing the entire string from the start
records = []
for fn in file_names:
    rec = dict()

    matches = capture_patterns.match(fn)
    ratio = ratio_pattern.search(fn)

    if matches is not None:
        rec.update(matches.groupdict())

    if ratio is not None:
        rec.update(ratio.groupdict())

    rec["file_name"] = fn
    records.append(rec)


## So we assign it to a dataframe for easy access
raw_frame = pd.DataFrame(records)
raw_frame.sort_values(by="t", inplace=True)  ## Sorting by unix epoch time, ascending

# print(raw_frame[["ratio", "metadata"]])


## I need to update the metadata column, to be a column of dict type.
## Isolating all the matches in the strings so as to prepare a
## list of tuple[str, float] types— which can be turned into a dictionary
capture_parameters = re.compile(r"(\w+[-m]?\d?p?)=([^_]*)_")


def parse_parameters(string):
    matches: list[tuple] = capture_parameters.findall(string)
    parameters = [(str(key), float(value)) for key, value in matches]
    return dict(parameters)


raw_frame.loc[:, "metadata"] = raw_frame["metadata"].map(parse_parameters)

no_ratio = raw_frame["ratio"].isna()
raw_frame.loc[no_ratio, "ratio"] = ""
# print(raw_frame["ratio"])


### Now we can mutate the dictionaries in memory and dynamically create columns
defined_metadata = raw_frame["metadata"].map(lambda d: len(d) > 0)
raw_frame = raw_frame[defined_metadata]
raw_frame["count"] = [dic.pop("num") for dic in raw_frame["metadata"]]
# raw_frame["count"] = raw_frame['metadata'].map(lambda d:
# raw_frame["count"].astype(np.int_)

# print(raw_frame["count"])
# exit()


## From the list, we then count the number of times each model shows up,
## assigning the number of each step accordingly.
## From which multi-index is generated for the dataframe
## so as to have an easier time choosing which file to plot
raw_frame["model_index"] = raw_frame.groupby("model").cumcount()
file_frame = raw_frame.set_index(["datasource", "model", "model_index"]).sort_index()

filtered_frame = file_frame.loc[data_source]

if args.model is not None:
    filtered_frame = filtered_frame.loc[args.model]


if args.filter is not None:
    filtered_frame = filtered_frame.loc[args.filter]

if args.opt:
    print(filtered_frame)
    exit()


if args.compare_plots:
    if args.index is None:
        model_choice = filtered_frame.index.values[-2]
    else:
        model_choice = args.index
    file_choice = filtered_frame.loc[model_choice : model_choice + 1].reset_index()
    data_file_1 = os.path.join(data_env, file_choice.loc[0, "file_name"])
    data_file_2 = os.path.join(data_env, file_choice.loc[1, "file_name"])
else:
    if args.index == None:
        model_choice = filtered_frame.index.values[-1]
    else:
        model_choice = args.index

    file_choice = filtered_frame.loc[model_choice]
    data_file_1 = os.path.join(data_env, str(file_choice["file_name"]))


if args.bin_count is not None:
    m: int = args.bin_count
else:
    m: int = file_choice.loc[0, "count"]


#######################
#   File name stuff   #
#######################

## The base name
file_name: str = data_source if args.model is None else args.model

## Subsequently adding stuff dependent on the resulting bool flags
if args.compare_plots:
    for item in file_choice.loc[0, "metadata"]:
        file_name += "_" + item
else:
    for item in file_choice["metadata"]:
        file_name += "_" + item


if not is_ode:
    if args.compare_plots:
        para_version = file_choice.loc[0, "ratio"]
    else:
        para_version = file_choice["ratio"]

    file_name += "_" + para_version


## Adding the epoch time to ensure the files don't overwrite eachother
file_name += "T"
if args.compare_plots:
    file_name += file_choice.loc[0, "t"]
    file_name += "W2"
else:
    file_name += file_choice["t"]
    file_name += "W1"

#######################

alpha = 0
beta = m

boxes = arange(alpha, beta + 1, 1, dtype=np.int_)
xlims = (alpha, beta)

curv_kwargs = {
    "linewidth": 4,
    "alpha": 0.7,
}

font_kwargs = {
    "fontsize": 12,
}

line_kwargs = {
    "linewidth": 3,
    "alpha": 0.8,
    "linestyle": "-.",
    "color": "black",
    "zorder": 11,
}
hist_kwargs = {
    "bins": boxes,
    "density": True,
    "edgecolor": "black",
    "alpha": 0.6,
    "align": "mid",
    # "normed": True,
}


walk_kwargs = {
    "linewidth": 1,
    "alpha": 0.7,
}


if is_ode:
    x_name = "c_1"
    y_name = "c_2"
else:
    x_name = "x"
    y_name = "y"


if data_source == "ssa":
    gillespies = ps.gillespie_plotters(
        curv_kwargs=curv_kwargs,
        font_kwargs=font_kwargs,
        hist_kwargs=hist_kwargs,
        line_kwargs=line_kwargs,
        walk_kwargs=walk_kwargs,
        x_name=x_name,
        y_name=y_name,
    )

    fig1 = figure(figsize=(5, 2.5))
    fig2 = figure(figsize=(5, 2.5))
    fig3 = figure(figsize=(5, 5))

    ax1 = fig1.add_subplot()
    ax2 = fig2.add_subplot()
    ax3 = fig3.add_subplot()

    figs = [fig1, fig2, fig3]
    axes = [ax1, ax2, ax3]

    if args.compare_plots:
        fig4 = figure(figsize=(5, 5))
        ax4 = fig4.add_subplot()

        figs.append(fig4)
        axes.append(ax4)

    if args.compare_plots:
        colors = ["b", "r"]
        files = [(ax1, data_file_1), (ax2, data_file_2)]
        for i, (ax, file) in enumerate(files):
            init_cond_string: str = "{}".format(file_choice.loc[i, "initcond"])
            numpy_data = np.load(file)
            time, states = numpy_data["time"], numpy_data["states"]
            gillespies.plot_walks(
                ax,
                time,
                states,
                colors[i],
                init_cond_string,
                args.include_starts,
                top=m,
            )
            gillespies.plot_hists(
                ax3,
                ax4,
                states,
                colors[i],
                init_cond_string,
            )

    else:
        init_cond_string: str = "{}".format(file_choice.loc["initcond"])
        numpy_data = np.load(data_file_1)
        time, states = numpy_data["time"], numpy_data["states"]
        colors = [(ax1, "r", x_name), (ax2, "b", y_name)]
        for i, (ax, color, name) in enumerate(colors):
            gillespies.plot_walk(
                ax,
                time,
                states[i, :],
                color,
                name,
                init_cond_string,
            )
        # gillespies.plot_walk(
        #     ax2,
        #     time,
        #     states,
        #     "b",
        #     y_name,
        #     init_cond_string,
        # )
        gillespies.plot_hist(
            ax3,
            states[0, :],
            "r",
            init_cond_string,
        )

    # file_name = "multiplot"

    # print(file_name)
    # exit()
    # ax1.hlines(anal_sol, 0, time_array_1[-1], colors=["r", "g"])
    # ax2.hlines(anal_sol, 0, time_array_2[-1], colors=["b", "g"])
    #
    # ax3.vlines(anal_sol[0], 0, 1 / m)
    # ax4.vlines(anal_sol[1], 0, 1 / m)

    if is_ode:
        ax3.set_title("Cell Population Densities")
        if args.compare_plots:
            ax4.set_title("Cell Population Densities")

    else:
        # para_match = re.compile(r"R(?P<ratio>b[<=>]n)R").search(file_choice["metadata"])
        # para_version = para_match.groupdict()["ratio"]  ## This is scuffed, but it works

        ax3.set_title("Distribution of Gene Copy Number\n" + para_version)
        if args.compare_plots:
            ax4.set_title("Distribution of Gene Copy Number\n" + para_version)

    # print(file_choice.loc[0, "t"])
    # exit()
    # ax3.hist(
    #     gillespie_results_2[0, :],
    #     **hist_kwargs,
    #     label="Start Condition [0,1]",
    #     color="b",
    # )plo

    ax3.set_xlim(xlims)
    ax3.set_ylim(bottom=0)

    ax3.set_ylabel("Fraction of States", fontsize=12)

    if args.compare_plots:
        ax3.set_xlabel(f"Distribution of the Count of ${x_name}$", fontsize=12)
        ax4.set_xlabel(f"Distribution of the Count of ${y_name}$", fontsize=12)
        ax4.set_xlim(xlims)
        ax4.set_ylim(bottom=0)
    else:
        ax3.set_xlabel("Distribution of cell states")

    # ax3.set_ylabel("Density", fontsize=12)
    # ax4.set_ylabel("Density", fontsize=12)

    for ax in axes:
        ax.legend(loc="upper right", fontsize=10)

    for fig in figs:
        fig.tight_layout()

# ax3.vlines(anal_sol, 0, 1, **line_kwargs, label="Analytical solution for $x$")

elif data_source == "phase":
    phaseies = ps.ode_plotters(
        curv_kwargs=curv_kwargs,
        font_kwargs=font_kwargs,
        hist_kwargs=hist_kwargs,
        line_kwargs=line_kwargs,
        walk_kwargs=walk_kwargs,
        x_name=x_name,
        y_name=y_name,
    )

    numpy_data = np.load(data_file_1)
    c1, c2 = numpy_data["c1"], numpy_data["c2"]
    dU, dV = numpy_data["dU"], numpy_data["dV"]

    speed = np.sqrt(dU**2 + dV**2)
    lw = np.log(speed / speed.max()) / 5
    stream_kwargs = {
        "density": 1.7,
        "linewidth": lw,
        "arrowstyle": "->",
        "color": "black",
        # "color":lw,
        # "norm":norm,
        # "cmap":"gist_heat_r",
        # "arrowsize":0,
        # "density":0.7,
        # "broken_streamlines":False,
    }

    fig, ax = plt.subplots(figsize=(6, 6))

    phaseies.plot_phase_space(ax, c1, c2, dU, dV, **stream_kwargs)
    ax.set_ylim(bottom=0)


## Taken from https://www.geeksforgeeks.org/save-multiple-matplotlib-figures-in-single-pdf-file-using-python/
def save_image(filename):
    # PdfPages is a wrapper around pdf
    # file so there is no clash and create
    # files with no error.
    filename += ".pdf"
    p = PdfPages(filename)

    # get_fignums Return list of existing
    # figure numbers
    fig_nums = plt.get_fignums()
    figs = [plt.figure(n) for n in fig_nums]

    # iterating over the numbers in list
    for fig in figs:
        # and saving the files
        fig.savefig(p, format="pdf")

    # close the object
    p.close()


if args.save:
    if is_ode:
        file_path = os.path.join(ode_dir, file_name)
    else:
        file_path = os.path.join(fpe_dir, file_name)

    save_image(file_path)

latest_file = os.path.join(fig_env, "latest_plot")
save_image(latest_file)

if args.show:
    show()

print("Done Plotting")

exit()
