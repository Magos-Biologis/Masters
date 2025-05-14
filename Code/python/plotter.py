import argparse
import os
import re

import numpy as np
import pandas as pd
import plottingstuff as ps
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
fpe_env = str(os.getenv("FPE_FIGURE_ENV"))
data_env = str(os.getenv("THESIS_DATA_PATH"))

ode_dir = os.path.join(fpe_env, "ode")
fpe_dir = os.path.join(fpe_env, "five_var")

# file_dir = ode_dir

## Parser
parser = argparse.ArgumentParser(
    prog="Plotter",
    description="A python3 script that plots the data generated elsewhere",
)

parser.add_argument(
    "model",
    # nargs=1,
    help="State which model to run",
    choices=[
        "2S",
        "2L",
        "5_2",
        "5_3",
        "ode_2",
        "ode_2_2",
        "ode_3",
        "ode_3_2",
        "ode_3_3",
        "ode_5_2",
        "ode_5_3",
    ],
    type=str,
)

parser.add_argument(
    "-st",
    "--steps",
    dest="steps",
    help="Number of Steps",
    type=int,
    default=int(1e5),
)
parser.add_argument(
    "-si",
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

parser.add_argument("-ks", "--parameters", dest="k", help="Test", type=float, default=1)


args = parser.parse_args()

# print(args.k1)

## Compiling the defaults and the choice of parameters


model: str = args.model
model_count = args.index


re_pattern: re.Pattern = re.compile(
    r"(?P<datasource>[^M]*)"
    r"M(?P<model>[^P]*)"  # Model specifier
    r"P(?P<metadata>[^S]*)"  # Parameters and such
    r"S(?P<steps>[^I]*)"  # Step count
    r"I(?P<initcond>\[[ 0-9]*\])C"  # Initial conditions
    r"T(?P<t>\d+)"
    r"\.(?P<filetype>\w*)$"
)


file_names = os.listdir(data_env)

records = []
for fn in file_names:
    ## So there are two methods, search, and match
    ## It seems that search is what works for what I want,
    ## match is more for if you are decomposing the entire string from the start
    matches = re_pattern.search(fn)
    if matches is not None:
        rec = matches.groupdict()
        rec["file_name"] = fn
        records.append(rec)

## So we assign it to a dataframe for easy access
raw_frame = pd.DataFrame(records)
raw_frame.sort_values(by="t", inplace=True)  ## Sorting by unix epoch time, ascending
# raw_frame["initcond"].astype()  ## Sorting by unix epoch time, ascending

## I need to update the metadata column, to be a column of dict type.
## We do this by updating the existing elements in the df
### First we seperate variable ratio metadata from non-that
find_ratios = lambda x: x[1] if len(x) > 1 else np.nan
raw_frame["ratio"] = raw_frame["metadata"].str.split("R").array.map(find_ratios)

## Using numpy's 'not' operator, ~, we can invert the truth and locate the relevant spots
## so as to remove the ratio from the other metadata
truths = ~raw_frame["ratio"].isna()
raw_frame.loc[truths, "metadata"] = (
    raw_frame.loc[truths, "metadata"].str.split("R").array.map(lambda x: x[0] + x[2])
)


### Making the remaining string in an array of elements for further processing
extract_count = lambda l: l[0].replace("num", "").replace("=", "")
remove_ends = lambda l: l[1:-1]

raw_frame.loc[:, "metadata"] = raw_frame["metadata"].str.split("_")
raw_frame["count"] = raw_frame["metadata"].map(extract_count).astype(int)
raw_frame.loc[:, "metadata"] = raw_frame["metadata"].map(remove_ends)


## From the list, we then count the number of times each model shows up,
## assigning the number of each step accordingly.
## From this, we make a multi-index for the dataframe so as to have an
## easier time choosing which file to plot
raw_frame["model_index"] = raw_frame.groupby("model").cumcount()
file_frame = raw_frame.set_index(["model", "model_index"]).sort_index()

model_frame = file_frame.loc[model]


if args.filter is not None:
    model_frame = model_frame.loc[args.filter]

if args.opt:
    print(model_frame)
    exit()


# print(file_options.tail(1).index.values)

if args.compare_plots:
    if args.index is None:
        model_choice = model_frame.index.values[-2]
    else:
        model_choice = args.index
    file_choice = model_frame.loc[model_choice : model_choice + 1].reset_index()
    data_file_1 = os.path.join(data_env, file_choice.loc[0, "file_name"])
    data_file_2 = os.path.join(data_env, file_choice.loc[1, "file_name"])
else:
    if args.index == None:
        model_choice = model_frame.index.values[-1]
    else:
        model_choice = args.index

    file_choice = model_frame.loc[model_choice]
    data_file_1 = os.path.join(data_env, str(file_choice["file_name"]))


check_ode = re.compile(r"ode").search(model)
if check_ode is None:
    is_ode = False
else:
    is_ode = True


if args.bin_count is not None:
    m: int = args.bin_count
else:
    m: int = file_choice.loc[0, "count"]


#######################
#   File name stuff   #
#######################

## The base name
file_name: str = model

## Subsequently adding stuff dependent on the resulting bool flags
if args.compare_plots:
    for item in file_choice.loc[0, "metadata"]:
        file_name += "_" + item
else:
    for item in file_choice["metadata"]:
        file_name += "_" + item

if not is_ode:
    if args.compare_plots:
        para_version = (
            file_choice.loc[0, "ratio"] if file_choice["ratio"].isna() is True else ""
        )
    else:
        para_version = file_choice["ratio"] if file_choice["ratio"].isna() is True else ""
    file_name += "_" + para_version[1:-1]


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
    init_cond_string: str = "{}".format(file_choice.loc[0, "initcond"])
    numpy_data = np.load(data_file_1)
    time, states = numpy_data["time"], numpy_data["states"]
    gillespies.plot_walks(
        ax1,
        time,
        states,
        "r",
        init_cond_string,
        args.include_starts,
    )
    gillespies.plot_hists(
        ax3,
        ax4,
        states,
        "b",
        init_cond_string,
    )

    numpy_data = np.load(data_file_2)
    time, states = numpy_data["time"], numpy_data["states"]
    gillespies.plot_walks(
        ax2,
        time,
        states,
        "b",
        init_cond_string,
        args.include_starts,
    )
    gillespies.plot_hists(
        ax3,
        ax4,
        states,
        "r",
        init_cond_string,
    )

else:
    init_cond_string: str = "{}".format(file_choice.loc["initcond"])
    numpy_data = np.load(data_file_1)
    time, states = numpy_data["time"], numpy_data["states"]
    gillespies.plot_walk(
        ax1,
        time,
        states,
        "r",
        x_name,
        init_cond_string,
    )
    gillespies.plot_walk(
        ax2,
        time,
        states,
        "b",
        y_name,
        init_cond_string,
    )
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


# ax3.vlines(anal_sol, 0, 1, **line_kwargs, label="Analytical solution for $x$")


for ax in axes:
    ax.legend(loc="upper right", fontsize=10)

for fig in figs:
    fig.tight_layout()


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

if args.show:
    show()

print("Done Plotting")

exit()
