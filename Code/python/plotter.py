#!../.venv/bin/python3
import argparse
import os
import re

import numpy as np
import pandas as pd
import plottingstuff as ps
from matplotlib.backends.backend_pdf import PdfPages
from pylab import *

# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
fig_env = str(os.getenv("THESIS_FIGURE_PATH"))
fpe_env = str(os.getenv("FPE_FIGURE_ENV"))
data_env = str(os.getenv("THESIS_DATA_PATH"))

ode_dir = os.path.join(fpe_env, "ode")
fpe_dir = os.path.join(fpe_env, "five_var")

# file_dir = ode_dir

## Parser
parser = argparse.ArgumentParser(
    prog="Plotter",
    description="A python3 script that plots the data generated elsewhere",
    prefix_chars="-+",
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
    "model",
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


bool_args = parser.add_argument_group(
    "Boolean arguments",
    "Arguments that exist as a boolean flag",
)
bool_args.add_argument(
    "-o",
    "--options",
    dest="opts",
    help="List options",
    action="store_true",
)
bool_args.add_argument(
    "-sh",
    "--show",
    dest="show",
    help="Show Plot",
    action="store_true",
)
bool_args.add_argument(
    "-ns",
    "--no-save",
    dest="save",
    help="Don't save plot",
    action="store_false",
)
bool_args.add_argument(
    "-nb",
    "--not-both",
    dest="both",
    help="Don't plot both x and y",
    action="store_false",
)
bool_args.add_argument(
    "-cp",
    "--compare-plots",
    dest="compare_plots",
    help="Compare Plots",
    action="store_true",
)
bool_args.add_argument(
    "-is",
    "--include-starts",
    dest="include_starts",
    help="Include starts in legend",
    action="store_true",
)
bool_args.add_argument(
    "-fp",
    "---fixed-points",
    dest="plot_fixedpoints",
    help="Include the fixed points in the histogram and walks",
    action="store_true",
)


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
    "-i",
    "--index",
    dest="index",
    help="Which index",
    type=int,
)


## Taken from my parameter flag in the sim script
filter_captures = re.compile(r"(\w[^=]*)=\s?([^ ,]*)")
kwarg_string_parser = re.compile(r"(\w[^=]*)=\s?([^ ,]*)")


def parse_filters(string):
    matches: list[tuple] = filter_captures.findall(string)
    parameters = [
        (str(key).replace("-", "_").replace("'", "p").replace(" ", ""), float(value))
        for key, value in matches
    ]
    return parameters


def parse_kwarg_string(string):
    matches: list[tuple] = kwarg_string_parser.findall(string)
    parameters = [(str(key).replace(" ", ""), value) for key, value in matches]
    return dict(parameters)


parser.add_argument(
    "-f",
    "--filter",
    dest="filter",
    help='Provide a list of parameters to filter by, "par = __ "',
    type=parse_filters,
)


kwarg_parse = parser.add_argument_group(
    "Kwarg",
    "Arguments for modifying the kwargs for various aspects of the plots",
    prefix_chars="+",
    argument_default=dict(),
)


kwarg_parse.add_argument("+c", "++curve", dest="kwarg_curve", type=parse_kwarg_string)
kwarg_parse.add_argument("+f", "++font", dest="kwarg_font", type=parse_kwarg_string)
kwarg_parse.add_argument("+h", "++hist", dest="kwarg_hist", type=parse_kwarg_string)
kwarg_parse.add_argument("+l", "++line", dest="kwarg_line", type=parse_kwarg_string)
kwarg_parse.add_argument("+w", "++walk", dest="kwarg_walk", type=parse_kwarg_string)
kwarg_parse.add_argument("+s", "++stream", dest="kwarg_stream", type=parse_kwarg_string)
kwarg_parse.add_argument("+p", "++plot", dest="kwarg_plot", type=parse_kwarg_string)


args = parser.parse_args()


# print(args.filter)
# exit()

## Compiling the defaults and the choice of parameters


data_source: str = args.source
model: str = args.model
model_count: int = args.index

is_ode = False
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


## So there are two methods, search, and match
## It seems that search is what works for what I want,
## match is more for if you are decomposing the entire string from the start
records = []
file_names = os.listdir(data_env)
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


## I need to update the metadata column, to be a column of dict type.
## Isolating all the matches in the strings so as to prepare a
## list of tuple[str, float] types— which can be turned into a dictionary
## For that I just use the same technique as for the filter flag type
metadata_captures = re.compile(r"(\w[^=]*)=([^_]*)_")


def parse_metadata(string):
    matches: list[tuple] = metadata_captures.findall(string)
    parameters = [(str(key).replace("-", "_"), float(value)) for key, value in matches]
    return dict(parameters)


raw_frame.loc[:, "metadata"] = raw_frame["metadata"].map(parse_metadata)


# print(raw_frame.loc[:, "metadata"])
# exit()

no_ratio = raw_frame["ratio"].isna()
raw_frame.loc[no_ratio, "ratio"] = ""


### Now we can mutate the dictionaries in memory and dynamically create columns
defined_metadata = raw_frame["metadata"].map(lambda d: len(d) > 0)
raw_frame = raw_frame.loc[defined_metadata]

ssa_index = raw_frame.loc[raw_frame["datasource"] == "ssa"].index
raw_frame.loc[ssa_index, "count"] = [
    dic.pop("num") for dic in raw_frame.loc[ssa_index, "metadata"]
]


## From the list, we then count the number of times each model shows up,
## assigning the number of each step accordingly.
## From which multi-index is generated for the dataframe
## so as to have an easier time choosing which file to plot
raw_frame["model_index"] = raw_frame.groupby("model").cumcount()
raw_frame = raw_frame.set_index(["datasource", "model", "model_index"]).sort_index()

parameters = raw_frame["metadata"].apply(pd.Series)
file_frame = pd.concat([raw_frame, parameters], axis=1)

filtered_frame = file_frame.loc[data_source]

if args.model is not None:
    filtered_frame = filtered_frame.loc[args.model]


if args.filter is not None:
    for str, val in args.filter:
        filtered_frame = filtered_frame.loc[filtered_frame[str] == val]

if args.opts:
    print("{}".format(args.model))
    print(filtered_frame[["initcond", "count", "steps", "ratio", "t"]])
    exit()


if args.compare_plots:
    if args.index is None:
        model_choice = filtered_frame.index.values[-2]
    else:
        model_choice = args.index
    file_choice = filtered_frame.loc[model_choice : model_choice + 1].reset_index()
else:
    if args.index == None:
        model_choice = filtered_frame.index.values[-1]
    else:
        model_choice = args.index

    file_choice = filtered_frame.loc[model_choice]


if args.bin_count is not None:
    m: int = args.bin_count
else:
    m: int = file_choice.loc[0, "count"]


#######################
#   File name stuff   #
#######################

## The base name
file_name: str = "{}:{}:".format(data_source, model).replace("_", "-")

para_version: str = ""
if not is_ode:
    if args.compare_plots:
        para_version = r"${}$".format(file_choice.loc[0, "ratio"])
    else:
        para_version = r"${}$".format(file_choice["ratio"])

    # if para_version != "":
    if model == "5_2":
        file_name += "{}:".format(para_version).replace("$", "")


## Adding the epoch time to ensure the files don't overwrite eachother
if args.compare_plots:
    epoch = file_choice.loc[0, "t"]
    data_set_quantity = 2
else:
    epoch = file_choice["t"]
    data_set_quantity = 1

file_name += "I{}:".format(model_choice)
file_name += "C{}:".format(data_set_quantity)
file_name += "T{}".format(epoch)

#######################

alpha = 0
beta = m

boxes = arange(alpha, beta + 1, 1, dtype=np.int_)
xlims = (alpha, beta)

curv_kwargs: dict = {
    "linewidth": 4,
    "alpha": 0.7,
}
font_kwargs: dict = {
    "fontsize": 12,
}
line_kwargs: dict = {
    "linewidth": 3,
    "alpha": 0.8,
    "linestyle": "-.",
    "color": "black",
    "zorder": 11,
}
hist_kwargs: dict = {
    "bins": boxes,
    "density": True,
    "edgecolor": "black",
    "alpha": 0.6,
    "align": "mid",
    # "normed": True,
}
walk_kwargs: dict = {
    "linewidth": 1,
    "alpha": 0.7,
}
plot_kwargs: dict = {
    "top": m,
}


curv_kwargs.update(args.kwarg_curve)
font_kwargs.update(args.kwarg_font)
line_kwargs.update(args.kwarg_line)
hist_kwargs.update(args.kwarg_hist)
walk_kwargs.update(args.kwarg_walk)
plot_kwargs.update(args.kwarg_plot)

if is_ode:
    x_name = r"$c_1$"
    y_name = r"$c_2$"
else:
    x_name = r"$x$"
    y_name = r"$y$"


names = [x_name, y_name]
colors = ["b", "r"]

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

    fig1 = plt.figure(figsize=(5, 2.5))
    fig2 = plt.figure(figsize=(5, 2.5))
    fig3 = plt.figure(figsize=(5, 5))
    fig4 = plt.figure(figsize=(5, 5))

    ax1 = fig1.add_subplot()
    ax2 = fig2.add_subplot()
    ax3 = fig3.add_subplot()
    ax4 = fig4.add_subplot()

    # figs = [fig1, fig2, fig3]
    figs = [plt.figure(i) for i in plt.get_fignums()]
    all_axes = [ax1, ax2, ax3]
    walk_axes = [ax1, ax2]
    hist_axes = [ax3]

    if args.compare_plots:
        figs.append(fig4)
        all_axes.append(ax4)
        hist_axes.append(ax4)

    if args.compare_plots:
        for i, ax in enumerate(walk_axes):
            init_cond_string: str = "{}".format(file_choice.loc[i, "initcond"])
            file_path = os.path.join(data_env, file_choice.loc[i, "file_name"])

            numpy_data = np.load(file_path)
            time = numpy_data["time"]
            states = numpy_data["states"]
            # print(time)
            # continue

            gillespies.plot_steps(
                ax=ax,
                time=time,
                results=states,
                color=colors[i],
                xstart=init_cond_string,
                plot_starts=args.include_starts,
                plot_kwargs=plot_kwargs,
            )

            gillespies.plot_hists(
                hist_axes[i],
                states,
                colors,
                init_cond_string,
            )

    else:
        init_cond_string: str = "{}".format(file_choice.loc["initcond"])
        file_path = os.path.join(data_env, file_choice.loc["file_name"])

        numpy_data = np.load(file_path)
        time, states = numpy_data["time"], numpy_data["states"]

        for i, ax in enumerate(walk_axes):
            plot_kwargs.update(label=names[i], color=colors[i])

            gillespies.plot_walk(
                ax,
                time=time,
                steps=states[i, :],
                plot_kwargs=plot_kwargs,
                xstart=init_cond_string,
                plot_starts=args.include_starts,
            )

            y_max = max(states[i, :].max(), 100)
            ax.set_yticks([y for y in np.linspace(0, y_max + 1, 5, dtype=np.int_)])
            ax.set_xlim(left=0)
            ax.set_ylim(bottom=0, top=y_max)

        for i, ax in enumerate(hist_axes):
            gillespies.plot_hist(
                ax,
                results=states[i, :],
                color=colors[i],
                label=names[i],
            )

    if is_ode:
        for ax in hist_axes:
            ax.set_title(
                "Cell Population Densities",
                fontdict=font_kwargs,
            )
    else:
        for ax in hist_axes:
            ax.set_title(
                "Distribution of Gene Copy Number\n {}".format(para_version),
                fontdict=font_kwargs,
            )

    for i, ax in enumerate(hist_axes):
        ax.set_xlim(xlims)
        ax.set_ylim(bottom=0)
        ax.set_ylabel(
            "Fraction of States",
            fontdict=font_kwargs,
        )
        if args.compare_plots:
            ax.set_xlabel(
                "Counts",
                fontdict=font_kwargs,
            )
        else:
            ax.set_xlabel(
                "Cell Counts",
                fontdict=font_kwargs,
            )

    # ax3.set_ylabel("Density", fontsize=12)
    # ax4.set_ylabel("Density", fontsize=12)

    for ax in all_axes:
        ax.legend(loc="upper right", fontsize=10)

    for fig in figs:
        fig.tight_layout()

    if args.plot_fixedpoints:
        inputs = [*filtered_frame.loc[model_choice, "metadata"].items()]
        input_dict = dict([(string.replace("-", "_"), val) for string, val in inputs])
        if args.compare_plots:
            gillespies.plot_walk_fixed(ax1, model, "x", xmax=1e10, parameters=input_dict)
        else:
            gillespies.plot_walk_fixed(ax1, model, "x", xmax=1e10, parameters=input_dict)
            gillespies.plot_walk_fixed(ax2, model, "y", xmax=1e10, parameters=input_dict)
            gillespies.plot_hist_fixed(ax3, model, "x", ymax=1, parameters=input_dict)

        # print(filtered_frame.loc[model_choice, "file_name"])
        # print(filtered_frame.loc[model_choice, "metadata"])
        # print(input_dict)

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

    file_path = os.path.join(data_env, file_choice.loc["file_name"])
    numpy_data = np.load(file_path)

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
        # "broken_streamlines":False,
    }
    stream_kwargs.update(args.kwarg_stream)

    fig1 = figure(figsize=(6, 6))
    ax1 = fig1.add_subplot()

    phaseies.plot_phase_space(ax1, c1, c2, dU, dV, **stream_kwargs)

    ax1.set_xlim(left=0, right=100)
    ax1.set_ylim(bottom=0, top=100)


figs = [plt.figure(i) for i in plt.get_fignums()]
for fig in figs:
    fig.tight_layout()


## Taken from https://www.geeksforgeeks.org/save-multiple-matplotlib-figures-in-single-pdf-file-using-python/
def save_image(filename):
    # PdfPages is a wrapper around pdf
    # file so there is no clash and create
    # files with no error.
    filename += ".pdf"

    ## We use a "with" statement to keep the system self contained
    with PdfPages(filename) as p:
        # get_fignums Return list of existing
        # figure numbers
        fig_nums = plt.get_fignums()
        figs = [plt.figure(n) for n in fig_nums]

        # iterating over the numbers in list
        for fig in figs:
            # and saving the files
            fig.savefig(p, format="pdf")

    # # close the object
    # p.close()


latest_file = os.path.join(fig_env, "latest_plot")
if args.save:
    if is_ode:
        file_path = os.path.join(ode_dir, file_name)
    else:
        file_path = os.path.join(fpe_dir, file_name)

    save_image(file_path)


if len(figs) > 1:
    save_image(latest_file)
else:
    plt.figure(1).savefig(latest_file + ".png", format="png")


if args.show:
    show()

print("Done Plotting/Saving")

exit()
