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
phs_env = str(os.getenv("PHS_FIGURE_ENV"))
ode_env = str(os.getenv("ODE_FIGURE_ENV"))

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
        "ode",
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
        "ode_7_2",
        "ode_8_3",
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
    "--fixed-points",
    dest="plot_fixedpoints",
    help="Include the fixed points in the histogram and walks",
    action="store_true",
)
bool_args.add_argument(
    "-op",
    "--on-phase",
    dest="plot_on_phase",
    help="Plots the stochastic system on the associated phase diagram",
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
# kwarg_string_parser = re.compile(r"(\w[^=]*)=\s?([^ ,]*)")


def parse_filters(string):
    matches: list[tuple] = filter_captures.findall(string)
    parameters = [
        (str(key).replace("-", "_").replace("'", "p").replace(" ", ""), float(value))
        for key, value in matches
    ]
    return parameters


# @njit
def string_or_float(input: str) -> str | float:
    if input.lower() == "true" or input.lower() == "false":
        return bool(input)
    else:
        try:
            return float(input)
        except:
            return str(input)


def parse_kwarg_string(string):
    matches: list[tuple] = filter_captures.findall(string)
    parameters = [
        (str(key).replace(" ", ""), string_or_float(value)) for key, value in matches
    ]
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
kwarg_parse.add_argument(
    "+rc", "++rcparams", dest="kwarg_rcparams", type=parse_kwarg_string
)


args = parser.parse_args()


## Compiling the defaults and the choice of parameters


data_source: str = args.source
model: str = args.model
model_count: int = args.index


is_ode = False
check_ode = re.compile(r"ode").search(model)
if check_ode is not None:
    is_ode = True

if data_source == "phase" or data_source == "ode":
    is_ode = True


ratio_pattern: re.Pattern = re.compile(r"R(?P<ratio>b.n)R")
capture_patterns: re.Pattern = re.compile(
    r"(?P<datasource>^[^M]*)"  # phase, ssa, etc.
    r"M(?P<model>[^P]*)"  # Model specifier
    r"P(?P<rawmetadata>[^SR]*)"  # Parameters and
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


raw_frame["metadata"] = raw_frame["rawmetadata"].map(parse_metadata)


no_ratio = raw_frame["ratio"].isna()
raw_frame.loc[no_ratio, "ratio"] = ""


### Now we can mutate the dictionaries in memory and dynamically create columns
defined_metadata = raw_frame["metadata"].map(lambda d: len(d) > 0)
raw_frame = raw_frame.loc[defined_metadata]

## That dict.pop has a default for missing entries is a literal godsend
ssa_index = raw_frame.loc[raw_frame["datasource"] == "ssa"].index
raw_frame["count"] = [d.pop("num", np.nan) for d in raw_frame.loc[:, "metadata"]]


## From the list, we then count the number of times each model shows up,
## assigning the number of each step accordingly.
## From which multi-index is generated for the dataframe
## so as to have an easier time choosing which file to plot
raw_frame["model_index"] = raw_frame.groupby("model").cumcount()
file_frame = raw_frame.set_index(["datasource", "model", "model_index"]).sort_index()

parameters = file_frame["metadata"].apply(pd.Series)
file_frame = pd.concat([file_frame, parameters], axis=1)
sourced_frame = (
    file_frame.loc[(data_source, model)]
    .dropna(axis="columns", how="all")
    .reset_index(drop=True)
)


## Just brute forces a sorting method, which is designed to fail if the list
## gets too small. Returns the next best in that case
def filter_frames(filters: list[tuple[str, float]], df: pd.DataFrame) -> pd.DataFrame:
    for string, value in filters:
        try:
            df_temp = df.loc[df[string] == value]
            if df_temp.empty:
                continue
            df = df_temp
        except:
            assert type(df) is pd.DataFrame, "Not a dataframe"
            continue

    return df


if args.filter is not None:
    sourced_frame = filter_frames(args.filter, sourced_frame).reset_index()


if args.opts:
    if data_source == "ssa":
        print("{}".format(args.model))
        print(sourced_frame[["initcond", "count", "steps", "ratio", "t"]])
        exit()
    elif is_ode:
        print("{}".format(args.model))
        print(sourced_frame[["metadata"]])
        exit()


## When choosing which index
if args.index is None:
    if args.compare_plots:
        m_index = sourced_frame.index.values[-2]
    else:
        m_index = sourced_frame.index.values[-1]
else:
    m_index = args.index

## For if it is a comparison
## We try to set as many things up here to make it easier in the future
### Placing things in lists so we can use indices to generalize
initial_list = []
path_list = []
if args.compare_plots:
    data_set_quantity = 2
    file_choice = sourced_frame.loc[m_index : m_index + 1].reset_index(drop=True)
    path_list.extend([path for path in file_choice.pop("file_name")])
    initial_list.extend([init for init in file_choice.pop("initcond")])
    data_frame = file_choice.loc[0]
else:
    data_set_quantity = 1
    file_choice = sourced_frame.loc[m_index]
    path_list.append(file_choice.pop("file_name"))
    initial_list.append(file_choice.pop("initcond"))
    data_frame = file_choice


#######################
#   File name stuff   #
#######################

## The base name
file_name: str = "{}:{}:".format(data_source, model).replace("_", "-")
if model == "5_2":
    para_version = data_frame["ratio"]
    file_name += "{}:".format(para_version).replace("$", "")


## Adding the epoch time to ensure the files don't overwrite eachother
epoch = data_frame["t"]
if data_source != "phase":
    c_data = data_set_quantity
else:
    c_data = file_choice["rawmetadata"]

# Format rapidfire
file_name += "I{}:".format(m_index)
if args.plot_fixedpoints:
    file_name += "C{}".format(c_data)
    file_name += "{}:".format("fp")
else:
    file_name += "C{}:".format(c_data)
file_name += "T{}".format(epoch)

#######################

## For custom plotting counts
if args.bin_count is not None:
    m: int = args.bin_count
else:
    m: int = data_frame["count"]

alpha = 0 - 1
beta = m + 1  # file_choice["count"]

boxes = arange(alpha, beta + 1, 1, dtype=np.int_)
xlims = (alpha, beta)
ylims = (alpha, beta)

##########################
#   Default Plot Stuff   #
##########################


curv_kwargs: dict = {
    "linewidth": 4,
    "alpha": 0.6,
    "zorder": 3,
}
font_kwargs: dict = {
    "fontsize": 12,
}
line_kwargs: dict = {
    "linewidth": 1,
    "alpha": 0.5,
    "linestyle": "-.",
    "zorder": 5,
}
hist_kwargs: dict = {
    "align": "left",
    "alpha": 0.4,
    "bins": boxes,
    "density": True,
    "histtype": "stepfilled",
    "linewidth": 1.25,
    "edgecolor": "k",
}
walk_kwargs: dict = {
    "linewidth": 1,
    "alpha": 0.6,
}
plot_kwargs: dict = {
    "top": m,
}
stream_kwargs: dict = {
    "density": 1.7,
    "color": "k",
    "arrowstyle": "->",
}
rcparams_kwargs: dict = {
    "axes.labelsize": 20,
    "axes.titleweight": "bold",
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    # "axes.titlecolor": "white",
    # "xtick.labelcolor": "white",
    # "ytick.labelcolor": "white",
    # "savefig.facecolor": "#c0c0ca",
}


if is_ode:
    plot_kwargs.update(x_name=r"$c_1$")
    plot_kwargs.update(y_name=r"$c_2$")
    plot_kwargs.update(z_name=r"$b$")
    plot_kwargs.update(w_name=r"$m$")
else:
    plot_kwargs.update(x_name=r"$x$")
    plot_kwargs.update(y_name=r"$y$")
    plot_kwargs.update(z_name=r"$n$")
    plot_kwargs.update(w_name=r"$b$")


if is_ode:
    if data_source == "ode":
        plot_kwargs.update(z_name=r"$m$")
        plot_kwargs.update(w_name=r"$m$")
    else:
        if model == "ode_8_3":
            plot_kwargs.update(z_name=r"$b$")
        else:
            plot_kwargs.update(z_name=r"$n$")
else:
    plot_kwargs.update(z_name=r"$n$")
    plot_kwargs.update(w_name=r"$b$")


##########################

curv_kwargs.update(args.kwarg_curve)
font_kwargs.update(args.kwarg_font)
hist_kwargs.update(args.kwarg_hist)
line_kwargs.update(args.kwarg_line)
walk_kwargs.update(args.kwarg_walk)
plot_kwargs.update(args.kwarg_plot)
stream_kwargs.update(args.kwarg_stream)

rcparams_kwargs.update(args.kwarg_rcparams)
plt.rcParams.update(**rcparams_kwargs)

legend_font_size = font_kwargs.pop("legend_fontsize", 11)

x_name = plot_kwargs.pop("x_name")
y_name = plot_kwargs.pop("y_name")
z_name = plot_kwargs.pop("z_name")
w_name = plot_kwargs.pop("w_name")

x_color = plot_kwargs.pop("x_color", "r")
y_color = plot_kwargs.pop("y_color", "b")
z_color = plot_kwargs.pop("z_color", "g")
# w_color = plot_kwargs.pop("w_color", "k")


names = [x_name, y_name, z_name, w_name]
colors = [x_color, y_color, z_color]

#############################################################################
#                                                                           #
#   I've thought about it, realistically this isn't super memory intense,   #
#   so I'm just gonna load both plot classes into memory and not really     #
#   care one way or the other                                               #
#                                                                           #
#############################################################################

plot_class_dict = {
    "curvargs": curv_kwargs,
    "fontargs": font_kwargs,
    "histargs": hist_kwargs,
    "lineargs": line_kwargs,
    "walkargs": walk_kwargs,
    "x_name": x_name,
    "y_name": y_name,
    "n_name": z_name,
    "b_name": w_name,
}

gillespies = ps.gillespiePlotters(**plot_class_dict)
odeies = ps.odePlotters(stream_kwargs=stream_kwargs, **plot_class_dict)

#############################################################################

fig1 = plt.figure(num=1, figsize=(5, 2.5))
fig3 = plt.figure(num=3, figsize=(5, 5))

if args.compare_plots:
    fig2 = plt.figure(num=2, figsize=(5, 2.5))
    fig4 = plt.figure(num=4, figsize=(5, 5))

figs = [plt.figure(i) for i in plt.get_fignums()]
axes = [fig.add_subplot() for fig in figs]

## By partitioning the axes list where the first half is the hists, and the
## second half is the walks, we can just use the data_set_quantity to split it
walk_axes = axes[:data_set_quantity]
hist_axes = axes[data_set_quantity:]


## Trying to use a 'switch case' to make the plotting hopefully smoother
## the idea is that we can seperate the source of the data here, and plot it after?
numpy_datas = []
match data_source:
    case "ssa":
        if args.compare_plots:
            for i, ax in enumerate(walk_axes):
                init_cond_string: str = "{}".format(initial_list[i])
                file_path = os.path.join(data_env, path_list[i])
                numpy_data = np.load(file_path)
                numpy_datas.append(numpy_data)

                time = numpy_data["time"]
                states = numpy_data["states"]

                plot_kwargs.update(label=names[i], color=colors[i])

                gillespies.plot_steps(
                    ax=ax,
                    time=time,
                    results=states,
                    xstart=init_cond_string,
                    plot_starts=args.include_starts,
                    plot_kwargs=plot_kwargs,
                )

                gillespies.plot_hists(
                    hist_axes[i],
                    states,
                )

        else:
            # init_cond_string: str = "{}".format(data_frame["initcond"])
            file_path = os.path.join(data_env, path_list[0])

            numpy_data = np.load(file_path)
            time, states = numpy_data["time"], numpy_data["states"]

            for i, ax in enumerate(walk_axes):
                plot_kwargs.update(label=names[i], color=colors[i])

                gillespies.plot_walk(
                    ax,
                    time=time,
                    steps=states[i, :],
                    plot_kwargs=plot_kwargs,
                    plot_starts=args.include_starts,
                )

                gillespies.plot_hist(
                    hist_axes[i], states[i, :], label=names[i], color=colors[i]
                )

                y_max = min(states[i, :].max(), 100)
                ax.set_yticks([y for y in np.linspace(0, y_max + 1, 5, dtype=np.int_)])
                ax.set_xlim(left=0)
                ax.set_ylim(bottom=0, top=y_max)
        #
        #     if model in ["2L"]:
        #         # for i, ax in enumerate(hist_axes):
        #         x_temp = states[0, :]
        #         y_temp = states[1, :]
        #         hist_axes[0].hist2d(
        #             x_temp,
        #             y_temp,
        #             bins=boxes,
        #             density=True,
        #             label="Distribution of chemical species",
        #         )
        #     else:
        #         for i, ax in enumerate(hist_axes):
        #             gillespies.plot_hist(
        #                 ax,
        #                 results=states[i, :],
        #                 color=colors[i],
        #                 label=names[i],
        #             )
        # break
        # plt.show()
        # exit()
        if is_ode:
            for ax in hist_axes:
                ax.set_title(
                    "Cell Population Densities",
                    fontdict=font_kwargs,
                )
        else:
            for ax in hist_axes:
                if model == "5_2":
                    ax.set_title(
                        "Distribution of Gene Copy Number\n"
                        + "${}$".format(para_version.replace("$", "")),
                        fontdict=font_kwargs,
                    )
                else:
                    ax.set_title(
                        "Distribution of Gene Copy Number",
                        fontdict=font_kwargs,
                    )

        for i, ax in enumerate(hist_axes):
            ax.set_xlim(xlims)
            ax.set_ylim(bottom=0)
            ax.set_ylabel(
                "Fraction of States",
                fontdict=font_kwargs,
            )
            ax.set_xlabel(
                "Counts",
                fontdict=font_kwargs,
            )

        # ax3.set_ylabel("Density", fontsize=12)
        # ax4.set_ylabel("Density", fontsize=12)
        # for ax in hist_axes:
        #     lines = ax.get_lines()
        #     for ln in lines:
        #         ln.set_alpha(1)
        #         ln.set_linewidth(10)

        for ax in axes:
            ax.legend(loc="upper right", fontsize=legend_font_size)

        for fig in figs:
            fig.tight_layout()

        if args.plot_fixedpoints:
            inputs = data_frame["metadata"]
            # input_dict = dict([(string.replace("-", "_"), val) for string, val in inputs])

            if args.compare_plots:
                for i, ax in enumerate(walk_axes):
                    odeies.plot_fixed("ode", ax, inputs, xmax=100, ymax=100)
                    gillespies.plot_walk_fixed(
                        ax,
                        model,
                        "x",
                        xmax=1e10,
                        parameters=inputs,
                        color=x_color,
                    )
                    gillespies.plot_walk_fixed(
                        ax,
                        model,
                        "y",
                        xmax=1e10,
                        parameters=inputs,
                        color=y_color,
                    )

                for i, ax in enumerate(hist_axes):
                    odeies.plot_fixed("ode", ax, inputs, axis="x", xmax=100, ymax=100)
                    gillespies.plot_hist_fixed(
                        ax,
                        model,
                        "x",
                        ymax=1,
                        parameters=inputs,
                        color=x_color,
                    )
                    gillespies.plot_hist_fixed(
                        ax,
                        model,
                        "y",
                        ymax=1,
                        parameters=inputs,
                        color=y_color,
                    )
            else:
                gillespies.plot_walk_fixed(
                    walk_axes[0], model, "x", xmax=1e10, parameters=inputs
                )
                gillespies.plot_walk_fixed(
                    walk_axes[1], model, "y", xmax=1e10, parameters=inputs
                )

                gillespies.plot_hist_fixed(
                    hist_axes[0], model, "x", ymax=1, parameters=inputs
                )

        if args.plot_on_phase:
            plt.close("all")

            fig1 = plt.figure(1, figsize=(6, 6))
            ax1 = fig1.add_subplot()

            phase_frame = file_frame.loc[("phase", "jesper")]
            data_frame["metadata"]["n1"] *= 2
            data_frame["metadata"]["n2"] *= 2
            phase_filters: list[tuple[str, float]] = [*data_frame["metadata"].items()]

            filtered_phase = filter_frames(phase_filters, phase_frame)
            phase_name = filtered_phase.reset_index().loc[0, "file_name"]
            phase_path = os.path.join(data_env, phase_name)
            phase_data = np.load(phase_path)

            c1, c2 = phase_data["c1"], phase_data["c2"]
            dU, dV = phase_data["dU"], phase_data["dV"]

            odeies.plot_phase_space(ax1, c1, c2, dU, dV)
            odeies.plot_trajectories(ax1, states[0, :], states[1, :])

            ax1.set_xlim(0, 100)
            ax1.set_ylim(0, 100)

    case "phase":
        plt.close("all")

        file_path = os.path.join(data_env, path_list[0])
        numpy_data = np.load(file_path)

        c1, c2 = numpy_data["c1"], numpy_data["c2"]
        dU, dV = numpy_data["dU"], numpy_data["dV"]

        c1_null = numpy_data["c1_nullcline"]
        c2_null = numpy_data["c2_nullcline"]

        fig1 = plt.figure(1, figsize=(6, 6))
        ax1 = fig1.add_subplot()

        # fig2 = plt.figure(figsize=(6, 6))
        # ax2 = fig2.add_subplot()

        odeies.plot_phase_space(ax1, c1, c2, dU, dV)
        odeies.plot_nullclines(ax1, *c2_null, *c1_null)

        # level_set = numpy_data["level"]
        # odeies._plot_phase_curve(ax1, *c1_null, label="$c_2$ Nullcline", color="b")
        # odeies._plot_phase_curve(ax1, *level_set, label="Levelset", color="g")

        if args.plot_fixedpoints:
            parameters: dict[str, float] = data_frame["metadata"]
            odeies.plot_fixed(data_source, ax1, parameters, xmax=100, ymax=100)

        ax1.set_xlim(left=0, right=100)
        ax1.set_ylim(bottom=0, top=100)

        ax1.legend(loc="upper right", fontsize=legend_font_size)

        ax1.set_xlabel("Number of $c_1$")
        ax1.set_ylabel("Number of $c_2$")

    case "ode":
        plt.close("all")

        file_path = os.path.join(data_env, path_list[0])
        numpy_data = np.load(file_path)

        time = numpy_data["time"]
        solutions = numpy_data["solutions"]

        fig1 = plt.figure(1, figsize=(6, 4))
        ax1 = fig1.add_subplot()

        # fig2 = plt.figure(figsize=(6, 6))
        # ax2 = fig2.add_subplot()

        odeies.plot_curves(ax1, time, solutions)

        if args.plot_fixedpoints:
            parameters: dict[str, float] = data_frame["metadata"]
            odeies.plot_fixed(
                data_source, ax1, parameters, total_vars=2, xmax=time[-1], ymax=100
            )

        ax1.set_xlim(left=0, right=time[-1])
        ax1.set_ylim(bottom=0)

        ax1.legend(loc="upper right", fontsize=legend_font_size)

        ax1.set_xlabel("Time")
        ax1.set_ylabel("Value of Variable")
    case _:
        print("Fucked it")
        exit()


figs = [plt.figure(i) for i in plt.get_fignums()]
for fig in figs:
    fig.tight_layout()

if model == "2L":
    hist_axes[0].set_title("2D histogram of chemical species distribution")


## Taken from https://www.geeksforgeeks.org/save-multiple-matplotlib-figures-in-single-pdf-file-using-python/
# PdfPages is a wrapper around pdf file so there is no clash
# and create files with no error.
def save_image(filename):
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


latest_file = os.path.join(fig_env, "latest_plot")
if args.save:
    ## Python Switch Case
    match data_source:
        case "ssa":
            if is_ode:
                file_path = os.path.join(ode_dir, file_name)
            else:
                file_path = os.path.join(fpe_dir, file_name)
        case "phase":
            file_path = os.path.join(phs_env, file_name)
        case "ode":
            file_path = os.path.join(ode_env, file_name)
        case _:
            file_path = os.path.join(data_env, file_name)

    ## Just making sure the shit is sorted
    save_image(file_path)

if args.show or not args.save:
    show()

# print("TESTING MODE")
# exit()


print("Done Plotting/Saving")

exit()
