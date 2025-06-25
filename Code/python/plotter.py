#!../.venv/bin/python3
import argparse
import json
import os
import re

import matplotlib as mpl
import numpy as np
import pandas as pd
import plottingstuff as ps
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# from pylab import *

# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
fig_env = str(os.getenv("THESIS_FIGURE_PATH"))
fpe_env = str(os.getenv("FPE_FIGURE_ENV"))
phs_env = str(os.getenv("PHS_FIGURE_ENV"))
ode_env = str(os.getenv("ODE_FIGURE_ENV"))

data_env = str(os.getenv("THESIS_DATA_PATH"))


presentation_fig_env = os.path.join(
    str(os.getenv("THESIS_REPORT")), "..", "Presentation", "images"
)

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
        "anal",
    ],
    type=str,
)

parser.add_argument(
    "model",
    help="State which model to run",
    # choices=[
    #     "jesper",
    #     "ode",
    #     "2S",
    #     "2L",
    #     "5_2",
    #     "5_3",
    #     "5_2_fixed",
    #     "5_3_fixed",
    #     "ode_2",
    #     "ode_2_2",
    #     "ode_3",
    #     "ode_3_2",
    #     "ode_3_2_alt",
    #     "ode_3_3",
    #     "ode_5_2",
    #     "ode_5_3",
    #     "ode_7_2",
    #     "ode_8_3",
    # ],
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
    # default=int(100),
)
parser.add_argument(
    "-b",
    "--bins",
    dest="bin_count",
    help="Number of bins",
    type=int,
    # default=int(100),
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
check_ode = re.compile(r"(ODE|ode)").search(model)
if check_ode is not None:
    is_ode = True

if data_source == "phase" or data_source == "ode":
    is_ode = True


ratio_pattern: re.Pattern = re.compile(r"R(?P<ratio>b.n)R")
capture_patterns: re.Pattern = re.compile(
    r"(?P<datasource>^[^M]*)"  # phase, ssa, etc.
    r"M(?P<model>[^L]*)"  # Model specifier
    r"L(?P<language_source>[^TP]*)"
    r"[TP](?P<serial>\d+)"
    r"\.(?P<filetype>\w*)$"
)


## Now that I have refactored the metadata storage, the organization process
## becomes a lot less complicated
records = []
file_names = os.listdir(data_env)
for fn in file_names:
    rec = dict()

    matches = capture_patterns.match(fn)

    if matches is not None:
        rec.update(matches.groupdict())

    rec["file_name"] = fn
    records.append(rec)


## So we assign it to a dataframe for easy access
raw_frame = pd.DataFrame(records)
raw_frame.sort_values(by="serial", inplace=True)  ## Sorting by unix epoch time, ascending

## From the list, we then count the number of times each model shows up,
## assigning the number of each step accordingly.
## From which multi-index is generated for the dataframe
## so as to have an easier time choosing which file to plot
file_frame = raw_frame.set_index(["datasource", "model"])
file_frame["model_index"] = file_frame.groupby("model").cumcount()


## Now we select only the entries with the entries we care about
sourced_frame = file_frame.loc[(data_source, model)]
sourced_frame = sourced_frame.set_index(["model_index"]).sort_index()


# print(sourced_frame.loc[0])
# exit()

# ## Just brute forces a sorting method, which is designed to fail if the list
# ## gets too small. Returns the next best in that case
# def filter_frames(filters: list[tuple[str, float]], df: pd.DataFrame) -> pd.DataFrame:
#     for string, value in filters:
#         try:
#             df_temp = df.loc[df[string] == value]
#             if df_temp.empty:
#                 continue
#             df = df_temp
#         except:
#             assert type(df) is pd.DataFrame, "Not a dataframe"
#             continue
#
#     return df


if args.filter is not None:
    print("Not updated yet")  # TODO: this thing
    exit()
    sourced_frame = filter_frames(args.filter, sourced_frame).reset_index()


if args.opts:
    match data_source:
        case "ssa":
            print("{}".format(args.model))
            # print(sourced_frame[["initcond", "count", "steps", "ratio", "t"]])
            exit()
        case "ode":
            print("{}".format(args.model))
            # print(sourced_frame[["metadata"]])
            exit()
        case "phase":
            print("{}".format(args.model))
            # print(sourced_frame[["m0", "k1", "k2"]])
            exit()


## When choosing which index
if args.index is None:
    if args.compare_plots:
        m_index = sourced_frame.index.values[-2]
    else:
        m_index = sourced_frame.index.values[-1]
else:
    m_index = args.index


## if should have a fixed_point plotted
if args.plot_fixedpoints:
    set_fixed_points = True
else:
    set_fixed_points = False


## For if it is a comparison
## We try to set as many things up here to make it easier in the future
### Placing things in lists so we can use indices to generalize
initial_list = []
path_list = []
if args.compare_plots:
    data_set_quantity = 2
    file_choice = sourced_frame.loc[m_index : m_index + 1].reset_index(drop=True)
    path_list.extend([path for path in file_choice.pop("file_name")])
    # initial_list.extend([init for init in file_choice.pop("initcond")])
    data_frame = file_choice.loc[0]
else:
    data_set_quantity = 1
    file_choice = sourced_frame.loc[m_index]
    path_list.append(file_choice.pop("file_name"))
    # initial_list.append(file_choice.pop("initcond"))
    data_frame = file_choice


# print(data_frame)
# exit()

#######################
#   File name stuff   #
#######################

## The base name
file_name: str = "{}:{}:".format(data_source, model).replace("_", "-")
if model == "5_2":
    para_version = data_frame["ratio"]
    file_name += "{}:".format(para_version).replace("$", "")


## Adding the epoch time to ensure the files don't overwrite eachother
epoch = data_frame["serial"]
if data_source != "phase":
    c_data = data_set_quantity
else:
    c_data = file_choice["rawmetadata"]

# Format rapidfire
file_name += "I{}:".format(m_index)
file_name += "C{}".format(c_data)
if set_fixed_points:
    file_name += "{}:".format("fp")
file_name += "T{}".format(epoch)


#############################################################################

## Trying to use a 'switch case' to make the plotting hopefully smoother
## the idea is that we can seperate the source of the data here, and plot it after?
numpy_datas = []
metadatas = []
for i, path in enumerate(path_list):
    file_path = os.path.join(data_env, path)
    numpy_data = np.load(file_path)
    numpy_datas.append(numpy_data)
    ## As the metadata is encoded in utf-8 for julia compatibility, we need
    ## to decode it here for the json module to properly read it
    metadatas.append(json.loads(numpy_data["metadata"].tobytes()))

#######################

# print(metadatas[0])
# exit()

## For custom plotting counts
if data_source == "ssa":
    if args.bin_count is not None:
        m: int = args.bin_count
    else:
        m: int = metadatas[0]["number_of_particles"]
else:
    m = 100

alpha = 0 - 1
beta = m + 1  # file_choice["count"]

boxes = np.arange(alpha, beta + 1, 1, dtype=np.int_)
xlims = (alpha, beta)
ylims = (alpha, beta)

#######################


##########################
#   Default Plot Stuff   #
##########################


curv_kwargs: dict = {
    "linewidth": 4,
    "alpha": 0.6,
    "zorder": 3,
}
font_kwargs: dict = {
    "fontsize": 16,
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
    "savefig.facecolor": "#c0c0ca",
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
        if model == "ODE8Par3Var":
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
w_color = plot_kwargs.pop("w_color", "m")


names = [x_name, y_name, z_name, w_name]
colors = [x_color, y_color, z_color, w_color]

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

# Define the colors for the custom colormap
# colors = [(0, 0, 1), (1, 1, 0), (1, 0, 0)]  # Blue, Yellow, Red
# cmap_name = "custom_cmap"

# [ (0.267004, 0.004874, 0.329415, 1.0) ] # Lowest color of viridis
# [ (0.9463847846200787, 0.9656862745098039, 0.9656862667892201, 1.0) ] # Lowest of 25 parts bone_r


# Create the custom colormap
custom_cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["darkolivegreen", "snow"])


pre_plasma = plt.cm.plasma
pre_bone = plt.cm.bone
pre_color_1 = pre_plasma(np.linspace(0, 1, 128))
pre_color_2 = pre_bone(np.linspace(0, 1, 128))

stacked_colors = np.vstack((pre_color_1, pre_color_2))
custom_cmap = mpl.colors.LinearSegmentedColormap.from_list("custom_cmap", stacked_colors)

gillespies = ps.gillespiePlotters(**plot_class_dict)
odeies = ps.odePlotters(stream_kwargs=stream_kwargs, **plot_class_dict)

## Using the switch case to 'prime' the plotting space so to speak
match data_source:
    case "anal":
        fig1 = plt.figure(num=1, figsize=(5, 2.5))

    case "ssa":
        fig1 = plt.figure(num=1, figsize=(5, 2.5))
        fig3 = plt.figure(num=3, figsize=(6, 5))

        if args.compare_plots:
            fig2 = plt.figure(num=2, figsize=(5, 2.5))
            fig4 = plt.figure(num=4, figsize=(5, 5))

    case "ode":
        fig1 = plt.figure(num=1, figsize=(5, 2.5))

        if args.compare_plots:
            fig2 = plt.figure(num=2, figsize=(5, 2.5))

    case "phase":
        fig1 = plt.figure(num=1, figsize=(5, 5))

## Programmatically accessing the figs and axes
figs = [plt.figure(i) for i in plt.get_fignums()]
axes = [fig.add_subplot() for fig in figs]


## By partitioning the axes list where the first half is the hists, and the
## second half is the walks, we can just use the data_set_quantity to split it
walk_axes = axes[:data_set_quantity]
hist_axes = axes[data_set_quantity:]


### The actual plotting portion
match data_source:
    case "anal":
        pass

    case "ssa":
        for i, data in enumerate(numpy_datas):
            walk_axis = walk_axes[i]
            hist_axis = hist_axes[i]

            # print(args.language_source)
            # print(typeof(args.language_source))
            # exit()
            time = data["time"]
            if data_frame["language_source"] == "python":
                states = data["states"]
            elif data_frame["language_source"] == "julia":
                states = data["states"].T
            else:
                print("Fuck")
                exit()

            # print(states.shape)
            # print(time.shape)
            # plt.plot(time, states[i, :])
            # plt.show()
            # exit()

            try:
                metadata = metadatas[i]
            except:
                init_cond_string: str = "{}".format("[100 0]")
            else:
                init_cond_string: str = "{}".format(metadata["initial_condition"])

            plot_kwargs.update(label=names[i], color=colors[i])
            gillespies.plot_walks(
                ax=walk_axis,
                time=time,
                results=states,
                xstart=init_cond_string,
                plot_starts=args.include_starts,
                plot_kwargs=plot_kwargs,
            )

            # gillespies.plot_hists(
            #     ax=hist_axis,
            #     results=states,
            # )
            # color=colors[i],
            # label=names[i],

            # print(metadata)
            # exit()

            walk_axis.set_xlabel("Time", fontdict=font_kwargs)
            walk_axis.set_ylabel("Count", fontdict=font_kwargs)

            walk_axis.set_xlim(left=0)  # , right=30000)
            walk_axis.set_ylim(bottom=0, top=m)
            # ax.set_yticks([y for y in range(y_min, y_max + 1, axis_step)])

            if is_ode:
                for ax in hist_axes:
                    ax.set_title(
                        "Cell Population Densities",
                        fontdict=font_kwargs,
                    )
            else:
                for ax in hist_axes:
                    if model == "Novel5Par2Var":
                        ax.set_title(
                            "Distribution of Gene Copy Number\n"
                            + "${}$".format(para_version.replace("$", "")),
                            fontdict=font_kwargs,
                        )
                    else:
                        ax.set_title(
                            "Distribution of Gene Copy Number to Resistence Protein",
                            fontdict=font_kwargs,
                        )

            # print(metadata)
            # exit()

            file_name += ":2dHist"
            hist_axis.grid(False)
            hist2d_boxes = np.arange(alpha, beta + 10, 1)
            hist_axis.hist2d(
                states[0, :],
                states[1, :],
                bins=(hist2d_boxes - 0.5, hist2d_boxes - 0.5),
                density=True,
                norm=mpl.colors.LogNorm(),
                cmap=mpl.cm.bone_r,
                # cmap=mpl.cm.plasma,
                label="Heatmap of $c_1$ and $c_2$",
            )
            hist_axis.set_xlabel(
                "Count of {}".format(x_name),
                fontdict=font_kwargs,
            )
            hist_axis.set_ylabel(
                "Count of {}".format(y_name),
                fontdict=font_kwargs,
            )

            # axe =
            # import matplotlib.pyplot as plt
            # my_cmap = plt.cm.jet
            # my_cmap.set_under("w", 1)
            # ...
            # fig3.colorbar(im, ax=hist_axes[i], label="Fraction of Samples")
            # phase_frame = file_frame.loc[("phase", "jesper")]
            # # data_frame["metadata"]["n1"] *= 2
            # # data_frame["metadata"]["n2"] *= 2
            # phase_filters: list[tuple[str, float]] = [*data_frame["metadata"].items()]
            # filtered_phase = filter_frames(phase_filters, phase_frame)
            # phase_name = filtered_phase.reset_index().loc[0, "file_name"]
            # phase_path = os.path.join(data_env, phase_name)
            # phase_data = np.load(phase_path)
            # c1, c2 = phase_data["c1"], phase_data["c2"]
            # dU, dV = phase_data["dU"], phase_data["dV"]
            # c1_null = phase_data["c1_nullcline"]
            # c2_null = phase_data["c2_nullcline"]
            # level = phase_data["level"]
            # odeies.plot_phase_space(
            #     hist_axes[i],
            #     c1,
            #     c2,
            #     dU,
            #     dV,
            #     linewidth=0.75,
            #     color="k",
            # )
            # # odeies.plot_trajectories(hist_axes[i], *level / 2, label="Levelset curve")
            # # odeies.plot_trajectories(hist_axes[i], *c1_null)
            # # odeies.plot_trajectories(hist_axes[i], *c2_null)
            # axe = hist_axes[i].contourf(*plot_space, output, 20)
            # fig3.colorbar(axe)
            # hist_axes[i].set_xlim(0, beta)
            # hist_axes[i].set_ylim(0, beta)
            # file_name += ":2dhist"
            # file_name += ":compared"
            # file_name += ":twice"
            # plt.show()
            # exit()
            # shared_ax = hist_axes[i].twinx()
            # temp_k1 = metadata["k1"]
            # temp_k2 = metadata["k2"]
            # temp_n = data_frame["count"]

            # break
            # from gillespie import analytical as dga
            #
            # # two_part_anal = dga.SimpleTwoChemical(k1=temp_k1, k2=temp_k2, n=m)
            # two_part_anal = dga.MultivariateTwoChemical(k1=temp_k1, k2=temp_k2, n=m)
            #
            # x_domain = np.linspace(0, 1, 500, dtype=np.float64)
            # # y_range = two_part_anal.stationary(x_domain)
            # x_domain = np.linspace(0, 1, 500, dtype=np.float64)
            # y_range = two_part_anal.stationary(x_domain, x_domain)
            # hist_axes[i].contour(y_range)
            # # shared_ax.yaxis.tick_right()
            # # shared_ax.plot(x_domain * m, y_range)
            #
            # hist_axes[i].set_xlim(left=alpha, right=beta)
            # hist_axes[i].set_ylim(bottom=alpha, top=beta)
            # shared_ax.set_ylim(bottom=0)
            # shared_ax.hlines([1, 0.1, 0.01], 0, m)
            # shared_ax.set_ylim(bottom=0, top=y_range.max())
            # plt.show()
            # exit()

            # for i, ax in enumerate(walk_axes):
            #     plot_kwargs.update(label=names[i], color=colors[i])
            #
            #     gillespies.plot_walk(
            #         ax,
            #         time=time,
            #         steps=states[i, :],
            #         plot_kwargs=plot_kwargs,
            #         plot_starts=args.include_starts,
            #     )
            #
            #     gillespies.plot_hist(
            #         hist_axes[i], states[i, :], label=names[i], color=colors[i]
            #     )
            #
            #     y_max = min(states[i, :].max(), 100)
            #     ax.set_yticks([y for y in np.linspace(0, y_max + 1, 5, dtype=np.int_)])
            #     ax.set_xlim(left=0)
            #     ax.set_ylim(bottom=0, top=y_max)

    #
    #     for i, ax in enumerate(hist_axes):
    #         ax.set_xlim(xlims)
    #         ax.set_ylim(bottom=0)
    #         ax.set_ylabel(
    #             "Fraction of States",
    #             fontdict=font_kwargs,
    #         )
    #         ax.set_xlabel(
    #             "Counts",
    #             fontdict=font_kwargs,
    #         )
    #     for fig in figs:
    #         fig.tight_layout()
    #
    #         inputs = data_frame["metadata"]
    #         # input_dict = dict([(string.replace("-", "_"), val) for string, val in inputs])
    #
    #         if args.compare_plots:
    #             for i, ax in enumerate(walk_axes):
    #                 odeies.plot_fixed("ode", ax, inputs, xmax=100, ymax=100)
    #                 gillespies.plot_walk_fixed(
    #                     ax,
    #                     model,
    #                     "x",
    #                     xmax=1e10,
    #                     parameters=inputs,
    #                     color=x_color,
    #                 )
    #                 gillespies.plot_walk_fixed(
    #                     ax,
    #                     model,
    #                     "y",
    #                     xmax=1e10,
    #                     parameters=inputs,
    #                     color=y_color,
    #                 )
    #
    #             for i, ax in enumerate(hist_axes):
    #                 odeies.plot_fixed("ode", ax, inputs, axis="x", xmax=100, ymax=100)
    #                 gillespies.plot_hist_fixed(
    #                     ax,
    #                     model,
    #                     "x",
    #                     ymax=1,
    #                     parameters=inputs,
    #                     color=x_color,
    #                 )
    #                 gillespies.plot_hist_fixed(
    #                     ax,
    #                     model,
    #                     "y",
    #                     ymax=1,
    #                     parameters=inputs,
    #                     color=y_color,
    #                 )
    #         else:
    #             gillespies.plot_walk_fixed(
    #                 walk_axes[0], model, "x", xmax=1e10, parameters=inputs
    #             )
    #             gillespies.plot_walk_fixed(
    #                 walk_axes[1], model, "y", xmax=1e10, parameters=inputs
    #             )
    #
    #             gillespies.plot_hist_fixed(
    #                 hist_axes[0], model, "x", ymax=1, parameters=inputs
    #             )
    #
    #     if args.plot_on_phase:
    #         plt.close("all")
    #
    #         fig1 = plt.figure(1, figsize=(6, 6))
    #         ax1 = fig1.add_subplot()
    #
    #         phase_frame = file_frame.loc[("phase", "jesper")]
    #         data_frame["metadata"]["n1"] *= 2
    #         data_frame["metadata"]["n2"] *= 2
    #         phase_filters: list[tuple[str, float]] = [*data_frame["metadata"].items()]
    #
    #         filtered_phase = filter_frames(phase_filters, phase_frame)
    #         phase_name = filtered_phase.reset_index().loc[0, "file_name"]
    #         phase_path = os.path.join(data_env, phase_name)
    #         phase_data = np.load(phase_path)
    #
    #         c1, c2 = phase_data["c1"], phase_data["c2"]
    #         dU, dV = phase_data["dU"], phase_data["dV"]
    #
    #         odeies.plot_phase_space(ax1, c1, c2, dU, dV)
    #         odeies.plot_trajectories(ax1, states[0, :], states[1, :])
    #
    #         ax1.set_xlim(0, 100)
    #         ax1.set_ylim(0, 100)
    #
    #
    case "ode":
        for i, data in enumerate(numpy_datas):
            time = data["time"]
            solutions = data["solutions"]

            axis = walk_axes[i]

            fig1 = plt.figure(1, figsize=(6, 4))
            ax1 = fig1.add_subplot()

            odeies.plot_curves(axis, time, solutions)

            if args.plot_fixedpoints:
                parameters: dict[str, float] = metadatas[0]
                odeies.plot_fixed(
                    data_source,
                    axis,
                    parameters,
                    total_vars=2,
                    xmax=time[-1],
                    ymax=100,
                )

            ax1.set_xlim(left=0, right=time[-1])
            ax1.set_ylim(bottom=0)

            ax1.legend(loc="upper right", fontsize=legend_font_size)

            ax1.set_xlabel("Time")
            ax1.set_ylabel("Value of Variable")

    case "phase":
        data = numpy_datas[0]

        c1, c2 = data["c1"], data["c2"]
        dU, dV = data["dU"], data["dV"]

        c1_null = data["c1_nullcline"]
        c2_null = data["c2_nullcline"]

        odeies.plot_phase_space(ax1, c1, c2, dU, dV)
        odeies.plot_nullclines(ax1, *c2_null, *c1_null)

        if args.plot_fixedpoints:
            try:
                metadata = metadatas[i]
            except:
                pass
            else:
                parameters: dict[str, float] = metadata["metadata"]
                odeies.plot_fixed(data_source, ax1, parameters, xmax=100, ymax=100)

        ax1.set_xlim(left=0, right=100)
        ax1.set_ylim(bottom=0, top=100)

        ax1.legend(loc="upper right", fontsize=legend_font_size)

        ax1.set_xlabel("Number of $c_1$")
        ax1.set_ylabel("Number of $c_2$")

        # file_path = os.path.join(data_env, path_list[0])
        # numpy_data = np.load(file_path)
        #
        # fig1 = plt.figure(1, figsize=(6, 6))
        # ax1 = fig1.add_subplot()

        # fig2 = plt.figure(figsize=(6, 6))
        # ax2 = fig2.add_subplot()

        # level_set = numpy_data["level"]
        # odeies._plot_phase_curve(ax1, *c1_null, label="$c_2$ Nullcline", color="b")
        # odeies._plot_phase_curve(ax1, *level_set, label="Levelset", color="g")

    case _:
        print("Fucked it")
        exit()

# if model == "2S":
#     hist_axes[0].set_title("2D histogram of chemical species distribution")
#
#

for ax in walk_axes:
    ax.legend(loc="upper right", fontsize=legend_font_size)

figs = [plt.figure(i) for i in plt.get_fignums()]
for fig in figs:
    fig.tight_layout()


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
            print("what")
            exit()
            pass

    file_path = os.path.join(presentation_fig_env, file_name)

    ## Just making sure the shit is sorted
    save_image(file_path)

if args.show or not args.save:
    plt.show()
    print("Done Plotting")
    exit()

# print("TESTING MODE")
# exit()


print("Saved or something")

exit()
