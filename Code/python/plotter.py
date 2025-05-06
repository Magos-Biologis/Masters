import os
import re
from pprint import pprint as pp

from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

import myodestuff as mpp

import numpy as np
import pandas as pd

## regex


# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_env = str(os.getenv("FPE_FIGURE_ENV"))
data_env = str(os.getenv("THESIS_DATA_PATH"))

ode_dir = os.path.join(figure_env, "ode")
fpe_dir = os.path.join(figure_env, "fpe")

file_dir = ode_dir

# file_name = "multipage_2var_fpe_ssa"


# file_name = "simulated_fpe"


# is_ode = True
model: str = "ode_3_2"
# model: None = None


var_count = 3

m = 100
# margin = 0.5

alpha = 0
beta = m

re_pattern = re.compile(
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
    m = re_pattern.search(fn)
    if m is not None:
        rec = m.groupdict()
        rec["file_name"] = fn
        records.append(rec)

## So we assign it to a dataframe for easy access
raw_frame = pd.DataFrame(records)

## I need to update the metadata column, to be a column of dict type.
## We do this by updating the existing elements in the df
### First we seperate variable ratio metadata from non-that
raw_frame["ratio"] = (
    raw_frame["metadata"]
    .str.split("R")
    .array.map(lambda x: x[1] if len(x) > 1 else np.nan)
)

## Using numpy's 'not' operator, ~, we can invert the truth and locate the relevant spots
## so as to remove the ratio from the other metadata
truths = ~raw_frame["ratio"].isna()
raw_frame.loc[truths, "metadata"] = (
    raw_frame.loc[truths, "metadata"].str.split("R").array.map(lambda x: x[0] + x[2])
)
### Making the remaining string in an array of elements for further processing
raw_frame.loc[:, "metadata"] = raw_frame["metadata"].str.split("_")
raw_frame["size"] = raw_frame["metadata"].map(lambda l: l[0].replace("num", ""))
raw_frame.loc[:, "metadata"] = raw_frame["metadata"].map(lambda l: l[1:-1])


## From the list, we then count the number of times each model shows up,
## assigning the number of each step accordingly.
## From this, we make a multi-index for the dataframe so as to have an
## easier time choosing which file to plot
raw_frame["count_index"] = raw_frame.groupby("model").cumcount()
file_frame = raw_frame.set_index(["model", "count_index"]).sort_index()


file_options = file_frame.loc[model]
file_choice = file_options.loc[1:2].reset_index()

# pp(file_choice["initcond"])
# exit()

# print(file_choice)
# exit()

data_file_1 = os.path.join(data_env, file_choice.loc[0, "file_name"])
data_file_2 = os.path.join(data_env, file_choice.loc[1, "file_name"])


check_ode = re.compile(r"ode").search(model)  # file_choice["model"])
if check_ode is None:
    is_ode = False
else:
    is_ode = True


# pp(test["time"])
# pp(states)
# exit()

boxes = arange(alpha, beta + 1, 1, dtype=np.int_)
xlims = (alpha, beta)
hist_kwargs = {
    "bins": boxes,
    "density": True,
    "edgecolor": "black",
    "alpha": 0.6,
    "align": "mid",
    # "normed": True,
}

curve_kwargs = {
    "linewidth": 4,
    "alpha": 0.7,
}

walk_kwargs = {
    "linewidth": 1,
    "alpha": 0.7,
}

line_kwargs = {
    "linewidth": 3,
    "alpha": 0.8,
    "linestyle": "-.",
    "color": "black",
    "zorder": 11,
}


fig1 = figure(figsize=(5, 2.5))
fig2 = figure(figsize=(5, 2.5))
fig3 = figure(figsize=(5, 5))
fig4 = figure(figsize=(5, 5))

ax1 = fig1.add_subplot()
ax2 = fig2.add_subplot()
ax3 = fig3.add_subplot()
ax4 = fig4.add_subplot()

figs = [fig1, fig2, fig3, fig4]
axes = [ax1, ax2, ax3, ax4]


if is_ode:
    # parameters = mpp.parameter_class(2, m_0, ks, ns, qs, ws)
    # init_conds1 = np.array([1, 0, 0])
    # model1 = mpp.ODEModel((0, time[-1]), parameters, init_conds1)
    # anal_sol = model1.roots()

    x_name = "c_1"
    y_name = "c_2"
else:
    x_name = "x"
    y_name = "y"


def plot_walk(
    ax,
    time: ndarray[tuple[int], dtype[float64]],
    results: ndarray[tuple[int, int], dtype[int_]],
    color: str,
    xstart: str = "[m,0]",
) -> None:
    ax.step(
        time,
        results[0, :],
        color=color,
        label=f"Walk of ${x_name}$ with start {xstart}",
        **walk_kwargs,
    )
    ax.step(
        time,
        results[1, :],
        color="g",
        label=f"Walk of ${y_name}$ with start {xstart}",
        **walk_kwargs,
    )

    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)

    ax.set_xlim(left=0)

    ax.set_yticks([0, 30, 60, 90, 120, 150])
    ax.set_ylim(bottom=0, top=beta)

    ax3.hist(
        results[0, :],
        **hist_kwargs,
        label=f"Start Condition {xstart}",
        color=color,
    )
    ax4.hist(
        results[1, :],
        **hist_kwargs,
        label=f"Start Condition {xstart}",
        color=color,
    )


numpy_data = np.load(data_file_1)
time, states = numpy_data["time"], numpy_data["states"]
plot_walk(ax1, time, states, "r", f"{file_choice.loc[0, 'initcond']}")

numpy_data = np.load(data_file_2)
time, states = numpy_data["time"], numpy_data["states"]
plot_walk(ax2, time, states, "b", f"{file_choice.loc[1, 'initcond']}")


file_name = "multiplot"
for item in file_choice.loc[0, "metadata"]:
    file_name += "_" + item

# print(file_name)
# exit()
if is_ode:
    # ax1.hlines(anal_sol, 0, time_array_1[-1], colors=["r", "g"])
    # ax2.hlines(anal_sol, 0, time_array_2[-1], colors=["b", "g"])
    #
    # ax3.vlines(anal_sol[0], 0, 1 / m)
    # ax4.vlines(anal_sol[1], 0, 1 / m)

    ax3.set_title("Cell Population Densities")
    ax4.set_title("Cell Population Densities")

    file_path = os.path.join(file_dir, file_name)

else:
    # para_match = re.compile(r"R(?P<ratio>b[<=>]n)R").search(file_choice["metadata"])
    # para_version = para_match.groupdict()["ratio"]  ## This is scuffed, but it works
    para_version = file_choice["ratio"]

    ax3.set_title("Distribution of Gene Copy Number\n" + para_version)
    ax4.set_title("Distribution of Gene Copy Number\n" + para_version)

    file_path = os.path.join(file_dir, file_name + "_" + para_version[1:-1])

# ax3.hist(
#     gillespie_results_2[0, :],
#     **hist_kwargs,
#     label="Start Condition [0,1]",
#     color="b",
# )plo

ax3.set_xlabel(f"Distribution of the Count of ${x_name}$", fontsize=12)
ax4.set_xlabel(f"Distribution of the Count of ${y_name}$", fontsize=12)

ax3.set_xlim(xlims)
ax3.set_ylim(bottom=0)
ax4.set_xlim(xlims)
ax4.set_ylim(bottom=0)


ax3.set_ylabel("Density", fontsize=12)
ax4.set_ylabel("Density", fontsize=12)

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


save_image(file_path)
# show()

# fig1.savefig(file_path + f"_x0_{init1[0]}_y0_{init1[1]}" + para_version + ".pdf")

exit()
