#!./.venv/bin/python
import os

import numba
from numba import njit

from pylab import *
from matplotlib.backends.backend_pdf import PdfPages


from scipy.stats import mode

# from old.gillespie import gillespie
import gillespie as dg


# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_env = str(os.getenv("FPE_FIGURE_ENV"))


# file_name = "simulated_fpe"
file_name = "multipage_2var_fpe_ssa"

print()

m = 100


b = 100
n = 50


## Simple if else so I can be lazy and not remember to add the naming conventions
if b == n:
    para_version = r"$b=n$"
elif b < n:
    para_version = r"$b<n$"
elif b > n:
    para_version = r"$b>n$"
else:
    para_version = "null"
    print("no numbers?")


# margin = 0.5
alpha = 0
beta = m

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

line_kwargs = {
    "linewidth": 3,
    "alpha": 0.8,
    "linestyle": "-.",
    "color": "black",
}


k = zeros(shape=(5, 2), dtype=float64)

k[0, 0] = k_1 = 1  # 0.55
k[0, 1] = k_m1 = 1  # 0.55


k[1, 0] = k_2 = 1
k[2, 0] = k_3 = 1
k[3, 0] = k_4 = 1
k[4, 0] = k_5 = 1


k[2, 0] *= b
k[3, 0] *= b
k[0, 0] *= n
k[1, 0] *= n

# print(k)
# exit()

x_fixed = -(k[2, 0] - k[0, 0]) / k[0, 1]

k = divide(k, m)


addons = [
    f"num{m}",
    f"n{n}",
    f"b{b}",
    f"k1{k_1}",
    f"k2{k_2}",
    f"k3{k_3}",
    f"k4{k_4}",
    f"k5{k_5}",
]

for addon in addons:
    file_name += "_" + addon


# print(file_name)
# exit()
# file_name += "_" + f"k1{k_1}"
# file_name += "_" + f"k2{k_2}"


# k_2 = 1
# file_name += "_" + f"w1{w_1}"
# file_name += "_" + f"w2{w_2}"
# file_name += "_" + f"n1{n_1}"
# file_name += "_" + f"n2{n_2}"


w_1 = 1
w_2 = 1
n_1 = 100
n_2 = 100


# k = array([k_1, k_2])
k_sols = array([k_1 / (k_1 + k_2), k_2 / (k_1 + k_2)])

# file_name += "_" + f"n{m}"


a_1 = m
a_2 = 0


aj_0 = array([k_1 * a_1, k_2 * a_2])


t_0 = 0
dt = 0


max_val = 0
if k_1 != k_2:
    if m == 1:
        max_val = log(2.022121436749997)
    elif m == 10:
        max_val = log(674763.2054860857)
    elif m == 100:
        max_val = log(6.367225127715182e43)
else:
    pass

x_array = linspace(0, 1, num=500, endpoint=False)[1:]


## Step function
@njit
def step_function(
    steps: int,
    x0: ndarray[tuple[int], dtype[float64]],
    v: list[ndarray[tuple[int], dtype[int_]]],
) -> tuple[ndarray[tuple[int], dtype[float64]], ndarray[tuple[int, int], dtype[int_]]]:
    ## v_j

    ## Other
    gillespie_results = empty((len(x0), steps), dtype=int_)
    time_array = empty(steps, dtype=float64)

    gillespie_results[:, 0] = x0

    for i in range(1, steps):
        x: ndarray[tuple[int], dtype[int_]] = gillespie_results[:, i - 1]
        a_j: ndarray[tuple[int], dtype[float64]] = dg.aj_5_2(x, k)
        j, dt = dg.ssa_event(a_j)

        if j == -1:
            break

        gillespie_results[:, i] = add(gillespie_results[:, i - 1], v[j])
        time_array[i] = time_array[i - 1] + dt

    final_step: int = i
    return time_array[:final_step], gillespie_results[:, :final_step]


# time_array, gillespie_results = step_function(100000, x_0)

# steps = 100
steps = 100_000
# steps = 1_000_000
# steps = 10_000_000
# steps = 100_000_000
# steps = 1_000_000_000


init1 = array([m, 0, m])
init2 = array([0, m, m])

init1 = array([m, 0])
init2 = array([1, 0])
# init2 = array([0, m])

init_conds = [init1, init2]

time_array_1, gillespie_results_1 = step_function(steps, init1, dg.transitions["vj_5_2"])
time_array_2, gillespie_results_2 = step_function(steps, init2, dg.transitions["vj_5_2"])


fig1 = figure(figsize=(5, 2.5))
fig2 = figure(figsize=(5, 2.5))
fig3 = figure(figsize=(5, 5))

ax1 = fig1.add_subplot()
ax2 = fig2.add_subplot()
ax3 = fig3.add_subplot()

figs = [fig1, fig2, fig3]
axes = [ax1, ax2, ax3]


def plot_walk(
    ax,
    time: ndarray[tuple[int], dtype[float64]],
    results: ndarray[tuple[int, int], dtype[int_]],
    color: str,
    xstart: str = "[m,0]",
) -> None:
    ax.step(time, results[0, :], color=color, label=f"Walk of $x$ with start {xstart}")
    ax.step(time, results[1, :], color="g", label=f"Walk of $y$ with start {xstart}")

    ax.set_xlabel("Time", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)

    ax.set_xlim(left=0)

    ax.set_yticks([0, 30, 60, 90, 120, 150])
    ax.set_ylim(bottom=0, top=150)


plot_walk(ax1, time_array_1, gillespie_results_1, "r")
plot_walk(ax2, time_array_2, gillespie_results_2, "b", "[1,0]")


ax3.hist(
    gillespie_results_1[0, :],
    **hist_kwargs,
    label="Start Condition [m,0]",
    color="r",
)
ax3.hist(
    gillespie_results_2[0, :],
    **hist_kwargs,
    label="Start Condition [1,0]",
    color="b",
)

ax3.set_xlim(xlims)
ax3.set_ylim(bottom=0)

ax3.set_title("Distribution of Gene Copy Number\n" + para_version)
ax3.set_xlabel("Count", fontsize=12)
ax3.set_ylabel("Density", fontsize=12)

ax3.vlines([x_fixed], 0, 1, **line_kwargs, label="Analytical solution for $x$")


for ax in axes:
    ax.legend(loc="upper right", fontsize=10)

for fig in figs:
    fig.tight_layout()

file_path = os.path.join(figure_env, "five_var", file_name + "_" + para_version[1:-1])


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
show()

# fig1.savefig(file_path + f"_x0_{init1[0]}_y0_{init1[1]}" + para_version + ".pdf")

exit()
