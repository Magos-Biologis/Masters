#!../.venv/bin/python3
import os

##
# os.sched_setaffinity(0, {3})

# import numba
import numpy as np
from numba import njit

from pylab import *
from matplotlib.backends.backend_pdf import PdfPages


from scipy.stats import mode

# from old.gillespie import gillespie
import gillespie as dg
from gillespie import propensities as dgp
from gillespie import steppers as dgs

import myodestuff as mpp
# from myPyPlotting import ODEModel


# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_env = str(os.getenv("FPE_FIGURE_ENV"))
file_dir = os.path.join(figure_env, "five_var")


# file_name = "simulated_fpe"
file_name = "multipage_2var_fpe_ssa"

print()

m = 100
# margin = 0.5

alpha = 0
beta = m


b = 800
b = 1
n = 100

is_ode = False

ks = np.empty(shape=2)
ns = np.empty(shape=2)
ws = np.empty(shape=2)
qs = np.empty(shape=2)

m_0 = 0

ks[0] = k_1 = 1
ks[1] = k_2 = 1
ns[0] = n_1 = 100
ns[1] = n_2 = 80
ws[0] = w_1 = 0.15
ws[1] = w_2 = 0.15

qs[0] = q_1 = 0.85
qs[1] = q_2 = 0.85


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


init1 = array([m / 2 - 1, 1, m / 2], dtype=int_)[0:2]
init2 = array([1, m / 2 - 1, m / 2], dtype=int_)[0:2]


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


if not is_ode:
    k = zeros(shape=(5, 2), dtype=float64)

    k[0, 0] = k_1 = 1  # 0.55
    k[0, 1] = k_m1 = 1  # 0.55

    k[1, 0] = k_2 = 1
    k[2, 0] = k_3 = 1
    k[3, 0] = k_4 = 1
    k[4, 0] = k_5 = 1

    k[0:2, 0] *= n
    k[2:4, 0] *= b

    # print(k)
    # exit()

    anal_sol = -(k[2, 0] - k[0, 0]) / k[0, 1]

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

    k_sols = array([k_1 / (k_1 + k_2), k_2 / (k_1 + k_2)])

else:
    k = zeros(shape=(4, 2), dtype=float64)

    k[0, :] = ks[:]
    k[0, :] *= u

    k[1, :] = ks
    k[2, :] = ks / ns
    k[3, :] = ws

    addons = [
        f"num{m}",
        f"u{n}",
        f"k1{k_1}",
        f"k2{k_2}",
        f"n1{n_1}",
        f"n2{n_2}",
        f"w1{w_1}",
        f"w2{w_2}",
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


# k = array([k_1, k_2])

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


# time_array, gillespie_results = step_function(100000, x_0)

# steps = 100
# steps = 100_000
# steps = 1_000_000
# steps = 10_000_000
steps = 100_000_000
# steps = 1_000_000_000


init_conds = [init1, init2]

import time

# time_array_1, gillespie_results_1 = steppy.step_function(steps)

model: str = "5_2"


# transitions = dgp.transitions["vj_" + model]

steppy1 = dgs.ssa_stepper(model, init1, k)  # transitions, k)
steppy2 = dgs.ssa_stepper(model, init2, k)  # , transitions

print()

t0_1 = time.time()
time_array_1, gillespie_results_1 = steppy1.step_function(steps)
t1_1 = time.time()

print("Stepper one done \n\tTime taken: ", t1_1 - t0_1)

t0_2 = time.time()
time_array_2, gillespie_results_2 = steppy2.step_function(steps)
t1_2 = time.time()

print("Stepper two done \n\tTime taken: ", t1_2 - t0_2)


exit()


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
    parameters = mpp.parameter_class(2, m_0, ks, ns, qs, ws)
    init_conds1 = np.array([1, 0, 0])
    model1 = mpp.ODEModel((0, time_array_1[-1]), parameters, init_conds1)
    anal_sol = model1.roots()


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
        label=f"Walk of $x$ with start {xstart}",
        **walk_kwargs,
    )
    ax.step(
        time,
        results[1, :],
        color="g",
        label=f"Walk of $y$ with start {xstart}",
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


# if len(init2) == 2:
#     if array_equal(init2, array([1, 0])):
#         start_string = "[1,0]"
#     elif array_equal(init2, array([0, m])):
#         start_string = "[0,m]"
# elif len(init2) == 3:
#     if array_equal(init2, array([1, 0, m])):
#         start_string = "[1,0,m]"
#     elif array_equal(init2, array([0, m, m])):
#         start_string = "[0,m,m]"
#     elif array_equal(init2, array([1, m-1, m])):
#         start_string = "[1,m,m]"
# else:
#     start_string = "undefined"

plot_walk(ax1, time_array_1, gillespie_results_1, "r", f"{init1}")
plot_walk(ax2, time_array_2, gillespie_results_2, "b", f"{init2}")

if is_ode:
    ax1.hlines(anal_sol, 0, time_array_1[-1], colors=["r", "g"])
    ax2.hlines(anal_sol, 0, time_array_2[-1], colors=["b", "g"])

    ax3.vlines(anal_sol[0], 0, 1 / m)
    ax4.vlines(anal_sol[1], 0, 1 / m)

# ax3.hist(
#     gillespie_results_2[0, :],
#     **hist_kwargs,
#     label="Start Condition [0,1]",
#     color="b",
# )

ax3.set_xlim(xlims)
ax3.set_ylim(bottom=0)
ax4.set_xlim(xlims)
ax4.set_ylim(bottom=0)

ax3.set_title("Distribution of Gene Copy Number\n" + para_version)
ax3.set_xlabel("Distribution of the Count of $x$", fontsize=12)
ax3.set_ylabel("Density", fontsize=12)

ax4.set_title("Distribution of Gene Copy Number\n" + para_version)
ax4.set_xlabel("Distribution of the Count of $y$", fontsize=12)
ax4.set_ylabel("Density", fontsize=12)

# ax3.vlines(anal_sol, 0, 1, **line_kwargs, label="Analytical solution for $x$")


for ax in axes:
    ax.legend(loc="upper right", fontsize=10)

for fig in figs:
    fig.tight_layout()

file_path = os.path.join(file_dir, file_name + "_" + para_version[1:-1])


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


# save_image(file_path)
show()

# fig1.savefig(file_path + f"_x0_{init1[0]}_y0_{init1[1]}" + para_version + ".pdf")

exit()
