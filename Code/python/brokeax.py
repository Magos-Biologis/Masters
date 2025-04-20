#!./.venv/bin/python
import os
from pprint import pprint as pp
from sys import exit

import numpy as np
from numpy import polynomial as npp


from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter
# import sympy as sy

from myPyPlotting.ODEModel import ODEModel
from myPyPlotting.parameter_class import parameter_class

print("")


figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_path = os.path.join(figure_env, "ode")

file_name = "ode_solution"


### Variables
# Growth Constant, same because same cell
k = np.zeros(2)
n = np.zeros(2)
w = np.zeros(2)
q = np.zeros(2)

k[0] = 2
k[1] = 2

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = 100
n[1] = 100

w[0] = 0.035
w[1] = 0.015

w[0] = 0.015
n[1] = 90

q[0] = 0.999
q[1] = 0.8

m_0 = 0.0

omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]

k_1 = k[0] / n[0]
k_2 = k[1] / n[1]

c1_min = w[1] / k_1
c2_min = w[0] / k_2

c1_max = omega_1 / k_1
c2_max = omega_2 / k_2

alpha = 0
c1_0 = n[0] - alpha
c2_0 = alpha

filename_addendum = (
    "_m0"
    + f"{m_0}"
    + "_k1"
    + f"{k[0]}"
    + "_k2"
    + f"{k[1]}"
    + "_n1"
    + f"{n[0]}"
    + "_n2"
    + f"{n[1]}"
    + "_w1"
    + f"{w[0]}"
    + "_w2"
    + f"{w[1]}"
)

file_name += filename_addendum


dt = 0.01
# t_end = 25
t_end = 50
# t_end = 150
t_array = np.arange(0, t_end, dt)
# sol1 = np.zeros((len(t_array), 3))


parameters = parameter_class(2, m_0, k, n, q, w)
init_conds1 = np.array([c1_0, c2_0, m_0])

### Integration

model1 = ODEModel((0, t_end), parameters, init_conds1)
t_array, sol1 = model1.integrate()
solutions = model1.roots()

cq = sum(n) / len(n)


color_list = ["b", "r", "g"]
name_list = [
    "Non-resistant Type",
    "Resistant Type",
    "Antibiotic",
]
dt_name_list = [
    r"$c_1$ fixed point at $c_1^*$",  # {c1_fixed}",
    r"$c_2$ fixed point at $c_2^*$",  # {c2_fixed}",
    # f"$m$ fixed point at {m_fixed}",
]


## Plotting

fig, axes = plt.subplots(2, 1, sharex=True)

plt.rcParams.update(
    {
        "axes.labelsize": 20,
        "axes.titleweight": "bold",
        # "axes.titlecolor": "white",
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        # "xtick.labelcolor": "white",
        # "ytick.labelcolor": "white",
        # "savefig.facecolor": "#c0c0ca",
    }
)


plt.style.use("bmh")
# axes[0].set_facecolor("#c0c0ca")
# axes[1].set_facecolor("#c0c0ca")
#
# axes[0].tick_params("x", color="#c0c0ca")


for i, curve in enumerate(sol1):
    # break
    if i == 2:
        break
    color = color_list[i]
    curve_name = name_list[i]
    dt_curve_name = dt_name_list[i]
    # zero_line = dt_zeroes[i]
    zero_line = solutions[i]

    for ax in axes:
        ax.hlines(
            zero_line,
            0,
            t_end,
            label=dt_curve_name,
            color=color,
            linestyle="dashed",
            linewidth=1,
        )

        ax.plot(t_array, curve, label=curve_name, color=color)


# for ax in axes:


# print()
# exit()

# total = sol1.T[0] + sol1.T[1]
# percent_total = np.divide(sol1.T[0], n[0]) + np.divide(sol1.T[1], n[1])

fixed_point_tick_labels = [
    r"$c_1^*$",
    r"$c_2^*$",
]

axes[0].set_ylim(91, 100)
axes[1].set_ylim(0, 9)

top_y_ticks = np.append([100, 95], solutions[0])
bot_y_ticks = np.append([0, 5], solutions[1])

axes[0].set_yticks(top_y_ticks)
axes[1].set_yticks(bot_y_ticks)


def top_format(val, pos):
    if val == solutions[0]:
        return fixed_point_tick_labels[0]
    else:
        return int(np.round(val, 3))


def bot_format(val, pos):
    if val == solutions[1]:
        return fixed_point_tick_labels[1]
    else:
        return int(np.round(val, 3))


axes[0].yaxis.set_major_formatter(FuncFormatter(top_format))
axes[1].yaxis.set_major_formatter(FuncFormatter(bot_format))


for ax in axes:
    ax.set_xlim(0, t_end)

    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")


axes[1].set_xlabel("Time")
#
axes[0].spines.bottom.set_visible(False)
axes[0].xaxis.tick_top()
axes[0].tick_params("x", colors="white")

axes[1].spines.top.set_visible(False)


d = 0.5  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(
    marker=[(-1, -d), (1, d)],
    markersize=12,
    linestyle="none",
    color="k",
    mec="k",
    mew=1,
    clip_on=False,
)

axes[0].plot([0, 1], [0, 0], transform=axes[0].transAxes, **kwargs)
axes[1].plot([0, 1], [1, 1], transform=axes[1].transAxes, **kwargs)

axes[0].legend(loc="upper right")
plt.tight_layout()

file_path = os.path.join(figure_path, file_name)
# plt.savefig(file_path + ".pdf")
plt.show()
