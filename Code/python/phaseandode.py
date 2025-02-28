#!./.venv/bin/python
import os
from pprint import pprint as pp
from re import I
from sys import exit

from matplotlib import pyplot as plt
from matplotlib.figure import figaspect
from matplotlib.lines import lineStyles
from matplotlib.pyplot import FuncFormatter, xlim


import numpy as np
from numpy import polynomial as npp

# from numpy.polynomial import Polynomial
from numpy import sum


# import sympy as sy

print("")


figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
ode_path = os.path.join(figure_env, "ode")

# figure_path = "./figs"


### Variables
# Growth Constant, same because same cell
k = np.zeros(2)
n = np.zeros(2)
w = np.zeros(2)
q = np.zeros(2)

k[0] = 1
k[1] = 1

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = 90
n[1] = 50

w[0] = 0.035
w[1] = 0.015


# w[0] = 0.015
# w[1] = 0.14
# n[0] = 55
# n[1] = 85

q[0] = 0.999
q[1] = qm = 0.8


alpha = 1
c1_0 = n[0] - alpha
c2_0 = alpha
m_0 = 15.0
# m_0 = 0.0

addendum = (
    "_m0"
    + f"{m_0}"
    + "_n1"
    + f"{n[0]}"
    + "_n2"
    + f"{n[1]}"
    + "_w1"
    + f"{w[0]}"
    + "_w2"
    + f"{w[1]}"
)

phase_figure_file = os.path.join(figure_env, "nomed_phaseplane")
ode_figure_file = os.path.join(ode_path, "ode_solution")

phase_figure_file += addendum
ode_figure_file += addendum


pop_norm = np.sqrt(n[0] ** 2 + n[1] ** 2)
n1_wt = n[0] / pop_norm
n2_wt = n[1] / pop_norm

ct_1 = n[0] * ((k[0] - 1) / k[0])
ct_2 = n[1] * ((k[1] - 1) / k[1])

omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]

k_1 = k[0] / n[0]
k_2 = k[1] / n[1]

c1_min = w[1] / k_1
c2_min = w[0] / k_2

c1_max = omega_1 / k_1
c2_max = omega_2 / k_2

w_1 = w[0]
w_2 = w[1]

cq = sum(n) / len(n)

M = np.array(
    [
        [(omega_1 - k_1 * cq), w[1]],
        [w[0], (omega_2 - k_2 * cq)],
    ]
)


trace = M[0, 0] + M[1, 1]
det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]

eig1 = (-trace + np.sqrt(trace**2 - 4 * det)) / 2
eig2 = (-trace - np.sqrt(trace**2 - 4 * det)) / 2

char_M1 = M - eig1 * np.identity(2)
char_M2 = M - eig2 * np.identity(2)

state_vec = np.array([w[1] / sum(w), w[0] / sum(w)]) * n


## Polynomial


c1_coeffs = [
    n[1] * w[1],
    omega_1 - n[1] * (w[1] / n[0] + k_1),
    k_1 * ((n[1] / n[0]) - 1),
]
c2_coeffs = [
    n[0] * w[0],
    omega_2 - n[0] * (w[0] / n[1] + k_2),
    k_2 * ((n[0] / n[1]) - 1),
]

c1_poly = npp.Polynomial(
    coef=c1_coeffs,
    symbol="_c1",
)
c2_poly = npp.Polynomial(
    coef=c2_coeffs,
    symbol="_c2",
)


c1_roots = c1_poly.roots()
c2_roots = c2_poly.roots()

c1_root = [root for root in c1_roots if c1_min < root < c1_max]
c2_root = [root for root in c2_roots if c2_min < root < c2_max]


c1s = np.arange(0.1, 125, 0.5)
c2s = np.arange(0.1, 125, 0.5)

c1, c2 = np.meshgrid(c1s, c2s)

m = np.linspace(0.1, m_0 + k[1], len(c1s))


def dc1(c1, c2, m):
    return c1 * (omega_1 - k_1 * c1) + c2 * (w_2 - k_1 * c1) - qm * m * c1


def dc2(c1, c2, m):
    return c2 * (omega_2 - k_2 * c2) + c1 * (w_1 - k_2 * c2)


def dm(c1, c2, m):
    return m_0 - qm * m * c2


dc1_U = dc1(c1, c2, m)
dc2_V = dc2(c1, c2, m)
# dm_V  = (omega_2 - k_2 * (c1 + c2)) * c2 + w[0] * c1


dt = 0.01
t_end = 30
t_array = np.arange(0, t_end, dt)

parameters = [k, w, q, n]


cq = sum(n) / len(n)


def step_function(xs, dt) -> np.ndarray:
    c_1, c_2, m = xs
    ct = c_1 + c_2
    # ct = cq

    c1_dot = (k[0] * (1 - ct / n[0]) - w_1) * c_1 + w_2 * c_2 - q[0] * m * c_1
    c2_dot = (k[1] * (1 - ct / n[1]) - w_2) * c_2 + w_1 * c_1
    m_dot = -q[1] * m * c_2 + m_0

    return dt * np.array([c1_dot, c2_dot, m_dot])


def integrate(initial, t_array):
    steps = np.zeros((len(t_array), 3))
    steps[0, :] = initial
    for j, _ in enumerate(t_array):
        steps[j, :] += steps[j - 1, :]
        steps[j, :] += step_function(steps[j - 1, :], dt)  # print(steps[j, :])

    return t_array, steps


init_conds1 = [c1_0, c2_0, m_0]
_, outcome1 = integrate(init_conds1, t_array)

init_conds2 = [20, 30, m_0]
_, outcome2 = integrate(init_conds2, t_array)

speed = np.sqrt(dc1_U**2 + dc2_V**2)
lw = 5 * speed / speed.max()
lw = 1


plt.rcParams.update(
    {
        "axes.labelsize": 20,
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        "axes.titleweight": "bold",
    }
)


### Plotting

x_lims = 0, 100
y_lims = 0, 100

colors = ["blue", "red", "green"]
name_list = [r"$c_1$", r"$c_2$", r"m"]
solutions = [c1_root, c2_root]


w, h = figaspect(1 / 3)
fig, axes = plt.subplots(1, 2, figsize=(w, h))


for i, curve in enumerate(outcome1.T):
    # break
    # if i == 2:
    #     break
    color = colors[i]
    curve_name = name_list[i]

    # dt_curve_name = dt_name_list[i]
    # zero_line = solutions[i]

    # axes[1].hlines(
    #     zero_line,
    #     0,
    #     t_end,
    #     color=color,
    #     linestyle="dashed",
    #     linewidth=1,
    # )

    axes[1].plot(t_array, curve, label=curve_name, color=color)


total = outcome1.T[0] + outcome1.T[1]
percent_total = np.divide(outcome1.T[0], n[0]) + np.divide(outcome1.T[1], n[1])

fixed_point_tick_labels = [
    "\t" + r"$c_1^*$",
    r"$c_2^*$",
]

axes[1].set_ylim(0, np.max(n))
axes[1].set_xlim(0, t_end)

# new_y_ticks = np.append(axes[1].get_yticks(), solutions)
new_y_ticks = axes[1].get_yticks()
# new_y_ticks = np.append([0, 25, 50, 75, 100], solutions)
axes[1].set_yticks(new_y_ticks)


def fixed_point_format(val, pos):
    if val == solutions[0]:
        return fixed_point_tick_labels[0]
    elif val == solutions[1]:
        return fixed_point_tick_labels[1]
    else:
        return int(np.round(val, 3))


axes[1].yaxis.set_major_formatter(FuncFormatter(fixed_point_format))

axes[1].yaxis.tick_right()
axes[1].yaxis.set_label_position("right")

axes[1].set_xlabel("Time")
axes[1].set_ylabel("Count")

axes[1].hlines([m_0], *x_lims, color=colors[2], linestyle="--", alpha=0.7, linewidth=0.5)

axes[1].legend(loc="upper right")


axes[0].streamplot(
    c1,
    c2,
    dc1_U,
    dc2_V,
    density=1.5,
    color="gray",
    linewidth=lw,
    arrowstyle="->",
    # l=0.7,
    # cmap="viridis",
    # arrowsize=0,
    # broken_streamlines=False,
)


axes[0].set_xlabel(r"$c_1$")
axes[0].set_ylabel(r"$c_2$")

axes[0].set_xlim(*x_lims)
axes[0].set_ylim(*y_lims)


# c1f = omega_1 / k_1
# c2f_b = (omega_1 - k_1 * c1f) / (w[1] - c1f)


sol_alpha = 0.6


def example_plot(ax):
    ax.plot(
        outcome1[:, 0],
        outcome1[:, 1],
        color=colors[0],
        label=r"Init conds, $c_1="
        + f"{int(init_conds1[0])}"
        + "$, $c_2="
        + f"{int(init_conds1[1])}"
        + "$",
        alpha=sol_alpha,
    )
    ax.scatter(outcome1[0, 0], outcome1[0, 1], color=colors[0], alpha=sol_alpha)

    ax.plot(
        outcome2[:, 0],
        outcome2[:, 1],
        color=colors[1],
        label=r"Init conds, $c_1="
        + f"{int(init_conds2[0])}"
        + "$, $c_2="
        + f"{int(init_conds2[1])}"
        + "$",
        alpha=sol_alpha,
    )
    ax.scatter(outcome2[0, 0], outcome2[0, 1], color=colors[1], alpha=sol_alpha)


## Nullcline

# axes[0].plot([0, n[1] + m_0], [n[0] - m_0, 0])


# def c1_sol(c1) -> float:
#     num = -c1 * (omega_1 - k_1 * c1)
#     den = w_2 - k_1 * c1
#     return num / den


def c1_sol(c1, m) -> float:
    num = c1 * (qm * m - (omega_1 - k_1 * c1))
    den = w_2 - k_1 * c1
    return num / den


def c2_sol(c2) -> float:
    num = -c2 * (omega_2 - k_2 * c2)
    den = w_1 - k_2 * c2
    return num / den


def m_sol(c2) -> float:
    num = m_0
    den = qm * c2
    return num / den


# def c2_sol1(m) -> float:
#     num = m_0
#     den = qm * m
#     return num / den


#
def m_sol2(c1, c2) -> float:
    num = c1 * (omega_1 - k_1 * c1) + c2 * (w_2 - k_1 * c1)
    den = q[1]
    return num / den


null_alpha = 0.3
null_width = 1


axes[0].plot(
    c2_sol(c2s),
    c2s,
    label=r"$c_1$ Nullcline",
    color=colors[0],
    alpha=null_alpha,
    linewidth=null_width,
)
axes[0].plot(
    c1s,
    c1_sol(c1s, m),
    label=r"$c_2$ Nullcline",
    color=colors[1],
    alpha=null_alpha,
    linewidth=null_width,
)
# axes[0].plot(
#     m,
#     c2_sol1(m),
#     label=r"$c_2$ Nullcline",
#     # color=colors[1],
#     # alpha=null_alpha,
#     # linewidth=null_width,
# )


def level_c2(c1):
    result = n[1] * (1 - c1 / n[0])
    return result


# axes[0].plot(c1s, level_c2(c1s))
axes[0].plot(outcome1.T[0], outcome1.T[1], color=colors[0])
axes[0].plot(outcome2.T[0], outcome2.T[1], color=colors[1])


# c1_sol_array = np.arange(c1_min + 0.01, 100, 3)
# c2_sol_array = np.arange(c2_min + 0.01, 100, 3)
# t_sol_arry = np.arange(1, 100, 0.5)
# for i, c10 in enumerate(c1_sol_array):
#     for j, c20 in enumerate(c2_sol_array):
#         init_conds = [c10, c20, m_0]
#         _, outcome = integrate(init_conds, t_sol_arry)
#         axes[0].plot(outcome.T[0], outcome.T[1], alpha=0.2)


def plot_null(ax):
    ax.vlines(
        c1_root,
        *x_lims,
        color=colors[0],
        linestyle="--",
        alpha=null_alpha,
        zorder=11,
    )
    ax.hlines(
        c2_root,
        *y_lims,
        color=colors[1],
        linestyle="--",
        alpha=null_alpha,
        zorder=11,
    )

    ax.legend(framealpha=1, loc="upper center")


c2_test = n[1] * (1 - c1s / n[0])

origin = np.array([[0, 0], [0, 0]])
extremes = np.array([[0, n[0]], [n[1], 0]])

points = [[], []]
margin = 1 / c1.shape[0]
for i in range(c1.shape[0]):
    for j in range(c1.shape[1]):
        c1_norm = c1[i, j] / n[0]
        c2_norm = c2[i, j] / n[1]

        c_sum = c1_norm + c2_norm
        # print(c_sum)

        if 1 - margin <= c_sum <= 1 + margin:
            points[0].append(c1[i, j])
            points[1].append(c2[i, j])
        else:
            continue


# ax1.yaxis.set_ticklabels()

# solutions =


# new_x_ticks = np.append(axes[0].get_xticks(), c1_root)
# new_y_ticks = np.append(axes[0].get_yticks(), c2_root)

new_x_ticks = axes[0].get_xticks()
new_y_ticks = axes[0].get_yticks()

axes[0].set_xticks(new_x_ticks)
axes[0].set_yticks(new_y_ticks)


fixed_point_tick_labels = [r"$c_1^*$", r"$c_2^*$"]
# fixed_point_tick_labels = [r"$c_1^* =" + f"{round(c1_root[0], 2)}" + r"$", r"$c_2^* =" + f"{round(c2_root[0], 2)}" + r"$",]

# plot_null(axes[0])


def x_format(val, pos):
    if val == c1_root:
        return fixed_point_tick_labels[0]
    else:
        return int(val)


def y_format(val, pos):
    if val == c2_root:
        return fixed_point_tick_labels[1]
    else:
        return int(val)


# print((outcome1.T[0][-1] + outcome1.T[1][-1]) / n[0])
# print((outcome1.T[0][-1] + outcome1.T[1][-1]) / n[1])
# print(outcome1.T[0][-1] / n[0] + outcome1.T[1][-1] / n[1])
# print(outcome1.T[0][-1] / n[1] + outcome1.T[1][-1] / n[0])

# print(c2_root)


axes[0].xaxis.set_major_formatter(FuncFormatter(x_format))
axes[0].yaxis.set_major_formatter(FuncFormatter(y_format))


# plt.tight_layout()


# axes[0].savefig(phase_figure_file + ".pdf", format="pdf")
# axes[1].savefig(ode_figure_file + ".pdf", format="pdf")
plt.show()
# axes[0].imshow()
