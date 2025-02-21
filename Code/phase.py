#!./.venv/bin/python
import os
from pprint import pprint as pp
from re import I
from sys import exit

from matplotlib import pyplot as plt
import numpy as np
from numpy import sum

# import sympy as sy

print("")


figure_path = str(os.getenv("THESIS_FIGURE_PATH"))
# figure_path = "./figs"


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

w[0] = 0.015
w[1] = 0.035


n[1] = 90
w[1] = 0.015


q[0] = 0.999
q[1] = qm = 0.8


figure_file = os.path.join(
    figure_path,
    "nomed_phaseplane"
    + "_n1"
    + f"{n[0]}"
    + "_n2"
    + f"{n[1]}"
    + "_w1"
    + f"{w[0]}"
    + "_w2"
    + f"{w[1]}",
)


pop_norm = np.sqrt(n[0] ** 2 + n[1] ** 2)
n1_wt = n[0] / pop_norm
n2_wt = n[1] / pop_norm

ct_1 = n[0] * ((k[0] - 1) / k[0])
ct_2 = n[1] * ((k[1] - 1) / k[1])

omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]

k_1 = k[0] / n[0]
k_2 = k[1] / n[1]

# w1w2 = w[0] / w[1]
# w2w1 = w[1] / w[0]

alpha = 1
c1_0 = n[0] - alpha
c2_0 = alpha
m_0 = 0.0

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


## Plotting

x_lims = 0, 100
y_lims = 0, 100

# x_lims = 0, 1
# y_lims = 0, 1

c1s = np.arange(0.1, 125, 0.5)
c2s = np.arange(0.1, 125, 0.5)
# ct = c1s + c2s
# ct = n[0]

c1, c2 = np.meshgrid(c1s, c2s)


def dc1(c1, c2):
    return (omega_1 - k_1 * (c1 + c2)) * c1 + w[1] * c2


def dc2(c1, c2):
    return (omega_2 - k_2 * (c1 + c2)) * c2 + w[0] * c1


# ratio = np.divide(dc1, dc2)
# slope1 = -w[1] / (omega_1 - k_1 * n[0])

# slope1 = c1s * ((w[1]) / (omega_1 - k_1 * cq))
# slope2 = c2s * ((w[0]) / (omega_2 - k_2 * cq))

dc1_c1 = w[1]
dc1_c2 = -(omega_1 - k_1 * cq)

dc2_c1 = -(omega_2 - k_2 * cq)
dc2_c2 = w[0]

f1 = np.array([dc1_c1, dc1_c2])
f2 = np.array([dc2_c1, dc2_c2])


slope1 = dc1_c1 / dc1_c2
slope2 = dc2_c1 / dc2_c2

# average_slope = np.divide(slope1, slope2)

dc1_U = (omega_1 - k_1 * (c1 + c2)) * c1 + w[1] * c2
dc2_V = (omega_2 - k_2 * (c1 + c2)) * c2 + w[0] * c1


dt = 0.01
t_end = 250
t_array = np.arange(0, t_end, dt)

parameters = [k, w, q, n]


cq = sum(n) / len(n)


def step_function(xs, dt) -> np.ndarray:
    c_1, c_2, m = xs
    ct = c_1 + c_2
    # ct = cq

    c1_dot = (k[0] * (1 - ct / n[0]) - w[0]) * c_1 + w[1] * c_2 - q[0] * m * c_1
    c2_dot = (k[1] * (1 - ct / n[1]) - w[1]) * c_2 + w[0] * c_1
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
lw = 4 * speed / speed.max()


fig, ax = plt.subplots()


ax.streamplot(
    c1,
    c2,
    dc1_U,
    dc2_V,
    density=1.5,
    cmap="viridis",
    linewidth=lw,
    # arrowsize=0,
    arrowstyle="->",
    # broken_streamlines=False,
)


ax.set_xlabel(r"$c_1$")
ax.set_ylabel(r"$c_2$")

ax.set_xlim(*x_lims)
ax.set_ylim(*y_lims)

c1f = omega_1 / k_1
c2f_b = (omega_1 - k_1 * c1f) / (w[1] - c1f)
# c1f = -omega_1/k_1

# ax.axline(f1, f2, color="orange")

# ax.axline((n[0], 0), (0, n[1]), color="blue", alpha=0.25)
# ax.axline(f2, f1, color="cyan")
# ax.plot(c1s, c1s * slope1)
# ax.plot(c2s, c2s * slope2)


# ratio = np.divide(steps[0, :], steps[1, :])
# ax.plot(t_array, ratio, label="progression")

sol_alpha = 0.6
colors = ["red", "blue"]


def example_plot():
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


# example_plot()


c2_slope = (-omega_1 + k_1 * c1s) / (1 - k_1)
# c2_slope = (-k_1) / (1 - k_1)

# ax.axline((0, omega_1 / (1 - k_1)), slope=c2_slope, color="k")
ax.plot(c2_slope, c2s)
dc1_f = dc1(*f1)
dc2_f = dc2(*f2)
dcf = np.array([dc1_f, dc2_f])
pd = dcf / sum(dcf) * 100

# dcf = np.array([dc1_f, dc2_f])
# p0 = (dcf / sum(dcf)) * n
# p2 = sum

## Nullcline


def c1_sol(c1) -> float:
    num = -(c1 * (omega_1 - k_1 * c1))
    den = w[1] - k_1 * c1
    return num / den


def c2_sol(c2) -> float:
    num = -(c2 * (omega_2 - k_2 * c2))
    den = w[0] - k_2 * c2
    return num / den


c1_min = w[1] / k_1
c2_min = w[0] / k_2

c1_sol_array = np.arange(c1_min + 0.01, 100, 0.2)
c2_sol_array = np.arange(c2_min + 0.01, 100, 0.2)

null_alpha = 0.3
null_width = 3
# colors = ["k", "k"]

ax.plot(
    c1_sol_array,
    c1_sol(c1_sol_array),
    label=r"$c_1$ Nullcline",
    color=colors[0],
    alpha=null_alpha,
    linewidth=null_width,
)
ax.plot(
    c2_sol(c2_sol_array),
    c2_sol_array,
    label=r"$c_2$ Nullcline",
    color=colors[1],
    alpha=null_alpha,
    linewidth=null_width,
)

c2_null = (omega_2 * w[1] - omega_1 * omega_2) / (k_1 * (w[1] - omega_1))

# ax.scatter(c1_sol(t_array[1:]), c2_sol(t_array[1:]))
# ax.scatter(c1_sol(n[0]), c2_sol(n[1]))


# ax.vlines([0], *x_lims, color=colors[0])
# ax.hlines([0], *y_lims, color=colors[1])

# c1_min = n[0] * w[1] / k[0]
# c2_min = n[1] * w[0] / k[1]


# ax.vlines([c1_min], *x_lims, color=colors[0], linestyles="--", alpha=null_alpha)
# ax.hlines([c2_min], *y_lims, color=colors[1], linestyles="--", alpha=null_alpha)

print(f1, "\n", f2)
print(pd)
print("\n")
print(c1_min, "  ", c2_min)


plt.legend()
plt.tight_layout()

plt.savefig(figure_file + ".pdf", format="pdf")
plt.show()
