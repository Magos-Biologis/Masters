#!./.venv/bin/python
import os
from pprint import pprint as pp
from sys import exit

import numpy as np
from matplotlib import pyplot as plt
# import sympy as sy

# Growth Constant, same because same cell
k = 2
k_1 = k_2 = k


# Population cap (purely aesthetic if n₁ = n₂)
n_1 = 100
# n_2 = n_1
n_2 = 99

pop_mean = (n_1 + n_2) / 2
n1_wt = n_1 / pop_mean
n2_wt = n_2 / pop_mean

w_1 = 0.01
w_2 = 0.02

# transition_matrix = sy.Matrix([[1 - w_1, w_1], [w_2, 1 - w_2]])
# init_condition = sy.Matrix([[n_1], [0]])
#
# print(transition_matrix)
# exit()

q_1 = 0.999
q_m = 0.8

c_10 = n_1
c_20 = 0.0
m_0 = 0.00

c_t0 = c_10 + c_20

omega_1 = k_1 - w_1
omega_2 = k_2 - w_2

ku1 = k_1 / n_1
ku2 = k_2 / n_2

w1w2 = w_1 / w_2
w2w1 = w_2 / w_1

# eta_1 = omega_1 - ku1
# eta_2 = omega_2 - ku2

# trace = (i_0 * q_1 - omega_2 - omega_1)
# dete  = (omega_1 - i_0 * q_1) * omega_2 - w_1 * w_2

# disc  = trace**2 - dete

dt = 0.01
t_end = 2000

t_array = np.arange(0, t_end, dt)

steps = np.zeros((len(t_array), 3))
steps[0, :] = n_1, 0, 0
# steps[0, :] = 0, n_2, 0


parameters = [
    k_1,
    k_2,
    w_1,
    w_2,
    q_1,
    q_m,
    n_1,
    n_2,
]


def step_function(xs, dt) -> np.ndarray:
    c_1, c_2, m = xs
    ct = c_1 + c_2

    c1_dot = c_1 * (k_1 * (1 - ct / n_1) - w_1) + w_2 * c_2 - q_1 * m * c_1
    c2_dot = c_2 * (k_2 * (1 - ct / n_2) - w_2) + w_1 * c_1
    m_dot = -q_m * m * c_2 + m_0

    return dt * np.array([c1_dot, c2_dot, m_dot])


for j, _ in enumerate(t_array):
    steps[j, :] += steps[j - 1, :]
    steps[j, :] += step_function(steps[j - 1, :], dt)  # print(steps[j, :])


# ax.plot(steps.T[0], steps.T[1], label="c1/c2")


# c1_fixed = (((q_1 * m_0) / q_m) * (eta_2 / (ku2 - w_1))) / (
#     eta_1 + ((w_2 - ku1) * (ku2 - w_1)) / eta_2
# )
# c2_fixed = 1
# c1_fixed = 1 - w_1
# c2_fixed = w_1
# step_size_quiver = 0.03
# c1_start_pts = np.arange(0, 1, step_size_quiver) * n_1
# c2_start_pts = np.arange(0, 1, step_size_quiver) * n_2
#
# c1_ax, c2_ax = np.meshgrid(c1_start_pts, c2_start_pts)
#
# m_ax = 0.0
# ct_ax = c1_ax + c2_ax
# dc1 = c1_ax * (k_1 * (1 - ct_ax / n_1) - w_1) + w_2 * c2_ax - q_1 * m_ax * c1_ax
# dc2 = c2_ax * (k_2 * (1 - ct_ax / n_2) - w_2) + w_1 * c1_ax


# dc1, dc2 = np.gradient(c1_c2_gradnt)
c10 = n_1  # n_1 * w_2 / w_1  # n_1
c20 = 0  # n_2 * w_1 / w_2  # n_1

eig1 = 1
eig2 = (1 - w_1 - w_2) ** 2
det_of_matrix_2 = (((1 - w_1) - eig1) * ((1 - w_2) - eig2)) - w_1 * w_2

w1_frac = w_1 / (w_1 + w_2)
w2_frac = w_2 / (w_1 + w_2)
n1_frac = n_1 / (n_1 + n_2)
n2_frac = n_2 / (n_1 + n_2)

pop_ratio = n_2 / n_1
pop_norm2 = n_2 / np.abs((n_1 + n_2) / 2)

state_vec1 = c10 * ((1 - w_1) ** 2 + w_1 * w_2) + c20 * w_2 * (2 - (w_1 + w_2))
state_vec2 = c10 * w_1 * (2 - (w_1 + w_2)) + c20 * ((1 - w_2) ** 2 + w_1 * w_2)
# state_vec1 = n_1 * ((1 - w_1) ** 2 + w_1 * w_2) + n_2 * w_2 * (2 - (w_1 + w_2))
# state_vec2 = n_1 * w_1 * (2 - (w_1 + w_2)) + n_2 * ((1 - w_2) ** 2 + w_1 * w_2)

stavec1 = n_1 * ((1 - omega_1 + k_1) ** 2 + w_1 * w_2)
stavec2 = n_1 * ((1 - omega_1 + k_1) * w_1 + w_1 * (1 - omega_2 + k_2 * (n_1 / n_2)))

normed_transition = np.sqrt(state_vec1**2 + state_vec2**2)

if m_0 == 0:
    # c1_fixed = w_1**2 - w_2**2 - 2 * w_1 + 2 * w_2 + 1
    # c2_fixed = -(w_1**2) + 2 * w_1 + (-1 + w_2) ** 2
    # c1_fixed = det_of_matrix**2
    # c2_fixed = 1 - c1_fixed
    c1_fixed = stavec1
    c2_fixed = stavec2
    # c1_fixed *= n_1
    # c2_fixed *= n_2
    m_fixed = 0
else:
    c1_fixed = (1 - w_1) ** 2 + w_1 * w_2
    c2_fixed = (1 - w_1) * w_1 + w_1 * (1 - w_2)
    m_fixed = m_0 / q_m

dt_zeroes = [c1_fixed, c2_fixed, m_fixed]

color_list = ["b", "r", "g"]
name_list = [
    "Non-resistant Type",
    "Resistant Type",
    "Antibiotic",
]
dt_name_list = [
    f"$c_1$ fixed point at {c1_fixed}",
    f"$c_2$ fixed point at {c2_fixed}",
    f"$m$ fixed point at {m_fixed}",
]

fig, ax = plt.subplots()

# Defining color
# n_c = -((w_1 + w_2) ** (-1))
# colorings = np.sqrt(((dc1 - n_c) / 2) * 2 + ((dc2 - n_c) / 2) * 2)

# Creating plot
# ax.quiver(c1_ax, c2_ax, dc1, dc2)
# ax.xaxis.set_ticks([])
# ax.yaxis.set_ticks([])
# ax.set_aspect("equal")
#
# plt.show()

for i, curve in enumerate(steps.T):
    color = color_list[i]
    curve_name = name_list[i]
    dt_curve_name = dt_name_list[i]

    ax.plot(t_array, curve, label=curve_name, color=color)
    ax.hlines(
        dt_zeroes,
        0,
        t_end,
        label=dt_curve_name,
        color=color,
        linestyle="dashed",
        linewidth=1,
    )

# plt.yscale("power")
plt.legend()
plt.show()


# plt.savefig()
