#!./.venv/bin/python
import os
from pprint import pprint as pp
from sys import exit

import numpy as np
from matplotlib import pyplot as plt
# import sympy as sy


figure_path = str(os.getenv("THESIS_FIGURE_PATH"))


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


w[0] = 0.01
w[1] = 0.03

q[0] = 0.999
q[1] = 0.8


pop_norm = np.sqrt(n[0] ** 2 + n[1] ** 2)
n1_wt = n[0] / pop_norm
n2_wt = n[1] / pop_norm

# transition_matrix = sy.Matrix([[1 - w[0], w[0]], [w[1], 1 - w[1]]])
# init_condition = sy.Matrix([[n_1], [0]])
#
# print(transition_matrix)
# exit()

ct_1 = n[0] * ((k[0] - 1) / k[0])
ct_2 = n[1] * ((k[1] - 1) / k[1])

# c_t0 = c_10 + c_20

omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]


k_1 = k[0] / n[0]
k_2 = k[1] / n[1]

w1w2 = w[0] / w[1]
w2w1 = w[1] / w[0]

alpha = 1
c_10 = n[0] - alpha
c_20 = alpha
m_0 = 0


# eta_1 = omega_1 - ku1
# eta_2 = omega_2 - ku2
# trace = (i_0 * q_1 - omega_2 - omega_1)
# dete  = (omega_1 - i_0 * q_1) * omega_2 - w[0] * w[1]
# disc  = trace**2 - dete

dt = 0.01
t_end = 250
t_array = np.arange(0, t_end, dt)
steps = np.zeros((len(t_array), 3))

parameters = [k, w, q, n]

c_tt = 2 * (n[0] * n[1]) / (n[0] + n[1])
# c_tt = 2 * (1 / n[0] + 1 / n[1])
# c_tt = n[0]


def step_function(xs, dt) -> np.ndarray:
    c_1, c_2, m = xs
    ct = c_1 + c_2
    # ct = c_tt

    c1_dot = (k[0] * (1 - ct / n[0]) - w[0]) * c_1 + w[1] * c_2 - q[0] * m * c_1
    c2_dot = (k[1] * (1 - ct / n[1]) - w[1]) * c_2 + w[0] * c_1
    m_dot = -q[1] * m * c_2 + m_0

    return dt * np.array([c1_dot, c2_dot, m_dot])


c10 = n[0]  # n_1 * w[1] / w[0]  # n_1
c20 = n[1]  # n_2 * w[0] / w[1]  # n_1

w1_frac = w[0] / (w[0] + w[1])
w2_frac = w[1] / (w[0] + w[1])
n1_frac = n[0] / (n[0] + n[1])
n2_frac = n[1] / (n[0] + n[1])

pop_ratio = n[1] / n[0]
pop_norm2 = n[1] / np.abs((n[0] + n[1]) / 2)

stavec1 = n[0] * ((1 - omega_1 + k[0]) ** 2 + w[0] * w[1])
stavec2 = n[0] * (
    (1 - omega_1 + k[0]) * w[0] + w[0] * (1 - omega_2 + k[1] * (n[0] / n[1]))
)


print("")
# if m_0 == 0:
# m_trace = (omega_1 - ku1 * ct_1) + (omega_2 - ku2 * ct_2)
# m_det = (omega_1 - ku1 * ct_1) * (omega_2 - ku2 * ct_2) - w[0] * w[1]
# m_disc = m_det**2 - 4 * m_trace
# c1_fixed = (m_trace + np.sqrt(m_disc)) / 2
# c2_fixed = (m_trace - np.sqrt(m_disc)) / 2
# c1_fixed *= n_1
# c2_fixed *= n_2
# if m_det <= 0:
#     print("Determinant less than zero")
# elif m_det >= 0:
#     print("Determinant greater than zero")
# if m_trace <= 0:
#     print("Trace less than zero")
# elif m_trace >= 0:
#     print("Trace greater than zero")
# if m_disc == 0:
#     print("Discriminant less than zero")
# elif m_disc <= 0:
#     print("Discriminant less than zero")
# elif m_disc >= 0:
#     print("Discriminant greater than zero")
# print(m_disc)
# c1_fixed = w[1] * n_1
# c2_fixed = n_2 / (1 - (w[0] + w[1]))
# c2_fixed = w[0] / (w[0] + w[1])

# c1_fixed = w[1] / (w[0] + w[1])
# c1_fixed = (w[1] / w[0]) * ((n[0] * n[1]) / (n[0] + n[1]))

c_T = 2 * (n[0] ** (-1) + n[1] ** (-1))
n1n2 = n[0] * n[1] / (n[0] + n[1])


def quadratic(a, b, c):
    disc = np.sqrt(b**2 - 4 * a * c)
    out1 = (-b + disc) / (2 * a)
    out2 = (-b - disc) / (2 * a)

    return out1, out2


a = 1
b = 2 * k[0] * c_T - omega_1 - omega_2
c = (
    c_T**2 * k[0] ** 2
    + 2 * (-omega_1 / 2 - omega_2 / 2) * k[0] * c_T
    - w[0] * w[1]
    + omega_1 * omega_2
)
quad = np.zeros(2)
quad = quadratic(a, b, c)

# c1_fixed = w[1] / quad[0]
# c2_fixed = w[1] / quad[1]


# c1_fixed = quad1
# c2_fixed = quad2
# c1_fixed = ((q[0] * m_0) - (w[1] * n[0])) / (omega_1 - (k_1 * n[0]) - w[1])
# c2_fixed = -(w[0] * n[0]) / (omega_2 - k_2 * n[0] - w[0])

c1_fixed = n[0]
c2_fixed = (1 - (w[0] / w[1])) * n[1]

# -(w[1] / w[0]) * n[1]
# ((w[0] * w[1]) / ((2 * w[0] + w[1]) * (2 * w[1] + w[0]))) * n[0]
# c2_fixed = w1_frac * c_tt

# c1_fixed = w2_frac * n[0]
# c2_fixed = w1_frac * n[1]
m_fixed = m_0 / q[1]

# c1_fixed = ((q_1) - (w[1] * n[0])) / (omega_1 - k_1 * n[0] - w[1])
# c2_fixed = -(w[0] * c_T) / (omega_2 - k_2 * c_T - w[0])

# steps[0, :] = c_10 - c2_fixed, c2_fixed, m_0
steps[0, :] = c_10 - c_20, c_20, m_0
for j, _ in enumerate(t_array):
    steps[j, :] += steps[j - 1, :]
    steps[j, :] += step_function(steps[j - 1, :], dt)  # print(steps[j, :])

dt_zeroes = [c1_fixed, c2_fixed, m_fixed]

color_list = ["b", "r", "g"]
name_list = [
    "Non-resistant Type",
    "Resistant Type",
    "Antibiotic",
]
dt_name_list = [
    r"$c_1$ fixed point at ${w_2 n_1 \over w_1 + w_2}$",  # {c1_fixed}",
    r"$c_2$ fixed point at ${w_1 n_2 \over w_1 + w_2}$",  # {c2_fixed}",
    f"$m$ fixed point at {m_fixed}",
]


## Plotting

fig, ax = plt.subplots()

# ax.plot(steps.T[0], steps.T[1])
# ax.plot(t_array, (steps.T[0] / steps.T[1]))
# ax.hlines(c2_fixed, 0, t_end)

for i, curve in enumerate(steps.T):
    # break
    if i == 2:
        break
    color = color_list[i]
    curve_name = name_list[i]
    dt_curve_name = dt_name_list[i]
    zero_line = dt_zeroes[i]

    ax.plot(t_array, curve, label=curve_name, color=color)
    ax.hlines(
        zero_line,
        0,
        t_end,
        label=dt_curve_name,
        color=color,
        linestyle="dashed",
        linewidth=1,
    )

ax.set_xlabel("Time")
ax.set_ylabel("Percentage")

ax.set_xlim(0, t_end)
# ax.set_ylim(0, 100)

figure_file = os.path.join(figure_path, "plot.pdf")
plt.legend()
# plt.savefig(figure_file)
plt.show()
