#!./.venv/bin/python
import os
from pprint import pprint as pp
from sys import exit

import numpy as np
from numpy import polynomial as npp


from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter
# import sympy as sy

print("")


figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_path = os.path.join(figure_env, "ode")


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

# n[1] = 90
# w[1] = 0.015

q[0] = 0.999
q[1] = 0.8

m_0 = 0

omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]

k_1 = k[0] / n[0]
k_2 = k[1] / n[1]

c1_min = w[1] / k_1
c2_min = w[0] / k_2

alpha = 0
c1_0 = n[0] - alpha
c2_0 = alpha

filename = (
    "ode_solution"
    + "_m0"
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


dt = 0.01
t_end = 25
# t_end = 150
t_array = np.arange(0, t_end, dt)
sol1 = np.zeros((len(t_array), 3))

parameters = [k, w, q, n]

cq = sum(n) / len(n)

### Integration


def step_function(xs, dt) -> np.ndarray:
    c1, c2, m = xs
    ct = c1 + c2

    c1_dot = (k[0] * (1 - ct / n[0]) - w[0]) * c1 + w[1] * c2 - q[0] * m * c1
    c2_dot = (k[1] * (1 - ct / n[1]) - w[1]) * c2 + w[0] * c1
    m_dot = -q[1] * m * c2 + m_0

    return dt * np.array([c1_dot, c2_dot, m_dot])


def integrate(initial, t_array):
    steps = np.zeros((len(t_array), 3))
    steps[0, :] = initial
    for j, _ in enumerate(t_array):
        steps[j, :] += steps[j - 1, :]
        steps[j, :] += step_function(steps[j - 1, :], dt)  # print(steps[j, :])

    return t_array, steps


init_conds1 = [c1_0, c2_0, m_0]
_, sol1 = integrate(init_conds1, t_array)


M = np.array(
    [
        [(omega_1 - k_1), w[0]],
        [w[1], (omega_2 - k_2)],
    ]
)


probability_sol = np.zeros(2)

solution_norm = M[0, 1] + M[1, 0]
probability_sol[0] = M[1, 0] / solution_norm
probability_sol[1] = M[0, 1] / solution_norm
probability_ticks = w / sum(w)


# print(solutions)


c1_fixed = n[0]
c2_fixed = (1 - (w[0] / w[1])) * n[1]
m_fixed = m_0 / q[1]


dt_zeroes = [c1_fixed, c2_fixed, m_fixed]

# total = steps[:, 0] + steps[:, 1]
# steps[:, 0] = steps[:, 0] / total
# steps[:, 1] = steps[:, 1] / total

# print(2 * ((n[0] * n[1]) / (n[0] + n[1])))
# print(total)
# exit()


color_list = ["r", "b", "g"]
name_list = [
    "Non-resistant Type",
    "Resistant Type",
    "Antibiotic",
]
dt_name_list = [
    r"$c_1$ fixed point at $c_1^*$",  # {c1_fixed}",
    r"$c_2$ fixed point at $c_2^*$",  # {c2_fixed}",
    f"$m$ fixed point at {m_fixed}",
]

## Polynomial

coeffs = [
    n[1] * w[1],
    omega_1 - n[1] * (w[1] / n[0] + k_1),
    k_1 * ((n[1] / n[0]) - 1),
]
poly = npp.Polynomial(
    coef=coeffs,
    symbol="_c1",
)
roots = poly.roots()

c1_root = [root for root in roots if root > c1_min]
c2_root = (1 - c1_root / n[0]) * n[1]

solutions = np.array([c1_root, c2_root])

# print(sum(solutions) / n[0] + sum(solutions) / n[1])
# exit()

## Plotting

fig, ax = plt.subplots()


for i, curve in enumerate(sol1.T):
    # break
    if i == 2:
        break
    color = color_list[i]
    curve_name = name_list[i]
    dt_curve_name = dt_name_list[i]
    # zero_line = dt_zeroes[i]
    zero_line = solutions[i]

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

ax.set_xlim(0, t_end)
ax.set_ylim(0, 100)
# ax.set_ylim(0, 100)

# print()
# exit()

total = sol1.T[0] + sol1.T[1]
percent_total = np.divide(sol1.T[0], n[0]) + np.divide(sol1.T[1], n[1])

fixed_point_tick_labels = [
    r"$c_1^*$",
    r"$c_2^*$",
]


new_y_ticks = np.append(ax.get_yticks(), solutions)
ax.set_yticks(new_y_ticks)


def fixed_point_format(val, pos):
    if val == solutions[0]:
        return fixed_point_tick_labels[0]
    elif val == solutions[1]:
        return fixed_point_tick_labels[1]
    else:
        return int(np.round(val, 3))


ax.yaxis.set_major_formatter(FuncFormatter(fixed_point_format))

ax.set_xlabel("Time")
ax.set_ylabel("Count")

plt.legend(loc="center right")
plt.tight_layout()

file_path = os.path.join(figure_path, filename)
plt.savefig(file_path + ".pdf")
plt.show()
