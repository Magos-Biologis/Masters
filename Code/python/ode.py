#!./.venv/bin/python
import os
from pprint import pprint as pp
from sys import exit

import numpy as np
from numpy import polynomial as npp


from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter
# import sympy as sy

from myPyPlotting import parameter_class, ODEModel

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

k[:] = 0.2

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = 100.0
n[1] = 90.0

w[0] = 0.015
w[1] = 0.015

# w[1] = 0.14
# n[0] = 55
# n[1] = 85

# q[0] = 0.999
# q[1] = 0.8

q[:] = 0
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

# parameters = [k, w, q, n]
parameters = parameter_class(2, k, n, q, w)

cq = sum(n) / len(n)

### Integration


def step_function(xs) -> np.ndarray:
    step = np.empty_like(xs)

    ct = sum(xs[0:2])
    c1, c2, m = xs

    step[0] = (k[0] * (1 - ct / n[0])) * c1 - w[0] * c1 + w[1] * c2 - q[0] * m * c1
    step[1] = (k[1] * (1 - ct / n[1])) * c2 - w[1] * c2 + w[0] * c1
    step[2] = -q[1] * m * c2 + m_0

    return step


def integrate(initial, t_array):
    steps = np.zeros((len(t_array), 3))
    steps[0, :] = initial
    for j, _ in enumerate(t_array):
        steps[j, :] += steps[j - 1, :]
        steps[j, :] += dt * step_function(steps[j - 1, :])  # print(steps[j, :])

    return t_array, steps


init_conds1 = [c1_0, c2_0, m_0]
# _, sol1 = integrate(init_conds1, t_array)

integrator = ODEModel((0, t_end), parameters, init_conds1)
t_array, sol1 = integrator.plot_system()


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


# c1_fixed = n[0]
# c2_fixed = (1 - (w[0] / w[1])) * n[1]
# m_fixed = m_0 / q[1]
# dt_zeroes = [c1_fixed, c2_fixed, m_fixed]

# total = steps[:, 0] + steps[:, 1]
# steps[:, 0] = steps[:, 0] / total
# steps[:, 1] = steps[:, 1] / total

# print(2 * ((n[0] * n[1]) / (n[0] + n[1])))
# print(total)
# exit()


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

solutions = np.array([c1_root, c2_root])

# print(sum(solutions) / n[0] + sum(solutions) / n[1])

# print(sol1.T[0][-1] + sol1.T[1][-1])
# exit()


## Plotting

fig, ax = plt.subplots()

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
# ax.set_facecolor("#c0c0ca")


for i, curve in enumerate(sol1.T):
    # break
    if i == 2:
        break
    color = color_list[i]
    curve_name = name_list[i]
    # dt_curve_name = dt_name_list[i]
    # zero_line = dt_zeroes[i]
    zero_line = solutions[i]

    ax.hlines(
        zero_line,
        0,
        t_end,
        # label=dt_curve_name,
        color=color,
        linestyle="dashed",
        linewidth=1,
    )

    ax.plot(t_array, curve, label=curve_name, color=color)


total = sol1.T[0] + sol1.T[1]
percent_total = np.divide(sol1.T[0], n[0]) + np.divide(sol1.T[1], n[1])

fixed_point_tick_labels = [
    "\t" + r"$c_1^*$",
    r"$c_2^*$",
]

ax.set_ylim(0, np.max(n))
ax.set_xlim(0, t_end)

new_y_ticks = np.append(ax.get_yticks(), solutions)
# new_y_ticks = np.append([0, 25, 50, 75, 100], solutions)
ax.set_yticks(new_y_ticks)


def fixed_point_format(val, pos):
    if val == solutions[0]:
        return fixed_point_tick_labels[0]
    elif val == solutions[1]:
        return fixed_point_tick_labels[1]
    else:
        return int(np.round(val, 3))


ax.yaxis.set_major_formatter(FuncFormatter(fixed_point_format))

ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")

ax.set_xlabel("Time")
ax.set_ylabel("Count")

ax.legend(loc="center right")
plt.tight_layout()

file_path = os.path.join(figure_path, file_name)
# plt.savefig(file_path + ".pdf")

# print(sol1[-1])
plt.show()
