import os
from parameter_class import parameter_class
from ODEModel import ODEModel


import collections

from sys import exit
from pprint import pp

import numpy as np
from numpy import polynomial as npp

from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter


# pass


# class ODEModel:
#     """Represents an ODE model that uses ODEParameters."""
#
#     def __init__(self, params):
#         self.params = params
#
#     def rhs(self, t, y):
#         # Logistic growth: dy/dt = α*y - β*y².
#         # (It’s like life: growth tempered by too much competition.)
#         α = self.params.alpha
#         β = self.params.beta
#         return α * y - β * y**2


# test2 = ODEModel(t_span, parameters)
# test = parameter_types(m=2)
# test2 = parameter_class(test)
print()
# pp(test)
# pp(test2.show_parameters())
# print()
# # print(test2.p)
# # print()
# print(test2)

### Variables
# Growth Constant, same because same cell
k = np.zeros(2)
n = np.zeros(2)
w = np.zeros(2)
q = np.zeros(2)

k[:] = 0.25

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

dt = 0.01
# t_end = 25
# t_end = 50
t_end = 150
t_array = np.arange(0, t_end, dt)
# sol1 = np.zeros((len(t_array), 3))

# parameters = [k, w, q, n]
parameters = parameter_class(2, m_0, k, n, q, w)

init_conds1 = [c1_0, c2_0, m_0]
integrator = ODEModel((0, t_end), parameters, init_conds1)
t_array, sol1 = integrator.integrate()


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

for i, curve in enumerate(sol1):
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

ax.set_ylim(bottom=0, top=n.max())
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

# ax.legend(loc="center right")

ax.legend()
plt.tight_layout()

# file_path = os.path.join(figure_path, file_name)
# plt.savefig(file_path + ".pdf")

# print(sol1[-1])
plt.show()

exit()


@dataclass
class base:
    def __init__(
        self,
        figure_name: str,
        subdirectory: str | None = None,
        figure_enviroment: str | None = None
        if type(os.getenv("THESIS_FIGURE_PATH")) is None
        else os.getenv("THESIS_FIGURE_PATH"),
    ) -> None:
        self.figure_name = figure_name
        self.subdirectory = subdirectory
        self.figure_env = figure_enviroment
