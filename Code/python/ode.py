# !./.venv/bin/python
import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter

# import sympy as sy
# import myodestuff
from myodestuff import ODEModel, ODEParameters

print("")


# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
# figure_path = os.path.join(figure_env, "ode")

figure_path = str(os.getenv("ODE_FIGURE_ENV"))

file_name = "ode_solution"

### Variables
# Growth Constant, same because same cell
n = np.zeros(2, dtype=np.float64)
k = np.zeros(2, dtype=np.float64)
w = np.zeros(2, dtype=np.float64)
q = np.zeros(2, dtype=np.float64)

k[:] = 0.2

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = 100
n[1] = 90

w[0] = 0.15
w[1] = 0.015

# w[1] = 0.14
# n[0] = 55
# n[1] = 85

# q[0] = 0.999
# q[1] = 0.8

q[:] = 1
m_0 = 0.0


alpha = 0
c1_0 = n[0] - alpha
c2_0 = alpha


parameter_dict = {
    "m0": m_0,
    "k1": k[0],
    "k2": k[1],
    "n1": n[0],
    "n2": n[1],
    "w1": w[0],
    "w2": w[1],
    "q1": q[0],
    "q2": q[1],
}


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
# t_end = 50
t_end = 150
t_array = np.arange(0, t_end, dt)
# sol1 = np.zeros((len(t_array), 3))

init_conds1 = np.array([c1_0, c2_0, m_0])
# print(len(init_conds1))
# exit()
#
parameters = ODEParameters(**parameter_dict)
model = ODEModel(parameters, (0, t_end), init_conds1)


# print()

# line1, line2 = model.roots()
#
# plt.scatter(line1[0, :], line1[1, :], label="level set")
# plt.scatter(line2[0, :], line2[1, :], label="nullcline")
#
# plt.legend()
# # plt.plot(line2[0, :], line2[1, :])
# # plt.plot(*line2)
#
#
# plt.show()


# exit()

### Integration


# parameters1 = parameter_class(*parameters)

t_array, sol1 = model.integrate()
solutions = model.roots()

# print(solutions)
# exit()

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
        alpha=0,
    )

    ax.plot(t_array, curve, label=curve_name, color=color)


# ax.plot(solutions[0], solutions[1], label="hi")

# total = sol1.T[0] + sol1.T[1]
# percent_total = np.divide(sol1.T[0], n[0]) + np.divide(sol1.T[1], n[1])

fixed_point_tick_labels = [
    "\t" + r"$c_1^*$",
    r"$c_2^*$",
]

ax.set_ylim(0, n.max())
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

ax.legend(loc="upper right")
plt.tight_layout()

file_path = os.path.join(figure_path, file_name)
# plt.savefig(file_path + ".pdf")

# print(sol1[-1])
plt.show()
