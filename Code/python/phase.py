#!./.venv/bin/python
import os

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter

# from numpy.polynomial import Polynomial
from myodestuff import ODEModel, parameter_class

figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
phase_path = os.path.join(figure_env, "phase")


file_name = "phaseplane"
# figure_path = "./figs"


three_d = False

### Variables
# Growth Constant, same because same cell
k = np.zeros(2)
k1 = np.zeros(2)

n = np.zeros(2)
n1 = np.zeros(2)
n2 = np.zeros(2)

w = np.zeros(2)
q = np.zeros(2)

k[0] = 1.2
k[1] = 0.2

k1[0] = 0.6
k1[1] = 1.2

# k[0] = 2
# k[1] = 2

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = 100
n[1] = 80

# n1[0] = 100
# n1[1] = 100

n2[0] = 100
n2[1] = 100

w[0] = 0.015
w[1] = 0.015

# w[1] = 0.035

# w[1] = 0.14
# n[0] = 55
# n[1] = 85


q[1] = qm = 0.8
q[:] = 0.999


omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]

k_1 = k[0] / n[0]
k_2 = k[1] / n[1]

c1_min = w[1] / k_1
c2_min = w[0] / k_2

c1_max = omega_1 / k_1
c2_max = omega_2 / k_2


alpha = 1
c1_0 = n[0] - alpha
c2_0 = alpha
m_0: int = 0

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
t_end = 250
t_array = np.arange(0, t_end, dt)

parameters1 = parameter_class(m=2, m_0=m_0, k=k, n=n, q=q, w=w)
# parameters2 = parameter_class(m=2, m_0=m_0, k=k1, n=n1, q=q, w=w)
# parameters3 = parameter_class(m=2, m_0=0, k=k, n=n2, q=q, w=w)


init_conds1 = np.array([c1_0, c2_0, m_0])
init_conds2 = np.array([c2_0, c1_0, 0])
init_conds3 = np.array([c1_0, c2_0, m_0])
# init_conds3 = np.array([100, 100, 100])

# if three_d:
# else:
#     init_conds1 = np.array([c1_0, c2_0,])
#     init_conds2 = np.array([c2_0, c1_0])
#     init_conds3 = np.array([100, 100])

model1 = ODEModel((0, t_end), parameters1, init_conds1)
# model2 = ODEModel((0, t_end), parameters2, init_conds2)
# model3 = ODEModel((0, t_end), parameters3, init_conds3)


### Integration

# parameters1 = parameter_class(*parameters)

t_array, sol1 = model1.integrate()
# t_array, sol2 = model2.integrate()
# t_array, sol3 = model3.integrate()

c1_root, c2_root = model1.roots()

## Plotting

x_lims = 0, 100
y_lims = 0, 100

# x_lims = 0, 1
# y_lims = 0, 1


c1s = np.arange(0.1, 125, 0.5)
c2s = np.arange(0.1, 125, 0.5)
ms = np.arange(0.0, m_0 + 1, 0.05)


c1, c2, m = np.meshgrid(c1s, c2s, ms)
c1, c2 = np.meshgrid(c1s, c2s)


# print(c1[0, :, 0])
# exit()


def vector_space(c1, c2, m=0):
    dc1 = (k[0] * (1 - (c1 + c2) / n[0]) - w[0]) * c1 + w[1] * c2 - q[0] * m * c1
    dc2 = (k[1] * (1 - (c1 + c2) / n[1]) - w[1]) * c2 + w[0] * c1
    dm = m_0 - q[1] * m * c2

    return (dc1, dc2, dm)


# vector_field = vector_space(c1, c2, m)
vector_field = vector_space(c1, c2)

dU = vector_field[0]
dV = vector_field[1]
dW = vector_field[2]


# speed = np.sqrt(dU[:, :, 0] ** 2 + dV[:, :, 0] ** 2)
speed = np.sqrt(dU**2 + dV**2)
lw = 4 * speed / speed.max()

stream_kwargs = {
    "density": 1.7,
    "linewidth": lw,
    "arrowstyle": "->",
    "color": "black",
    # "color":lw,
    # "norm":norm,
    # "cmap":"gist_heat_r",
    # "arrowsize":0,
    # "density":0.7,
    # "broken_streamlines":False,
}


resolution = c1.shape[0]
# c1s = np.arange(c1_min, c1_max, 0.5)
# c2s = np.arange(c2_min, c2_max, 0.5)
c1range = np.linspace(c1_min + 0.001, c1_max - 0.001, resolution)
c2range = np.linspace(c2_min + 0.001, c2_max - 0.001, resolution)


plt.rcParams.update(
    {
        "axes.labelsize": 20,
        "axes.titleweight": "bold",
        "xtick.labelsize": 15,
        "ytick.labelsize": 15,
        # "axes.titlecolor": "white",
        # "xtick.labelcolor": "white",
        # "ytick.labelcolor": "white",
        # "savefig.facecolor": "#c0c0ca",
    }
)

null_kwargs = {
    "linewidth": 1.2,
    "alpha": 0.8,
}

fixed_kwargs = {
    "linestyle": "--",
    "linewidth": 1.2,
    "alpha": 0.8,
    "zorder": 11,
}


# norm = mpl.colors.Normalize(vmin=lw.min(), vmax=lw.max())
# norm = mpl.colors.LogNorm(vmin=lw.min(), vmax=lw.max())

### Plotting

plt.style.use("bmh")

if three_d:
    fig = plt.figure()
    ax = plt.axes(projection="3d")

    ax.set_zlabel(
        r"$m$",
    )
else:
    fig, ax = plt.subplots()

    # ax.streamplot(c1[:, :, 0], c2[:, :, 0], dU[:, :, 0], dV[:, :, 0], **stream_kwargs)
    ax.streamplot(c1, c2, dU, dV, **stream_kwargs)


# ax.set_facecolor("#c0c0ca")
ax.set_facecolor("white")

ax.tick_params(color="white")
# ax.spines["all"].set_color("white")

ax.set_xlabel(
    r"$c_1$",
)
ax.set_ylabel(
    r"$c_2$",
)

ax.set_xlim(*x_lims)
ax.set_ylim(*y_lims)


c1f = omega_1 / k_1
c2f_b = (omega_1 - k_1 * c1f) / (w[1] - c1f)


sol_alpha = 0.0
colors = ["b", "r"]


def example_plot():
    ax.plot(
        sol1[:, 0],
        sol1[:, 1],
        color=colors[0],
        label=r"Init conds, $c_1="
        + f"{int(init_conds1[0])}"
        + "$, $c_2="
        + f"{int(init_conds1[1])}"
        + "$",
        alpha=sol_alpha,
    )
    ax.scatter(sol1[0, 0], sol1[0, 1], color=colors[0], alpha=sol_alpha)

    ax.plot(
        sol2[:, 0],
        sol2[:, 1],
        color=colors[1],
        label=r"Init conds, $c_1="
        + f"{int(init_conds2[0])}"
        + "$, $c_2="
        + f"{int(init_conds2[1])}"
        + "$",
        alpha=sol_alpha,
    )
    ax.scatter(sol2[0, 0], sol2[0, 1], color=colors[1], alpha=sol_alpha)


## Nullcline


def c1_sol(c1, c2):
    num = q[0] * (m_0 / (q[1] * c2)) * c1 - (c1 * (omega_1 - k_1 * c1))
    den = w[1] - k_1 * c1
    return num / den


def c2_sol(c2):
    num = -(c2 * (omega_2 - k_2 * c2))
    den = w[0] - k_2 * c2
    return num / den


dx = 0.01

# c1_sol_array = np.linspace(c1_min + dx, 100, resolution)
# c2_sol_array = np.linspace(c2_min + dx, 100, resolution)
# m_sol_array = np.linspace(0.0, m_0 + dx, resolution, endpoint=True)


# colors = ["k", "k"]


def plot_null():
    c1_col = "blue"
    c2_col = "red"

    ax.plot(
        # c2_sol(c2)[0, :, 0],
        # c2[:, 0, 0],
        c2_sol(c2range),
        c2range,
        label=r"$c_1$ Nullcline",
        color=c1_col,
        **null_kwargs,
    )
    ax.plot(
        # c1[0, :, 0],
        # c1_sol(c1, m)[:, 0, 0],
        c1range,
        c1_sol(c1range, c2range),
        label=r"$c_2$ Nullcline",
        color=c2_col,
        **null_kwargs,
    )

    ax.vlines(
        c1_root,
        *x_lims,
        color=c1_col,
        **fixed_kwargs,
    )
    ax.hlines(
        c2_root,
        *y_lims,
        color=c2_col,
        **fixed_kwargs,
    )

    # ax.plot(
    #     c1_sol_array,
    #     (1 - c1_sol_array * k[0] / n[0]) * n[1] / k[1],
    #     zorder=12,
    #     label=r"$c_2 = \frac{n_2}{k_2} \left(1 - \frac{k_1 c_1}{n_1}\right)$",
    #     color="indigo",
    # )
    # ax.plot(
    #     c1_sol_array,
    #     (1 - c1_sol_array / (n[0] * k[0])) * k[0] * n[1],
    #     zorder=12,
    #     label=r"$c_2 = k_2 n_2 \left(1 - \frac{c_1}{k_1 n_1}\right)$",
    #     color="firebrick",
    # )

    # ax.plot([0, n[0]], [n[1], 0], label="ratio", color="hotpink", alpha=0.5)


# def pop_curve():
#     c1 = n[0] - m_0 * q[0]
#     c2 = n[1]
#
#     return [c1, 0], [0, c2]


# c2_test = n[1] * (1 - c1s / n[0])

origin = np.array([[0, 0], [0, 0]])
extremes = np.array([[0, n[0]], [n[1], 0]])

# points = [[], []]
# margin = 1 / c1.shape[0]
# for i in range(c1.shape[0]):
#     for j in range(c1.shape[1]):
#         c1_norm = c1[i, j] / n[0]
#         c2_norm = c2[i, j] / n[1]
#
#         c_sum = c1_norm + c2_norm
#         # print(c_sum)
#
#         if 1 - margin <= c_sum <= 1 + margin:
#             points[0].append(c1[i, j])
#             points[1].append(c2[i, j])
#         else:
#             continue


# ax1.yaxis.set_ticklabels()

# solutions =


new_x_ticks = np.append(ax.get_xticks(), c1_root)
new_y_ticks = np.append(ax.get_yticks(), c2_root)

ax.set_xticks(new_x_ticks)
ax.set_yticks(new_y_ticks)


fixed_point_tick_labels = [r"$c_1^*$", r"$c_2^*$"]
# fixed_point_tick_labels = [r"$c_1^* =" + f"{round(c1_root[0], 2)}" + r"$", r"$c_2^* =" + f"{round(c2_root[0], 2)}" + r"$",]

#


# ax.scatter(*init_conds1)
def plot_3d(sol, ax):
    ax.plot3D(sol[0, :], sol[1, :], sol[2, :], label="Trajectory")
    ax.scatter3D(sol[0, 0], sol[1, 0], sol[2, 0], color="k")


if three_d:
    plot_3d(sol1, ax)
    plot_3d(sol2, ax)
    plot_3d(sol3, ax)
else:
    plot_null()


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


# print(c2_root)


ax.xaxis.set_major_formatter(FuncFormatter(x_format))
ax.yaxis.set_major_formatter(FuncFormatter(y_format))


# plt.tight_layout()


# print(sol1[:, -1])
# print(sol2[:, -1])
# print(sol3[:, -1])


ax.legend(framealpha=1, loc="upper right")
figure_file = os.path.join(phase_path, file_name)

plt.savefig(figure_file + ".pdf", format="pdf")
# plt.show()
