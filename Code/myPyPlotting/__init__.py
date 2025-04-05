# from myPyPlotting import *

import os
import itertools
import collections
from dataclasses import dataclass, field

from sys import exit
from pprint import pp

import numpy as np
from numpy import polynomial as npp
from numpy._core.multiarray import dtype

from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter


@dataclass
class parameter_class:
    m: int

    k: float | np.ndarray[tuple[int], np.dtype[np.float64]]
    n: float | np.ndarray[tuple[int], np.dtype[np.float64]]
    q: float | np.ndarray[tuple[int], np.dtype[np.float64]]
    w: float | np.ndarray[tuple[int], np.dtype[np.float64]]

    W: None | np.ndarray[tuple[int, int], np.dtype[np.float64]] = None
    K: None | np.ndarray[tuple[int, int], np.dtype[np.float64]] = None

    ## The dataclass equivilant of initialization functions
    def __post_init__(self):
        self.k = self.__size_verification(self.k)
        self.n = self.__size_verification(self.n)
        self.q = self.__size_verification(self.q)
        self.w = self.__size_verification(self.w)

        self.W, self.K = self.__set_matrices()

    def __size_verification(self, parameter) -> np.ndarray:
        output = (
            np.ones(shape=self.m, dtype=np.float64)[:] * parameter
            if type(parameter) is not np.ndarray[tuple[self.m], np.dtype[np.float64]]
            else parameter
        )
        return output

    def __set_matrices(self) -> tuple:  # np.ndarray:
        assert type(self.k) is np.ndarray
        assert type(self.n) is np.ndarray
        ran = range(self.m)
        # [tuple[parameters.m], np.dtype[np.float64]]
        # [tuple[parameters.m], np.dtype[np.float64]]

        ones_no_diag = np.ones((self.m, self.m)) - np.diag(np.ones(self.m))

        w_temp = self.w
        for i in range(self.m - 1):
            w_temp = np.append(w_temp, self.w)
        w_reshaped = np.reshape(w_temp, (self.m, self.m))
        indices = np.diag_indices(self.m)
        for i, j in zip(*indices):
            w_reshaped[i, j] = self.k[i]

        matrix_W = self.W if type(self.W) is None else w_reshaped

        matrix_K = (
            self.K
            if type(self.W) is None
            else np.array([[self.k[i] / self.n[j] for j in ran] for i in ran])
        )

        return matrix_W, matrix_K


class ODEModel:
    def __init__(
        self,
        t_range: tuple[int, int],
        parameters: parameter_class | tuple[int, float, float, float, float],
        initial_condition: list[float]
        | np.ndarray[tuple[int], np.dtype[np.float64]]
        | None = None,
        **kwargs,
    ):
        assert type(parameters) is parameter_class, "Incorrect format of parameters"

        self.p: parameter_class = parameters

        self.m: int = self.p.m
        self.ran: range = range(len(initial_condition))
        self.init_cond = initial_condition

        self.n_max = kwargs.get("n_max", self.p.n)
        self.t_span: tuple[int, int] = t_range
        self.dt: np.float64 = kwargs.get("dt", 0.01)
        self.t_array: np.ndarray = np.arange(*self.t_span, self.dt)

    def __step_function(self, xs) -> np.ndarray:
        assert type(self.p) is parameter_class, "Incorrect format of parameters"
        step = np.zeros_like(xs)

        ct = sum(xs)

        return step

    def __integrate_model(self):
        # steps = np.zeros((self.m, len(self.t_array)))
        steps = np.zeros((len(self.init_cond), len(self.t_array)))
        steps[:, 0] = self.init_cond

        for j, _ in enumerate(self.t_array):
            steps[:, j] += steps[:, j - 1]
            steps[:, j] += self.dt * self.__normal_model(steps[:, j - 1])

        return self.t_array, steps

    def __normal_model(
        self, xs: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        # assert (
        #     all(type(self.p.k), type(self.p.n), type(self.p.q), type(self.p.w))
        #     is np.ndarray[tuple[int], np.dtype[np.float64]]
        # )

        assert len(xs) == 3

        step = np.empty_like(xs)

        c1, c2, m = xs
        ct = sum(xs[0:2])

        step[0] = (
            (self.p.k[0] * (1 - ct / self.p.n[0])) * c1
            - self.p.w[0] * c1
            + self.p.w[1] * c2
            - self.p.q[0] * m * c1
        )
        step[1] = (
            (self.p.k[1] * (1 - ct / self.p.n[1])) * c2
            - self.p.w[1] * c2
            + self.p.w[0] * c1
        )
        step[2] = -self.p.q[1] * m * c2 + m_0

        return step

    def __general_model(self, xs: np.ndarray):
        step = np.zeros_like(xs)
        ct = sum(xs)

        for i, j in itertools.product(self.ran, repeat=2):
            if i != j:
                step[i] += self.p.W[j, i] * xs[j] - self.p.W[i, j] * xs[i]
            else:
                step[i] += self.p.k[i] * xs[i] - (self.p.k[i] / self.p.n[i]) * ct * xs[i]

        return step

    def system(self):
        pass

    def __general_no_med(self):
        pass

    def show_parameters(self):
        return self.p

    def plot_system(self):
        return self.__integrate_model()
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


# test = parameter_types(1, 2, 3, 4, 5, 6)
# k = np.array([[6.1, 6.3], [5.6, 5.9]])
k = np.array([6.1, 6.3])
n = np.array([100, 90])
# test = (3, k, 3, 4, 5)
test = (3, 1, 3, 4, 5)
parameters = parameter_class(*test)

t_end = 50
t_span = (0, t_end)


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

dt = 0.01
# t_end = 25
# t_end = 50
t_end = 150
t_array = np.arange(0, t_end, dt)
# sol1 = np.zeros((len(t_array), 3))

# parameters = [k, w, q, n]
parameters = parameter_class(2, k, n, q, w)

init_conds1 = [c1_0, c2_0, m_0]
integrator = ODEModel((0, t_end), parameters, init_conds1)
t_array, sol1 = integrator.plot_system()


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

ax.set_ylim(bottom=0)
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
