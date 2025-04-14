#!./.venv/bin/python
import os

import numba
from numba import njit
from pylab import *

# from old.gillespie import gillespie
import gillespie as dg


# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_env = str(os.getenv("FPE_FIGURE_ENV"))


file_name = "simulated_fpe"

print()

m = 100


boxes = arange(0, m + 1, 1, dtype=np.int_)


# print(boxes.shape)
# exit()

margin = 0.5
x_lims = (0 - margin, m + margin)

# boxes = range(n)
# print(boxes)
# exit()
# boxes = n+1

hist_kwargs = {
    "bins": boxes,
    "density": True,
    "edgecolor": "black",
    "alpha": 0.7,
    "align": "mid",
    # "normed": True,
}

curve_kwargs = {
    "linewidth": 4,
    "alpha": 0.7,
}

line_kwargs = {
    "linewidth": 3,
    "alpha": 0.6,
    "linestyle": "dashed",
    "color": "black",
}


k_1 = 1
file_name += "_" + f"k1{k_1}"
k_1 /= m


k_2 = 3
# k_2 = 1
file_name += "_" + f"k2{k_2}"
k_2 /= m

w_1 = 1
file_name += "_" + f"w1{w_1}"
w_1 /= m

w_2 = 1
file_name += "_" + f"w2{w_2}"
w_2 /= m


n_1 = 100
file_name += "_" + f"n1{n_1}"
n_1 /= m

n_2 = 100
file_name += "_" + f"n2{n_2}"
n_2 /= m


k = array([k_1 / (k_1 + k_2), k_2 / (k_1 + k_2)])

file_name += "_" + f"n{m}"


a_1 = m
a_2 = 0


aj_0 = array([k_1 * a_1, k_2 * a_2])


# def aj(x):
#     a, b = x
#     return array([k_1 * a, k_2 * b])


## S_j
s_j = []

## x_0
x_0 = array([a_1, a_2])

# print(v_vec)
# exit()

t_0 = 0
dt = 0


def B(x, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)

    n = kwargs.get("nt", 1)

    return (1 / n) * (k2 + (k1 - k2) * x)


def C(x, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)

    # if k1 - k2 != 0:
    num1 = 2 * k1
    num2 = (k1 - k2) * x - k2 * log(1 + ((k1 - k2) * x) / k2)
    den = (k1 - k2) ** 2

    prod = num1 * num2 * (den ** (-1))

    result = x - prod
    return array(result)


def expon(x, **kwargs):
    n = kwargs.get("nt", 1)
    max_value = kwargs.get("max_value", 0)

    exponent_input = 2 * n * C(x, **kwargs)

    return exp(subtract(exponent_input, max_value))


# @njit
def p_s(x, **kwargs):
    top = expon(x, **kwargs)
    bot = B(x, **kwargs)
    return top / bot


steps = 500
x_array = linspace(0, 1, num=steps, endpoint=False)[1:]

kwarg_dict = {"k1": k_1, "k2": k_2, "nt": m}


@njit
def aj(
    x: ndarray[tuple[int], dtype[float64 | int_]],
) -> ndarray[tuple[int], dtype[float64]]:
    a_1 = k_1 * x[0]
    a_2 = k_2 * x[1]

    # a_m1 =  k_1 * x[0] ** 2
    # a_m2 =  k_2 * x[1] ** 2

    return array([a_1, a_2], dtype=float64)


# @njit
# def aj(
#     x: ndarray[tuple[int], dtype[float64 | int_]],
# ) -> ndarray[tuple[int], dtype[float64]]:
#     a_1 = k_1 * x[0]
#     # a_m1 = (k_1 / n_1) * x[0] ** 2
#     a_m1 = k_1 * x[0] ** 2
#     # a_m1 = (k_1 / n_1) * x[0] * (x[0] + x[1])
#     # a_m1 = k_1 / (1 + (x[0] / n_1))
#     a_2 = w_1 * x[0]
#     a_3 = k_2 * x[1]
#     # a_m3 = (k_2 / n_2) * x[1] ** 2
#     a_m3 = k_2 * x[1] ** 2
#     # a_m3 = (k_2 / n_2) * x[1] * (x[0] + x[1])
#     # a_m3 = k_2 / (1 + (x[1] / n_2))
#     a_4 = w_2 * x[1]
#     return array([a_1, a_m1, a_2, a_3, a_m3, a_4], dtype=float64)


## Step function
@njit
def step_function(
    steps: int, x0: ndarray[tuple[int], dtype[float64]]
) -> tuple[ndarray[tuple[int], dtype[float64]], ndarray[tuple[int], dtype[int_]]]:
    ## v_j
    v1 = array([-1, 1], dtype=int_)
    vm1 = array([1, -1], dtype=int_)
    v2 = array([1, -1], dtype=int_)
    vm2 = array([-1, 1], dtype=int_)
    v = [v1, vm1, v2, vm2]

    # v1 = array([1, 0], dtype=int_)
    # vm1 = array([-1, 0], dtype=int_)
    # v2 = array([-1, 0], dtype=int_)
    # v3 = array([0, 1], dtype=int_)
    # vm3 = array([0, -1], dtype=int_)
    # v4 = array([1, -1], dtype=int_)
    # v = [v1, vm1, v2, v3, vm3, v4]

    ## Other
    gillespie_results = empty((2, steps), dtype=int_)
    time_array = empty(steps, dtype=float64)

    gillespie_results[:, 0] = x0
    x = x0.astype(float64)

    broke_loop = False
    for i in range(1, steps + 1):
        a_j = aj(gillespie_results[:, i - 1])
        j, dt = dg.ssa_event(a_j)

        if j == -1:
            broke_loop = True
            final_step: int = i
            break

        gillespie_results[:, i] = add(gillespie_results[:, i - 1], v[j])
        time_array[i] = time_array[i - 1] + dt

    if not broke_loop:
        return time_array, gillespie_results
    else:
        return time_array[:final_step], gillespie_results[:, :final_step]


# time_array, gillespie_results = step_function(100000, x_0)

steps = 1000000

time_array_1, gillespie_results_1 = step_function(steps, array([m, 0]))
time_array_2, gillespie_results_2 = step_function(steps, array([0, m]))
# time_array_2, gillespie_results_2 = step_function(steps, array([0, n]))


if m == 1:
    max_val = log(2.022121436749997)
elif m == 10:
    max_val = log(674763.2054860857)
elif m == 100:
    max_val = log(6.367225127715182e43)
else:
    max_val = 0

analytical_results = p_s(x_array, **kwarg_dict, max_value=max_val)


fig, ax = subplots(2, 2, figsize=(10, 10))

ax[0, 0].step(time_array_1, gillespie_results_1[0, :])
ax[0, 0].step(time_array_1, gillespie_results_1[1, :])

ax[0, 0].set_ylim(bottom=0)
ax[0, 0].set_xlim(left=0)

ax[0, 1].step(time_array_2, gillespie_results_2[0, :])
ax[0, 1].step(time_array_2, gillespie_results_2[1, :])

ax[0, 1].set_ylim(bottom=0)
ax[0, 1].set_xlim(left=0)

ax[1, 0].hist(gillespie_results_1[0, :], **hist_kwargs)
# ax[1, 0].hist(gillespie_results_1[1, :], **hist_kwargs)
ax[1, 0].set_ylim(bottom=0)
ax[1, 0].set_xlim(*x_lims)

ax[1, 1].hist(gillespie_results_2[0, :], **hist_kwargs)
# ax[1, 1].hist(gillespie_results_2[1, :], **hist_kwargs)
ax[1, 1].set_ylim(bottom=0)
ax[1, 1].set_xlim(*x_lims)

hist, edges = histogram(gillespie_results_1[0, :], bins=boxes, density=True)
scale = hist.max()
scaled_result = multiply(analytical_results, scale)


show()

# print(scale)
# exit()
# ax[1, 0].plot(x_array * n, scaled_result, **curve_kwargs)
# ax[1, 1].plot(x_array * n, scaled_result, **curve_kwargs)


# fig1, ax1 = subplots()
# fig2, ax2 = subplots()
#
# ax1.hist(gillespie_results_1[0, :], **hist_kwargs)
#
# for i, axes in enumerate([ax1, ax2]):
#     axes.set_xlim(0, n)
#     axes.set_ylim(0, hist.max() + hist.max() / 10)
#
#     axes.plot(x_array * n, scaled_result, **curve_kwargs)
#     axes.vlines([n * k[1]], 0, 1, **line_kwargs)
#
# ax2.hist(gillespie_results_2[0, :], **hist_kwargs)
#
#
# file_path = os.path.join(figure_env, file_name)
#
# show()

# print(file_path)
# exit()

# exit()
# fig1.savefig(file_path + "_x0_n_0" + ".pdf")
# fig2.savefig(file_path + "_x0_1_n-1" + ".pdf")
