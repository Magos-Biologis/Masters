#!./.venv/bin/python
import os

import numba
from numba import njit
from pylab import *


figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_path = os.path.join(figure_env, "fpe")

# file_name = "ode_solution"

file_name = "simulated_fpe"

print()


b = 100


alpha = 0
beta = 1


n = 100

x_lims = (0, n)

boxes = arange(0, n + 1, 1, dtype=np.int_)

# boxes = range(n)
# print(boxes)
# exit()
# boxes = n+1

hist_kwargs = {
    "bins": boxes,
    "density": True,
    "edgecolor": "black",
    "alpha": 0.7,
    "align": "left",
    # "normed": True,
}

curve_kwargs = {
    "linewidth": 5,
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
k_1 /= n


k_2 = 1
file_name += "_" + f"k2{k_2}"
k_2 /= n


k = array([k_1 / (k_1 + k_2), k_2 / (k_1 + k_2)])

file_name += "_" + f"n{n}"


a_1 = n
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

    num1 = 2 * k1
    num2 = (k1 - k2) * x - k2 * log(1 + ((k1 - k2) * x) / k2)
    den = (k1 - k2) ** 2
    result = x - (num1 * num2) / den
    return result


def expon(x, **kwargs):
    n = kwargs.get("nt", 1)
    max_value = kwargs.get("max_value", 0)

    return exp(2 * n * C(x, **kwargs) - max_value)


# @njit
def p_s(x, **kwargs):
    calced = kwargs.get("log_results", 1)

    top = expon(x, **kwargs)
    bot = B(x, **kwargs)
    return top / bot


steps = 500
x_array = linspace(alpha, beta, num=steps, endpoint=False)[1:]

kwarg_dict = {"k1": k_1, "k2": k_2, "nt": n}

# exponential_array = log(2 * kwarg_dict["nt"] * C(x_array, **kwarg_dict))
# exponential_array = log(C(x_array, **kwarg_dict))

# plot(x_array, results)
# print(rand())


@njit
def gillespie(
    aj: ndarray[tuple[int], dtype[float64]] | list[float],
) -> tuple[int, float]:
    # assert type(aj) is ndarray
    j: int = -1

    if aj.min() <= 0.0:
        tau: float = 0
        return j, tau

    # np.random.seed(seed)
    r = rand(2)
    while r[0] == 0.0:
        r[0] = rand()

    a_0: float = aj.sum()
    tau: float = log(1 / r[0]) / a_0
    ra_0: float = r[1] * a_0

    for k, _ in enumerate(aj):
        j: int = k
        s_j: float = aj[: k + 1].sum()
        if s_j >= ra_0:
            break

    return j, tau


# from old.gillespie import gillespie


@njit
def aj(
    x: ndarray[tuple[int], dtype[float64 | int_]],
) -> ndarray[tuple[int], dtype[float64]]:
    return array([k_1 * x[0], k_2 * x[1]])


## Step function
@njit
def step_function(
    steps: int, x0: ndarray[tuple[int], dtype[float64]]
) -> tuple[ndarray[tuple[int], dtype[float64]], ndarray[tuple[int], dtype[int_]]]:
    ## v_j
    v1 = array([-1, 1], dtype=int_)
    v2 = array([1, -1], dtype=int_)
    v = [v1, v2]

    zero_vector = zeros_like(v1)

    ## Other
    gillespie_results = empty((2, steps), dtype=int_)
    time_array = empty(steps, dtype=float64)

    gillespie_results[:, 0] = x0
    x = x0.astype(float64)

    for i in range(1, steps + 1):
        a_j = aj(gillespie_results[:, i - 1])
        j, dt = gillespie(a_j)

        if j != -1:
            x = add(x, v[j])  ## Using numpy for speed I guess

        else:  ## When j is -1, we need to check to bounds
            if x[0] == 0:
                x = add(x, v[1])
            elif x[1] == 0:
                x = add(x, v[0])

        gillespie_results[:, i] = x
        time_array[i] = time_array[i - 1] + dt

    return time_array, gillespie_results


# time_array, gillespie_results = step_function(100000, x_0)

steps = 1000000

time_array_1, gillespie_results_1 = step_function(steps, array([n - 1, 1]))
# time_array_2, gillespie_results_2 = step_function(steps, array([0, n]))
time_array_2, gillespie_results_2 = step_function(steps, array([1, n - 1]))


if n == 1:
    max_val = log(2.022121436749997)
elif n == 10:
    max_val = log(674763.2054860857)
elif n == 100:
    max_val = log(6.367225127715182e43)
else:
    max_val = 0

# val = 10**n
# log(val)


fig, ax = subplots(2, 2, figsize=(10, 10))
# fig1, ax1 = subplots()
# fig2, ax2 = subplots()

# print(boxes)
# exit()

ax[0, 0].step(time_array_1, gillespie_results_1[0, :])
ax[0, 0].step(time_array_1, gillespie_results_1[1, :])

ax[0, 0].set_ylim(bottom=0)
ax[0, 0].set_xlim(left=0)

ax[0, 1].step(time_array_2, gillespie_results_2[0, :])
ax[0, 1].step(time_array_2, gillespie_results_2[1, :])

ax[0, 1].set_ylim(bottom=0)
ax[0, 1].set_xlim(left=0)

# gillespie_hist11, gillespie_edges11 = histogram(gillespie_results_1[0, :])
# gillespie_hist12, gillespie_edges12 = histogram(gillespie_results_1[1, :])
#
# gillespie_hist21, gillespie_edges21 = histogram(gillespie_results_2[0, :])
# gillespie_hist22, gillespie_edges22 = histogram(gillespie_results_2[1, :])


ax[1, 0].hist(gillespie_results_1[0, :], **hist_kwargs)
ax[1, 0].hist(gillespie_results_1[1, :], **hist_kwargs)
ax[1, 0].set_ylim(bottom=0)
ax[1, 0].set_xlim(*x_lims)

ax[1, 1].hist(gillespie_results_2[0, :], **hist_kwargs)
ax[1, 1].hist(gillespie_results_2[1, :], **hist_kwargs)
ax[1, 1].set_ylim(bottom=0)
ax[1, 1].set_xlim(*x_lims)


# ax1.hist(gillespie_results_1[0, :], **hist_kwargs)

# for i, axes in enumerate([ax1, ax2]):
#     axes.plot(x_array * n, hist.max() * results, **curve_kwargs)
#     axes.set_xlim(0, n)
#     axes.vlines([n * k[1]], 0, 1, **line_kwargs)
#     axes.set_ylim(0, hist.max() + hist.max() / 10)
#
# ax2.hist(gillespie_results_2[0, :], **hist_kwargs)


file_path = os.path.join(figure_path, file_name)


# fig1.show()
# fig2.show()

show()

exit()
fig1.savefig(file_path + "_x0_n_0" + ".pdf")
fig2.savefig(file_path + "_x0_1_n-1" + ".pdf")
