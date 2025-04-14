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

m: int = 75
n = 1000
# n /= m

# b = 0
b = 150
# b = 1000

b /= m


alpha = 0
beta = 1


boxes = arange(0, m + 1, 1)
# boxes = n+1

hist_kwargs = {
    "bins": boxes,
    "density": True,
    "edgecolor": "black",
    # "normed": True,
    # "align": "mid",
}

curve_kwargs = {
    "linewidth": 5,
    "alpha": 0.7,
}

bayes_kwargs = {
    "linewidth": 1,
    "alpha": 0.7,
}

line_kwargs = {
    "linewidth": 3,
    "alpha": 0.6,
    "linestyle": "dashed",
    "color": "black",
}

k_1 = 0.15
k_m1 = 0.15
file_name += "_" + f"k1{k_1}"


k_2 = 1
file_name += "_" + f"k2{k_2}"

k_3 = 1
file_name += "_" + f"k3{k_3}"

k_4 = 1
file_name += "_" + f"k4{k_4}"

k_5 = 0.5
file_name += "_" + f"k5{k_5}"

k = divide(array([k_1, k_m1, k_2, k_3, k_4, k_5]), m)

file_name += "_" + f"n{m}"


# a_1 = m
# a_2 = 0
#
#
# aj_0 = array([k_1 * a_1, k_2 * a_2])


@njit
def gillespie(
    # a0: ndarray[tuple[int], dtype[float64]] | list[float],
    aj: ndarray[tuple[int], dtype[float64]],  # | list[float],
    # seed: int = 1984,
) -> tuple[int, float]:
    # assert type(aj) is ndarray

    a_0: float = aj.sum()
    j: int = -1

    if a_0 <= 0.0:
        tau: float = 0
        return j, tau

    # np.random.seed(seed)
    r: ndarray[tuple[int], dtype[float64]] = rand(2)
    while r[0] == 0:
        r[0] = rand()

    tau: float = log(1 / r[0]) / a_0

    for k, _ in enumerate(aj):
        j: int = k
        s_j: float = aj[: k + 1].sum()
        ra_0: float = r[1] * a_0

        if s_j > ra_0:
            break

    return j, tau


# from old.gillespie import gillespie


@njit
def aj(xs: ndarray[tuple[int], dtype[float64]]) -> ndarray[tuple[int], dtype[float64]]:
    x = xs[0]  # - 1
    y = xs[1]

    a_1 = k[0] * n * x
    a_m1 = k[1] * x**2
    a_2 = k[2] * n * x
    a_3 = k[3] * b * x
    a_4 = k[4] * b * y

    return array([a_1, a_m1, a_2, a_3, a_4])


@njit
def aj(xs: ndarray[tuple[int], dtype[float64]]) -> ndarray[tuple[int], dtype[float64]]:
    x = xs[0]  # - 1
    y = xs[1]
    n = xs[2]

    a_1 = k[0] * n * x
    a_m1 = k[1] * x**2
    a_2 = k[2] * n * x
    a_3 = k[3] * b * x
    a_4 = k[4] * b * y
    a_5 = k[5] * y

    return array([a_1, a_m1, a_2, a_3, a_4, a_5])


## Step function
@njit
def step_function(
    steps: int, x0: ndarray[tuple[int], dtype[float64]]
) -> tuple[ndarray[tuple[int], dtype[float64]], ndarray[tuple[int], dtype[int_]]]:
    ## v_j

    # v1 = array([1, 0], dtype=int_)
    # vm1 = array([-1, 0], dtype=int_)
    # v2 = array([0, 1], dtype=int_)
    # v3 = array([-1, 0], dtype=int_)
    # v4 = array([0, -1], dtype=int_)

    v1 = array([1, 0, -1], dtype=int_)
    vm1 = array([-1, 0, 1], dtype=int_)
    v2 = array([0, 1, -1], dtype=int_)
    v3 = array([-1, 0, 1], dtype=int_)
    v4 = array([0, -1, 1], dtype=int_)
    v5 = array([0, -1, 1], dtype=int_)

    v = [v1, vm1, v2, v3, v4, v5]

    ## Other
    gillespie_results = empty((len(x0), steps), dtype=int_)
    time_array = empty(steps, dtype=float64)

    gillespie_results[:, 0] = x0
    # try:
    # except KeyError:
    #     print("bitch")
    #     pass
    # x = x0.astype(float64)

    broke_loop = False
    for i in range(1, steps + 1):
        a_j = aj(gillespie_results[:, i - 1])
        j, dt = gillespie(a_j)

        if j == -1:
            broke_loop = True
            final_step: int = i
            break

        if gillespie_results[0, i - 1] <= 0:
            continue
        else:
            gillespie_results[:, i] = add(gillespie_results[:, i - 1], v[j])
            time_array[i] = time_array[i - 1] + dt

    if not broke_loop:
        return time_array, gillespie_results
    else:
        return time_array[:final_step], gillespie_results[:, :final_step]


# time_array, gillespie_results = step_function(100000, x_0)

steps = 1_000_000_000

# time_array_2, gillespie_results_2 = step_function(steps, array([0, n]))
# time_array_1, gillespie_results_1 = step_function(steps, array([m, 0], dtype=float64))
# time_array_2, gillespie_results_2 = step_function(steps, array([1, 0], dtype=float64))

time_array_1, gillespie_results_1 = step_function(steps, array([m, 0, n], dtype=float64))
time_array_2, gillespie_results_2 = step_function(
    steps, array([1, 0, n - 1], dtype=float64)
)

times = [time_array_1, time_array_2]
gillespies = [gillespie_results_1, gillespie_results_2]

fig, ax = subplots(2, 2, figsize=(12, 10))


ran = range(len(gillespie_results_1[:, 0]))

initial_cond_1 = array([m, 0, n], dtype=float64)
initial_cond_2 = array([1, 0, n], dtype=float64)

conds = [initial_cond_1, initial_cond_2]

cabber_pillers = [ax[0, 0], ax[0, 1]]
for i, aks in enumerate(cabber_pillers):
    for piller in gillespies[i]:
        aks.step(times[i], piller, **bayes_kwargs)

    # Axis modifications after the plotting, to ensure it doesn't get overwritten
    aks.set_ylim(bottom=0)
    aks.set_xlim(left=0)

# for i in gillespie_results_1:
#     print(i)
# exit()
#

# ax[0, 0].step(time_array_1, gillespie_results_1[0, :], **bayes_kwargs)
# ax[0, 1].step(time_array_2, gillespie_results_2[0, :], **bayes_kwargs)
# ax[0, 0].step(time_array_1, gillespie_results_1[1, :], **bayes_kwargs)
# ax[0, 1].step(time_array_2, gillespie_results_2[1, :], **bayes_kwargs)
# ax[0, 0].step(time_array_1, gillespie_results_1[2, :], **bayes_kwargs)
# ax[0, 1].step(time_array_2, gillespie_results_2[2, :], **bayes_kwargs)

# fig1, ax1 = subplots()
# ax1.plot(gillespie_results_1[0, :], gillespie_results_1[2, :])


ax[1, 0].hist(gillespie_results_1[0, :], **hist_kwargs)
# ax[1, 0].hist(gillespie_results_1[1, :], **hist_kwargs)
ax[1, 0].set_ylim(bottom=0)
ax[1, 0].set_xlim(0, m)


ax[1, 1].hist(gillespie_results_2[0, :], **hist_kwargs)
# ax[1, 1].hist(gillespie_results_2[1, :], **hist_kwargs)
ax[1, 1].set_ylim(bottom=0)
ax[1, 1].set_xlim(0, m)


# hist, edges = histogram(gillespie_results_1[0, :], density=True, bins=boxes)
# ax[1, 0].plot(x_array * m, hist.max() * results)
# ax[1, 1].plot(x_array * m, hist.max() * results)
# ax[1, 0].plot(x_array * n, results / n)

# ax[1, 1].hist(gillespie_results_1[0, :], bins=boxes, density=True)

# ax[0, 1].hist2d(gillespie_results_1[0, :], gillespie_results_1[1, :], bins=(m, m))


show()
