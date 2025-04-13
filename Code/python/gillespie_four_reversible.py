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

m = 100

n = 1
b = 20000000000000

n /= m
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

k_5 = 1
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
    aj: ndarray[tuple[int], dtype[float64]] | list[float],
    # seed: int = 1984,
) -> tuple[int, float]:
    # assert type(aj) is ndarray

    if aj.min() <= 0.0:
        j: int = -1
        tau: float = 0
        return j, tau

    # np.random.seed(seed)
    r = rand(2)
    while r[0] == 0:
        r[0] = rand()

    a_0: float = aj.sum()

    j: int = -1
    tau: float = log(1 / r[0]) / a_0

    for k, _ in enumerate(aj):
        s_j: float = aj[: k + 1].sum()
        ra_0: float = r[1] * a_0

        if s_j > ra_0:
            j: int = k
            break

    return j, tau


# from old.gillespie import gillespie


@njit
def aj(x: ndarray[tuple[int], dtype[float64]]) -> ndarray[tuple[int], dtype[float64]]:
    x, y = x

    return array(
        [
            k[0] * n * x,
            k[1] * x**2,
            k[2] * n * x,
            k[3] * b * x,
            k[4] * b * y,
            k[5] * b * y,
        ]
    )


## Step function
@njit
def step_function(
    steps: int, x0: ndarray[tuple[int], dtype[float64]]
) -> tuple[ndarray[tuple[int], dtype[float64]], ndarray[tuple[int], dtype[int16]]]:
    ## v_j
    v1 = array([1, 0], dtype=int16)
    v2 = array([0, 1], dtype=int16)
    v3 = array([-1, 0], dtype=int16)
    v4 = array([0, -1], dtype=int16)
    v5 = array([0, 0], dtype=int16)

    vm1 = array([-1, 0], dtype=int16)

    v = [v1, vm1, v2, v3, v4, v5]

    ## Other
    gillespie_results = empty((2, steps), dtype=int16)
    time_array = empty(steps, dtype=float64)
    x = x0.astype(float64)
    gillespie_results[:, 0] = x

    # while j != -1:
    for i in range(steps - 1):
        ai = aj(x)
        j, dt = gillespie(ai)
        if j != -1:
            x[:] = add(x, v[j])  ## Using numpy for speed I guess
        else:
            x[:] = x

        if x[0] < 1:
            x[0] = 1

        time_array[i + 1] = time_array[i] + dt
        gillespie_results[:, i + 1] = x

    return time_array, gillespie_results


# time_array, gillespie_results = step_function(100000, x_0)

steps = 2000000

time_array_1, gillespie_results_1 = step_function(steps, array([m - 1, 1]))
# time_array_2, gillespie_results_2 = step_function(steps, array([0, n]))
time_array_2, gillespie_results_2 = step_function(steps, array([1, m - 1]))


fig, ax = subplots(2, 2, figsize=(12, 10))

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
