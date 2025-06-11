#!./.venv/bin/python
import os

import numba
from numba import njit
from pylab import *

import gillespie as dg

figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_path = os.path.join(figure_env, "fpe")

# file_name = "ode_solution"

file_name = "simulated_fpe"

print()

m: int = 100
n = 1500
# n /= m

b = 0
# b = 150
# b = 1000
# b = 10000

b /= m + n


alpha = 0
beta = 1


boxes = arange(0, m + 1, 1)
# boxes = n+1

hist_kwargs = {
    "bins": boxes,
    "density": True,
    "edgecolor": "black",
    "alpha": 0.7,
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


k = zeros(shape=(5, 2), dtype=float64)


k[0, 0] = k_1 = 1  # 0.55
k[0, 1] = k_m1 = 1  # 0.55
file_name += "_" + f"k1{k_1}"


k[1, 0] = k_2 = 1
file_name += "_" + f"k2{k_2}"

k[2, 0] = k_3 = 1
file_name += "_" + f"k3{k_3}"

k[3, 0] = k_4 = 1
file_name += "_" + f"k4{k_4}"

k[4, 0] = k_5 = 0.5
file_name += "_" + f"k5{k_5}"


# print(k)
# exit()

k = divide(k, m + n)
file_name += "_" + f"n{m}"


# @njit
# def aj(xs: ndarray[tuple[int], dtype[float64]]) -> ndarray[tuple[int], dtype[float64]]:
#     x = xs[0]  # - 1
#     y = xs[1]
#
#     a_1 = k[0] * n * x
#     a_m1 = k[1] * x**2
#     a_2 = k[2] * n * x
#     a_3 = k[3] * b * x
#     a_4 = k[4] * b * y
#
#     return array([a_1, a_m1, a_2, a_3, a_4])


@njit
def aj(xs: ndarray[tuple[int], dtype[int_]]) -> ndarray[tuple[int], dtype[float64]]:
    x, y, n = xs

    a_1 = k[0, 0] * n * x
    a_m1 = k[0, 1] * x**2
    a_2 = k[1, 0] * n * x
    a_3 = k[2, 0] * b * x
    a_4 = k[3, 0] * b * y
    a_5 = k[4, 0] * y

    return array([a_1, a_m1, a_2, a_3, a_4, a_5])


## Step function
@njit
def step_function(
    steps: int,
    x0: ndarray[tuple[int], dtype[float64]],
    v: list[ndarray[tuple[int], dtype[int_]]],
) -> tuple[ndarray[tuple[int], dtype[float64]], ndarray[tuple[int], dtype[int_]]]:
    ## v_j

    ## Other
    gillespie_results = empty((len(x0), steps), dtype=int_)
    time_array = empty(steps, dtype=float64)

    gillespie_results[:, 0] = x0

    for i in range(1, steps):
        x: ndarray[tuple[int], dtype[int_]] = gillespie_results[:, i - 1]
        a_j = aj(x)
        j, dt = dg.ssa_event(a_j)

        if j == -1:
            break

        gillespie_results[:, i] = add(gillespie_results[:, i - 1], v[j])
        time_array[i] = time_array[i - 1] + dt

    final_step: int = i
    return time_array[:final_step], gillespie_results[:, :final_step]


# time_array, gillespie_results = step_function(100000, x_0)

steps = 1_000_000
# steps = 10_000_000
# steps = 100_000_000
# steps = 1_000_000_000


v1 = array([1, 0], dtype=int_)
vm1 = array([-1, 0], dtype=int_)
v2 = array([0, 1], dtype=int_)
v3 = array([-1, 0], dtype=int_)
v4 = array([0, -1], dtype=int_)
vj_4 = [v1, vm1, v2, v3, v4]


v1 = array([1, 0, -1], dtype=int_)
vm1 = array([-1, 0, 1], dtype=int_)
v2 = array([0, 1, -1], dtype=int_)
v3 = array([-1, 0, 1], dtype=int_)
v4 = array([0, -1, 1], dtype=int_)
v5 = array([0, -1, 1], dtype=int_)
vj_5 = [v1, vm1, v2, v3, v4, v5]


init1 = array([m, 0, n], dtype=float64)
init2 = array([1, 0, n + (m - 1)], dtype=float64)


time_array_1, gillespie_results_1 = step_function(steps, init1, vj_5)
time_array_2, gillespie_results_2 = step_function(steps, init2, vj_5)

times = [time_array_1, time_array_2]
gillespies = [gillespie_results_1, gillespie_results_2]

fig, ax = subplots(2, 2, figsize=(12, 10))


ran = range(len(gillespie_results_1[:, 0]))


cabber_pillers = [ax[0, 0], ax[0, 1]]
for i, aks in enumerate(cabber_pillers):
    for piller in gillespies[i]:
        aks.step(times[i], piller, **bayes_kwargs)

    # Axis modifications after the plotting, to ensure it doesn't get overwritten
    aks.set_ylim(bottom=0)
    aks.set_xlim(left=0)

# for aks in zip(cabber_pillers):
#     for piller in gillespies[0]:
#         aks.step(times[0], piller, **bayes_kwargs)
#
#     # Axis modifications after the plotting, to ensure it doesn't get overwritten
#     aks.set_ylim(bottom=0)
#     aks.set_xlim(left=0)

ant_hills = [ax[1, 0], ax[1, 1]]

ax[1, 0].hist(gillespie_results_1[0, :], **hist_kwargs, label="GCN")
ax[1, 1].hist(gillespie_results_2[0, :], **hist_kwargs, label="GCN")

# ax[1, 0].hist(gillespie_results_1[1, :], **hist_kwargs, label="Antibodies")
# ax[1, 1].hist(gillespie_results_2[1, :], **hist_kwargs, label="Antibodies")

for aks in ant_hills:
    aks.set_ylim(bottom=0)
    aks.set_xlim(left=0, right=m)

# hist, edges = histogram(gillespie_results_1[0, :], density=True, bins=boxes)
# ax[1, 0].plot(x_array * m, hist.max() * results)
# ax[1, 1].plot(x_array * m, hist.max() * results)
# ax[1, 0].plot(x_array * n, results / n)

# ax[1, 1].hist(gillespie_results_1[0, :], bins=boxes, density=True)

# ax[0, 1].hist2d(gillespie_results_1[0, :], gillespie_results_1[1, :], bins=(m, m))


show()
