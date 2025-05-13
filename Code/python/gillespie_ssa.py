#!../.venv/bin/python3
import os
import re

##
# os.sched_setaffinity(0, {3})
# import numba
import numpy as np

# from old.gillespie import gillespie
from gillespie import steppers as dgs
from pylab import *

# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_env = str(os.getenv("FPE_FIGURE_ENV"))
data_env = str(os.getenv("THESIS_DATA_PATH"))

# file_dir = os.path.join(figure_env, "five_var")
# file_name = "multipage_2var_fpe_ssa"


# file_name = "simulated_fpe"


model: str = "5_3"
var_count = 3


# steps = 100
# steps = 100_000
steps = 1_000_000
# steps = 10_000_000
# steps = 100_000_000
# steps = 1_000_000_000

m = 100
# margin = 0.5

alpha = 0
beta = m

inx1, iny1 = m - 1, 1
inx2, iny2 = 1, m - 1

init1 = array([inx1, iny1, m], dtype=int_)[0:var_count]
init2 = array([inx2, iny2, m], dtype=int_)[0:var_count]


b = 100
n = 100
# b = 1

u = 1

k = zeros(shape=(5, 2), dtype=float64)

ks = np.empty(shape=2)
ns = np.empty(shape=2)
ws = np.empty(shape=2)
qs = np.empty(shape=2)

m_0 = 0

ks[0] = k_1 = 0.7
ks[1] = k_2 = 0.6

ns[0] = n_1 = 100
ns[1] = n_2 = 90

ws[0] = w_1 = 0.15
ws[1] = w_2 = 0.15

qs[0] = q_1 = 0.85
qs[1] = q_2 = 0.85


file_name = "ssa"
file_name += f"M{model}"
check_ode = re.compile(r"ode").search(model)

if check_ode is None:
    is_ode = False
else:
    is_ode = True

addons = []
addons.append(f"num={init1.sum()}")
if not is_ode:
    k[0, 0] = k_1 = 1  # 0.55
    k[0, 1] = k_m1 = 1  # 0.55

    k[1, 0] = k_2 = 1
    k[2, 0] = k_3 = 1
    k[3, 0] = k_4 = 1
    k[4, 0] = k_5 = 1

    if var_count == 2:
        k[0:2, 0] *= n
        k[2:4, 0] *= b

        addons.append(f"n={n}")

    if var_count == 3:
        k[2:4, 0] *= b

    # print(k)
    # exit()

    # anal_sol = -(k[2, 0] - k[0, 0]) / k[0, 1]

    addons.append(
        [
            f"b={b}",
            f"k1={k_1}",
            f"km1={k_m1}",
            f"k2={k_2}",
            # f"km2{k_m2}",a
            f"k3={k_3}",
            # f"km3{k_m3}",
            f"k4={k_4}",
            # f"km4{k_m4}",
            f"k5={k_5}",
            # f"km5{k_m5}",
        ]
    )

    k_sols = array([k_1 / (k_1 + k_2), k_2 / (k_1 + k_2)])

else:
    # k[0, :] *= u

    k[0, :] = ks
    k[1, :] = ks * u
    k[2, :] = ks / ns
    k[3, :] = ws

    k[4, :] = ns

    addons.append(
        [
            f"num={m}",
            f"u={u}",
            f"k1={k_1}",
            f"k2={k_2}",
            f"n1={n_1}",
            f"n2={n_2}",
            f"w1={w_1}",
            f"w2={w_2}",
        ]
    )


file_name += "P"
for addon in addons:
    file_name += addon + "_"

if not is_ode and (var_count == 2):
    if b == n:
        para_version = r"$b=n$"
    elif b < n:
        para_version = r"$b<n$"
    elif b > n:
        para_version = r"$b>n$"
    else:
        para_version = "null"
        print("no numbers?")

    ## Simple if else so I can be lazy and not remember to add the naming conventions
    file_name += f"R{para_version[1:-1]}R"


# print(file_name)
# exit()
# file_name += "_" + f"k1{k_1}"
# file_name += "_" + f"k2{k_2}"


# k_2 = 1
# file_name += "_" + f"w1{w_1}"
# file_name += "_" + f"w2{w_2}"
# file_name += "_" + f"n1{n_1}"
# file_name += "_" + f"n2{n_2}"


# k = array([k_1, k_2])

# file_name += "_" + f"n{m}"


a_1 = m
a_2 = 0


aj_0 = array([k_1 * a_1, k_2 * a_2])


t_0 = 0
dt = 0


max_val = 0
if k_1 != k_2:
    if m == 1:
        max_val = log(2.022121436749997)
    elif m == 10:
        max_val = log(674763.2054860857)
    elif m == 100:
        max_val = log(6.367225127715182e43)
else:
    pass

x_array = linspace(0, 1, num=500, endpoint=False)[1:]


## Step function


# time_array, gillespie_results = step_function(100000, x_0)


file_name += "S{:.0e}S".format(steps)

# print(file_name)
# exit()

init_conds = [init1, init2]

import time

# time_array_1, gillespie_results_1 = steppy.step_function(steps)


# transitions = dgp.transitions["vj_" + model]

steppy1 = dgs.ssa_stepper(model, init1, k)  # transitions, k)
steppy2 = dgs.ssa_stepper(model, init2, k)  # , transitions


date_time = time.time()

t0_1 = time.time()
time_results_1, gillespie_results_1 = steppy1.step_function(steps)
t1_1 = time.time()

print("Stepper one done \n\tTime taken: ", t1_1 - t0_1)

t0_2 = time.time()
time_results_2, gillespie_results_2 = steppy2.step_function(steps)
t1_2 = time.time()


print("Stepper two done \n\tTime taken: ", t1_2 - t0_2)


np.savez(
    os.path.join(data_env, file_name + f"I{init1}C" + f"T{round(t1_1)}"),
    time=time_results_1,
    states=gillespie_results_1,
)

np.savez(
    os.path.join(data_env, file_name + f"I{init2}C" + f"T{round(t1_2)}"),
    time=time_results_2,
    states=gillespie_results_2,
)

exit()
