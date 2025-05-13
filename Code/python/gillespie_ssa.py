#!../.venv/bin/python3
import argparse
import os
import re
import time

##
# os.sched_setaffinity(0, {3})
# import numba
import numpy as np

# from old.gillespie import gillespie
from gillespie import steppers as dgs
from pylab import *

figure_env = str(os.getenv("FPE_FIGURE_ENV"))
data_env = str(os.getenv("THESIS_DATA_PATH"))

## Adding a cmdline parser so I can change the options without opening the file
## This is very much a brute force, no real elegance to be had
parser = argparse.ArgumentParser(
    prog="Gillespie Stepper",
    description="A python3 script that runs the gillespie SSA algorithm for various models I've made",
)

parser.add_argument(
    "model",
    # nargs=1,
    help="State which model to run",
    choices=[
        "2S",
        "2L",
        "5_2",
        "5_3",
        "ode_2",
        "ode_2_2",
        "ode_3",
        "ode_3_2",
        "ode_5_2",
        "ode_5_3",
    ],
    type=str,
)


### As the type can be a function, we'll just assert it as an integer
def make_int(input: str) -> int:
    return int(input)


parser.add_argument(
    "-st",
    "--steps",
    dest="steps",
    help="Number of Steps",
    type=make_int,
    default=10_000,
)
parser.add_argument(
    "-si",
    "--size",
    dest="size",
    help="Size of System",
    type=int,
    default=int(100),
)

parser.add_argument(
    "-ic",
    "--initial_conds",
    nargs="*",
    dest="initial_conds",
    help="Initial Conditions",
    type=int,
)
# default=[99, 1, 100],

parser.add_argument("-ks", "--parameters", dest="k", help="Test", type=float, default=1)


args = parser.parse_args()

# print(args.k1)

## Compiling the defaults and the choice of parameters


model: str = args.model
step_count: int = args.steps

set_of_3 = [
    "5_3",
    "ode_3",
    "ode_3_2",
    "ode_5_3",
]

set_of_2 = [
    "2S",
    "2L",
    "5_2",
    "ode_2",
    "ode_2_2",
    "ode_5_2",
]

# print(set(model))
# print(set(model).intersection(set(set_of_2)))
# exit()

# print(model)

if model in set_of_2:
    var_count = 2
elif model in set_of_3:
    var_count = 3
else:
    print("What?")
    exit()


m: int = args.size
inx, iny = m - 1, 1

initial = args.initial_conds
if initial is None:
    initial = [inx, iny, 0]

initial = np.array(initial, dtype=int_)[0:var_count]

# print(initial)
# exit()
# initial = np.array([*args.init], dtype=int_)[0:var_count]


alpha = 0
beta = m

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


file_name: str = "ssa"
file_name += f"M{model}"
check_ode = re.compile(r"ode").search(model)

if check_ode is None:
    is_ode = False
else:
    is_ode = True

addons = []
addons.append(f"num={initial.sum()}")
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

    addons.extend(
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

    addons.extend(
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


file_name += "S{:.0e}S".format(step_count)
steppy1 = dgs.ssa_stepper(model, initial, k)  # transitions, k)

# steppy2 = dgs.ssa_stepper(model, init2, k)  # , transitions


date_time = time.time()

t0_1 = time.time()
time_results_1, gillespie_results_1 = steppy1.step_function(step_count)
t1_1 = time.time()

print("Stepper done \n\tTime taken: ", t1_1 - t0_1)

# t0_2 = time.time()
# time_results_2, gillespie_results_2 = steppy2.step_function(step_count)
# t1_2 = time.time()
#
#
# print("Stepper two done \n\tTime taken: ", t1_2 - t0_2)


np.savez(
    os.path.join(data_env, file_name + f"I{initial}C" + f"T{round(t1_1)}"),
    time=time_results_1,
    states=gillespie_results_1,
)

# np.savez(
#     os.path.join(data_env, file_name + f"I{init2}C" + f"T{round(t1_2)}"),
#     time=time_results_2,
#     states=gillespie_results_2,
# )

exit()
