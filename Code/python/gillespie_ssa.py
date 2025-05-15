#!../.venv/bin/python3
import argparse
import os
import re
import time

# from old.gillespie import gillespie
# from gillespie import steppers as dgs
import gillespie as dgs

##
# os.sched_setaffinity(0, {3})
# import numba
import numpy as np
from gillespie.parameter_class import parameter_default
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
        "ode_3_2_alt",
        "ode_3_3",
        "ode_5_2",
        "ode_5_3",
    ],
    type=str,
)

parser.add_argument(
    "-ns",
    "--no-save",
    dest="save",
    action="store_false",
)


### As the type can be a function, we'll just assert it as an integer
def make_int(input: str) -> int:
    return int(input)


parser.add_argument(
    "-s",
    "--steps",
    dest="steps",
    help="Number of Steps",
    type=make_int,
    default=10_000,
)
parser.add_argument(
    "-n",
    "--size",
    dest="size",
    help="Total units in system",
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


def intify(string) -> list[int]:
    return [int(val) for val in string.split()]


def floatify(string) -> list[float]:
    return [float(val) for val in string.split()]


def parse_parameters(string):
    matches: list[tuple] = re.findall(r"(\w+[-m]?\d*'?)=(\S*)", string)
    parameters = [
        (
            str(key).replace("-", "_").replace("m", "_").replace("'", "p"),
            np.float64(value),
        )
        for key, value in matches
    ]
    return dict(parameters)


parser.add_argument(
    "-p",
    "--parameters",
    dest="parameters",
    help="Takes a string of equalities, and creates a dict from it",
    type=parse_parameters,
    default=parameter_default,
)


# parser.add_argument(
#     "-ks",
#     "--k-params",
#     dest="ks",
#     help="Test",
#     type=floatify,
# )
# parser.add_argument(
#     "-ws",
#     "--w-params",
#     dest="ws",
#     help="Test",
#     type=floatify,
#     default=[0.15, 0.15],
# )
# parser.add_argument(
#     "-ns",
#     "--n-params",
#     dest="ns",
#     help="Test",
#     type=intify,
#     default=[100, 90],
# )


args = parser.parse_args()


# print(args.parameters)
# exit()


## Compiling the defaults and the choice of parameters
model: str = args.model
step_count: int = args.steps


set_of_3 = [
    "5_3",
    "ode_3",
    "ode_3_3",
    "ode_5_3",
]

set_of_2 = [
    "2S",
    "2L",
    "5_2",
    "ode_2",
    "ode_2_2",
    "ode_3_2",
    "ode_3_2_alt",
    "ode_5_2",
]

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


if parameter_default != args.parameters:
    for key in args.parameters:
        parameter_default[key] = args.parameters[key]


parameters = dgs.ParameterClass(**parameter_default)


ks = np.array([parameters.k1, parameters.k2])
qs = np.array([parameters.q1, parameters.q2])

ns = np.array([parameters.n1, parameters.n2])
ws = np.array([parameters.w1, parameters.w2])

m_0 = parameters.m0


file_name: str = "ssa"
file_name += f"M{model}"
check_ode = re.compile(r"ode").search(model)

if check_ode is None:
    is_ode = False
else:
    is_ode = True

addons = []
addons.append(f"num={initial.sum()}")
if var_count == 2:
    addons.append("n={}".format(parameters.n))

alpha = 0
beta = m

if not is_ode:
    addons.extend(
        [
            "b={}".format(parameters.b),
            "k1={}".format(parameters.k1),
            "k-1={}".format(parameters.k_1),
            "k2={}".format(parameters.k2),
            # "km2{}".format(parameters.k_2),a
            "k3={}".format(parameters.k3),
            # "km3{}".format(parameters.k_3),
            "k4={}".format(parameters.k4),
            # "km4{}".format(parameters.k_4),
            "k5={}".format(parameters.k5),
            # "km5{}".format(parameters.k_5),
        ]
    )

    # k_sols = array([k_1 / (k_1 + k_2), k_2 / (k_1 + k_2)])

else:
    addons.extend(
        [
            "k1={}".format(parameters.k1),
            "k2={}".format(parameters.k2),
            "n1={}".format(parameters.n1),
            "n2={}".format(parameters.n2),
            "w1={}".format(parameters.w1),
            "w2={}".format(parameters.w2),
        ]
    )


file_name += "P"
for addon in addons:
    file_name += addon + "_"

if not is_ode and (var_count == 2):
    if parameters.b == parameters.n:
        para_version = r"$b=n$"
    elif parameters.b < parameters.n:
        para_version = r"$b<n$"
    elif parameters.b > parameters.n:
        para_version = r"$b>n$"
    else:
        para_version = "null"
        print("no numbers?")

    ## Simple if else so I can be lazy and not remember to add the naming conventions
    file_name += f"R{para_version[1:-1]}R"


t_0 = 0
dt = 0


max_val = 0
if parameters.k1 != parameters.k2:
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
steppy1 = dgs.ssa_stepper(model, initial, parameters)  # transitions, k)


date_time = time.time()

t0_1 = time.time()
time_results_1, gillespie_results_1 = steppy1.step_function(step_count)
t1_1 = time.time()

print("Stepper done \n\tTime taken: ", t1_1 - t0_1)

# print("Test Mode")
# exit()

# t0_2 = time.time()
# time_results_2, gillespie_results_2 = steppy2.step_function(step_count)
# t1_2 = time.time()
#
#
# print("Stepper two done \n\tTime taken: ", t1_2 - t0_2)


final_name1 = file_name + f"I{initial}C" + f"T{round(t1_1)}"
full_file_path1 = os.path.join(data_env, final_name1)

if save:
    np.savez(
        full_file_path1,
        time=time_results_1,
        states=gillespie_results_1,
    )

    print('saved as "{}"'.format(final_name1))
else:
    print('file name is "{}"'.format(final_name1))

# np.savez(
#     os.path.join(data_env, file_name + f"I{init2}C" + f"T{round(t1_2)}"),
#     time=time_results_2,
#     states=gillespie_results_2,
# )

exit()
