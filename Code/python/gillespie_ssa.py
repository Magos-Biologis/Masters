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
        "ode_7_2",
        "ode_8_3",
    ],
    type=str,
)


parser.add_argument(
    "-r",
    "--repeats",
    dest="repeats",
    type=int,
    default=1,
)

bool_args = parser.add_argument_group(
    "Boolean arguments",
    "Arguments that exist as a boolean flag",
)
bool_args.add_argument(
    "-ns",
    "--no-save",
    dest="save",
    action="store_false",
)
bool_args.add_argument(
    "--test-parameters",
    dest="test",
    help="Breaks the file to print out the parameter flag",
    action="store_true",
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


parameter_parser = re.compile(r"(\w+[-m]?\d*'?)\s?=\s?([^ ,]*)")


def parse_parameters(string):
    matches: list[tuple] = parameter_parser.findall(string)
    parameters = [
        (
            str(key).replace("-", "_").replace("m", "_").replace("'", "p"),
            float(value),
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


args = parser.parse_args()


if args.test:
    print(args.parameters)
    exit()

# exit()

## Compiling the defaults and the choice of parameters
model: str = args.model
step_count: int = args.steps


set_of_3 = [
    "5_3",
    "ode_3",
    "ode_3_3",
    "ode_5_3",
    "ode_8_3",
]

set_of_2 = [
    "2S",
    "2L",
    "5_2",
    "ode_2",
    "ode_2_2",
    "ode_3_2",
    "ode_3_2_alt",
    "ode_7_2",
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

# initial = args.initial_conds
if args.initial_conds is None:
    initial = np.array([inx, iny, 0, 0, 0], dtype=int_)[:var_count]
else:
    initial = np.array([*args.initial_conds, 0, 0], dtype=int_)[:var_count]


# print(args.initial_conds)
# print(initial)
# exit()


if parameter_default != args.parameters:
    parameter_default.update(args.parameters)

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

alpha = 0
beta = m

if not is_ode:
    addons.extend(
        [
            "n={}".format(parameters.n),
            "b={}".format(parameters.b),
            "k1={}".format(parameters.k1),
            "k-1={}".format(parameters.k_1),
            "k2={}".format(parameters.k2),
            "k-2={}".format(parameters.k_2),
            "k3={}".format(parameters.k3),
            "k-3={}".format(parameters.k_3),
            "k4={}".format(parameters.k4),
            "k-4={}".format(parameters.k_4),
            "k5={}".format(parameters.k5),
            "k-5={}".format(parameters.k_5),
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
            "m0={}".format(parameters.m0),
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


## Step function


steppy1 = dgs.ssa_stepper(model, initial, parameters)  # transitions, k)


for _ in range(args.repeats):
    date_time = time.time()
    t0 = time.time()
    time_results, state_results = steppy1.step_function(step_count)
    t1 = time.time()

    print("Stepper done \n\tTime taken: ", t1 - t0)

    save_name = (
        file_name
        + "S{:.0e}S".format(len(time_results))
        + "I{}C".format(initial)
        + "T{}".format(t1).replace(".", "")
    )
    full_file_path = os.path.join(data_env, save_name)

    if args.save:
        np.savez(
            full_file_path,
            time=time_results,
            states=state_results,
        )

        print('saved as "{}"'.format(save_name))
    else:
        print('file name is "{}"'.format(save_name))

if args.repeats > 1:
    print("Done repeating")
# np.savez(
#     os.path.join(data_env, file_name + f"I{init2}C" + f"T{round(t1_2)}"),
#     time=time_results_2,
#     states=gillespie_results_2,
# )

exit()
