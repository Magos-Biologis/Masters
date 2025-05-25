# !./.venv/bin/python
import argparse
import os
import re

import numpy as np

# import sympy as sy
# import myodestuff
from myodestuff import ODEModel, ODEParameters

data_env = str(os.getenv("THESIS_DATA_PATH"))

parameter_default = {
    "k1": 1,
    "k2": 1,
    "n1": 100,
    "n2": 90,
    "w1": 0.015,
    "w2": 0.015,
    "q1": 0.85,
    "q2": 0.85,
    "m0": 0,
}

parser = argparse.ArgumentParser(
    prog="Gillespie Stepper",
    description="A python3 script that runs the gillespie SSA algorithm for various models I've made",
)


def parse_parameters(string):
    matches: list[tuple] = re.findall(r"(\w[^=]*)=\s?([^, ]*)", string)
    parameters = [
        (
            str(key).replace("-", "_").replace(" ", "").replace("'", "p"),
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

parser.add_argument(
    "-ic",
    "--initial_conds",
    nargs="*",
    dest="initial_conds",
    help="Initial Conditions",
    type=float,
)

args = parser.parse_args()


file_name: str = "ode"

model: str = "jesper"
# figure_path = "./figs"


three_d = False

### Variables
# Growth Constant, same because same cell
k = np.zeros(2)
n = np.zeros(2)
w = np.zeros(2)
q = np.zeros(2)

if parameter_default != args.parameters:
    parameter_default.update(args.parameters)


m_0 = parameter_default["m0"]

k[0] = parameter_default["k1"]
k[1] = parameter_default["k2"]

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = parameter_default["n1"]
n[1] = parameter_default["n2"]


w[0] = parameter_default["w1"]
w[1] = parameter_default["w2"]


q[0] = parameter_default["q1"]
q[1] = qm = parameter_default["q2"]


# "num={}".format(0),
filename_addendum = [
    "m0={}".format(m_0),
    "k1={}".format(k[0]),
    "k2={}".format(k[1]),
    "n1={}".format(n[0]),
    "n2={}".format(n[1]),
    "w1={}".format(w[0]),
    "w2={}".format(w[1]),
    "q1={}".format(q[0]),
    "q2={}".format(q[1]),
]


dt = 0.01
t_end = parameter_default.pop("t_end", 100)
# t_array = np.arange(0, t_end, dt)


if args.initial_conds is None:
    initial = np.array([n[0], 0, 0], dtype=float)
else:
    initial = np.array([*args.initial_conds, 0, 0], dtype=float)[:3]


parameters = ODEParameters(**parameter_default)
ode_model = ODEModel(parameters=parameters, t_range=(0, t_end), initial_condition=initial)


### Integration
t_array, sol = ode_model.integrate()

# roots = ode_model.roots(use_old=True)

# exit()
import time

t = time.time()

file_name += "M" + model
file_name += "P"
for entry in filename_addendum:
    file_name += entry + "_"
file_name += "S"
file_name += "IC"
file_name += "T" + "{}".format(t).replace(".", "")

full_file_path = os.path.join(data_env, file_name)

np.savez(
    full_file_path,
    time=t_array,
    solutions=sol,
)

print("data saved at {}".format(file_name + ".npz"))

exit()
