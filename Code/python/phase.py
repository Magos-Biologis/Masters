#!../.venv/bin/python3
import argparse
import json
import os
import re
from copy import deepcopy as dc

import numpy as np

# from numpy.polynomial import Polynomial
from myodestuff import ODEModel, ODEParameters

figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
phase_path = os.path.join(figure_env, "phase")

data_env = str(os.getenv("THESIS_DATA_PATH"))


import argparse

parameter_default = {
    "k1": 1,
    "k2": 1,
    "n1": 100,
    "n2": 90,
    "w1": 0.15,
    "w2": 0.15,
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

args = parser.parse_args()


file_name: str = "phase"

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

metadata_dict = dc(parameter_default)

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


omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]

k_1 = k[0] / n[0]
k_2 = k[1] / n[1]

c1_min = w[1] / k_1
c2_min = w[0] / k_2

c1_max = omega_1 / k_1
c2_max = omega_2 / k_2


alpha = 1
c1_0 = n[0] - alpha
c2_0 = alpha

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


file_name += "M" + model
file_namestring = (
    "{}".format(m_0)
    + "{}".format(k[0])
    + "{}".format(k[1])
    + "{}".format(n[0])
    + "{}".format(n[1])
    + "{}".format(w[0])
    + "{}".format(w[1])
    + "{}".format(q[0])
    + "{}".format(q[1])
).replace(".", "")
file_name += "P{}".format(file_namestring)


# print(file_namestring)
# exit()


# file_name += "T" + "{}".format(0)


dt = 0.01
t_end = 250
t_array = np.arange(0, t_end, dt)

parameters = ODEParameters(**parameter_default)
# parameters2 = parameter_class(m=2, m_0=m_0, k=k1, n=n1, q=q, w=w)
# parameters3 = parameter_class(m=2, m_0=0, k=k, n=n2, q=q, w=w)


init_conds = np.array([c1_0, c2_0, m_0])
init_conds2 = np.array([c2_0, c1_0, 0])
init_conds3 = np.array([c1_0, c2_0, m_0])
# init_conds3 = np.array([100, 100, 100])

# if three_d:
# else:
#     init_conds1 = np.array([c1_0, c2_0,])
#     init_conds2 = np.array([c2_0, c1_0])
#     init_conds3 = np.array([100, 100])

ode_model = ODEModel(parameters, t_range=(0, t_end), initial_condition=init_conds)


### Integration

# parameters1 = parameter_class(*parameters)

t_array, sol1 = ode_model.integrate()
# t_array, sol2 = model2.integrate()
# t_array, sol3 = model3.integrate()

c1_root, c2_root, m_root = ode_model.roots()

## Plotting

x_lims = 0, 100
y_lims = 0, 100

# x_lims = 0, 1
# y_lims = 0, 1


c1s = np.arange(0.1, 125, 0.5)
c2s = np.arange(0.1, 125, 0.5)
ms = np.arange(0.0, m_0 + 1, 0.05)


c1, c2, m = np.meshgrid(c1s, c2s, ms)
c1, c2 = np.meshgrid(c1s, c2s)


# print(c1[0, :, 0])
# exit()


def vector_space(c1, c2):
    m = m_0 / (c2 * q[1])

    dc1 = (k[0] * (1 - (c1 + c2) / n[0]) - w[0]) * c1 + w[1] * c2 - q[0] * m * c1
    dc2 = (k[1] * (1 - (c1 + c2) / n[1]) - w[1]) * c2 + w[0] * c1

    return (dc1, dc2)


# vector_field = vector_space(c1, c2, m)
vector_field = vector_space(c1, c2)

dU = vector_field[0]
dV = vector_field[1]

# final_name = file_name
full_file_path = os.path.join(data_env, file_name)

# print(c1.shape, c2.shape)
# print(dU.shape, dV.shape)
# exit()

# q[0] * (m_0 / (q[1] * c2)) * c1 -


margin = 0.01
c1s = np.linspace(c1_min, c1_max, 201, endpoint=False)[1:]
c2s = np.linspace(c2_min, c2_max, 201, endpoint=False)[1:]


def c1_sol(c1):
    num = -c1 * (omega_1 - k_1 * c1)
    den = w[1] - k_1 * c1
    return num / den


def c2_sol(c2):
    num = -c2 * (omega_2 - k_2 * c2)
    den = w[0] - k_2 * c2
    return num / den


c1_nullcline = np.array([c2_sol(c2s), c2s])
c2_nullcline = np.array([c1s, c1_sol(c1s)])

#
import time

t = time.time()

metadata_dict.update(
    {
        "data_source": "phase space",
        "model": "simple compartment",
        "date": t,
        "t_end": t_end,
    }
)

#
# c1_nullcline = ode_model._find_c1_nullcline(mesh=c_mesh)
# c2_nullcline = ode_model._find_c2_nullcline(c2s)

c_mesh = np.meshgrid(c1s, c2s)
cs_levelset = ode_model._find_level_set(mesh=c_mesh)


# from .metadata_json import Numpy2Native


class Numpy2Native(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


metadata_json = json.dumps(metadata_dict, cls=Numpy2Native, ensure_ascii=False).encode(
    "utf-8"
)

np.savez(
    full_file_path,
    c1=c1,
    c2=c2,
    dU=dU,
    dV=dV,
    c1_nullcline=c1_nullcline,
    c2_nullcline=c2_nullcline,
    level=cs_levelset,
    metadata=metadata_json,
)
exit()
