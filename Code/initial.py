#!./.venv/bin/python
import os
from sys import exit

from pprint import pprint as pp

import numpy as np
from matplotlib import pyplot as plt


dt = 0.001
t_end = 100

t_array = np.arange(0, t_end, dt) 
steps   = []


# Growth Constant, same because same cell
k = k_1 = k_2 = 2

# Population cap (purely aesthetic)
n = n_1 = n_2 = 1

w_1 = 0.01
w_2 = 0.01

q_1 = 0.999
q_2 = 0.8


c_10 = 1.0 * n
c_20 = 0.0 * n
i_0 = 0.02

c_t0 = c_10 + c_20

init_cond = np.array([n, 0, 0])
steps.append(init_cond)

omega_1 = k_1*(1-c_t0/n)-w_1
omega_2 = k_2*(1-c_t0/n)-w_2

#trace = (i_0 * q_1 - omega_2 - omega_1)
#dete  = (omega_1 - i_0 * q_1) * omega_2 - w_1 * w_2

#disc  = trace**2 - dete

c1_fixed = -((w_1 / k_1 - 1) * n_2 + c_20)


parameters = [
    k_1, k_2,
    w_1, w_2,
    q_1, q_2,
    n_1, n_2,
]

#dd_x = lambda c1: lambda c2: lambda I: k_1*c1 
def d_c1(c1, c2, I): 
    ct   = c1 + c2
    diff = c1 * ( k_1 * ( 1 - ct/n_1 ) - w_1 ) + w_2 * c2 - q_1 * I * c1
    return diff

def d_c2(c1, c2, I): 
    ct   = c1 + c2
    diff = c2 * ( k_2 * ( 1 - ct/n_2 ) - w_2 ) + w_1 * c1
    return diff

def d_i(c1, c2, I): 
    diff = - q_2 * I * c2 + i_0
    return diff


def dynamical_system(xs):

    c_1, c_2, i = xs

    dx = np.array([
        d_c1(c_1, c_2, i),
        d_c2(c_1, c_2, i),
        d_i(c_1, c_2, i)
        ])

    return dx




for j, _ in enumerate(t_array[0:-1]):
    current = steps[j]
    difference = dynamical_system(steps[j])

    value = current + difference * dt
    steps.append(value)


names = ["Non-resistant Type", "Resistant Type", "Antibiotic", ]

fig, ax = plt.subplots()

ax.plot(t_array, steps, label=names)
#ax.hlines(c1_fixed, 0, t_end)

#[((-dete + np.sqrt(dete))/2),((-dete - np.sqrt(dete))/2)]

plt.legend()
plt.show()


#plt.savefig()
