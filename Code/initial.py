#!./.venv/bin/python
import os
from sys import exit

from pprint import pprint as pp

import numpy as np
from matplotlib import pyplot as plt


dt = 0.001


# Growth Constant, same because same cell
k = k_1 = k_2 = 2

# Population cap (purely aesthetic)
n = n_1 = n_2 = 1

w_1 = 0.01
w_2 = 0.01

q_1 = 0.999
q_2 = 0.8


c_10 = 0.9 * n
c_20 = 0.1 * n
i_0 = 0.3


parameters = [
    k_1,
    k_2,
    w_1,
    w_2,
    q_1,
    q_2,
    n_1,
    n_2,
]

def dynamical_system(xs):

    c_1, c_2, i = xs
    c_t         = c_1 + c_2

    c_1_dot = c_1 * ( k_1 * ( 1 - c_t/n_1 ) - w_1 ) + w_2 * c_2 - q_1 * i * c_1
    c_2_dot = c_2 * ( k_1 * ( 1 - c_t/n_2 ) - w_2 ) + w_1 * c_1

    i_dot = - q_2 * i * c_2

    dx = np.array([ c_1_dot, c_2_dot, i_dot ])

    return dx


t_array = np.arange(0, 100, dt) 
steps   = []

init_cond = np.array([c_10, c_20, i_0])
steps.append(init_cond)


for i, _ in enumerate(t_array[0:-1]):
    current = steps[i]
    next    = dynamical_system(steps[i])

    value   = current + next * dt
    steps.append(value)


plt.plot(t_array, steps, label=["Non-resistant Type", "Resistant Type", "Antibiotic", ])
plt.legend()
plt.show()


#plt.savefig()
