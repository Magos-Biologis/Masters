#!./.venv/bin/python
import os

import numba
from numba import njit
from pylab import *


from gillespie import gillespie


m = 100

k_1 = 1
# file_name += "_" + f"k1{k_1}"
k_1 /= m


k_2 = 3
# k_2 = 1
# file_name += "_" + f"k2{k_2}"
k_2 /= m


@njit
def aj(
    x: ndarray[tuple[int], dtype[float64 | int_]],
) -> ndarray[tuple[int], dtype[float64]]:
    a_1 = k_1 * x[0]
    # a_m1 =  k_1 * x[0] ** 2
    a_2 = k_2 * x[1]
    # a_m2 =  k_2 * x[1] ** 2
    return array([a_1, a_2], dtype=float64)


v1 = array([-1, 1], dtype=int_)
v2 = array([1, -1], dtype=int_)

# vm1 = array([1, -1], dtype=int_)
# vm2 = array([-1, 1], dtype=int_)

v = [v1, v2]


steps = 1000
ssa = gillespie(array([m, 0]), aj, v)

time_array, gillespie_results = ssa.generate()


plot(time_array, gillespie_results)
show()
