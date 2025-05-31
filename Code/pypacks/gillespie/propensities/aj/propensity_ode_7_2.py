import numpy as np
from gillespie.stochastic_parameter_class import ParameterClass
from numba import njit

vj = [
    np.array([-1, 1], dtype=np.int_),
    np.array([1, -1], dtype=np.int_),
    np.array([1, 0], dtype=np.int_),
    np.array([-1, 0], dtype=np.int_),
    np.array([0, 1], dtype=np.int_),
    np.array([0, -1], dtype=np.int_),
    np.array([-1, 0], dtype=np.int_),
    np.array([0, -1], dtype=np.int_),
    np.array([-1, 0], dtype=np.int_),
    np.array([0, 0], dtype=np.int_),
]


@njit
def main(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    p: ParameterClass,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a logistically restricted system that looks like the ode system
    It is presumed that the substitutions have been made.
    """
    c1, c2 = cs

    a_1 = p.w1 * c1
    a_m1 = p.w2 * c2

    a_2 = p.k1p * c1
    a_m2 = p.k1 * c1 * c1

    a_3 = p.k2p * c2
    a_m3 = p.k2 * c2 * c2

    a_4 = p.k1 * c1 * c2
    a_5 = p.k2 * c1 * c2

    a_6 = p.k1 * c1 * c2
    a_7 = p.k2 * c1 * c2

    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3, a_4, a_5, a_6, a_7])
