import numpy as np
from numba import njit

from gillespie.parameter_class import ParameterClass


@njit
def main(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    p: ParameterClass,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a logistically restricted system that looks like the ode system
    """
    c1, c2, n = cs
    ct = c1 + c2

    a_1 = p.k1 * n * c1
    a_m1 = p.k1 * c1 * ct

    a_2 = p.k2 * n * c2
    a_m2 = p.k2 * c2 * ct

    a_3 = p.w1 * c1
    a_m3 = p.w2 * c2
    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])
