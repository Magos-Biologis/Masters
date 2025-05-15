import numpy as np
from numba import njit

from gillespie.parameter_class import ParameterClass


@njit
def aj_ode_3_2_alt(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    params: ParameterClass,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    # k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    """
    For a logistically restricted system that looks like the ode system
    """
    c1, c2 = cs

    k1 = params.k1
    k2 = params.k2

    w1 = params.w1
    w2 = params.w2

    a_1 = w1 * c1
    a_m1 = w2 * c2

    a_2 = k1 * params.n * c1
    a_m2 = k1 * c1 * c1

    a_3 = k2 * params.n * c2
    a_m3 = k2 * c2 * c2

    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])
