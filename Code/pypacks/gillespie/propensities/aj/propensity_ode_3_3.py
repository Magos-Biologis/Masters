import numpy as np
from numba import njit


@njit
def aj_ode_3_3(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a logistically restricted system that looks like the ode system
    """
    c1, c2, n = cs
    # ct = c1 + c2

    k1, k2 = k[0, :]
    k1n1, k2n2 = k[2, :]
    w1, w2 = k[3, :]

    # n1, n2 = k[4, :]
    # k1n1, k2n2 = k1 / n1, k2 / n2

    a_1 = w1 * c1
    a_m1 = w2 * c2

    a_2 = k1 * c1 * n
    a_m2 = k1n1 * c1 * c2

    a_3 = k2 * c2 * n
    a_m3 = k2n2 * c2 * c2

    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])
