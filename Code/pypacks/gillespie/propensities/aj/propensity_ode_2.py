import numpy as np
from numba import njit


@njit
def aj_ode_2(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a logistically restricted system that looks like the ode system
    It is presumed that the substitutions have been made.
    """
    c1, c2 = cs
    ct = c1 + c2

    k1, k2 = k[0, :]
    k1p, k2p = k[1, :]
    k1n1, k2n2 = k[2, :]
    w1, w2 = k[3, :]

    # omega1 = k1 - w1
    # omega2 = k2 - w2

    a_1 = k1p * c1
    a_m1 = k1n1 * c1 * c1

    a_2 = k2p * c2
    a_m2 = k2n2 * c2 * c2

    a_3 = w1 * c1
    a_m3 = w2 * c2
    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])
