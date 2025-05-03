import numpy as np
from numba import njit


@njit
def aj_ode_2(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    # ks: np.ndarray[tuple[int], np.dtype[np.float64]],
    # ns: np.ndarray[tuple[int], np.dtype[np.float64]],
    # ws: np.ndarray[tuple[int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a logistically restricted system that looks like the ode system
    """
    c1, c2 = cs
    ct = c1 + c2
    # k1n, k2n = k[0, :]
    k1, k2 = k[1, :]
    n1, n2 = k[2, :]
    w1, w2 = k[3, :]

    k1n = k1 / n1
    k2n = k2 / n2

    omega1 = k1 - w1
    omega2 = k2 - w2

    a_1 = omega1 * c1
    a_m1 = k1n * c1 * c1

    a_2 = omega2 * c2
    a_m2 = k2n * c2 * c2

    a_3 = w1 * c1
    a_m3 = w2 * c2
    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])
