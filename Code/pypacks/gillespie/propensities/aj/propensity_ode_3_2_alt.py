import numpy as np
from numba import njit


@njit
def aj_ode_3_2_alt(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    params: dict[str, np.float64],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    # k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    """
    For a logistically restricted system that looks like the ode system
    """
    c1, c2 = cs

    k1: np.float64 = params["k1"]
    k2: np.float64 = params["k2"]

    k1n1: np.float64 = params["k1n1"]
    k2n2: np.float64 = params["k2n2"]

    w1: np.float64 = params["w1"]
    w2: np.float64 = params["w2"]

    # n1, n2 = k[4, :]

    a_1 = w1 * c1
    a_m1 = w2 * c2

    a_2 = k1n1 * c1
    a_m2 = k1 * c1 * c1

    a_3 = k2n2 * c2
    a_m3 = k2 * c2 * c2

    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])
