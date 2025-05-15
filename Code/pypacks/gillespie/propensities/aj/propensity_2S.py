import numpy as np
from numba import njit

from gillespie.parameter_class import ParameterClass


@njit
def aj_2S(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    p: ParameterClass,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For the most simple two part system
    """
    x, y = xs

    a_1 = p.k1 * x
    a_2 = p.k_1 * y

    return np.array([a_1, a_2])
