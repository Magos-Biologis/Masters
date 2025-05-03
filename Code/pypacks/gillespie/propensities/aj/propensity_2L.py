import numpy as np
from numba import njit


@njit
def aj_2L(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a simple logistically restricted system
    """
    x, y = xs
    a_1 = k[0, 0] * x
    a_2 = k[0, 1] * y * y
    return np.array([a_1, a_2])
