import numpy as np
from numba import njit


@njit
def aj_5_3(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    Don't forget that that function assumes the substitutions
    k₃' := k₃ * b and k₄' := k₄ * b are inside of the reaction matrix k
    """
    x, y, n = xs
    a_1 = k[0, 0] * n * x
    a_m1 = k[0, 1] * x * x
    a_2 = k[1, 0] * n * x
    a_3 = k[2, 0] * x
    a_4 = k[3, 0] * y
    a_5 = k[4, 0] * y
    return np.array([a_1, a_m1, a_2, a_3, a_4, a_5])
