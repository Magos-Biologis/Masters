import numpy as np
from numba import njit


@njit
def aj_5_2(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    Don't forget that that function assumes the substitutions
    k₁' := nk₁, k₂' := nk₂ , k₃' := bk₃, and k₄' := bk₄ are inside
    of the reaction matrix k
    """
    x: np.int_ = xs[0]
    y: np.int_ = xs[1]
    a_1: np.float64 = k[0, 0] * x
    a_m1: np.float64 = k[0, 1] * x * x
    a_2: np.float64 = k[1, 0] * x
    a_3: np.float64 = k[2, 0] * x
    a_4: np.float64 = k[3, 0] * y
    a_5: np.float64 = k[4, 0] * y
    return np.array([a_1, a_m1, a_2, a_3, a_4, a_5])
