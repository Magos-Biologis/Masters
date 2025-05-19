import numpy as np
from numba import njit

from gillespie.parameter_class import ParameterClass


@njit
def main(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    p: ParameterClass,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    Don't forget that that function assumes the substitutions
    k₃' := k₃ * b and k₄' := k₄ * b are inside of the reaction matrix k
    """
    x, y, n = xs

    a_1 = p.k1 * n * x
    a_m1 = p.k_1 * x * x

    a_2 = p.k1 * n * x
    a_3 = p.k2 * x
    a_4 = p.k3p * y
    a_5 = p.k4p * y

    return np.array([a_1, a_m1, a_2, a_3, a_4, a_5])
