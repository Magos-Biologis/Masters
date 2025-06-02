import numpy as np
from gillespie.stochastic_parameter_class import ParameterClass
from numba import njit

vj = [
    np.array([1, 0, -1], dtype=np.int_),
    np.array([-1, 0, 1], dtype=np.int_),
    np.array([0, 1, -1], dtype=np.int_),
    np.array([-1, 0, 1], dtype=np.int_),
    np.array([0, -1, 1], dtype=np.int_),
    np.array([0, -1, 1], dtype=np.int_),
]


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

    a_1 = p.k1 * n * (x - 1)
    a_m1 = p.k_1 * (x - 1) ** 2

    a_2 = p.k1 * n * (x - 1)
    a_3 = p.k2 * (x - 1)
    a_4 = p.k3p * y
    a_5 = p.k4p * y

    return np.array([a_1, a_m1, a_2, a_3, a_4, a_5])
