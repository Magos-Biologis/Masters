import numpy as np
from gillespie.stochastic_parameter_class import ParameterClass
from numba import njit

vj = [
    np.array([1, 0], dtype=np.int_),
    np.array([-1, 0], dtype=np.int_),
    np.array([0, 1], dtype=np.int_),
    np.array([-1, 0], dtype=np.int_),
    np.array([0, -1], dtype=np.int_),
    np.array([0, -1], dtype=np.int_),
]


@njit
def main(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    p: ParameterClass,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    Don't forget that that function assumes the substitutions
    k₁' := nk₁, k₂' := nk₂ , k₃' := bk₃, and k₄' := bk₄ are inside
    of the reaction matrix k
    """
    x, y = xs

    a_1 = p.k1p * x
    a_m1 = p.k_1 * x * x

    a_2 = p.k2p * x
    a_3 = p.k3p * x
    a_4 = p.k4p * y
    a_5 = p.k5 * y

    return np.array([a_1, a_m1, a_2, a_3, a_4, a_5])
