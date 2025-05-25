import numpy as np
from gillespie.parameter_class import ParameterClass
from numba import njit

vj = [
    np.array([-1, 1], dtype=np.int_),
    np.array([1, -1], dtype=np.int_),
]


@njit
def main(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    p: ParameterClass,
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a simple logistically restricted system
    """
    x, y = xs
    a_1 = p.k1p * x
    a_2 = p.k_1 * y * y

    return np.array([a_1, a_2])
