import numpy as np
from numba import njit


@njit
def aj_2S(
    xs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For the most simple two part system
    """
    x, y = xs
    a_1 = k[0, 0] * x
    a_2 = k[0, 1] * y
    return np.array([a_1, a_2])


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


@njit
def ode_2_aj(
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


@njit
def ode_2_2_aj(
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
    kn1, kn2 = k[2, :]
    w1, w2 = k[3, :]

    omega1 = k1 - w1
    omega2 = k2 - w2

    a_1 = omega1 * c1
    a_m1 = kn1 * c1 * c1

    a_2 = omega2 * c2
    a_m2 = kn2 * c2 * c2

    a_3 = w1 * c1
    a_m3 = w2 * c2
    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])


@njit
def ode_3_aj(
    cs: np.ndarray[tuple[int], np.dtype[np.int_]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    """
    For a logistically restricted system that looks like the ode system
    """
    c1, c2, u = cs
    ct = c1 + c2
    k1, k2 = k[1, :]
    n1, n2 = k[2, :]
    w1, w2 = k[3, :]

    a_1 = k1 * u * c1
    a_m1 = (k1 / n1) * c1 * ct

    a_2 = k2 * u * c2
    a_m2 = (k2 / n2) * c2 * ct

    a_3 = w1 * c1
    a_m3 = w2 * c2
    return np.array([a_1, a_m1, a_2, a_m2, a_3, a_m3])
