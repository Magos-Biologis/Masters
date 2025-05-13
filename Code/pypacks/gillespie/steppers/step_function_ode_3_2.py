import numpy as np
from numba import njit

from gillespie import ssa_event
from gillespie.propensities import aj_ode_3_2 as propensities


@njit
def main(
    steps: int | np.int_,
    x0: np.ndarray[tuple[int], np.dtype[np.int_]],
    v: list[np.ndarray[tuple[int], np.dtype[np.int_]]],
    k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
) -> tuple[
    np.ndarray[tuple[int], np.dtype[np.float64]],
    np.ndarray[tuple[int, int], np.dtype[np.int_]],
]:
    m: int = int(sum(x0))
    scaled_k: np.ndarray[tuple[int, int], np.dtype[np.float64]] = np.divide(k, m)

    ## Other
    gillespie_results = np.empty((len(x0), steps), dtype=np.int_)
    time_array = np.empty(steps, dtype=np.float64)

    gillespie_results[:, 0] = x0

    for i in range(1, steps):
        x: np.ndarray[tuple[int], np.dtype[np.int_]] = gillespie_results[:, i - 1]
        a_j: np.ndarray[tuple[int], np.dtype[np.float64]] = propensities(x, scaled_k)
        j, dt = ssa_event(a_j)

        if j == -1:
            break

        gillespie_results[:, i] = np.add(gillespie_results[:, i - 1], v[j])
        time_array[i] = time_array[i - 1] + dt

    final_step: int = i
    return time_array[:final_step], gillespie_results[:, :final_step]
