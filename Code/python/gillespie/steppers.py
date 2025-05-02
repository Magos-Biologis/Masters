import numpy as np
from numba import njit

# from gillespie import ssa_event
import gillespie as dg
import gillespie.propensities as dgp


class ssa_stepper:
    steps: int = 10_000

    vj: dict[str, list[np.ndarray[tuple[int], np.dtype[np.int_]]]] = dgp.transitions

    def __init__(
        self,
        x0: np.ndarray[tuple[int], np.dtype[np.int_]],
        k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        # v: list[np.ndarray[tuple[int], np.dtype[np.int_]]],
    ) -> None:
        self.x0 = x0
        self.k = k

    def step_function(
        self,
        # x0: np.ndarray[tuple[int], np.dtype[np.int_]],
        # v: list[np.ndarray[tuple[int], np.dtype[np.int_]]],
        # k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        variation: str,
        steps: int | np.int_ = steps,
    ) -> tuple[
        np.ndarray[tuple[int], np.dtype[np.float64]],
        np.ndarray[tuple[int, int], np.dtype[np.int_]],
    ]:
        return step_function_5_3(variation, steps, self.x0, self.vj[variation], self.k)


@njit
def step_function_5_3(
    variation: str,
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
    if variation == "5_2":
        function = dgp.aj_5_2
    else:
        function = dgp.aj_5_3

    for i in range(1, steps):
        x: np.ndarray[tuple[int], np.dtype[np.int_]] = gillespie_results[:, i - 1]
        a_j: np.ndarray[tuple[int], np.dtype[np.float64]] = function(x, scaled_k)
        j, dt = dg.ssa_event(a_j)

        if j == -1:
            break

        gillespie_results[:, i] = np.add(gillespie_results[:, i - 1], v[j])
        time_array[i] = time_array[i - 1] + dt

    final_step: int = i
    # print(final_step)
    return time_array[:final_step], gillespie_results[:, :final_step]
