#
## Stealing Jesper's code and making it into Python
##
# import numpy as np
#
# from numba import njit
# from dataclasses import dataclass


from gillespie.analytical import *
# from gillespie.propensities import *
# import gillespie.steppers as steppers
# from gillespie import *


# spec = [
#     ("x0", nb.float64[:]),
#     # ("propensity_function", nb.types.FunctionType(nb.types.float64[:](nb.types.int_[:]))),
#     ("vj", nb.int_[:]),
#     ("steps", int),
# ]
#
#
# @jitclass(spec)
class ssa:
    def __init__(
        self,
        x0: np.ndarray[tuple[int], np.dtype[np.float64]],
        vj: list[np.ndarray[tuple[int], np.dtype[np.int_]]],
        # propensity_tuple: tuple[
        #     function,
        #     list[np.ndarray[tuple[int], np.dtype[np.int_]]],
        # ],
        steps: int = 1000,
    ) -> None:
        self.x0 = x0
        self.v = vj

        self.steps = steps
        # self.props, self.v = propensity_tuple

        self.time = np.empty(shape=self.steps, dtype=np.float64)
        self.results = np.empty(shape=(2, self.steps), dtype=np.int_)

    def generate(
        self,
        aj,
    ) -> tuple[
        np.ndarray[tuple[int, int], np.dtype[np.float64]],
        np.ndarray[tuple[int], np.dtype[np.int_]],
    ]:
        self.results[:, 0] = self.x0
        x = self.x0.astype(np.float64)

        for i in range(1, self.steps):
            a_j = aj(self.results[:, i - 1])
            j, dt = ssa_event(a_j)

            if j == -1:
                break

            self.results[:, i] = np.add(self.results[:, i - 1], self.v[j])
            self.time[i] = self.time[i - 1] + dt

        final_step: int = i

        return self.time[:final_step], self.results[:, :final_step]


@njit
def ssa_event(
    aj: np.ndarray[tuple[int], np.dtype[np.float64]],
) -> tuple[int, np.float64]:
    """
    Gillespie event determination

    Input
    propensities: The probability of an event, eg. A+B → A+A  propensity is k*A*B

    Output
    eventId: The index of event as defined by the propensities
    dt: Time step

    Example
    >> # For reaction scheme  A+B → P and A → P using mass action
    >> k1 = 0.01; k2=1.0; B = 30; A = 80; t = 0.0;
    >> propensities = [k1*A*B, k2*A];  [eventid dt] = gillespie(prop);
    """

    r: np.ndarray[tuple[int], np.dtype[np.float64]] = np.random.rand(2)
    while r[0] == 0.0:
        r[0] = np.random.rand()

    j: int = -1
    a_0: np.float64 = aj.sum()

    if a_0 <= 0:
        tau: np.float64 = np.float64(0)
        print("womp womp")
        return j, tau

    tau: np.float64 = np.log(1 / r[0]) / a_0
    ra_0: np.float64 = r[1] * a_0

    for n, _ in enumerate(aj):
        j += 1  ## Cause its funny, why not

        s_j: np.float64 = aj[0 : n + 1].sum()

        if s_j >= ra_0:
            break

    return j, tau


# @dataclass(frozen=True)
# class transitions(object):
#     vj_4: list[np.ndarray] = [
#         np.array([1, 0], dtype=np.int_),
#         np.array([-1, 0], dtype=np.int_),
#         np.array([0, 1], dtype=np.int_),
#         np.array([-1, 0], dtype=np.int_),
#         np.array([0, -1], dtype=np.int_),
#     ]

# v1 =
# vm1 =
# v2 =
# v3 =
# v4 =

# def __post_init__(self) -> None:
#     pass
# self.vj_4 = [self.v1, self.vm1, self.v2, self.v3, self.v4]


# trans = transitions

# print()
# print(trans["vj_4"])
