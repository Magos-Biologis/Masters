##
##
##

from importlib import import_module

import numpy as np
from numba import njit

from .analytical import simple_two_system
from .parameter_class import ParameterClass  # , parameter_default
from .propensities import transitions


@njit
def ssa_event(
    aj: np.ndarray[tuple[int], np.dtype[np.float64]],
) -> tuple[int, np.float64]:
    """
    Stealing Jesper's code and making it into Python

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


class ssa_stepper:
    step_module: str = "gillespie.steppers"
    steps: int = 10_000

    vj: dict[str, list[np.ndarray[tuple[int], np.dtype[np.int_]]]] = transitions

    def __init__(
        self,
        variation: str,
        x0: np.ndarray[tuple[int], np.dtype[np.int_]],
        params: ParameterClass,
    ) -> None:
        self.sub_module = self.step_module + ".step_function_" + variation
        self.state_changes = self.vj["vj_" + variation]

        self.x0 = x0
        self.params = params

    def step_function(
        self,
        steps: int | np.int_ = steps,
    ) -> tuple[
        np.ndarray[tuple[int], np.dtype[np.float64]],
        np.ndarray[tuple[int, int], np.dtype[np.int_]],
    ]:
        modu = import_module(self.sub_module)
        step_function = getattr(modu, "main")

        return step_function(steps, self.x0, self.state_changes, self.params)


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
