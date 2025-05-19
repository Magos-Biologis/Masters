##
##
##

from importlib import import_module

import numpy as np
from numba import njit

from .analytical import simple_two_system
from .parameter_class import ParameterClass  # , parameter_default
from .propensities import transitions

# from .step_function import main as stepper


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


@njit
def step_function(
    steps: int | np.int_,
    x0: np.ndarray[tuple[int], np.dtype[np.int_]],
    v: list[np.ndarray[tuple[int], np.dtype[np.int_]]],
    k: ParameterClass,
    propensities,
) -> tuple[
    np.ndarray[tuple[int], np.dtype[np.float64]],
    np.ndarray[tuple[int, int], np.dtype[np.int_]],
]:
    m: int = x0.sum()
    k.divide(m)

    gillespie_results = np.empty((len(x0), steps), dtype=np.int_)
    time_array = np.empty(steps, dtype=np.float64)

    gillespie_results[:, 0] = x0

    for i in range(1, steps):
        x: np.ndarray[tuple[int], np.dtype[np.int_]] = gillespie_results[:, i - 1]
        a_j: np.ndarray[tuple[int], np.dtype[np.float64]] = propensities(x, k)
        j, dt = ssa_event(a_j)

        if j == -1:
            break

        gillespie_results[:, i] = np.add(gillespie_results[:, i - 1], v[j])
        time_array[i] = time_array[i - 1] + dt

    final_step: int = i
    return time_array[:final_step], gillespie_results[:, :final_step]


class ssa_stepper:
    step_module: str = "gillespie.steppers"
    propensity_module: str = "gillespie.propensities.aj"
    steps: int = 10_000

    vj: dict[str, list[np.ndarray[tuple[int], np.dtype[np.int_]]]] = transitions

    def __init__(
        self,
        model: str,
        x0: np.ndarray[tuple[int], np.dtype[np.int_]],
        params: ParameterClass,
    ) -> None:
        """
        model <- expects the name of the model as str
        x0 <- expects initial condition as int array
        params <- expects jit class named ParameterClass
        """

        self.sub_step_module = self.step_module + ".step_function_" + model
        self.sub_propensity_module = self.propensity_module + ".propensity_" + model
        self.state_changes = self.vj["vj_" + model]

        self.x0 = x0
        self.params = params

    def step_function(
        self,
        steps: int | np.int_ = steps,
    ) -> tuple[
        np.ndarray[tuple[int], np.dtype[np.float64]],
        np.ndarray[tuple[int, int], np.dtype[np.int_]],
    ]:
        modu = import_module(self.sub_propensity_module)
        propensity = getattr(modu, "main")

        return step_function(
            steps,
            self.x0,
            self.state_changes,
            self.params,
            propensity,
        )
