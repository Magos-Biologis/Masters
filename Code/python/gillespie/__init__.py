##
## Stealing Jesper's code and making it into Python
##
import numpy as np

import numba as nb
from numba import typed as nbt
from numba.experimental import jitclass


spec = [
    ("x0", nb.float64[:]),
    ("propensity_function", nb.types.FunctionType(nb.types.float64[:](nb.types.int_[:]))),
    ("vj", nb.int_[:]),
    ("steps", nb.int_),
]


@jitclass(spec)
class gillespie:
    def __init__(
        self,
        x0,  #: np.ndarray[tuple[int], np.dtype[np.float64]],
        propensity_function,  #: nb.FunctionType,  # np.ndarray[tuple[int], np.dtype[np.float64]],
        vj,
        steps=1000,
    ) -> None:
        self.x0 = x0

        # self.props: nb.FunctionType = self._propensity(propensity_function)
        # self.v = self._transitions(v)

        self.props = propensity_function
        self.v = vj

        self.steps = steps

        self.time = np.empty(self.steps, dtype=np.float64)
        self.results = np.empty((2, self.steps), dtype=np.int_)

    def generate(
        self,
    ) -> tuple[
        np.ndarray[tuple[int], np.dtype[np.float64]],
        np.ndarray[tuple[int], np.dtype[np.int_]],
    ]:
        self.results[:, 0] = self.x0
        x = self.x0.astype(np.float64)

        broke_loop: bool = False
        for i in range(1, self.steps + 1):
            # a_j: np.ndarray[tuple[int], np.dtype[np.float64]] = self.props(
            #     self.results[:, i - 1]
            # )
            a_j = self.props(self.results[:, i - 1])
            j, dt = self._ssa(a_j)

            if j == -1:
                broke_loop: bool = True
                final_step: np.int_ = i
                break

            self.results[:, i] = np.add(self.results[:, i - 1], self.v[j])
            self.time[i] = self.time[i - 1] + dt

        if broke_loop:
            return self.time[:final_step], self.results[:, :final_step]

        return self.time, self.results

    # def _propensity(self, func: function) -> nb.FunctionType:
    #     return nb.FunctionType(func)
    #
    # def _transitions(self, v: list) -> nbt.NumbaList:
    #     return nbt.NumbaList(v)

    def _ssa(self, aj: np.ndarray[tuple[int], np.dtype[np.float64]]) -> tuple[int, float]:
        """
        Gillespie event determination

        Input
        propensities: The probability of an event, eg. A+B->A+A  propensity is k*A*B

        Output
        eventId: The index of event as defined by the propensities
        dt: Time step

        Example
        >> # For reaction scheme  A+B->P and A->P using mass action
        >> k1 = 0.01; k2=1.0; B = 30; A = 80; t = 0.0;
        >> propensities = [k1*A*B, k2*A];  [eventid dt] = gillespie(prop);
        """

        r: np.ndarray[tuple[int], np.dtype[np.float64]] = np.random.rand(2)
        while r[0] == 0.0:
            r[0] = np.random.rand()

        j: int = -1
        a_0: float = aj.sum()

        if a_0 <= 0:
            tau: float = 0
            return j, tau

        tau: float = np.log(1 / r[0]) / a_0
        ra_0: float = r[1] * a_0

        for n, _ in enumerate(aj):
            j += 1  ## Cause its funny, why not

            s_j: float = aj[0 : n + 1].sum()

            if s_j >= ra_0:
                break

        return j, tau
