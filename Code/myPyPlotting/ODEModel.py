import os
import itertools

from parameter_class import parameter_class

import numpy as np
from numpy import polynomial as npp


class ODEModel:
    def __init__(
        self,
        t_range: tuple[int, int],
        parameters: parameter_class | tuple[int, float, float, float, float],
        initial_condition: list[float]
        | np.ndarray[tuple[int], np.dtype[np.float64]]
        | None = None,
        **kwargs,
    ):
        assert type(parameters) is parameter_class, "Incorrect format of parameters"

        self.p: parameter_class = parameters

        self.m: int = self.p.m
        self.ran: range = (
            range(len(initial_condition))
            if type(initial_condition) is not None
            else range(self.p.m)
        )

        self.n_max = kwargs.get("n_max", np.max(self.p.n))

        self.init_cond = (
            initial_condition
            if type(initial_condition) is not None
            else [self.n_max] + self.m * [0]
        )

        self.t_span: tuple[int, int] = t_range
        self.dt: np.float64 = kwargs.get("dt", 0.01)
        self.t_array: np.ndarray = np.arange(*self.t_span, self.dt)

    def __step_function(self, xs) -> np.ndarray:
        assert type(self.p) is parameter_class, "Incorrect format of parameters"
        step = np.zeros_like(xs)

        ct = sum(xs)

        return step

    def __integrate_base_model(self):
        assert type(self.ran) is range
        assert len(self.init_cond) == 3

        # steps = np.zeros((self.m, len(self.t_array)))
        steps = np.zeros((3, len(self.t_array)))
        steps[:, 0] = self.init_cond

        for j, _ in enumerate(self.t_array):
            steps[:, j] += steps[:, j - 1]
            steps[:, j] += self.dt * self.__normal_model(steps[:, j - 1])

        return self.t_array, steps

    def __integrate_general_model(self):
        assert type(self.ran) is range

        # steps = np.zeros((self.m, len(self.t_array)))
        steps = np.zeros((self.m, len(self.t_array)))
        steps[:, 0] = self.init_cond

        for j, _ in enumerate(self.t_array):
            steps[:, j] += steps[:, j - 1]
            steps[:, j] += self.dt * self.__general_model(steps[:, j - 1])

        return self.t_array, steps

    def __normal_model(
        self, xs: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        # assert (
        #     all(type(self.p.k), type(self.p.n), type(self.p.q), type(self.p.w))
        #     is np.ndarray[tuple[int], np.dtype[np.float64]]
        # )

        assert len(xs) == 3

        step = np.empty_like(xs)

        c1, c2, m = xs
        ct = sum(xs[0:2])

        step[0] = (
            (self.p.k[0] * (1 - ct / self.p.n[0])) * c1
            - self.p.w[0] * c1
            + self.p.w[1] * c2
            - self.p.q[0] * m * c1
        )
        step[1] = (
            (self.p.k[1] * (1 - ct / self.p.n[1])) * c2
            - self.p.w[1] * c2
            + self.p.w[0] * c1
        )
        step[2] = -self.p.q[1] * m * c2 + self.p.m_0

        return step

    def __general_model(self, xs: np.ndarray):
        step = np.zeros_like(xs)
        ct = sum(xs)

        for i, j in itertools.product(self.ran, repeat=2):
            if i != j:
                step[i] += self.p.W[j, i] * xs[j] - self.p.W[i, j] * xs[i]
            else:
                step[i] += self.p.k[i] * xs[i] - (self.p.k[i] / self.p.n[i]) * ct * xs[i]

        return step

    def system(self):
        pass

    def __general_no_med(self):
        pass

    def show_parameters(self):
        return self.p

    def integrate(self, generalized: bool = False):
        assert type(generalized) is bool

        if generalized:
            return self.__integrate_general_model()
        else:
            return self.__integrate_base_model()

    # def plot_system(self):
    #     return self.__integrate_general_model()
