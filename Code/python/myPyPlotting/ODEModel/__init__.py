import os
import itertools

from myPyPlotting.parameter_class import parameter_class

import numpy as np
from numpy import polynomial as npp


class ODEModel:
    def __init__(
        self,
        t_range: tuple[int, int],
        parameters: parameter_class | tuple[int, float, float, float, float],
        initial_condition: np.ndarray[tuple[int], np.dtype[np.float64]] | None = None,
        **kwargs,
    ):
        assert type(parameters) is parameter_class, "Incorrect format of parameters"
        assert type(initial_condition) is np.ndarray, (
            "Initial conditions should be ndarray"
        )

        self.p: parameter_class = parameters

        self.m: int = self.p.m
        self.ran: range = (
            range(len(initial_condition))
            if initial_condition is not None
            else range(self.p.m)
        )

        self.n_max = kwargs.get("n_max", np.max(self.p.n))

        self.init_cond = (
            initial_condition
            if initial_condition is not None
            else [self.n_max] + self.m * [0]
        )

        self.t_span: tuple[int, int] = t_range
        self.dt: np.float64 = kwargs.get("dt", 0.01)
        self.t_array: np.ndarray = np.arange(*self.t_span, self.dt)

    def _step_function(self, xs) -> np.ndarray:
        assert type(self.p) is parameter_class, "Incorrect format of parameters"
        step = np.zeros_like(xs)

        ct = sum(xs)

        return step

    def _integrate_base_model(self):
        assert type(self.ran) is range

        assert type(self.init_cond) is np.ndarray
        assert len(self.init_cond) == 3

        # steps = np.zeros((self.m, len(self.t_array)))
        steps = np.zeros((3, len(self.t_array)))
        steps[:, 0] = self.init_cond

        for j, _ in enumerate(self.t_array):
            steps[:, j] += steps[:, j - 1]
            steps[:, j] += self.dt * self._normal_model(steps[:, j - 1])

        return self.t_array, steps

    def _integrate_general_model(self):
        assert type(self.ran) is range

        # steps = np.zeros((self.m, len(self.t_array)))
        steps = np.zeros((self.m, len(self.t_array)))
        steps[:, 0] = self.init_cond

        for j, _ in enumerate(self.t_array):
            steps[:, j] += steps[:, j - 1]
            steps[:, j] += self.dt * self._general_model(steps[:, j - 1])

        return self.t_array, steps

    def _normal_roots(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        # assert type(self.p) is parameter_class
        assert type(self.p.k) is np.ndarray
        assert type(self.p.n) is np.ndarray
        assert type(self.p.w) is np.ndarray
        assert type(self.p.q) is np.ndarray

        omega_1 = self.p.k[0] - self.p.w[0]
        omega_2 = self.p.k[1] - self.p.w[1]

        k_1 = self.p.k[0] / self.p.n[0]
        k_2 = self.p.k[1] / self.p.n[1]

        c1_min = self.p.w[1] / k_1
        c2_min = self.p.w[0] / k_2

        c1_max = omega_1 / k_1
        c2_max = omega_2 / k_2

        c1_coeffs = [
            self.p.n[1] * self.p.w[1],
            omega_1 - self.p.n[1] * (self.p.w[1] / self.p.n[0] + k_1),
            k_1 * ((self.p.n[1] / self.p.n[0]) - 1),
        ]
        c2_coeffs = [
            self.p.n[0] * self.p.w[0],
            omega_2 - self.p.n[0] * (self.p.w[0] / self.p.n[1] + k_2),
            k_2 * ((self.p.n[0] / self.p.n[1]) - 1),
        ]

        c1_poly = npp.Polynomial(
            coef=c1_coeffs,
            symbol="_c1",
        )
        c2_poly = npp.Polynomial(
            coef=c2_coeffs,
            symbol="_c2",
        )

        c1_roots = c1_poly.roots()
        c2_roots = c2_poly.roots()

        c1_root = [root for root in c1_roots if c1_min < root < c1_max]
        c2_root = [root for root in c2_roots if c2_min < root < c2_max]

        return np.array([c1_root, c2_root])

    def _normal_model(
        self, xs: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        # assert (
        #     all(type(self.p.k), type(self.p.n), type(self.p.q), type(self.p.w))
        #     is np.ndarray[tuple[int], np.dtype[np.float64]]
        # )

        assert type(self.p.k) is np.ndarray
        assert type(self.p.n) is np.ndarray
        assert type(self.p.w) is np.ndarray
        assert type(self.p.q) is np.ndarray

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

    def _general_model(self, xs: np.ndarray):
        assert type(self.p.k) is np.ndarray
        assert type(self.p.n) is np.ndarray
        assert type(self.p.q) is np.ndarray
        assert type(self.p.W) is np.ndarray

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

    def _general_no_med(self):
        pass

    def show_parameters(self):
        return self.p

    def integrate(self, generalized: bool = False):
        assert type(generalized) is bool

        if generalized:
            return self._integrate_general_model()
        else:
            return self._integrate_base_model()

    def roots(self, generalized: bool = False):
        assert type(generalized) is bool

        if generalized:
            return self._normal_roots()
        else:
            return self._normal_roots()

    # def plot_system(self):
    #     return self._integrate_general_model()
