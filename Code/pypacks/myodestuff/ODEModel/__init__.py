import itertools
import os
from typing import ParamSpecArgs

import numpy as np
from contourpy import contour_generator
from numpy import polynomial as npp

# from matplotlib import pyplot as plt
# from matplotlib import
from myodestuff.parameter_class import ODEParameters


class ODEModel:
    def __init__(
        self,
        parameters: ODEParameters
        | tuple[int, float, float, float, float, float]
        | dict[str, float],
        t_range: tuple[int, int] | None = None,
        initial_condition: np.ndarray[tuple[int], np.dtype[np.float64]] | None = None,
        **kwargs,
    ):
        # if type(parameters) is dict():
        try:
            parameters = ODEParameters(**parameters)
        except TypeError:
            assert type(parameters) is ODEParameters, "Incorrect format of parameters"

        if initial_condition is None:
            initial_condition = np.array([0, 0, 0])

        assert type(initial_condition) is np.ndarray, (
            "Initial conditions should be ndarray"
        )

        self.p: ODEParameters = parameters
        self.t_span: tuple[int, int] = t_range if t_range is not None else (0, 1)

        self.m: int = self.p.m
        self.ran: range = (
            range(len(initial_condition))
            if initial_condition is not None
            else range(self.p.m)
        )

        self.n_max = kwargs.pop("n_max", np.max(self.p.n))

        self.init_cond = (
            initial_condition
            if initial_condition is not None
            else [self.n_max] + self.m * [0]
        )

        self.dt: np.float64 = kwargs.pop("dt", np.float64(0.01))
        self.t_array: np.ndarray[tuple[int], np.dtype[np.float64]] = np.arange(
            *self.t_span, self.dt
        )

        self.var_count: int = len(self.init_cond)
        self.step_count: int = len(self.t_array)

    def _step(
        self, xs: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        return xs + self.dt * self._model(xs)

    def _eular_integrate(self, function) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        self._model = function

        assert type(self.ran) is range
        steps = np.empty((self.var_count, self.step_count), dtype=np.float64)
        steps[:, 0] = self.init_cond

        for i in range(1, self.step_count):
            steps[:, i] = self._step(steps[:, i - 1])

        return steps

    def _general_level_set_func(
        self, x, y
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        first_term = (x + y) * (self.p._k1 * x + self.p._k2 * y)
        second_term = self.p.k1 * x + self.p.k2 * y
        third_term = self.p.m0 * (self.p.q1 / self.p.q2) * (x / y)

        return first_term - second_term - third_term

    # third_term = 0

    def _find_level_set(
        self,
        grid_concentration: int,
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        x_domain = np.linspace(
            self.p.c1_min,
            self.p.c1_max,
            grid_concentration,
            endpoint=False,
        )
        y_domain = np.linspace(
            self.p.c2_min,
            self.p.c2_max,
            grid_concentration,
            endpoint=False,
        )

        x, y = np.meshgrid(x_domain, y_domain)
        z = self._general_level_set_func(x, y)
        contour_data = contour_generator(x, y, z)

        nullcline = contour_data.lines(0)[0]

        assert type(nullcline) is np.ndarray

        return nullcline.T

    def _find_c2_nullcline(
        self,
        c2: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        c1 = -((self.p.o2 - self.p._k2 * c2) * c2) / (self.p.w1 - self.p._k2 * c2)

        array = np.array([c1, c2])

        solution = array[:, c1 < self.p.c1_max]

        return solution

    def _get_intersection(self, **kwargs):
        dens = self.grid_resolution

        cs = self._find_level_set(dens)
        n1 = self._find_c2_nullcline(cs[1, :])

        cs_int, n1_int = np.intersect1d(cs[1, :], n1[1, :], return_indices=True)[1:]
        min_index = np.argmin(np.abs(cs[0, cs_int] - n1[0, n1_int]))

        roots = cs[:, min_index]

        return roots

    def _simple_level_set_roots(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
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

    def _medless_equal_normal_roots(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
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

    # def _medless_normal_roots(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
    #     # assert type(self.p) is parameter_class
    #     assert type(self.p.k) is np.ndarray
    #     assert type(self.p.n) is np.ndarray
    #     assert type(self.p.w) is np.ndarray
    #     assert type(self.p.q) is np.ndarray
    #
    #     omega_1 = self.p.k[0] - self.p.w[0]
    #     omega_2 = self.p.k[1] - self.p.w[1]
    #
    #     k_1 = self.p.k[0] / self.p.n[0]
    #     k_2 = self.p.k[1] / self.p.n[1]
    #
    #     c1_min = self.p.w[1] / k_1
    #     c2_min = self.p.w[0] / k_2
    #
    #     c1_max = omega_1 / k_1
    #     c2_max = omega_2 / k_2
    #
    #     c2_root = (omega_2 + self.p.w[0]) / 2 * self.p.k[1]
    #     c1_root = self.p.n[0] * (1 - c2_root / self.p.n[1])
    #
    #     return np.array([c1_root, c2_root])

    def _simple_three_part_compartment_model(
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

        dxs = np.empty_like(xs)

        c1, c2, m = xs
        ct = sum(xs[0:2])

        dxs[0] = (
            (self.p.k[0] * (1 - ct / self.p.n[0])) * c1
            - self.p.w[0] * c1
            + self.p.w[1] * c2
            - self.p.q[0] * m * c1
        )
        dxs[1] = (
            (self.p.k[1] * (1 - ct / self.p.n[1])) * c2
            - self.p.w[1] * c2
            + self.p.w[0] * c1
        )
        dxs[2] = -self.p.q[1] * m * c2 + self.p.m0

        return dxs

    def _generalized_compartment_model(self, xs: np.ndarray):
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

    def show_parameters(self):
        return self.p

    def integrate(self, generalized: bool = False):
        assert type(generalized) is bool

        if generalized:
            return self.t_array, self._eular_integrate(
                function=self._generalized_compartment_model
            )
        else:
            return self.t_array, self._eular_integrate(
                function=self._simple_three_part_compartment_model
            )

    def roots(
        self, generalized: bool = False, **kwargs
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        assert type(generalized) is bool
        self.grid_resolution = kwargs.pop("grid_density", 225)

        if generalized:
            return self._medless_equal_normal_roots()
        else:
            return self._get_intersection()
            # if self.p.m0 == 0:
            # else:
            # return self._get_intersection()
            # return self._simple_level_set_roots()

            # assert type(self.p.k) is np.ndarray
            # assert type(self.p.n) is np.ndarray
            # if self.p.k[0] == self.p.k[1] or self.p.n[0] == self.p.n[1]:
            #     return self._medless_equal_normal_roots()
            # else:
            #     return self._no_treat_level_set()

    def nullclines(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        return np.empty(shape=(2, 2))

    # def plot_system(self):
    #     return self._integrate_general_model()
