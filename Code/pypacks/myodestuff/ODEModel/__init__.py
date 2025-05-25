import itertools
import os
from typing import ParamSpecArgs

import numpy as np
from contourpy import contour_generator
from matplotlib import pyplot as plt
from numpy import polynomial as npp

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
        first_term = self.p.k1 * x + self.p.k2 * y
        second_term = (x + y) * (self.p._k1 * x + self.p._k2 * y)
        third_term = self.p.m0 * (self.p.q1 / self.p.q2) * x * y**-1

        return first_term - second_term - third_term

    # third_term = 0

    def _find_level_set(
        self,
        grid_concentration: int,
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        x_domain = np.linspace(
            self.p.c1_min + self.tolerance,
            self.p.c1_max,
            grid_concentration,
            endpoint=False,
        )
        y_domain = np.linspace(
            self.p.c2_min + self.tolerance,
            self.p.c2_max,
            grid_concentration,
            endpoint=False,
        )

        x, y = np.meshgrid(x_domain, y_domain)
        z = self._general_level_set_func(x, y)

        # plt.contour(x, y, z, [0])

        contour_data = contour_generator(x, y, z)
        nullcline = contour_data.lines(0)[0]

        assert type(nullcline) is np.ndarray
        return nullcline.T

    def _find_c1_nullcline(
        self,
        c1: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        c2 = -((self.p.o1 - self.p._k1 * c1) * c1) / (self.p.w2 - self.p._k1 * c1)

        solution = np.array([c1, c2])

        return solution

    def _find_c2_nullcline(
        self,
        c2: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        c1 = -((self.p.o2 - self.p._k2 * c2) * c2) / (self.p.w1 - self.p._k2 * c2)

        solution = np.array([c1, c2])

        return solution

    def _get_intersection(self, **kwargs) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        dens = self.grid_resolution
        self.tolerance = 2 ** (-4)

        cs = self._find_level_set(dens)
        # n1 = self._find_c1_nullcline(cs[0, :])
        n2 = self._find_c2_nullcline(cs[1, :])

        cs_int, n2_int = np.intersect1d(cs[1, :], n2[1, :], return_indices=True)[1:]
        min_index = np.argmin(np.abs(cs[0, cs_int] - n2[0, n2_int]))

        intersect = cs[:, min_index]

        # # intersect_idx = np.argwhere(np.diff(np.sign(cs - n2))).flatten()
        # plt.plot(*cs, color="r", label="level set")
        # # plt.plot(*n1, color="b", label="c1 nullcline")
        # plt.plot(*n2, color="g", label="c2 nullcline")
        # # plt.plot(cs[:, intersect_idx], color="k", label="overhere?")
        # plt.hlines([intersect[1]], 0, 100, color="r", alpha=0.2)
        # plt.vlines([intersect[0]], 0, 100, color="b", alpha=0.2)
        # plt.xlim(0, 100)
        # plt.ylim(0, 100)
        # plt.show()
        # exit()

        return intersect

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

    def _simple_three_part_compartment_model(
        self, xs: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
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
        old_method = kwargs.pop("use_old", False)

        if old_method:
            return self._medless_equal_normal_roots()

        if generalized:
            return self._medless_equal_normal_roots()
        else:
            c_roots = self._get_intersection()
            m_star = self.p.m0 * (self.p.q1 / self.p.q2) * (c_roots[0] / c_roots[1])
            return np.append(c_roots, m_star)

    def nullclines(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        return np.empty(shape=(2, 2))

    # def plot_system(self):
    #     return self._integrate_general_model()
