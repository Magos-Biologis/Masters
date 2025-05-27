import itertools
import os
from typing import ParamSpecArgs

import numpy as np
from contourpy import contour_generator
from matplotlib import pyplot as plt
from numpy import intersect1d
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
        third_term = self.p.m0 * (self.p.q1 / self.p.q2) * (x / y)

        return first_term - second_term - third_term

    # third_term = 0

    def _find_level_set(
        self, mesh=None
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        if mesh is None:
            x, y = self.mesh_manifold
        else:
            x, y = mesh
        z = self._general_level_set_func(x, y)

        # plt.contour(x, y, z, [0])

        contour_data = contour_generator(x, y, z)
        nullcline = contour_data.lines(0)[0]

        assert type(nullcline) is np.ndarray
        return nullcline.T

    def _find_c1_nullcline(
        self, mesh=None
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        if mesh is None:
            x, y = self.mesh_manifold
        else:
            x, y = mesh
        first_term = (self.p.q1 * x * self.p.m0) / (
            (self.p.w2 - self.p._k1 * x) * (self.p.q2 * y)
        )
        second_term = -((self.p.o1 - self.p._k1 * x) * x) / (self.p.w2 - self.p._k1 * x)

        z = first_term + second_term

        contour_data = contour_generator(x, y, z)
        nullcline = contour_data.lines(0)[0]

        assert type(nullcline) is np.ndarray
        return nullcline.T

    def _find_c2_nullcline(
        self,
        c2: np.ndarray[tuple[int, int], np.dtype[np.float64]],
    ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        c1 = -((self.p.o2 - self.p._k2 * c2) * c2) / (self.p.w1 - self.p._k2 * c2)

        # contour_data = contour_generator(x, y, z)
        # nullcline = contour_data.lines(0)[0]

        solution = np.array([c1, c2])
        less_than = solution[:, c1 < self.p.c1_max]
        greater_than = solution[:, c1 > self.p.c1_min]
        less_idx = np.where(np.intersect1d(c2, less_than[1, :]))[0]
        greater_idx = np.where(np.intersect1d(c2, greater_than[1, :]))[0]

        intersects = np.intersect1d(less_idx, greater_idx)

        # assert type(nullcline) is np.ndarray
        return solution[:, intersects]

    def _get_intersection(self, **kwargs) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        self.tolerance = 2 ** (-6)

        cs = self._find_level_set()
        n2 = self._find_c2_nullcline(cs[1, :])

        cs_idx, n2_idx = np.intersect1d(cs[1, :], n2[1, :], return_indices=True)[1:]

        narrowed_cs = cs[:, cs_idx]
        narrowed_n2 = n2[:, n2_idx]

        diff = narrowed_cs[0, :] - narrowed_n2[0, :]
        intersect_idx = np.where(np.diff(np.sign(diff)))[0]

        return cs[:, intersect_idx]

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

        old_method = kwargs.pop("use_old", False)
        self.grid_resolution = kwargs.pop("grid_density", 225)

        x_domain = np.linspace(
            self.p.c1_min,
            self.p.c1_max,
            self.grid_resolution,
            endpoint=False,
        )[1:]
        y_domain = np.linspace(
            self.p.c2_min,
            self.p.c2_max,
            self.grid_resolution,
            endpoint=False,
        )[1:]

        self.mesh_manifold = np.meshgrid(x_domain, y_domain)

        if old_method:
            cs = self._find_level_set()
            n2 = self._find_c2_nullcline(cs[1, :])

            cs_idx, n2_idx = np.intersect1d(cs[1, :], n2[1, :], return_indices=True)[1:]

            narrowed_cs = cs[:, cs_idx]
            narrowed_n2 = n2[:, n2_idx]

            intersect = self._medless_equal_normal_roots()

            plt.plot(*narrowed_cs, color="g", label="level set")
            # plt.plot(n1_x, n1_y, color="b", label="c1 nullcline")
            plt.plot(*narrowed_n2, color="r", label="c2 nullcline")
            plt.vlines([intersect[0]], 0, 100, color="b", alpha=0.2)
            plt.hlines([intersect[1]], 0, 100, color="r", alpha=0.2)
            plt.xlim(0, 100)
            plt.ylim(0, 100)
            plt.show()
            exit()
            return self._medless_equal_normal_roots()

        if generalized:
            return self._medless_equal_normal_roots()
        else:
            if self.p.m0 == 0 and (self.p.k1 == self.p.k2 or self.p.n1 == self.p.k2):
                return np.append(self._medless_equal_normal_roots(), 0)
            c_roots = self._get_intersection()
            m_star = self.p.m0 * (self.p.q1 * c_roots[0]) / (c_roots[1] * self.p.q2)
            return np.append(c_roots, m_star)

    def nullclines(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        return np.empty(shape=(2, 2))

    # def plot_system(self):
    #     return self._integrate_general_model()
