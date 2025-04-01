# from myPyPlotting import *

import os
import itertools
import collections
from dataclasses import dataclass

from sys import exit
from pprint import pp

import numpy as np
from numpy._core.multiarray import dtype


@dataclass
class parameter_class:
    m: int

    k: float | np.ndarray[tuple[int], np.dtype[np.float64]]
    n: float | np.ndarray[tuple[int], np.dtype[np.float64]]
    q: float | np.ndarray[tuple[int], np.dtype[np.float64]]
    w: float | np.ndarray[tuple[int], np.dtype[np.float64]]

    W: None | np.ndarray[tuple[int, int], np.dtype[np.float64]] = None
    K: None | np.ndarray[tuple[int, int], np.dtype[np.float64]] = None


class ODEModel:
    def __init__(
        self,
        parameters: parameter_class | tuple[int, float, float, float, float],
        # | type[tuple[Any, ...]]
        # | tuple[int, float, float, float, float],
        init_conds: np.ndarray | list = [],
        # m: int,
        # k: float | np.ndarray[tuple[int], np.dtype[np.float64]],
        # n: float | np.ndarray[tuple[int], np.dtype[np.float64]],
        # q: float | np.ndarray[tuple[int], np.dtype[np.float64]],
        # w: float | np.ndarray[tuple[int, int], np.dtype[np.float64]],
        # K: float | np.ndarray[tuple[int, int], np.dtype[np.float64]],
        # n_max: int = 100,
    ):
        if type(parameters) is not tuple:
            self.p = parameters
        else:
            self.p = parameter_class(*parameters)

        assert type(self.p) is parameter_class, "Incorrect format of parameters"

        ran = range(self.p.m)
        # ones_matrix = np.ones((self.p.m, self.p.m))

        self.p.k = (
            np.ones(shape=self.p.m, dtype=np.float64)[:] * self.p.k
            if type(self.p.k) is not np.ndarray[tuple[self.p.m], np.dtype[np.float64]]
            else self.p.k
        )
        self.p.n = (
            np.ones(shape=self.p.m, dtype=np.float64)[:] * self.p.n
            if type(self.p.n) is not np.ndarray[tuple[self.p.m], np.dtype[np.float64]]
            else self.p.n
        )
        self.p.q = (
            np.ones(shape=self.p.m, dtype=np.float64)[:] * self.p.q
            if type(self.p.q) is not np.ndarray[tuple[self.p.m], np.dtype[np.float64]]
            else self.p.q
        )

        ones_no_diag = np.ones((self.p.m, self.p.m)) - np.diag(np.ones(self.p.m))
        self.p.w = np.multiply(ones_no_diag, self.p.w)

        self.p.W = self.p.W if type(self.p.W) is None else self.p.w + np.diag(self.p.k)
        self.p.K = (
            self.p.K
            if type(self.p.W) is None
            else np.array([[self.p.k[i] / self.p.n[j] for j in ran] for i in ran])
        )

        def __set_matrices(self, K, W) -> None:  # np.ndarray:
            pass

        def __step_function(self, xs) -> np.ndarray:
            step = np.zeros_like(xs)
            ran = range(len(step))
            ct = sum(xs)

            for i, j in itertools.product(ran, repeat=2):
                if i != j:
                    step[i] += self.p.W[j, i] * xs[j] - self.p.W[i, j] * xs[i]
                else:
                    step[i] += k[i] * xs[i] - (k[i] / n[i]) * ct * xs[i]

            return step

        def __integrate_model(self):
            pass

        def system(self):
            pass

        def __general_no_med(self):
            pass


# class ODEModel:
#     """Represents an ODE model that uses ODEParameters."""
#
#     def __init__(self, params):
#         self.params = params
#
#     def rhs(self, t, y):
#         # Logistic growth: dy/dt = α*y - β*y².
#         # (It’s like life: growth tempered by too much competition.)
#         α = self.params.alpha
#         β = self.params.beta
#         return α * y - β * y**2


# test = parameter_types(1, 2, 3, 4, 5, 6)
# k = np.array([[6.1, 6.3], [5.6, 5.9]])
k = np.array([6.1, 6.3])
n = np.array([100, 90])
# test = (3, k, 3, 4, 5)
test = (3, 1, 3, 4, 5)
test = parameter_class(*test)
test2 = ODEModel(test)
# test = parameter_types(m=2)
# test2 = parameter_class(test)
print()
pp(test)
print()
# print(test2.p)
# print()
print(test2)

exit()


@dataclass
class base:
    def __init__(
        self,
        figure_name: str,
        subdirectory: str | None = None,
        figure_enviroment: str | None = None
        if type(os.getenv("THESIS_FIGURE_PATH")) is None
        else os.getenv("THESIS_FIGURE_PATH"),
    ) -> None:
        self.figure_name = figure_name
        self.subdirectory = subdirectory
        self.figure_env = figure_enviroment
