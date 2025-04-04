# from myPyPlotting import *

import os
import itertools
import collections
from dataclasses import dataclass, field

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

    ## The dataclass equivilant of initialization functions
    def __post_init__(self):
        self.k = self.__size_verification(self.k, self.m)
        self.n = self.__size_verification(self.n, self.m)
        self.q = self.__size_verification(self.q, self.m)
        self.w = self.__size_verification(self.w, self.m)

    def __size_verification(self, parameter, m: int) -> np.ndarray:
        output = (
            np.ones(shape=m, dtype=np.float64)[:] * parameter
            if type(parameter) is not np.ndarray[tuple[m], np.dtype[np.float64]]
            else parameter
        )
        return output

    def __set_matrices(self) -> tuple:  # np.ndarray:
        assert type(self.k) is np.ndarray
        assert type(self.n) is np.ndarray
        # [tuple[parameters.m], np.dtype[np.float64]]
        # [tuple[parameters.m], np.dtype[np.float64]]

        ones_no_diag = np.ones((self.m, self.m)) - np.diag(np.ones(self.m))

        w_temp = self.w
        for i in range(self.m - 1):
            w_temp = np.append(w_temp, self.w)
        w_reshaped = np.reshape(w_temp, (self.m, self.m))
        indices = np.diag_indices(self.m)
        for i, j in zip(*indices):
            w_reshaped[i, j] = self.k[i]

        # w_reshaped[i, j] = 0
        # print(w_reshaped)
        # # w_nodiag = np.delete(w_reshaped, indices)  # , 0)
        # # w_nodiag = np.insert(w_reshaped, indices)  # , 0)

        # matrix_W[indices[0], indices[1]]

        matrix_W = parameters.W if type(parameters.W) is None else w_reshaped
        matrix_K = (
            parameters.K
            if type(parameters.W) is None
            else np.array([[self.k[i] / self.n[j] for j in self.ran] for i in self.ran])
        )

        return matrix_W, matrix_K


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

        self.m = self.p.m
        self.ran = range(self.m)
        # ones_matrix = np.ones((self.p.m, self.p.m))

        # self.p.k = self.k = self.__size_verification(self.p.k, m)
        # self.p.n = self.n = self.__size_verification(self.p.n, m)
        # self.p.q = self.q = self.__size_verification(self.p.q, m)
        # self.p.w = self.w = self.__size_verification(self.p.w, m)

        # matrices = self.__set_matrices(self.p)
        # self.p.W = self.W = matrices[0]
        # self.p.K = self.K = matrices[1]

    def __step_function(self, xs, parameters) -> np.ndarray:
        step = np.zeros_like(xs)
        # ran = range(len(step))
        ct = sum(xs)

        for i, j in itertools.product(self.ran, repeat=2):
            if i != j:
                step[i] += self.W[j, i] * xs[j] - self.W[i, j] * xs[i]
            else:
                step[i] += self.k[i] * xs[i] - (self.k[i] / self.n[i]) * ct * xs[i]

        return step

    # def __test(self):
    #     return print("hello")

    def __integrate_model(self):
        pass

    def system(self):
        pass

    def __general_no_med(self):
        pass

    def show_parameters(self):
        return self.p


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
pp(test2.show_parameters())
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
