# from myPyPlotting import *

import os
import itertools


from sys import exit
from dataclasses import dataclass

import numpy as np


# import pprint


class parameter_types:
    m: int

    k: np.ndarray | float
    n: np.ndarray | float
    q: np.ndarray | float

    W: np.ndarray | float
    K: np.ndarray | float


class parameter_class(parameter_types):
    def __init__(
        self,
        m: int,
        k: float | list[float],
        q: float | list[float],
        w: float | list[float],
        n_max: int = 100,
    ):
        self.m = m
        ran = range(self.m)

        self.k = np.empty(self.m)
        self.n = np.empty(self.m)
        self.q = np.empty(self.m)

        self.W = np.empty((self.m, self.m))
        self.K = np.empty((self.m, self.m))

        self.k[:] = k
        self.q[:] = q
        self.n[:] = [n_max - i * 10 for i in ran]

        # n[:] =
        self.W[:] = w

        for i, j in itertools.product(ran, repeat=2):
            self.K[i, j] = self.k[i] / self.n[j]


test = parameter_class(1, 2, 3, 4, 5)
print(test)
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
