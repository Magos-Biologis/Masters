import os
from dataclasses import dataclass, field

import numpy as np


@dataclass
class parameter_class:
    m: int = 2
    k: float | np.ndarray[tuple[int], np.dtype[np.float64]] = 1
    n: float | np.ndarray[tuple[int], np.dtype[np.float64]] = 100
    q: float | np.ndarray[tuple[int], np.dtype[np.float64]] = 0.085
    w: float | np.ndarray[tuple[int], np.dtype[np.float64]] = 0.015
    generalize: bool = False
    m0: float = 0
    k1: float = 1
    k2: float = 1
    n1: float = 100
    n2: float = 100
    q1: float = 0.085
    q2: float = 0.085
    w1: float = 0.015
    w2: float = 0.015

    W: None | np.ndarray[tuple[int, int], np.dtype[np.float64]] = None
    K: None | np.ndarray[tuple[int, int], np.dtype[np.float64]] = None

    ## The dataclass equivilant of initialization functions
    def __post_init__(self):
        if self.generalize:
            self.k = self.__size_verification(self.k)
            self.n = self.__size_verification(self.n)
            self.q = self.__size_verification(self.q)
            self.w = self.__size_verification(self.w)
        else:
            self.k = np.array([self.k1, self.k2], dtype=float)
            self.n = np.array([self.n1, self.n2], dtype=float)
            self.q = np.array([self.q1, self.q2], dtype=float)
            self.w = np.array([self.w1, self.w2], dtype=float)

        self.W, self.K = self.__set_matrices()

    def __size_verification(self, parameter) -> np.ndarray:
        output = (
            np.ones(shape=self.m, dtype=np.float64)[:] * parameter
            if type(parameter) is not np.ndarray[tuple[self.m], np.dtype[np.float64]]
            else parameter
        )
        return output

    def __set_matrices(self) -> tuple:  # np.ndarray:
        assert type(self.k) is np.ndarray
        assert type(self.n) is np.ndarray
        ran = range(self.m)
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

        matrix_W = self.W if type(self.W) is None else w_reshaped

        matrix_K = (
            self.K
            if type(self.W) is None
            else np.array([[self.k[i] / self.n[j] for j in ran] for i in ran])
        )

        return matrix_W, matrix_K
