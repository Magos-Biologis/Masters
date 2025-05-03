import os

import numpy as np
from numba import njit


class simple_two_system:
    def __init__(self, k1: float | np.float64, k2: float | np.float64, n: int) -> None:
        self.n = np.int_(n)

        self.k1 = np.float64(k1)
        self.k2 = np.float64(k2)

        self.scale = np.log(2 * n)

    def _normalize_scale(self):
        pass

    def _b(
        self,
        x: np.ndarray[tuple[int], np.dtype[np.float64]],
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        result = self.k2 + (self.k1 - self.k2) * x
        return np.divide(result, self.n)

    def _equiv_c(
        self,
        x: np.ndarray[tuple[int], np.dtype[np.float64]],
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        return x**2

    def _noneq_c(
        self,
        x: np.ndarray[tuple[int], np.dtype[np.float64]],
    ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        num1: np.float64 = 2 * self.k1
        num2: np.ndarray[tuple[int], np.dtype[np.float64]] = (self.k1 - self.k2) * x
        num2 -= self.k2 * np.log(1 + ((self.k1 - self.k2) * x) / self.k2)

        den: np.float64 = (self.k1 - self.k2) ** 2

        return np.divide(np.multiply(num2, num1), den)

    def stationary(
        self,
        x: np.ndarray[tuple[int], np.dtype[np.float64]],
        n: int | np.int_ | None = None,
        scale: float | np.float64 | None = None,
    ) -> np.ndarray:
        _n: int | np.int_ = n if n is not None else self.n
        _scale: float | np.float64 = scale if scale is not None else self.scale

        c_x = np.empty_like(x, dtype=np.float64)

        if self.k1 != self.k2:
            c_x[:] = np.subtract(x, self._noneq_c(x))
        else:
            c_x[:] = np.subtract(x, self._equiv_c(x))

        top: np.ndarray[tuple[int], np.dtype[np.float64]] = np.exp(
            _n * np.subtract(c_x, _scale)
        )

        return np.divide(top, self._b(x), dtype=np.float64)
