import numpy as np
from contourpy import contour_generator

from .analytical_parameter_class import AnalParams


class TwoParameterNovel(AnalParams):
    def __init__(self, x_bounds, y_bounds, b=0, k_1=1, k3=1, k4=1, k5=1, *args, **kwargs):
        self.x_domain = self._domain_creation(x_bounds, **kwargs)
        self.y_domain = self._domain_creation(y_bounds, **kwargs)

        length = len(self.x_domain)
        abstol = kwargs.pop("atol", 1 / length)
        should_mask = kwargs.pop("mask", False)

        self.xx, self.yy = np.meshgrid(self.x_domain, self.y_domain, indexing="xy")
        super().__init__(*args, **kwargs)

        if should_mask:
            grid_mask = ~np.isclose(
                self.k2 * self.n * self.xx,
                (k4 * b + k5) * self.yy,
                atol=abstol,
            )
            self.x = np.ma.array(self.xx, mask=grid_mask)
            self.y = np.ma.array(self.yy, mask=grid_mask)
        else:
            self.x, self.y = np.meshgrid(self.x_domain, self.y_domain, indexing="xy")

        self.b = b

        self.k_1 = k_1
        self.k3 = k3
        self.k4 = k4
        self.k5 = k5

        self.k1p = self.k1 * self.n
        self.k2p = self.k2 * self.n
        self.k3p = self.k3 * self.b
        self.k4p = self.k4 * self.b

    def _something(self):
        temp = contour_generator()

    def _A(self):
        k1, k_1, k2, k3, k4, k5 = self.k1, self.k_1, self.k2, self.k3, self.k4, self.k5
        b, n = self.b, self.n
        x, y = self.x, self.y
        xA = (k1 * n - k_1 * x) * x - k3 * b * x
        yA = k2 * n * x - (k4 * b + k5) * y
        return np.array([xA, yA])

    def _B(self):
        k1, k_1, k2, k3, k4, k5 = self.k1, self.k_1, self.k2, self.k3, self.k4, self.k5
        b, n = self.b, self.n
        x, y = self.x, self.y
        xxB = (k1 * n + k_1 * x) * x + k3 * b * x
        yyB = k2 * n * x + (k4 * b + k5) * y
        return [[xxB, 0], [0, yyB]]

    def _potential_gradient(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        k1, k_1, k2, k3, k4, k5 = self.k1, self.k_1, self.k2, self.k3, self.k4, self.k5
        b, n = self.b, self.n
        x, y = self.x, self.y
        first = 2 * self._A()[0] - 2 * k_1 * x - (k1 * n + k3 * b)
        first /= self._B()[0][0]
        second = 2 * self._A()[1] - (k4 * b + k5)
        second /= self._B()[1][1]
        output = np.array([first, second])
        return output

    def gradient(self):
        output = self._potential_gradient()
        return (self.x, self.y), output
