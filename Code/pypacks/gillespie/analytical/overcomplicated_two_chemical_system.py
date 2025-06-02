import numpy as np
from contourpy import contour_generator

from .analytical_parameter_class import AnalParams


class MultivariateTwoChemical(AnalParams):
    def __init__(self, x_bounds, y_bounds, *args, **kwargs):
        self.x_domain = self._domain_creation(x_bounds, **kwargs)
        self.y_domain = self._domain_creation(y_bounds, **kwargs)
        density = kwargs.pop("density", 100)

        self.xx, self.yy = np.meshgrid(self.x_domain, self.y_domain, indexing="xy")

        length = len(self.x_domain)
        abstol = kwargs.pop("atol", 1 / length)
        should_mask = kwargs.pop("mask", False)

        if should_mask:
            diagonal_mask = ~np.isclose(self.xx + self.yy - 1, 0, atol=abstol)
            self.x = np.ma.array(self.xx, mask=diagonal_mask)
            self.y = np.ma.array(self.yy, mask=diagonal_mask)
        else:
            self.x, self.y = np.meshgrid(self.x_domain, self.y_domain, indexing="xy")

        super().__init__(*args, **kwargs)

    def _A(self):
        x = self.k2 * self.y - self.k1 * self.x
        y = -(self.k2 * self.y - self.k1 * self.x)
        return np.array([x, y])

    def _D(self):
        output = self.k2 * self.y + self.k1 * self.x
        return output

    def _sigma(self):
        z = self.x - self.y
        x = (self.n + z) / 2
        y = (self.n - z) / 2
        output = self.k2 * y + self.k1 * x
        return output

    def _integral_of_x(self):
        first_term = -self.x
        second_term = np.log(self._D()) * self.k2 * self.y / self.k1
        output = first_term + second_term
        return output

    def _integral_of_y(self):
        first_term = -self.y
        second_term = np.log(self._D()) * self.k1 * self.x / self.k2
        output = first_term + second_term
        return output

    def _potential_function(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        first_term = (self.k1) ** 2 * self.x + (self.k2) ** 2 * self.y
        first_term *= np.log(self._D())
        second_term = self.k1 * self.k2 * (self.x + self.y) / 2
        scale = (self.n * 4) / (self.k1 * self.k2)
        output = scale * (first_term - second_term)
        return output

    def _potential_gradient(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        k1, k2, n = self.k1, self.k2, self.n
        x, y = self.x, self.y
        sum_term = 2 * (k1 * x + k2 * y)
        divergence_term = k1 - k2
        first = 2 * n * self._A()[0] - divergence_term
        second = 2 * n * self._A()[1] + divergence_term
        output = np.array([first, second]) / sum_term
        return output

    # def _potential_function(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
    #     k1, k2, n = self.k1, self.k2, self.n
    #     x, y = self.x, self.y
    #     log_term = np.log(k1 * x + k2 * y)
    #     if k1 == k2:
    #         sum = np.log(x + y)
    #         fx = n * (2 * y * sum - x)
    #         fy = n * (2 * x * sum - y)
    #         output = fx + fy
    #     else:
    #         fx = (y * k2 / k1) * log_term
    #         fy = (x * k1 / k2) * log_term
    #         # standalone_terms =
    #         output = fx + fy  # + standalone_terms
    #     return output

    def get_equality(self):
        k1, k2 = self.k1, self.k2
        x = self.xx
        y = self.yy
        z = k1 * x + k2 * y - 1
        contour_data = contour_generator(self.x, self.y, z)
        nullcline = contour_data.lines(0)[0]

        assert type(nullcline) is np.ndarray
        return nullcline.T

    # x: np.ndarray[tuple[int], np.dtype[np.float64]],
    # y: np.ndarray[tuple[int], np.dtype[np.float64]],

    def stationary(self) -> tuple:
        potential = self._potential_function()

        output = np.exp(-potential / self._D())
        output /= self._D()

        return (self.x, self.y), output

    def gradient(self):
        output = self._potential_gradient()
        return (self.x, self.y), output
