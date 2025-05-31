import numpy as np
from contourpy import contour_generator

from .analytical_parameter_class import AnalParams


class MultivariateTwoChemical(AnalParams):
    def __init__(self, x_bounds, y_bounds, *args, **kwargs):
        self.x_domain = self._domain_creation(x_bounds, **kwargs)
        self.y_domain = self._domain_creation(y_bounds, **kwargs)

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

    def _F(self):
        _Fx = 2 * self.n * (self._A()[0]) / self._D()
        _Fy = 2 * self.n * (self._A()[1]) / self._D()

        return np.array([_Fx, _Fy])

    def _potential_function(self) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        first_term = (self.k1 * self.x) ** 2 + (self.k2 * self.y) ** 2
        first_term *= np.log(self._D())

        second_term = self.k1 * self.k2 * (self.x + self.y) / 2

        scale = (self.n * 4) / (self.k1 * self.k2)

        output = scale * (first_term - second_term)
        return output

    def get_equality(self):
        # plt.contour(x, y, z, [0])
        z = self.k1 * self.x + self.k2 * self.y - 1
        contour_data = contour_generator(self.x, self.y, z)
        nullcline = contour_data.lines(0)[0]

        assert type(nullcline) is np.ndarray
        return nullcline.T

    # x: np.ndarray[tuple[int], np.dtype[np.float64]],
    # y: np.ndarray[tuple[int], np.dtype[np.float64]],
    def _domain_creation(self, bounds: tuple[int, int], **kwargs):
        grid_density = kwargs.pop("density", 100)
        include_ends = kwargs.pop("endpoints", False)

        output = np.linspace(
            *bounds,
            grid_density + 1,
            endpoint=include_ends,
        )[1 - include_ends :]

        return output

    def stationary(
        self,
    ) -> tuple:
        self.grid = np.meshgrid(self.x_domain, self.y_domain)
        self.x, self.y = self.grid

        potential = self._potential_function()
        output = np.exp(-potential) / (self._D())  #

        return (self.x, self.y[::-1]), output


# @np.vectorize
# def func(x,y):
#     # some arbitrary function
#     return a[x] + a[y]
#     # print(output)


#
# def _potential_function_old(
#     self,
# ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
#     # first_term = np.log(self.x * self.k1 + self.y * self.k2) * self.x * self.k1
#     # second_term = self.k2 * (self.x - self.y) / 2
#
#     _D = self._D()
#     scale = 2 * self.n / self.k1
#
#     first_term = 2 * np.log(_D) * self.k1 * self.y - _D
#     second_term = self.k2**2 / _D
#
#     output = scale * (first_term + second_term)
#
#     return output
#
#
# def _potential_function(
#     self,
# ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
#     _D = self._D()
#     term = 2 * np.log(_D) * self.k2 * self.y / self.k1
#     output = term - self.x
#     return output
#
#
# def _potential_function_phat(
#     self,
# ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
#     # first_term = np.log(self.x * self.k1 + self.y * self.k2) * self.x * self.k1
#     # second_term = self.k2 * (self.x - self.y) / 2
#     _D = self._D()
#     first_term = (self.k1**2 * self.x + self.k2**2 * self.y) * np.log(_D)
#     second_term = -self.k1 * (self.k1 * self.x)
#     third_term = -self.k1 * self.k2 * (self.x + self.y) / 2
#
#     scale = (self.k2 * self.k1) ** -1
#
#     output = scale * (first_term + second_term + third_term)
#     return output
#
#
# def _potential_function_phat_again(
#     self,
# ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
#     # first_term = np.log(self.x * self.k1 + self.y * self.k2) * self.x * self.k1
#     # second_term = self.k2 * (self.x - self.y) / 2
#     _D = self._D()
#
#     first_term = (2 * np.log(_D) - 3) * self.k1**2 * self.x**2
#     second_term = (-2 * np.log(_D) - 1) * self.k2**2 * self.y**2
#
#     scale = (self.k2 * self.k1) ** -1
#
#     output = scale * (first_term + second_term)
#     return output
