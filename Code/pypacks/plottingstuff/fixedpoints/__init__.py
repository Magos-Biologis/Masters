import os

import numpy as np

from .ParameterClass import ParameterClass


class GillespieFixed(ParameterClass):
    def Novel_5_2(self) -> dict[str, np.dtype[np.float64]]:
        output = dict()

        output["x"] = x_fixed = -(self.k3p - self.k1p) / self.k_1
        output["y"] = y_fixed = -(self.k2p / self.k4p) * x_fixed

        return output

    def Novel_5_3(self) -> dict[str, np.dtype[np.float64]]:
        output = dict()

        output["n"] = n_fixed = -self.k3p / (2 * self.k1)
        output["x"] = x_fixed = (self.k3p - self.k1 * n_fixed) / self.k_1
        output["y"] = y_fixed = (self.k2 * n_fixed * x_fixed) / (self.k4p + self.k5)

        return output


class ODEFixed(ParameterClass):
    def Novel_5_3(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        n_fixed = -self.k3p / (2 * self.k1)
        x_fixed = (self.k3p - self.k1 * n_fixed) / self.k_1
        y_fixed = (self.k2 / (self.k4p + self.k5)) * n_fixed * x_fixed

        return np.array([x_fixed, y_fixed, n_fixed])
