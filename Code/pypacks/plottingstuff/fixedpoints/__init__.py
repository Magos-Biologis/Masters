import os

import numpy as np
from myodestuff.ODEModel import ODEModel
from myodestuff.parameter_class import parameter_class

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
    from myodestuff import ODEModel, parameter_class

    def ode_5_2(self) -> dict[str, np.dtype[np.float64]]:
        temp_dict = {
            "m0": self.m0,
            "k1": self.k1,
            "k2": self.k2,
            "w1": self.w1,
            "w2": self.w2,
            "q1": self.q1,
            "q2": self.q2,
            "n1": self.n1,
            "n2": self.n2,
        }

        params = parameter_class(**temp_dict)
        model = ODEModel(params)

        output = dict()
        sols = model.roots()
        output["x"] = sols[0]
        output["y"] = sols[1]

        return output

    def ode_3_2(self) -> dict[str, np.dtype[np.float64]]:
        temp_dict = {
            "m0": self.m0,
            "k1": self.k1,
            "k2": self.k2,
            "w1": self.w1,
            "w2": self.w2,
            "q1": self.q1,
            "q2": self.q2,
            "n1": self.n1,
            "n2": self.n2,
        }

        params = parameter_class(**temp_dict)
        model = ODEModel(params)

        output = dict()
        sols = model.roots()
        output["x"] = sols[0]
        output["y"] = sols[1]

        return output
