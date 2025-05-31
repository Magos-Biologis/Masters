import numpy as np

from .analytical_parameter_class import AnalParams


class SimpleTwoChemical(AnalParams):
    # print(self.n, self.k1, self.k2)

    def _b(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        result = (self.k2 + (self.k1 - self.k2) * self.x) * (self.n)
        return result

    def _equiv_c(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        return (1 - self.x) * self.x

    def _noneq_c(self) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        x = self.x
        in_log = (self.k2 + (self.k1 - self.k2) * x) / self.k2
        to_scale = (self.k1 - self.k2) * x - self.k2 * np.log(in_log)

        scalar = 2 * self.k1
        scalar /= (self.k1 - self.k2) ** 2

        output = x - scalar * to_scale
        return output

    def stationary(
        self,
        x: np.ndarray[tuple[int], np.dtype[np.float64]],
    ) -> np.ndarray:
        self.x = x

        if self.k1 == self.k2:
            c_x = self._equiv_c()
        else:
            c_x = self._noneq_c()

        top = 2 * self.n * c_x
        data_max = top.max()

        output = np.exp(top - data_max) / (self._b())

        return output

    # print(output)
