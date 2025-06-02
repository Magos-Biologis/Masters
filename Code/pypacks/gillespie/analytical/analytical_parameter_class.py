from dataclasses import dataclass

import numpy as np


@dataclass
class AnalParams:
    n: int | np.int_
    k1: float | np.float64
    k2: float | np.float64

    def _domain_creation(self, bounds: tuple[int, int], **kwargs):
        grid_density = kwargs.pop("density", 100)
        include_ends = kwargs.pop("endpoints", False)

        output = np.linspace(
            *bounds,
            grid_density + 1,
            endpoint=include_ends,
        )[1 - include_ends :]

        return output
