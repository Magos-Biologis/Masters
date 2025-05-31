from dataclasses import dataclass

import numpy as np


@dataclass
class AnalParams:
    n: int | np.int_
    k1: float | np.float64
    k2: float | np.float64
