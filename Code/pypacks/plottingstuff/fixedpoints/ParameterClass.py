from dataclasses import dataclass

from numpy import float64


@dataclass
class ParameterClass:
    n: float64 = float64(1)
    b: float64 = float64(1)

    k1: float64 = float64(1)
    k_1: float64 = float64(1)
    k2: float64 = float64(1)
    k_2: float64 = float64(1)
    k3: float64 = float64(1)
    k_3: float64 = float64(1)
    k4: float64 = float64(1)
    k_4: float64 = float64(1)
    k5: float64 = float64(1)
    k_5: float64 = float64(1)

    m0: float64 = float64(0.0)
    w1: float64 = float64(0.15)
    w2: float64 = float64(0.15)
    q1: float64 = float64(1)
    q2: float64 = float64(1)
    n1: float64 = float64(100)
    n2: float64 = float64(90)

    def __post_init__(self) -> None:
        self.k1p: float64 = self.k1 * self.n
        self.k2p: float64 = self.k2 * self.n

        self.k3p: float64 = self.k3 * self.b
        self.k4p: float64 = self.k4 * self.b
