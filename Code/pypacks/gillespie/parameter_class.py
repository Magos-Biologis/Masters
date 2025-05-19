from numba import float64
from numba.experimental import jitclass

spec = [
    ("n", float64),
    ("b", float64),
    ("k1", float64),
    ("k_1", float64),
    ("k2", float64),
    ("k_2", float64),
    ("k3", float64),
    ("k_3", float64),
    ("k4", float64),
    ("k_4", float64),
    ("k5", float64),
    ("k_5", float64),
    ("k6", float64),
    ("k_6", float64),
    ("k7", float64),
    ("k_7", float64),
    ("k8", float64),
    ("k_8", float64),
    ("m0", float64),
    ("w1", float64),
    ("w2", float64),
    ("q1", float64),
    ("q2", float64),
    ("n1", float64),
    ("n2", float64),
    ("k1p", float64),
    ("k2p", float64),
    ("k3p", float64),
    ("k4p", float64),
]


@jitclass(spec=spec)
class ParameterClass:
    def __init__(
        self,
        n,
        b,
        k1,
        k_1,
        k2,
        k_2,
        k3,
        k_3,
        k4,
        k_4,
        k5,
        k_5,
        k6,
        k_6,
        k7,
        k_7,
        k8,
        k_8,
        m0,
        w1,
        w2,
        q1,
        q2,
        n1,
        n2,
    ):
        # assign raw parameters
        self.n = n
        self.b = b
        self.k1 = k1
        self.k_1 = k_1
        self.k2 = k2
        self.k_2 = k_2
        self.k3 = k3
        self.k_3 = k_3
        self.k4 = k4
        self.k_4 = k_4
        self.k5 = k5
        self.k_5 = k_5
        self.k6 = k6
        self.k_6 = k_6
        self.k7 = k7
        self.k_7 = k_7
        self.k8 = k8
        self.k_8 = k_8

        self.m0 = m0
        self.w1 = w1
        self.w2 = w2
        self.q1 = q1
        self.q2 = q2
        self.n1 = n1
        self.n2 = n2

        self.k1p = self.k1 * self.n
        self.k2p = self.k2 * self.n
        self.k3p = self.k3 * self.b
        self.k4p = self.k4 * self.b

    def divide(self, m: float64):
        self.k1 = self.k1 / m
        self.k_1 = self.k_1 / m
        self.k2 = self.k2 / m
        self.k_2 = self.k_2 / m
        self.k3 = self.k3 / m
        self.k_3 = self.k_3 / m
        self.k4 = self.k4 / m
        self.k_4 = self.k_4 / m
        self.k5 = self.k5 / m
        self.k_5 = self.k_5 / m
        self.k6 = self.k6 / m
        self.k_6 = self.k_6 / m
        self.k7 = self.k7 / m
        self.k_7 = self.k_7 / m
        self.k8 = self.k8 / m
        self.k_8 = self.k_8 / m

        self.m0 = self.m0 / m
        self.w1 = self.w1 / m
        self.w2 = self.w2 / m
        self.q1 = self.q1 / m
        self.q2 = self.q2 / m
        # self.n1 = self.n1 / m
        # self.n2 = self.n2 / m

        self.k1p = self.k1p / m
        self.k2p = self.k2p / m
        self.k3p = self.k3p / m
        self.k4p = self.k4p / m

        return self


parameter_default = dict(
    [
        ("n", 100),
        ("b", 0),
        ("k1", 1),
        ("k_1", 1),
        ("k2", 1),
        ("k_2", 0),
        ("k3", 1),
        ("k_3", 0),
        ("k4", 1),
        ("k_4", 0),
        ("k5", 1),
        ("k_5", 0),
        ("k6", 1),
        ("k_6", 0),
        ("k7", 1),
        ("k_7", 0),
        ("k8", 1),
        ("k_8", 0),
        ("m0", 0),
        ("w1", 0.15),
        ("w2", 0.15),
        ("q1", 0.85),
        ("q2", 0.85),
        ("n1", 100),
        ("n2", 100),
    ]
)
