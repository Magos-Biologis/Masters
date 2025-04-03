#!./.venv/bin/python
from numba import njit
from pylab import *


def drift(x, y, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)
    k3 = kwargs.get("k3", 1)
    k4 = kwargs.get("k4", 1)
    k5 = kwargs.get("k5", 1)

    x_part = k1 * x + k3 * x * y - k4 * x
    y_part = k2 * x - k3 * x * y - k5 * y

    return [x_part, y_part]


def diffusion(x, y, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)
    k3 = kwargs.get("k3", 1)
    k4 = kwargs.get("k4", 1)
    k5 = kwargs.get("k5", 1)
    omega = kwargs.get("omega", 1)

    xx = k1 * x + k3 * x * y + k4 * x
    xy = -k3 * x * y
    yy = k2 * x + k3 * x * y + k5 * y

    return divide(array([[xx, xy], [xy, yy]]), omega)


k_1 = 1
k_2 = 1
k_3 = 1
k_4 = 1
k_5 = 1

omega = 1

print()
print()
