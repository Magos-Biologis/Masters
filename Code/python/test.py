#!./.venv/bin/python
from numba import njit
from pylab import *


b = 20


alpha = 1
beta = 2 * b


k_1 = 1
k_2 = 1


n = 10


def A(x):
    return k_2 - (k_1 + k_2) * x


def B(x):
    return (1 / n) * (k_2 + (k_1 - k_2) * x)


@njit
def stationary_p(x, b):
    top = exp(-2 * x) * (b + x) ** (4 * b - 1)
    return top / x


@njit
def p_s(x, b, k1, k2):
    top = exp(-2 * x) * (b + x) ** (4 * alpha - 1)
    return top / x


# def function(x):
#     output = empty_like(x)
#     for i, x in enumerate(x):
#         output[i] = B(x) ** (-1) * exp(2 * internal_results[i])
#     return output


x_array = np.linspace(alpha, beta, 1000)
# results = stationary_p(x_array, b)
results = p_s(x_array, b, k_1, k_2)

plot(x_array, results)
show()
