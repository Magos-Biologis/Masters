#!./.venv/bin/python
import enum
from pylab import *


k_1 = 1
k_2 = 1


n = 10


def A(x):
    return k_2 - (k_1 + k_2) * x


def B(x):
    return (1 / n) * (k_2 + (k_1 - k_2) * x)


# def exponent(y):
#     num = (
#         2 * k_1 * k_2 * log(abs((k_2 - k_1) * y - k_2))
#         + (k_2**2 - k_1**2) * y
#         - 2 * k_1 * k_2 * log(abs(k_2))
#         + 2 * k_1 * k_2 * log(k_2 - k_1)
#         - 2 * k_1 * log(k_1 - k_2) * k_2
#     )
#     den = (k_2 - k_1) ** 2
#
#     return num / den


# t_array = linspace(0, 1, 1000)
dx = 0.01
x_array = arange(0, 1, dx)


def step_function(x) -> ndarray:
    step = A(x) / B(x)

    return step


def integra(x):
    steps = zeros_like(x)
    steps[0] = 0
    for j, _ in enumerate(x_array):
        steps[j] += steps[j - 1]
        steps[j] += dx * step_function(steps[j - 1])  # print(steps[j, :])

    return steps


internal_results = integra(x_array)


def function(x):
    output = empty_like(x)
    for i, x in enumerate(x):
        output[i] = B(x) ** (-1) * exp(2 * internal_results[i])

    return output


results = function(x_array)

plot(x_array, results)
show()
