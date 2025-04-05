#!./.venv/bin/python
from numba import njit
from pylab import *


b = 100


alpha = 0
beta = 1


k_1 = 1
k_2 = 3


n = 1

# @njit
# def stationary_p(x, b):
#     top = exp(-2 * x) * (b + x) ** (4 * b - 1)
#     return top / x


def A(x):
    return k_2 - (k_1 + k_2) * x


def B(x, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)

    n = kwargs.get("nt", 1)

    return (1 / n) * (k2 + (k1 - k2) * x)


def C(x, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)

    num1 = 2 * k1
    num2 = (k1 - k2) * x - k2 * log(1 + ((k1 - k2) * x) / k2)
    den = (k1 - k2) ** 2
    result = x - (num1 * num2) / den
    return log(result)


def expon(x, **kwargs):
    n = kwargs.get("nt", 1)

    return exp(2 * n * C(x, **kwargs))


# @njit
def p_s(x, **kwargs):
    top = expon(x, **kwargs)
    bot = B(x, **kwargs)
    return top / bot


# def function(x):
#     output = empty_like(x)
#     for i, x in enumerate(x):
#         output[i] = B(x) ** (-1) * exp(2 * internal_results[i])
#     return output

kwarg_dict = {"k1": k_1, "k2": k_2, "nt": n}

# results = stationary_p(x_array, b)
# results = p_s(x_array, b, k_1, k_2)
# results = c(x_array, k1=k_1, k2=k_2)
# print(sum(results / max(results)))
# x_array = np.arange(alpha, beta, 0.001)

steps = 500
x_array = np.linspace(alpha, beta, num=steps, endpoint=False)[1:]

results = p_s(x_array, **kwarg_dict)

plot(x_array, results)

xlim(left=0, right=1)
ylim(bottom=0)

show()
