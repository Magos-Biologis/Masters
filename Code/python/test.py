#!./.venv/bin/python
from numba import njit
from pylab import *


b = 100


alpha = 0
beta = 1


k_1 = 1
k_2 = 3


n = 100


def A(x):
    return k_2 - (k_1 + k_2) * x


def B(x, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)

    n = kwargs.get("nt", 1)

    return (1 / n) * (k2 + (k1 - k2) * x)


@njit
def stationary_p(x, b):
    top = exp(-2 * x) * (b + x) ** (4 * b - 1)
    return top / x


def C(x, **kwargs):
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)

    num1 = 2 * k1
    num2 = (k1 - k2) * x - k2 * log(1 + ((k1 - k2) * x) / k2)
    den = (k1 - k2) ** 2
    result = x - (num1 * num2) / den
    return result


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

x_array = np.linspace(alpha, beta, 1000)
# results = stationary_p(x_array, b)
# results = p_s(x_array, b, k_1, k_2)

# results = c(x_array, k1=k_1, k2=k_2)
results = p_s(x_array, nt=n, k1=k_1, k2=k_2)
print(sum(results / max(results)))

plot(x_array, results)
show()
