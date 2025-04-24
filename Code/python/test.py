#!./.venv/bin/python
import numba
from numba import njit
from pylab import *

from scipy import stats as sps

print()  ## to make running in neovim easier
print()  ## to make running in neovim easier


b = 100


alpha = 0
beta = 1


n = 10

k_1 = 1
# k_1 /= n

k_2 = 1
# k_2 /= n


a_1 = n
a_2 = 0

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


# def C(x, **kwargs):
#     k1 = kwargs.get("k1", 1)
#     k2 = kwargs.get("k2", 1)
#
#     num1 = 2 * k1
#     num2 = (k1 - k2) * x - k2 * log(1 + ((k1 - k2) * x) / k2)
#     den = (k1 - k2) ** 2
#
#     prod = num1 * num2 * (den ** (-1))
#
#     result = x - (num1 * num2) / den
#
#     return result


def expon(x, **kwargs):
    n = kwargs.get("nt", 1)
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)
    max_value = kwargs.get("max_value", 0)

    num1 = 2 * k1
    num2 = (k1 - k2) * x - k2 * log(1 + ((k1 - k2) * x) / k2)
    den = (k1 - k2) ** 2

    return exp(2 * n * (x - (num1 * num2 / den)) - max_value)


def prop_expon(x, **kwargs):
    n = kwargs.get("nt", 1)
    k1 = kwargs.get("k1", 1)
    k2 = kwargs.get("k2", 1)
    max_value = kwargs.get("max_value", 0)

    if k1 != k2:
        num1 = 2 * k1
        num2 = (k1 - k2) * x - k2 * log(1 + ((k1 - k2) * x) / k2)
        den = (k1 - k2) ** 2
        res = (num1 * num2) / den
    else:
        res = (k1 / 2) * (x**2)

    return exp(2 * n * (x - res) - max_value)


# @njit
def p_s(x, **kwargs):
    top = expon(x, **kwargs)
    bot = B(x, **kwargs)
    return top / bot


def prop_p_s(x, **kwargs):
    top = prop_expon(x, **kwargs)
    bot = B(x, **kwargs)
    return top / bot


# def function(x):
#     output = empty_like(x)
#     for i, x in enumerate(x):
#         output[i] = B(x) ** (-1) * exp(2 * internal_results[i])
#     return output


# results = stationary_p(x_array, b)
# results = p_s(x_array, b, k_1, k_2)
# results = c(x_array, k1=k_1, k2=k_2)
# print(sum(results / max(results)))
# x_array = arange(alpha, beta, 0.001)


steps = 500
x_array = linspace(alpha, beta, num=steps, endpoint=False)[1:]

kwarg_dict = {"k1": k_1, "k2": k_2, "nt": n}

# exponential_array = log(2 * kwarg_dict["nt"] * C(x_array, **kwarg_dict))

# exponential_array = 2 * n * C(x_array, **kwarg_dict)
# expon_max = exp(exponential_array.max())


# results = p_s(x_array, **kwarg_dict)
results2 = prop_p_s(x_array, **kwarg_dict)


fig, ax = subplots()
# ax.plot(x_array, results)
ax.plot(x_array, results2)

# print(results.max()c


xlim(0, 1)
# ylim(bottom=0)


# yscale("log")
# plot(norm)
# plot(x_array, kde(x_array))
# plot(x_array, kde2(x_array))


show()

# diff_array = subtract(log(results), log(results2))
# greatest_diff = diff_array.max() - diff_array.min()

# print()
# print(diff_array)

# print("The largest difference in the functions, proportionally, is:", greatest_diff)
# print(diff_array)
# print(diff_array.max(), "\t", diff_array.min())

# def aj(x):
#     a, b = x
#     return array([k_1 * a, k_2 * b])
