#!../.venv/bin/python

import os

import numpy as np
import scipy as sy
import matplotlib.pyplot as plt

from numba import njit

from fipy import CellVariable, Grid1D, Viewer, DiffusionTerm


# from scipy.integrate import quad


figure_env: str | None = os.getenv("THESIS_FIGURE_PATH")
figure_dir: str | None = os.path.join(figure_env, "fpe")

file_name = "simple_fde"

k_11 = 1
k_21 = 1

k_12 = 1
k_22 = 1

n_T1 = 1
n_T2 = 10
n_T3 = 100

file_name += "_" + f"k1{k_11}" + "_" + f"k2{k_21}"
# file_name +=  "_" + f"n{n_T1}"

a = 0
b = 1


def A(x, k1, k2):
    return k2 - (k1 + k2) * x


def B(x, k1, k2, nT):
    return (k2 + (k1 - k2) * x) / nT


def integrand(x_prime, k1, k2, nT):
    return A(x_prime, k1, k2) / B(x_prime, k1, k2, nT)


@njit
def integrand_two(x, k1, k2, nT):
    num = nT * (
        2 * k1 * k2 * np.log(k2)
        - 2 * k1 * k2 * np.log(k2 + (k1 - k2) * x)
        + x * k1**2
        - x * k2**2
    )
    den = (k1 - k2) ** 2
    return num / den


def unnormalized_ps(x, k1, k2, nT):
    integral, _ = sy.integrate.quad(integrand, a, x, args=(k1, k2, nT))
    return np.exp(2 * integral) / B(x, k1, k2, nT)
    # return np.exp(2 * integrand_two(x, k1, k2, nT)) / B(x, k1, k2, nT)


# normalization_integral, _ = sy.integrate.quad(unnormalized_ps, a, b, args=(k_1, k_2, n_T))
# N = 1 / normalization_integral  # N ensures that the total probability is 1


# Define the normalized steady state distribution
def ps(x, k1, k2, nT):
    normalization_integral, _ = sy.integrate.quad(
        unnormalized_ps, a, b, args=(k1, k2, nT)
    )
    N = 1 / normalization_integral  # N ensures that the total probability is 1

    return N * unnormalized_ps(x, k1, k2, nT)


x_vals = np.linspace(a, b, 400)
ps_vals1 = np.array([ps(x, k_11, k_21, n_T1) for x in x_vals])
ps_vals2 = np.array([ps(x, k_11, k_21, n_T2) for x in x_vals])
ps_vals3 = np.array([ps(x, k_11, k_21, n_T3) for x in x_vals])

### ---------------------------------------------------------------------------
###
### ---------------------------------------------------------------------------

figsize = (8, 5)
fig, ax = plt.subplots()

plt.style.use("bmh")
ax.plot(x_vals, ps_vals1, label="$p_s(x)$ with $n_T = 1$", alpha=0.7, linestyle="dotted")
ax.plot(x_vals, ps_vals2, label="$p_s(x)$ with $n_T = 10$", alpha=0.7, linestyle="dashed")
ax.plot(x_vals, ps_vals3, label="$p_s(x)$ with $n_T = 100$", alpha=0.7)

ax.set_ylim(bottom=0)
ax.set_xlim(0, 1)
ax.set_xticks(np.linspace(0, 1, 11))

ax.set_xlabel("$x$")
ax.set_ylabel("$p_s(x)$")

ax.vlines(
    [k_21 / (k_11 + k_21)], ymin=0, ymax=100, linestyles="dashed", alpha=0.4, color="k"
)

plt.legend(loc="upper left")

# figure = os.path.join(figure_dir, file_name)
# plt.savefig(figure + ".pdf")

plt.show()
