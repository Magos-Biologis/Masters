import os

import numpy as np
import scipy as sy
import matplotlib.pyplot as plt

# from scipy.integrate import quad


figure_env: str | None = os.getenv("THESIS_FIGURE_PATH")
figure_dir = os.path.join(figure_env, "fpe")

file_name = "simple_fde"

k_1 = 1
k_2 = 1
n_T = 1000

file_name += "_" + f"k1{k_1}" + "_" + f"k2{k_2}" + "_" + f"n{n_T}"

a = 0
b = 1


def A(x, k1, k2):
    return k2 - (k1 + k2) * x


def B(x, k1, k2, nT):
    return (k2 + (k1 - k2) * x) / nT


def integrand(x_prime, k1, k2, nT):
    return A(x_prime, k1, k2) / B(x_prime, k1, k2, nT)


def unnormalized_ps(x):
    integral, _ = sy.integrate.quad(integrand, a, x)
    return np.exp(2 * integral) / B(x, k_1, k_2, n_T)


normalization_integral, _ = sy.integrate.quad(unnormalized_ps, a, b)
N = 1 / normalization_integral  # N ensures that the total probability is 1


# Define the normalized steady state distribution
def ps(x):
    return N * unnormalized_ps(x)


x_vals = np.linspace(a, b, 400)
ps_vals = np.array([ps(x) for x in x_vals])

# figsize=(8, 5)
fig, ax = plt.subplots()

plt.style.use("bmh")
ax.plot(x_vals, ps_vals, label="$p_s(x)$")

ax.set_ylim(bottom=0)
ax.set_xlim(0, 1)
ax.set_xticks(np.linspace(0, 1, 11))

ax.set_xlabel("$x$")
ax.set_ylabel("$p_s(x)$")

# plt.legend()

figure = os.path.join(figure_dir, file_name)

plt.savefig(figure + ".pdf")
# plt.show()
