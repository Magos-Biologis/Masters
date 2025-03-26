import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


k_1 = 1
k_2 = 1
n_T = 1


# --- Define your functions ---
# Feel free to change these functions. They represent your drift (A) and diffusion (B).
def A(x):
    # Example drift: a simple linear function (you can change this)
    # return 1.0 - x
    return k_2 - (k_1 + k_2) * x


def B(x):
    # Example diffusion: constant (you can change this)
    return (k_2 + (k_1 - k_2) * x) / n_T


# Lower limit for the integration
a = 0.0


# --- Define the integrand and the unnormalized steady state distribution ---
def integrand(x_prime):
    return A(x_prime) / B(x_prime)


def unnormalized_ps(x):
    # Compute the integral from a to x of A(x')/B(x')
    integral, _ = quad(integrand, a, x)
    # Note: We divide by B(x) as per the formula.
    return np.exp(2 * integral) / B(x)


# --- Compute the normalization constant ---
# Define the domain over which you want to analyze p_s(x)
x_min, x_max = 0.0, 1.0
normalization_integral, _ = quad(unnormalized_ps, x_min, x_max)
N = 1 / normalization_integral  # N ensures that the total probability is 1


# Define the normalized steady state distribution
def ps(x):
    return N * unnormalized_ps(x)


# --- Evaluate and plot p_s(x) ---
x_vals = np.linspace(x_min, x_max, 400)
ps_vals = np.array([ps(x) for x in x_vals])

plt.figure(figsize=(8, 5))
plt.plot(x_vals, ps_vals, label="$p_s(x)$")
plt.xlabel("x")
plt.ylabel("$p_s(x)$")
plt.title("Steady State Distribution $p_s(x)$")
plt.legend()
plt.grid(True)
plt.show()
