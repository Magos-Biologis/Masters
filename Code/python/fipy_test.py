#! ../.venv/bin/python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# --- PARAMETERS ---
# Reaction rate constants (you can change these as needed)
k1 = 0.8  # Replication: X -> 2X
k3_prime = 0.5  # Effective selective pressure (scaled by constant B)
k4_prime = 1.0  # Elimination rate constant (scaled by constant B)
k2 = 0.5  # (Not used in this one-dimensional reduction)
k5_prime = 0.5  # (Not used in this one-dimensional reduction)

y = 1.0  # Constant value for Y (treated as a parameter)
V = 100.0  # System volume

# Effective drift and diffusion parameters for x:
# Drift: A(x) = a * x, with a = k1 + k3'*y - k4'
# Diffusion: B(x) = D * x, with D = (k1 + k3'*y + k4')/V
a = k1 + k3_prime * y - k4_prime
D = (k1 + k3_prime * y + k4_prime) / V

print("Effective drift coefficient a =", a)
print("Effective diffusion coefficient D =", D)
if a >= 0:
    print(
        "Warning: a should be negative for a proper stationary distribution (decaying at infinity)."
    )


# --- STATIONARY DISTRIBUTION ---
# With reflecting boundaries, the stationary zero-current solution is obtained from:
# d/dx[B(x)P(x)] = 2 A(x) P(x)  ->  P(x) ∝ (1/x) exp((2a/D) x)
def P_unnorm(x):
    return (1.0 / x) * np.exp((2 * a / D) * x)


# Compute normalization constant over x in [1, ∞)
norm, err = quad(P_unnorm, 1, np.inf)
print("Normalization constant =", norm)


# Normalized stationary probability density
def P(x):
    return P_unnorm(x) / norm


# --- PLOTTING ---
# Define x values (from 1 to some cutoff; choose cutoff large enough so the tail is negligible)
x_vals = np.linspace(1, 10, 1000)
P_vals = np.array([P(x) for x in x_vals])

plt.figure(figsize=(8, 6))
plt.plot(x_vals, P_vals, label=r"Stationary $P(x)$")
plt.ylim(bottom=0)

plt.xlabel(r"$x$")
plt.ylabel(r"$P(x)$")
plt.title(
    "Stationary Solution of the FPE\n(Reflecting boundaries at $x=1$ and $x=\infty$)"
)
plt.legend()
# plt.grid(True)
plt.show()

### -------------------------------------------------------------------------

# from fipy import CellVariable, Grid1D, Viewer, DiffusionTerm
#
#
# nx = 200
#
# dx = 0.01
#
# L = nx * dx
#
# mesh = Grid1D(dx=dx, nx=nx)
# potential = CellVariable(mesh=mesh, name="potential", value=0.0)
# permittivity = 1
#
# electrons = CellVariable(mesh=mesh, name="e-")
#
# electrons.valence = -1
# charge = electrons * electrons.valence
#
# charge.name = "charge"
#
#
# potential.equation = DiffusionTerm(coeff=permittivity) + charge == 0
# potential.constrain(0.0, mesh.facesLeft)
#
# electrons.setValue(1.0)
# potential.equation.solve(var=potential)
# x = mesh.cellCenters[0]
#
# analytical = CellVariable(mesh=mesh, name="analytical solution", value=(x**2) / 2 - 2 * x)
#
# # from fipy import input
# x = mesh.cellCenters[0]
#
# electrons.setValue(0.0)
#
# electrons.setValue(1.0, where=x > L / 2.0)
# potential.equation.solve(var=potential)
# analytical.setValue(-x)
#
# analytical.setValue(((x - 1) ** 2) / 2 - x, where=x > L / 2)
#
# if __name__ == "__main__":
#     viewer = Viewer(vars=(charge, potential, analytical))
#     viewer.plot()
#     # input("Press any key to continue...")
#

### -------------------------------------------------------------------------
