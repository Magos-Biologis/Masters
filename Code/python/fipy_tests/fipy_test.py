#! ../.venv/bin/python3
from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer

m = Grid1D(nx=100, Lx=1.0)

v0 = CellVariable(mesh=m, hasOld=True, value=0.5)
v1 = CellVariable(mesh=m, hasOld=True, value=0.5)

v0.constrain(0, m.facesLeft)
v0.constrain(1, m.facesRight)
v1.constrain(1, m.facesLeft)
v1.constrain(0, m.facesRight)

eq0 = TransientTerm() == DiffusionTerm(coeff=0.01) - v1.faceGrad.divergence
eq1 = TransientTerm() == v0.faceGrad.divergence + DiffusionTerm(coeff=0.01)

vi = Viewer((v0, v1))

from builtins import range

for t in range(100):
    v0.updateOld()
    v1.updateOld()
    res0 = res1 = 1e100

    while max(res0, res1) > 0.1:
        res0 = eq0.sweep(var=v0, dt=1e-5)
        res1 = eq1.sweep(var=v1, dt=1e-5)

    if t % 10 == 0:
        vi.plot()

# The uncoupled method still works, but it can be advantageous to solve the two equations simultaneously. In this case, by coupling the equations, we can eliminate the explicit sources and dramatically increase the time steps:

v0.value = 0.5
v1.value = 0.5

eqn0 = TransientTerm(var=v0) == DiffusionTerm(0.01, var=v0) - DiffusionTerm(1, var=v1)
eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) + DiffusionTerm(0.01, var=v1)

eqn = eqn0 & eqn1

from builtins import range

for t in range(1):
    v0.updateOld()
    v1.updateOld()
    eqn.solve(dt=1.0e-3)
    vi.plot()

# It is also possible to pose the same equations in vector form:

v = CellVariable(mesh=m, hasOld=True, value=[[0.5], [0.5]], elementshape=(2,))

v.constrain([[0], [1]], m.facesLeft)
v.constrain([[1], [0]], m.facesRight)

eqn = TransientTerm([[1, 0], [0, 1]]) == DiffusionTerm([[[0.01, -1], [1, 0.01]]])

vi = Viewer((v[0], v[1]))

from builtins import range

for t in range(1):
    v.updateOld()
    eqn.solve(var=v, dt=1.0e-3)
    vi.plot()


### -------------------------------------------------------------------------

# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import quad
#
# # --- PARAMETERS ---
# # Reaction rate constants (you can change these as needed)
# k1 = 0.8  # Replication: X -> 2X
# k3_prime = 0.5  # Effective selective pressure (scaled by constant B)
# k4_prime = 1.0  # Elimination rate constant (scaled by constant B)
# k2 = 0.5  # (Not used in this one-dimensional reduction)
# k5_prime = 0.5  # (Not used in this one-dimensional reduction)
#
# y = 1.0  # Constant value for Y (treated as a parameter)
# V = 100.0  # System volume
#
# # Effective drift and diffusion parameters for x:
# # Drift: A(x) = a * x, with a = k1 + k3'*y - k4'
# # Diffusion: B(x) = D * x, with D = (k1 + k3'*y + k4')/V
# a = k1 + k3_prime * y - k4_prime
# D = (k1 + k3_prime * y + k4_prime) / V
#
# print("Effective drift coefficient a =", a)
# print("Effective diffusion coefficient D =", D)
# if a >= 0:
#     print(
#         "Warning: a should be negative for a proper stationary distribution (decaying at infinity)."
#     )
#
#
# # --- STATIONARY DISTRIBUTION ---
# # With reflecting boundaries, the stationary zero-current solution is obtained from:
# # d/dx[B(x)P(x)] = 2 A(x) P(x)  ->  P(x) ∝ (1/x) exp((2a/D) x)
# def P_unnorm(x):
#     return (1.0 / x) * np.exp((2 * a / D) * x)
#
#
# # Compute normalization constant over x in [1, ∞)
# norm, err = quad(P_unnorm, 1, np.inf)
# print("Normalization constant =", norm)
#
#
# # Normalized stationary probability density
# def P(x):
#     return P_unnorm(x) / norm
#
#
# # --- PLOTTING ---
# # Define x values (from 1 to some cutoff; choose cutoff large enough so the tail is negligible)
# x_vals = np.linspace(1, 10, 1000)
# P_vals = np.array([P(x) for x in x_vals])
#
# plt.figure(figsize=(8, 6))
# plt.plot(x_vals, P_vals, label=r"Stationary $P(x)$")
# plt.ylim(bottom=0)
#
# plt.xlabel(r"$x$")
# plt.ylabel(r"$P(x)$")
# plt.title(
#     "Stationary Solution of the FPE\n(Reflecting boundaries at $x=1$ and $x=\infty$)"
# )
# plt.legend()
# # plt.grid(True)
# plt.show()

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
