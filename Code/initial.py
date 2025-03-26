#!./.venv/bin/python
import os
from pprint import pprint as pp
from sys import exit

import numpy as np
from numpy import polynomial as npp


from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter
# import sympy as sy

print("")


figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_path = os.path.join(figure_env, "ode")

# figure_path = "./figs"


### Variables
# Growth Constant, same because same cell
k = np.zeros(2)
n = np.zeros(2)
w = np.zeros(2)
q = np.zeros(2)

k[:] = 1

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = 100
n[1] = 100

w[0] = 0.015
w[1] = 0.035

# n[1] = 90
# w[1] = 0.015

q[0] = 0.999
q[1] = 0.8

omega_1 = k[0] - w[0]
omega_2 = k[1] - w[1]

k_1 = k[0] / n[0]
k_2 = k[1] / n[1]


c1_min = w[1] / k_1
c2_min = w[0] / k_2


# transition_matrix = sy.Matrix([[1 - w[0], w[0]], [w[1], 1 - w[1]]])
# init_condition = sy.Matrix([[n_1], [0]])
#
# pop_norm = np.sqrt(n[0] ** 2 + n[1] ** 2)
# n1_wt = n[0] / pop_norm
# n2_wt = n[1] / pop_norm
# print(transition_matrix)
# exit()

# ct_1 = n[0] * ((k[0] - 1) / k[0])
# ct_2 = n[1] * ((k[1] - 1) / k[1])
#
# w1w2 = w[0] / w[1]
# w2w1 = w[1] / w[0]
#
# alpha = 0
# c1_0 = n[0] - alpha
# c2_0 = alpha
# m_0 = 0
#
#
# # eta_1 = omega_1 - ku1
# # eta_2 = omega_2 - ku2
# # trace = (i_0 * q_1 - omega_2 - omega_1)
# # dete  = (omega_1 - i_0 * q_1) * omega_2 - w[0] * w[1]
# # disc  = trace**2 - dete
#
# dt = 0.01
# t_end = 250
# t_array = np.arange(0, t_end, dt)
# steps = np.zeros((len(t_array), 3))
#
# parameters = [k, w, q, n]
#
#
# cq = sum(n) / len(n)
#
#
# def step_function(xs, dt) -> np.ndarray:
#     c_1, c_2, m = xs
#     ct = c_1 + c_2
#     # ct = cq
#
#     c1_dot = (k[0] * (1 - ct / n[0]) - w[0]) * c_1 + w[1] * c_2 - q[0] * m * c_1
#     c2_dot = (k[1] * (1 - ct / n[1]) - w[1]) * c_2 + w[0] * c_1
#     m_dot = -q[1] * m * c_2 + m_0
#
#     return dt * np.array([c1_dot, c2_dot, m_dot])
#
#
# steps[0, :] = c1_0, c2_0, m_0
# for j, _ in enumerate(t_array):
#     steps[j, :] += steps[j - 1, :]
#     steps[j, :] += step_function(steps[j - 1, :], dt)  # print(steps[j, :])
#
#
# M = np.array(
#     [
#         [(omega_1 - k_1), w[0]],
#         [w[1], (omega_2 - k_2)],
#     ]
# )
#
#
# probability_sol = np.zeros(2)
#
# solution_norm = M[0, 1] + M[1, 0]
# probability_sol[0] = M[1, 0] / solution_norm
# probability_sol[1] = M[0, 1] / solution_norm
# probability_ticks = w / sum(w)
#
#
# # print(solutions)
#
#
# c1_fixed = n[0]
# c2_fixed = (1 - (w[0] / w[1])) * n[1]
# m_fixed = m_0 / q[1]
#
#
# dt_zeroes = [c1_fixed, c2_fixed, m_fixed]
#
# # total = steps[:, 0] + steps[:, 1]
# # steps[:, 0] = steps[:, 0] / total
# # steps[:, 1] = steps[:, 1] / total
#
# # print(2 * ((n[0] * n[1]) / (n[0] + n[1])))
# # print(total)
# # exit()
#
#
# color_list = ["b", "r", "g"]
# name_list = [
#     "Non-resistant Type",
#     "Resistant Type",
#     "Antibiotic",
# ]
# dt_name_list = [
#     r"$c_1$ fixed point at $c_1^*$",  # {c1_fixed}",
#     r"$c_2$ fixed point at $c_2^*$",  # {c2_fixed}",
#     f"$m$ fixed point at {m_fixed}",
# ]
#
# ## Polynomial
#
# coeffs = [
#     n[1] * w[1],
#     omega_1 - n[1] * (w[1] / n[0] + k_1),
#     k_1 * ((n[1] / n[0]) - 1),
# ]
# poly = npp.Polynomial(
#     coef=coeffs,
#     symbol="_c1",
# )
# roots = poly.roots()
#
# c1_root = [root for root in roots if root > c1_min]
# c2_root = (1 - c1_root / n[0]) * n[1]
#
# solutions = np.array([c1_root, c2_root])
#
### Plotting

fig, ax = plt.subplots()

# # ax.plot(steps.T[0], steps.T[1])
# # ax.plot(t_array, (steps.T[0] / steps.T[1]))
# # ax.hlines(c2_fixed, 0, t_end)
#
# for i, curve in enumerate(steps.T):
#     # break
#     if i == 2:
#         break
#     color = color_list[i]
#     curve_name = name_list[i]
#     dt_curve_name = dt_name_list[i]
#     # zero_line = dt_zeroes[i]
#     zero_line = solutions[i]
#
#     ax.hlines(
#         zero_line,
#         0,
#         t_end,
#         label=dt_curve_name,
#         color=color,
#         linestyle="dashed",
#         linewidth=1,
#     )
#
#     ax.plot(t_array, curve, label=curve_name, color=color)
#
# ax.set_xlim(0, t_end)
# ax.set_ylim(0, 100)
# # ax.set_ylim(0, 100)
#
# # print()
# # exit()
#
# total = steps.T[0] + steps.T[1]
# percent_total = np.divide(steps.T[0], n[0]) + np.divide(steps.T[1], n[1])
#
# # pp(percent_total)
# # exit()
#
# # fixed_point_tick_labels = [
# #     r"$w_1 \over w_1 + w_2$" + f" = {np.round(solutions[0], 3)}",
# #     r"$w_2 \over w_1 + w_2$" + f" = {np.round(solutions[1], 3)}",
# # ]
#
# fixed_point_tick_labels = [
#     r"$c_1^*$",
#     r"$c_2^*$",
# ]
#
# # ax1.yaxis.set_ticklabels()
#
# new_y_ticks = np.append(ax.get_yticks(), solutions)
# ax.set_yticks(new_y_ticks)
#
#
# def fixed_point_format(val, pos):
#     if val == solutions[0]:
#         return fixed_point_tick_labels[0]
#     elif val == solutions[1]:
#         return fixed_point_tick_labels[1]
#     else:
#         return int(np.round(val, 3))
#
#
# ax.yaxis.set_major_formatter(FuncFormatter(fixed_point_format))

alpha0 = 1
a = 10


# def alpha(t):
#     return alpha0 * np.exp(-k[0] * t) + (k[1] * a / k[0]) * (1 - np.exp(-k[0] * t))
#
#
# import math
#
#
# def poisson(x, t):
#     return (np.exp(-alpha(t)) * alpha(t) ** x) / (math.factorial(x))
#
#
# # t_span = np.linspace(-50, -30, 1000)
# t_span = np.linspace(-20, 20, 100)
# x_span = np.arange(0, 2 * a, 1)

# ax.plot(x_span, poisson(x_span, 100))


# def t_pos(x):
#     return k[1] * a
#
#
# def t_neg(x):
#     return k[0] * x
#
#
# def P_s(x):
#     if x == 0:
#         return 1
#
#     num = t_pos(x - 1)
#     den = t_neg(x)
#     return (num / den) * P_s(x - 1)
#
#
# x_plot = [P_s(x) for x in x_span]
# ax.plot(x_span, x_plot)


ax.set_xlabel("$x$")
# ax.set_ylabel("Count")

# plt.legend()
plt.tight_layout()

figure_file = os.path.join(figure_path, "plot.pdf")
# plt.savefig(figure_file)
plt.show()
