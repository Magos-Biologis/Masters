#!./.venv/bin/python

from itertools import chain
import os
from pprint import pprint as pp
from sys import exit

import numpy as np
from numpy import max, sum
from numpy import random as npr

from matplotlib import pyplot as plt
from matplotlib.pyplot import FuncFormatter

# import sympy as sy

thesis_env = str(os.getenv("THESIS_FIGURE_PATH"))
figure_path = os.path.join(thesis_env, "markov")
# figure_path = "./figs"

print()

# Growth Constant, same because same cell
k = np.zeros(2)
n = np.zeros(2)
w = np.zeros(2)
q = np.zeros(2)

k[0] = 2
k[1] = 2

# Population cap (purely aesthetic if n₁ = n₂)
n[0] = 100
n[1] = 100


w[0] = 0.015
w[1] = 0.035

w[0] = 0.7
w[1] = 0.4

q[0] = 0.999
q[1] = 0.8

m_0 = 0

c_t = sum(n) / len(n)

k1 = k[0] / n[0]
k2 = k[1] / n[1]

om1 = k[0] - w[0]
om2 = k[1] - w[1]

# t_array = np.arange(0, t_end, dt)
# steps = np.zeros((len(t_array), 3))

var_name = [r"$x$", r"$y$"]

## Markov Matrix multiplication
# state_end = 150
state_end = 3
dt = 1

## Simulated Markov process
t_end = 1500
# t_end = 25000

a = 1
c1_0 = n[0] - a
c2_0 = a
m__0 = m_0

tran_mat3 = np.array(
    [
        # [k[0] * (1 - (1 / n[0]))-w[0], w[1], q[0] * m_0],
        # [w[0], k[1] * (1 - (1 / n[1]))-w[1], 0],
        [(om1 - k1), w[0], q[0] * m_0],
        [w[1], (om2 - k2), 0],
        [0, 0, 1],
    ]
)

tran_mat2 = np.array(
    [
        [(om1 - k1 * c_t), w[0]],
        [w[1], (om2 - k2 * c_t)],
    ]
)

tran_mat22 = np.array(
    [
        [(om1 - k1), w[0]],
        [w[1], (om2 - k2)],
    ]
)


det2 = tran_mat2[0, 0] * tran_mat2[1, 1] - tran_mat2[0, 1] * tran_mat2[1, 0]


def vector_norm(matrix):
    res_matrix = np.zeros_like(matrix)
    for n in range(matrix.shape[0]):
        intermediate = 0
        for m in range(matrix.shape[1]):
            intermediate += matrix[n, m]

        res_matrix[n, :] = matrix[n, :] / intermediate
    return res_matrix


# def step_function(xs, dt) -> np.ndarray:
#     c_1, c_2, m = xs
#
#     c1_dot = c_1 * tran_mat3[0, 0] * c_1 + tran_mat3[0, 1] * c_2 + tran_mat3[0, 2] * m
#     c2_dot = c_1 * tran_mat3[1, 0] * c_1 + tran_mat3[1, 1] * c_2 + tran_mat3[1, 2] * m
#     m_dot = c_1 * tran_mat3[2, 0] * c_1 + tran_mat3[2, 1] * c_2 + tran_mat3[2, 2] * m
#
#     return dt * np.array([c1_dot, c2_dot, m_dot])
# steps[0, :] = c1_0 - c2_0, c2_0, m__0
# for j, _ in enumerate(t_array):
#     steps[j, :] += steps[j - 1, :]
#     steps[j, :] += step_function(steps[j - 1, :], dt)  # print(steps[j, :])
# dt_zeroes = [c1_fixed, c2_fixed, m_fixed]

color_list = ["b", "r", "g"]
name_list = [
    "Non-resistant Type",
    "Resistant Type",
    "Antibiotic",
]

# dt_name_list = [
#     f"$c_1$ fixed point at {c1_fixed}",
#     f"$c_2$ fixed point at {c2_fixed}",
#     f"$m$ fixed point at {m_fixed}",
# ]


# row1_sum = np.abs(sum(tran_mat2[0, :]))
# row2_sum = np.abs(sum(tran_mat2[1, :]))

# print(tran_mat2, "   ", tran_mat2[0, :])
# exit()


def rownorm(matrix):
    row1_sum = np.abs(sum(matrix[0, :]))
    row2_sum = np.abs(sum(matrix[1, :]))

    matrix[0, :] = matrix[0, :] / row1_sum
    matrix[1, :] = matrix[1, :] / row2_sum

    return matrix


# tran_mat2 = rownorm(tran_mat2)
# tran_mat22 = rownorm(tran_mat22)


naive_mat = np.array([[1 - w[0], w[0]], [w[1], 1 - w[1]]])


M = naive_mat


# naive_mat_left_eigen = np.array(
#     [
#         [
#             naive_mat[1, 0] / (naive_mat[1, 0] + naive_mat[0, 1]),
#             naive_mat[0, 1] / (naive_mat[1, 0] + naive_mat[0, 1]),
#         ],
#         [1, -1],
#     ]
# )
# naive_norm_1 = sum([x**2 for x in naive_mat_left_eigen[0, :]])
# naive_norm_2 = sum([x**2 for x in naive_mat_left_eigen[1, :]])


# solutions[0] = w[1] / np.sum(w)
# solutions[1] = w[0] / np.sum(w)
# solutions[0] = c1_0 * (w[1] / (w[0] + w[1])) + c2_0
# solutions[1] = c1_0 * (w[0] / (w[0] + w[1])) - c2_0

# solutions[0] = w[0] / sum(w)  # - 1
# solutions[1] = w[1] / sum(w)  # + 1

solutions = np.zeros(2)

solution_norm = M[0, 1] + M[1, 0]
solutions[0] = M[0, 1] / solution_norm
solutions[1] = M[1, 0] / solution_norm

# print(M[0, 1])

c1_0 = 0.99
c2_0 = 1 - c1_0


# solutions *= n
# c1_0 *= n[0]
# c2_0 *= n[1]
# vector_norm = np.sum([x**2 for x in solution_vector])
# naive_solution
# pp(naive_mat_left_eigen)
# exit()
# pp(M)
# exit()

is_c1 = r"$c_1$"
is_c2 = r"$c_2$"
state = is_c1

c1_array = [c1_0]
c2_array = [c2_0]
t__array = [0]


### Matrix multiplication

state_vector = np.array([c1_0, c2_0])
states = []

t_axis = np.arange(0, state_end + 1, dt)
for t in t_axis:
    states.append(state_vector)
    state_vector = np.dot(state_vector, M)


fig, ax1 = plt.subplots()


# ax1.set_xticks(np.arange(0, state_end + dt, 1))
y_ticks = np.arange(0, 1.1, 0.25)
ax1.set_yticks(y_ticks)
ax1.set_yticklabels(y_ticks)

fixed_point_ticks = w / sum(w)
fixed_point_tick_labels = [r"$p_1 \over p_1 + p_2$", r"$p_2 \over p_1 + p_2$"]

new_y_ticks = np.append(ax1.get_yticks(), solutions)
# ax1.yaxis.set_ticklabels()
ax1.set_yticks(new_y_ticks)


def fixed_point_format(val, pos):
    if val == solutions[0]:
        return fixed_point_tick_labels[0]
    elif val == solutions[1]:
        return fixed_point_tick_labels[1]
    else:
        return val


ax1.yaxis.set_major_formatter(FuncFormatter(fixed_point_format))


ax1.set_xlim(left=0, right=state_end)
ax1.set_ylim(bottom=0, top=1)
# print(list(ax1.get_yticks()))


ax1.set_xlabel("Markov Chain Steps")
ax1.set_ylabel("Frequency")


ax1.hlines(
    solutions[1],
    0,
    state_end,
    color="blue",
    alpha=0.5,
    label=f"{var_name[0]} steady state",
    linestyle="dashed",
)
ax1.hlines(
    solutions[0],
    0,
    state_end,
    color="orange",
    alpha=0.5,
    label=f"{var_name[1]} steady state",
    linestyle="dashed",
)


ax1.plot(t_axis, states)

# plt.legend()
# plt.show()


### If else

tau = 0
npr.seed(69420)

u = np.zeros(3)
while t__array[tau] <= t_end:
    upsilon = tau + 1
    # print(t__array[tau])
    # if (tau + 1) % 10000 == 0:
    #     ax.scatter(t__array, c1_array)
    #     plt.pause(0.05)

    u = npr.uniform(0, 1, 3)
    while u[2] == 0:
        u = npr.uniform(0, 1, 3)

    c1_array.append(c1_array[tau])
    c2_array.append(c2_array[tau])
    t__array.append(t__array[tau])

    # np.append(c1_array, c1_array[tau])
    # np.append(c2_array, c2_array[tau])
    # np.append(t__array, t__array[tau])

    # scalar = (w[0] / w[1]) * c1_array[tau] + (w[0] / w[1]) * c2_array[tau]
    # time_modifier = -np.log(u[2]) / scalar

    # time_modifier = 0.002
    time_modifier = 1

    t__array[upsilon] += time_modifier

    if state == is_c1:
        if u[0] <= M[0, 0]:
            c1_array[upsilon] = c1_array[tau] + 1
            # c2_array[upsilon] = c2_array[tau] - 1
        else:
            # c1_array[upsilon] = c1_array[tau] - 1
            c2_array[upsilon] = c2_array[tau] + 1

            state = is_c2
    elif state == is_c2:
        if u[1] <= M[1, 1]:
            # c1_array[upsilon] = c1_array[tau] - 1
            c2_array[upsilon] = c2_array[tau] + 1
        else:
            c1_array[upsilon] = c1_array[tau] + 1
            # c2_array[upsilon] = c2_array[tau] - 1

            state = is_c1

    else:
        print("fuck you")
        break

    ## End of chain, increments at end of loop so arrays work properly
    tau += 1

# solution_vector = np.array([-w[0] / w[1], 1])
# vector_norm = sum([x**2 for x in solution_vector])
# vector_norm **= 1 / 2
# norm = np.sqrt(np.sum([x**2 for x in solution_vector]))
# solution_vector /= vector_norm


c1 = np.array(c1_array)
c2 = np.array(c2_array)
t = np.array(t__array)


chain_norm = np.add(c1, c2)

c1 = np.divide(c1, chain_norm)
c2 = np.divide(c2, chain_norm)

ax2 = ax1.twiny()
ax2.set_xlabel("Simulated Markov Steps")


ax2.plot(t, c1, label=var_name[0], color="blue")
ax2.plot(t, c2, label=var_name[1], color="orange")

# ax2.hlines(naive_norm_1, 0, t_end, color="blue")
# ax2.hlines(naive_norm_2, 0, t_end, color="orange")
# oldsol1 = -naive_mat_left_eigen[0, 0] * n[0]
# oldsol2 = naive_mat_left_eigen[0, 1] * n[1]

# ax2.hlines(solutions[1], 0, t_end, color="blue")
# ax2.hlines(solutions[0], 0, t_end, color="orange")

ax2.set_xlim(left=0, right=t_end)


# plt.xscale("log")

plt.legend()
plt.tight_layout()


file_path = os.path.join(figure_path, f"chain-{w[0]}-{w[1]}")
plt.savefig(file_path + ".pdf", format="pdf")


plt.show()

print("done")
exit()

for i in range(2):
    norm = sqrt(tran_mat2[i, 0] ** 2 + tran_mat2[i, 1] ** 2)
    # tran_mat2[i,:] *=

print(tran_mat2)
exit()

u = npr.random((3, 3))

a = 1
c1_0 = n[0] - a
c2_0 = a
m__0 = 0


populations = [c1_0, c2_0, m__0]


# norm = (om1 - k1 * c_t) * c1_0 + w[0] * c2_0
# norm = ((om1 - k1 * c_t) * c1_0 * w[0] + w[0] * c2_0 * (om1 - k1 * c_t)) / (
#     (om1 - k1 * c_t) * w[0]
# )

pp(norm)
# pp((c1_0 + c2_0) / (n[0] + n[1]))
exit()

t = 0
while t <= 100:
    t += 1
    if u[0, 0] <= tran_mat[0, 0]:
        pass
    elif tran_mat[0, 0] < u[0, 0] <= tran_mat[0, 1]:
        pass
    elif tran_mat[0, 1] < u[0, 0] <= tran_mat[0, 2]:
        pass


print(u)
