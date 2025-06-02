import os

from gillespie import analytical as dga
from pylab import *

data_env = str(os.getenv("FPE_FIGURE_ENV"))

domain = linspace(0, 1, 200)

xb = (0, 1)
yb = (0, 1)


kwags = {"density": 250}

system = dga.MultivariateTwoChemical(
    x_bounds=xb,
    y_bounds=yb,
    k1=1,
    k2=1,
    n=10,
    # mask=True,
)
# name = "ssa:2S:I7:C1:T1748705657984731:analytical"

# system = dga.MultivariateTwoChemical(
#     x_bounds=xb,
#     y_bounds=yb,
#     k1=1,
#     k2=2,
#     n=10,
#     **kwags,
#     # mask=True,
# )
# name = "ssa:2S:I8:C1:T17487132184568338:analytical"

# name += ":unmasked"

n = 10
b = 0
name = f"ssa:5_2:I?:C1:T?:analytical_n={n}_b={b}"
# system = dga.TwoParameterNovel(
#     x_bounds=xb,
#     y_bounds=yb,
#     k1=1,
#     k2=1,
#     n=n,
#     b=b,
#     # mask=True,
# )


grid, output = system.stationary()
x, y = grid
contours = output
# grid, gradient = system.gradient()
# x, y = grid
# u, v = gradient
# contours = np.sqrt(u**2 + v**2)
# equal_line = system.get_equality()


fig, ax = subplots(figsize=(6, 5))


# ax.plot(*equal_line)

axe = ax.contourf(
    x,
    y,  # [::-1],
    np.exp(contours),
    # np.exp(contours),
    # levels=256,  # 2 * n,
    # levels=10,
    cmap="bone_r",
    # norm=mpl.colors.LogNorm(),
)

# ax.set_facecolor((0.9463847846200787, 0.9656862745098039, 0.9656862667892201, 1.0))
# ax.set_facecolor((1.0, 1.0, 1.0, 1.0))
# ax.quiver(x, y, u, v)


# norm=mpl.colors.LogNorm(),
# ax.plot(*equal_line)
fig.colorbar(axe)

font_kwargs: dict = {
    "fontsize": 15,
}

ax.set_ylabel(
    # "Fraction of $c_2$",
    "Fraction of $y$",
    fontdict=font_kwargs,
)
ax.set_xlabel(
    # "Fraction of $c_1$",
    "Fraction of $x$",
    fontdict=font_kwargs,
)

ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
ax.set_xticks(ticks)
ax.set_yticks(ticks)

ax.set_xlim(0, 1.0)
ax.set_ylim(0, 1.0)

ax.set_title("Analytical Solution of the Distribution for Two Chemical Species")
# axe = ax.plot(*output)

# file_name = os.path.join(data_env, "five_var", name + ".pdf")
# file_name = os.path.join(data_env, "five_var", name + ".pdf")
file_name = os.path.join(data_env, name + ".pdf")

savefig(file_name, format="pdf")
show()
