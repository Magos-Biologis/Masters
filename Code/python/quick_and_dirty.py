import os

import numpy as np
from pylab import *

# data_env = str(os.getenv("FPE_FIGURE_ENV"))
fig_env = str(os.getenv("THESIS_FIGURE_PATH"))

# domain = linspace(0, 1, 200)
#
# xb = (0, 1)
# yb = (0, 1)


rcparams_kwargs: dict = {
    "axes.labelsize": 20,
    "axes.titleweight": "bold",
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "savefig.facecolor": "#c0c0ca",
}
# "axes.titlecolor": "white",
# "xtick.labelcolor": "white",
# "ytick.labelcolor": "white",

keys = rcparams_kwargs.keys()
values = rcparams_kwargs.values()

for k, v in zip(keys, values):
    matplotlib.rcParams[k] = v

fig, ax = subplots(figsize=(6, 5))

name = "Novel1VarData"


num_data = np.load("../julia/{}.npz".format(name))
x = num_data["x"]
logy = num_data["log(y)"]


f = np.exp(logy)

ax.plot(x, f)


# axe = ax.contourf(
#     np.exp(contours),
#     # np.exp(contours),
#     # levels=256,  # 2 * n,
#     # levels=10,
#     cmap="bone_r",
#     # norm=mpl.colors.LogNorm(),
# )
# ax.set_facecolor((0.9463847846200787, 0.9656862745098039, 0.9656862667892201, 1.0))
# ax.set_facecolor((1.0, 1.0, 1.0, 1.0))
# ax.quiver(x, y, u, v)
# norm=mpl.colors.LogNorm(),
# ax.plot(*equal_line)
# fig.colorbar(axe)


font_kwargs: dict = {
    "fontsize": 15,
}

ax.set_ylabel(
    # "Fraction of $c_2$",
    "Distribution",
    fontdict=font_kwargs,
)
ax.set_xlabel(
    # "Fraction of $c_1$",
    "Fraction of $x$",
    fontdict=font_kwargs,
)


# ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
# ax.set_xticks(ticks)
# ax.set_yticks(ticks)

ax.set_xlim(0 - 0.01, 1.0)
ax.set_ylim(0, 1.0)

ax.set_title("Analytical solution of 1 variable steady state")


# axe = ax.plot(*output)

# file_name = os.path.join(data_env, "five_var", name + ".pdf")
# file_name = os.path.join(data_env, "five_var", name + ".pdf")

file_name = os.path.join(fig_env, "{}.pdf".format(name))

savefig(file_name, format="pdf")
show()
