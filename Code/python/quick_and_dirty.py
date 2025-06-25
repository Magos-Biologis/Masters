import os

import numpy as np
from pylab import *

# data_env = str(os.getenv("FPE_FIGURE_ENV"))
fig_env = str(os.getenv("THESIS_FIGURE_PATH"))

# domain = linspace(0, 1, 200)
#
# xb = (0, 1)
# yb = (0, 1)


# rcparams_kwargs: dict = {
#     "axes.labelsize": 20,
#     "axes.titleweight": "bold",
#     "xtick.labelsize": 12,
#     "ytick.labelsize": 12,
#     # "figure.labelsize": 20,
#     "figure.titlesize": 20,
#     "savefig.facecolor": "#c0c0ca",
# }


matplotlib.rcParams["axes.titlesize"] = 18
matplotlib.rcParams["axes.labelsize"] = 16
matplotlib.rcParams["axes.titleweight"] = "bold"

matplotlib.rcParams["xtick.labelsize"] = 13
matplotlib.rcParams["ytick.labelsize"] = 13
matplotlib.rcParams["savefig.facecolor"] = "#c0c0ca"

matplotlib.rcParams["figure.titlesize"] = 20


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
)
ax.set_xlabel(
    # "Fraction of $c_1$",
    "Fraction of $x$",
)


# ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
# ax.set_xticks(ticks)
# ax.set_yticks(ticks)

ax.set_xlim(0 - 0.005, 1.0)
ax.set_ylim(0, 1.0)

ax.set_title("Analytical Solution of Steady State Distribution")


# axe = ax.plot(*output)

# file_name = os.path.join(data_env, "five_var", name + ".pdf")
# file_name = os.path.join(data_env, "five_var", name + ".pdf")

pres_path = os.path.join(fig_env, "..", "..", "Presentation", "images")
file_name = os.path.join(pres_path, "{}.pdf".format(name))

savefig(file_name, format="pdf")
show()
