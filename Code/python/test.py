import os

from gillespie import analytical as dga
from pylab import *

data_env = str(os.getenv("FPE_FIGURE_ENV"))

domain = linspace(0, 1, 200)

xb = (0, 1)
yb = (0, 1)

system = dga.MultivariateTwoChemical(x_bounds=xb, y_bounds=yb, k1=1, k2=2, n=10)


# grid, output = system.stationary(domain, domain)
grid, output = system.stationary()
cont = diag(output)[::-1]
output = diag(cont)

equal_line = system.get_equality()


scale = 1
plot_space = multiply(grid, scale)


fig, ax = subplots(figsize=(5, 5))
axe = ax.contourf(*plot_space, output, 25)
# ax.plot(*equal_line)
# fig.colorbar(axe)

font_kwargs: dict = {
    "fontsize": 12,
}

ax.set_ylabel(
    "Fraction of $y$",
    fontdict=font_kwargs,
)
ax.set_xlabel(
    "Fraction of $x$",
    fontdict=font_kwargs,
)

ticks = [0, 0.2, 0.4, 0.6, 0.8, 1]
ax.set_xticks(ticks)
ax.set_yticks(ticks)

# axe = ax.plot(*output)

# name = "ssa:2S:I7:C1:T1748705657984731:analytical"
name = "ssa:2S:I8:C1:T17487132184568338:analytical"
file_name = os.path.join(data_env, "five_var", name + ".pdf")

# savefig(file_name, format="pdf")
show()
