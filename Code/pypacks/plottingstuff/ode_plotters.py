import numpy as np
from matplotlib import axes as axe

from .plotting_class import plotting_class


class ode_plotters(plotting_class):
    def plot_curves(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        solutions: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        fixed_points: np.ndarray[tuple[int], np.dtype[np.float64]],
        color_list: list[str],
        name_list: list[str],
    ):
        for i, curve in enumerate(solutions):
            if i == 2:
                break

            color = color_list[i]
            curve_name = name_list[i]
            zero_line = fixed_points[i]

            ax.hlines(
                zero_line,
                0,
                time[-1],
                color=color,
                linestyle="dashed",
                **self.lineargs,
            )

            ax.plot(
                time,
                curve,
                **self.curvargs,
                label=curve_name,
                color=color,
            )

    def _fixed_point_format(self, val, pos):
        if val == solutions[0]:
            return fixed_point_tick_labels[0]
        elif val == solutions[1]:
            return fixed_point_tick_labels[1]
        else:
            return int(np.round(val, 3))

    def _plot_nullcline(
        self,
        ax: axe.Axes,
        x_points,
        y_points,
        **kwargs,
    ) -> None:
        ax.plot(
            x_points,
            y_points,
            **self.curvargs,
            **kwargs,
        )

    def plot_nullclines(
        self,
        ax: axe.Axes,
        x_range,
        y_range,
        dx_nullcline,
        dy_nullcline,
    ) -> None:
        self._plot_nullcline(
            ax,
            x_range,
            dy_nullcline(x_range),
            label=r"$c_2$ Nullcline",
            color=c2_col,
        )
        self._plot_nullcline(
            ax,
            dx_nullcline(y_range),
            y_range,
            label=r"$c_2$ Nullcline",
            color=c2_col,
        )

    def plot_phase_fixed_points(
        self,
        ax: axe.Axes,
    ) -> None:
        ax.vlines(
            c1_root,
            *x_lims,
            color=c1_col,
            **fixed_kwargs,
        )
        ax.hlines(
            c2_root,
            *y_lims,
            color=c2_col,
            **fixed_kwargs,
        )

    def plot_phase_space(
        self,
        ax: axe.Axes,
        c1,
        c2,
        dU,
        dV,
        stream_kwargs,
        **kwargs,
    ) -> None:
        ax.streamplot(c1, c2, dU, dV, **stream_kwargs)
