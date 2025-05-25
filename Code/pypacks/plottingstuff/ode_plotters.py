import numpy as np
from matplotlib import axes as axe
from matplotlib.ticker import FuncFormatter
from myodestuff import ODEModel

from .plotting_class import plotting_class


class odePlotters(plotting_class):
    def __init__(self, stream_kwargs: dict, *args, **kwargs):
        self.streamargs = stream_kwargs
        super().__init__(*args, **kwargs)
        # self.init_something(param2)

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

    def _plot_phase_curve(
        self,
        ax: axe.Axes,
        x_points: np.ndarray[tuple[int], np.dtype[np.float64]]
        | np.ndarray[tuple[int], np.dtype[np.int_]],
        y_points: np.ndarray[tuple[int], np.dtype[np.float64]]
        | np.ndarray[tuple[int], np.dtype[np.int_]],
        **kwargs,
    ) -> None:
        self.curvargs.update(kwargs)

        ax.plot(
            x_points,
            y_points,
            **self.curvargs,
        )

    def plot_trajectories(
        self,
        ax: axe.Axes,
        dx: np.ndarray[tuple[int], np.dtype[np.float64]]
        | np.ndarray[tuple[int], np.dtype[np.int_]],
        dy: np.ndarray[tuple[int], np.dtype[np.float64]]
        | np.ndarray[tuple[int], np.dtype[np.int_]],
        **kwargs,
    ) -> None:
        label = kwargs.pop("label", self.x_name)
        color = kwargs.pop("color", self.colors[0])
        self._plot_phase_curve(
            ax,
            dx,
            dy,
            label=label,
            color=color,
            **kwargs,
        )

    def plot_nullclines(
        self,
        ax: axe.Axes,
        dx_domain: np.ndarray[tuple[int], np.dtype[np.float64]],
        dx_range: np.ndarray[tuple[int], np.dtype[np.float64]],
        dy_domain: np.ndarray[tuple[int], np.dtype[np.float64]],
        dy_range: np.ndarray[tuple[int], np.dtype[np.float64]],
        **kwargs,
    ) -> None:
        self._plot_phase_curve(
            ax,
            dx_domain,
            dx_range,
            label="{} Nullcline".format(self.x_name),
            color=self.colors[0],
            **kwargs,
        )
        self._plot_phase_curve(
            ax,
            dy_domain,
            dy_range,
            label="{} Nullcline".format(self.y_name),
            color=self.colors[1],
            **kwargs,
        )

    def plot_phase_fixed_points(
        self,
        ax: axe.Axes,
        parameters: dict[str, float],
        **kwargs,
    ) -> None:
        col1 = kwargs.pop("color1", self.x_color)
        col2 = kwargs.pop("color2", self.y_color)

        model = ODEModel(parameters=parameters)
        sol1, sol2, sol3 = model.roots()

        new_x_ticks = np.array([sol1])
        new_y_ticks = np.array([sol2])
        ticks = np.arange(0, 101, 20)
        for tick in ticks:
            if np.abs(sol1 - tick) > 5:
                new_x_ticks = np.append(new_x_ticks, tick)

            if np.abs(sol2 - tick) > 5:
                new_y_ticks = np.append(new_y_ticks, tick)

        ax.set_xticks(new_x_ticks)
        ax.set_yticks(new_y_ticks)

        fixed_point_tick_labels = [
            "{}".format(self.x_name) + r"$^*$",
            "{}".format(self.y_name) + r"$^*$",
        ]

        def fixed_point_format(val, pos):
            if val == sol1:
                return fixed_point_tick_labels[0]
            elif val == sol2:
                return fixed_point_tick_labels[1]
            else:
                return int(np.round(val, 3))

        xmax = kwargs.pop("xmax", 100)
        ymax = kwargs.pop("ymax", 100)

        ax.xaxis.set_major_formatter(FuncFormatter(fixed_point_format))
        ax.yaxis.set_major_formatter(FuncFormatter(fixed_point_format))

        self._plot_vline(ax, sol1, ymax, color=col1, **kwargs)
        self._plot_hline(ax, sol2, xmax, color=col2, **kwargs)

    def plot_phase_space(
        self,
        ax: axe.Axes,
        x: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        y: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        u: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        v: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        **kwargs,
    ) -> None:
        scale = self.streamargs.pop("widthscalar", 6)
        speed = np.sqrt(u**2 + v**2)
        lw = speed / speed.max() * scale

        self.streamargs.update(kwargs)
        self.streamargs.update(linewidth=lw)

        ax.streamplot(x, y, u, v, **self.streamargs)

    def _plot_hline(
        self,
        ax: axe.Axes,
        fixed,
        max,
        **kwargs,
    ):
        self.lineargs.update(kwargs)
        ax.hlines([fixed], 0, max, **self.lineargs)

    def _plot_vline(
        self,
        ax: axe.Axes,
        fixed,
        max,
        **kwargs,
    ):
        self.lineargs.update(kwargs)
        ax.vlines([fixed], 0, max, **self.lineargs)
