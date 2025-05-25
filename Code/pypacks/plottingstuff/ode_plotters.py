import numpy as np
from matplotlib import axes as axe
from matplotlib.ticker import FuncFormatter
from myodestuff import ODEModel

from .plotting_class import PlottingData


class odePlotters(PlottingData):
    def __init__(self, stream_kwargs: dict, *args, **kwargs):
        self.streamargs = stream_kwargs
        super().__init__(*args, **kwargs)
        self.throw_in_lists()

    def plot_curves(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        solutions: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        limit_curves: int = 3,
        **kwargs,
    ):
        labels = kwargs.pop("labels", self.names)

        for i, curve in enumerate(solutions):
            if i == limit_curves:
                break
            kwargs.update(color=self.colors[i])
            kwargs.update(label=labels[i])

            self._plot_curve(
                ax,
                time,
                curve,
                **kwargs,
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

    def plot_fixed(
        self,
        plot_type: str,
        ax: axe.Axes,
        parameters: dict[str, float],
        total_vars: int = 2,
        **kwargs,
    ) -> None:
        model = ODEModel(parameters=parameters)
        fixed_points = model.roots()

        self.xmax = kwargs.pop("xmax", 100)
        self.ymax = kwargs.pop("ymax", 100)

        sol1, sol2, sol3 = fixed_points
        self.sol1 = sol1
        self.sol2 = sol2
        self.sol3 = sol3

        if plot_type == "ode":
            self._plot_ode_fixed_points(ax, fixed_points[:total_vars], **kwargs)
        elif plot_type == "phase":
            self._plot_phase_fixed_points(ax, fixed_points[:total_vars], **kwargs)
        else:
            print("wrong plottype")
            return None

    def _plot_ode_fixed_points(self, ax, fixed_points, axis: str = "y", **kwargs) -> None:
        if axis == "y":
            self._fixed_point_yaxis_setter(ax, fixed_points)

            for i, point in enumerate(fixed_points):
                kwargs.update(color=self.colors[i])
                self._plot_hline(ax, point, **kwargs)
        else:
            self._fixed_point_xaxis_setter(ax, fixed_points)

            for i, point in enumerate(fixed_points):
                kwargs.update(color=self.colors[i])
                self._plot_vline(ax, point, **kwargs)

    def _plot_phase_fixed_points(
        self,
        ax: axe.Axes,
        fixed_points,
        **kwargs,
    ) -> None:
        sol1, sol2, sol3 = fixed_points
        col1 = kwargs.pop("color1", self.x_color)
        col2 = kwargs.pop("color2", self.y_color)

        self._fixed_point_xaxis_setter(ax, sol1)
        self._fixed_point_yaxis_setter(ax, sol2)

        self._plot_vline(ax, sol1, color=col1, **kwargs)
        self._plot_hline(ax, sol2, color=col2, **kwargs)

    def _fixed_point_xaxis_setter(self, ax, fp) -> None:
        if np.array(fp, ndmin=1).shape[0] > 1:
            fixed_ticks = np.array(fp)
        else:
            fixed_ticks = np.array([fp])

        new_x_ticks = fixed_ticks
        ticks = np.arange(0, 101, 20)

        for tick in ticks:
            if all(np.abs(fixed_ticks - tick) > 5):
                new_x_ticks = np.append(new_x_ticks, tick)

        ax.set_xticks(new_x_ticks)

        def fixed_point_format(val, pos):
            if val in fixed_ticks:
                if val == self.sol1:
                    return "{}".format(self.x_name) + r"$^*$"
                elif val == self.sol2:
                    return "{}".format(self.y_name) + r"$^*$"
                elif val == self.sol3:
                    return "{}".format(self.n_name) + r"$^*$"
            else:
                return int(np.round(val, 3))

        ax.xaxis.set_major_formatter(FuncFormatter(fixed_point_format))

    def _fixed_point_yaxis_setter(self, ax, fp) -> None:
        if np.array(fp, ndmin=1).shape[0] > 1:
            fixed_ticks = np.array(fp)
        else:
            fixed_ticks = np.array([fp])

        new_y_ticks = fixed_ticks
        ticks = np.arange(0, 101, 20)

        for tick in ticks:
            if all(np.abs(fixed_ticks - tick) > 5):
                new_y_ticks = np.append(new_y_ticks, tick)

        ax.set_yticks(new_y_ticks)

        def fixed_point_format(val, pos):
            if val in fixed_ticks:
                if val == self.sol1:
                    return "{}".format(self.x_name) + r"$^*$"
                elif val == self.sol2:
                    return "{}".format(self.y_name) + r"$^*$"
                elif val == self.sol3:
                    return "{}".format(self.n_name) + r"$^*$"
            else:
                return int(np.round(val, 3))

        ax.yaxis.set_major_formatter(FuncFormatter(fixed_point_format))

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
        **kwargs,
    ):
        self.lineargs.update(kwargs)
        ax.hlines([fixed], 0, self.xmax, **self.lineargs)

    def _plot_vline(
        self,
        ax: axe.Axes,
        fixed,
        **kwargs,
    ):
        self.lineargs.update(kwargs)
        ax.vlines([fixed], 0, self.ymax, **self.lineargs)

    def _plot_curve(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        curve: np.ndarray[tuple[int], np.dtype[np.float64]],
        **kwargs,
    ):
        self.curvargs.update(kwargs)
        ax.plot(time, curve, **self.curvargs)
