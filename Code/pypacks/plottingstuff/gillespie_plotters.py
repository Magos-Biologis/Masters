import numpy as np
from matplotlib import axes as axe

from .plotting_class import plotting_class


class gillespie_plotters(plotting_class):
    def _plot_step(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        results: np.ndarray[tuple[int], np.dtype[np.int_]],
        color: str,
        **kwargs,
    ) -> None:
        ax.step(
            time,
            results,
            color=color,
            **self.walkargs,
            **kwargs,
        )

    def _plot_hist(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int], np.dtype[np.float64]],
        color: str | list[str],
        **kwargs,
    ) -> None:
        ax.hist(
            results,
            color=color,
            **self.histargs,
            **kwargs,
        )

    def plot_hist(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int], np.dtype[np.float64]],
        color: str,
        label: str = "[m,0]",
        **kwargs,
    ) -> None:
        # start_cond: str = kwargs.get("label", "Start Condition {}".format(label))

        self._plot_hist(
            ax,
            results=results,
            color=color,
            label=label,
            **kwargs,
        )

    def plot_hists(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        color: list[str],
        labels: str,
        **kwargs,
    ) -> None:
        self._plot_hist(
            ax,
            results[0, :],
            color=color[0],
            label="Distribution of {} states".format(self.x_name),
            **kwargs,
        )
        self._plot_hist(
            ax,
            results[1, :],
            color=color[1],
            label="Distribution of {} states".format(self.y_name),
            **kwargs,
        )

    def plot_walk(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        steps: np.ndarray[tuple[int], np.dtype[np.int_]],
        plot_kwargs: dict,
        xstart: str = "[m,0]",
        plot_starts: bool = False,
        **kwargs,
    ) -> None:
        y_min: int = plot_kwargs.get("bottom", 0)
        y_max: int = plot_kwargs.get("top", 100)
        axis_step: int = plot_kwargs.get("axis_step", 30)

        label = "Walk of {}".format(plot_kwargs.get("label", r"$x$"))
        color = "{}".format(plot_kwargs.get("color", "r"))

        if plot_starts:
            label += " with start {}".format(xstart)

        self._plot_step(
            ax,
            time,
            steps,
            color=color,
            label=label,
            **kwargs,
        )

        ax.set_xlabel("Time", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)

        ax.set_xlim(left=0)
        ax.set_yticks([y for y in range(y_min, y_max + 1, axis_step)])
        ax.set_ylim(bottom=y_min, top=y_max)

    def plot_steps(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        results: np.ndarray[tuple[int, int], np.dtype[np.int_]],
        color: str,
        plot_kwargs: dict,
        xstart: str = "[m,0]",
        plot_starts: bool = False,
        **kwargs,
    ) -> None:
        y_min: int = plot_kwargs.get("bottom", 0)
        y_max: int = plot_kwargs.get("top", 100)
        axis_step: int = plot_kwargs.get("axis_step", 30)

        x_label = "Walk of {}".format(self.x_name)
        y_label = "Walk of {}".format(self.y_name)

        if plot_starts:
            x_label += " with start {}".format(xstart)
            y_label += " with start {}".format(xstart)

        self._plot_step(
            ax,
            time,
            results[0, :],
            color=color,
            label=x_label,
            **kwargs,
        )
        self._plot_step(
            ax,
            time,
            results[1, :],
            color="g",
            label=y_label,
            **kwargs,
        )

        ax.set_xlabel("Time", fontdict=self.fontargs)
        ax.set_ylabel("Count", fontdict=self.fontargs)

        ax.set_xlim(left=0)
        ax.set_yticks([y for y in range(y_min, y_max + 1, axis_step)])
        ax.set_ylim(bottom=y_min, top=y_max)
