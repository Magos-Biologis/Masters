import numpy as np
from matplotlib import axes as axe

from .plotting_class import plotting_class


class gillespie_plotters(plotting_class):
    def _plot_hist(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int], np.dtype[np.float64]],
        color: str,
        xstart: str = "[m,0]",
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
        xstart: str = "[m,0]",
        **kwargs,
    ) -> None:
        start_cond: str = kwargs.get("label", "Start Condition {}".format(xstart))

        self._plot_hist(
            ax,
            results=results,
            color=color,
            label=start_cond,
            **kwargs,
        )

    def plot_hists(
        self,
        ax1: axe.Axes,
        ax2: axe.Axes,
        results: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        color: str,
        xstart: str = "[m,0]",
    ) -> None:
        self._plot_hist(
            ax1,
            results[0, :],
            color=color,
            label=f"Start Condition {xstart}",
        )
        self._plot_hist(
            ax2,
            results[0, :],
            color=color,
            label=f"Start Condition {xstart}",
        )

    def plot_walk(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        results: np.ndarray[tuple[int, int], np.dtype[np.int_]],
        color: str,
        label: str,
        xstart: str = "[m,0]",
        plot_starts: bool = False,
        **kwargs,
    ) -> None:
        y_min: int = kwargs.get("bottom", 0)
        y_max: int = kwargs.get("top", 100)
        axis_step: int = kwargs.get("axis_step", 30)

        x_label = "Walk of ${}$".format(label)

        if plot_starts:
            x_label += " with start {}".format(xstart)

        ax.step(
            time,
            results[0, :],
            color=color,
            label=x_label,
            **self.walkargs,
        )

        ax.set_xlabel("Time", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)

        ax.set_xlim(left=0)
        ax.set_yticks([y for y in range(y_min, y_max + 1, axis_step)])
        ax.set_ylim(bottom=y_min, top=y_max)

    def plot_walks(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        results: np.ndarray[tuple[int, int], np.dtype[np.int_]],
        color: str,
        xstart: str = "[m,0]",
        plot_starts: bool = False,
        **kwargs,
    ) -> None:
        y_min: int = kwargs.get("bottom", 0)
        y_max: int = kwargs.get("top", 100)
        axis_step: int = kwargs.get("axis_step", 30)

        x_label = "Walk of ${}$".format(self.x_name)
        y_label = "Walk of ${}$".format(self.y_name)

        if plot_starts:
            x_label += " with start {}".format(xstart)
            y_label += " with start {}".format(xstart)

        ax.step(
            time,
            results[0, :],
            color=color,
            label=x_label,
            **self.walkargs,
        )
        ax.step(
            time,
            results[1, :],
            color="g",
            label=y_label,
            **self.walkargs,
        )

        ax.set_xlabel("Time", fontsize=12)
        ax.set_ylabel("Count", fontsize=12)

        ax.set_xlim(left=0)
        ax.set_yticks([y for y in range(y_min, y_max + 1, axis_step)])
        ax.set_ylim(bottom=y_min, top=y_max)
