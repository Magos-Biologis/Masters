import os

import matplotlib as mpl
import numpy as np
from matplotlib import axes as axe
from matplotlib import pyplot as plt


class plotting_functions:
    def __init__(
        self,
        font_kwargs: dict,
        hist_kwargs: dict,
        line_kwargs: dict,
        walk_kwargs: dict,
        x_name: str = "x",
        y_name: str = "y",
    ) -> None:
        self.histargs: dict = hist_kwargs
        self.lineargs: dict = line_kwargs
        self.fontargs: dict = font_kwargs
        self.walkargs: dict = walk_kwargs

        self.x_name: str = x_name
        self.y_name: str = y_name

    def plot_hist(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int], np.dtype[np.float64]],
        color: str,
        xstart: str = "[m,0]",
    ) -> None:
        ax.hist(
            results,
            **self.histargs,
            label=f"Start Condition {xstart}",
            color=color,
        )
        # ax.hist(
        #     results[1, :],
        #     **self.histargs,
        #     label=f"Start Condition {xstart}",
        #     color=color,
        # )

    def plot_hists(
        self,
        ax1: axe.Axes,
        ax2: axe.Axes,
        results: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        color: str,
        xstart: str = "[m,0]",
    ) -> None:
        ax1.hist(
            results[0, :],
            **self.histargs,
            label=f"Start Condition {xstart}",
            color=color,
        )
        ax2.hist(
            results[1, :],
            **self.histargs,
            label=f"Start Condition {xstart}",
            color=color,
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
