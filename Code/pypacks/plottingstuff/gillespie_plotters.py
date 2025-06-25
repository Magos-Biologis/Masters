import numpy as np
from matplotlib import axes as axe

from .fixedpoints import GillespieFixed, ODEFixed
from .plotting_class import PlottingData


class gillespiePlotters(PlottingData):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.throw_in_lists()

    def _plot_step(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        results: np.ndarray[tuple[int], np.dtype[np.int_]],
        **kwargs,
    ) -> None:
        self.walkargs.update(kwargs)

        ax.step(
            time,
            results,
            **self.walkargs,
        )

    def _plot_hist(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int], np.dtype[np.float64]],
        **kwargs,
    ) -> None:
        self.histargs.update(kwargs)

        ax.hist(
            results,
            **self.histargs,
        )

    def plot_hist(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int], np.dtype[np.float64]],
        color: str,
        label: str | None = None,
        **kwargs,
    ) -> None:
        # start_cond: str = kwargs.get("label", "Start Condition {}".format(label))

        kwargs.update(facecolor=color)
        if label is not None:
            kwargs.update(label="{}".format(label))

        self._plot_hist(
            ax,
            results=results,
            **kwargs,
        )

    def plot_hists(
        self,
        ax: axe.Axes,
        results: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        **kwargs,
    ) -> None:
        x_label = "Distribution of {} states".format(self.x_name)
        y_label = "Distribution of {} states".format(self.y_name)
        n_label = "Distribution of {} states".format(self.n_name)

        labels = [x_label, y_label, n_label]
        patterns = ["/", "\\", "."]

        for i, result in enumerate(results):
            kwargs.update(
                facecolor=self.colors[i],
                label=labels[i],
            )
            self._plot_hist(
                ax,
                result,
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

        label = "Walk of {}".format(plot_kwargs.get("label", self.x_name))
        color = "{}".format(plot_kwargs.get("color", "r"))

        if plot_starts:
            label += " with start {}".format(xstart)

        kwargs.update(
            color=color,
            label=label,
        )
        self._plot_step(
            ax,
            time,
            steps,
            **kwargs,
        )

    def plot_walks(
        self,
        ax: axe.Axes,
        time: np.ndarray[tuple[int], np.dtype[np.float64]],
        results: np.ndarray[tuple[int, int], np.dtype[np.int_]],
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
        n_label = "Walk of {}".format(self.n_name)
        b_label = "Walk of {}".format(self.b_name)

        labels = [x_label, y_label, n_label, b_label]

        if plot_starts:
            x_label += " with start {}".format(xstart)
            y_label += " with start {}".format(xstart)

        for i, result in enumerate(results):
            kwargs.update(
                color=self.colors[i],
                label=labels[i],
            )
            self._plot_step(
                ax,
                time,
                result,
                **kwargs,
            )

    def _plot_hline(
        self,
        ax: axe.Axes,
        model: str,
        var: str,
        xmax: float,
        parameters: GillespieFixed | ODEFixed,
        **kwargs,
    ):
        # assert type(parameters) is GillespieFixed or , "You fucked the parameters somehow"

        if model in ["5_2", "5_3"]:
            if model == "5_2":
                fixed: dict[str, float] = parameters.Novel_5_2()
            elif model == "5_3":
                fixed: dict[str, float] = parameters.Novel_5_3()
        else:
            if model == "ode_3_2" or model == "ode_3_3":
                fixed: dict[str, float] = parameters.ode_3_2()
            elif model == "ode_5_2" or model == "ode_5_3" or model == "ode_8_3":
                fixed: dict[str, float] = parameters.ode_5_2()

        # print("hey")
        # print(fixed[var])
        self.lineargs.update(kwargs)
        ax.hlines([fixed[var]], 0, xmax, **self.lineargs)

    def _plot_vline(
        self,
        ax: axe.Axes,
        model: str,
        var: str,
        ymax: float,
        parameters: GillespieFixed | ODEFixed,
        **kwargs,
    ):
        if model in ["5_2", "5_3"]:
            if model == "5_2":
                fixed: dict[str, float] = parameters.Novel_5_2()
            elif model == "5_3":
                fixed: dict[str, float] = parameters.Novel_5_3()
        else:
            if model == "ode_3_2" or model == "ode_3_3":
                fixed: dict[str, float] = parameters.ode_3_2()
            elif model == "ode_5_2" or model == "ode_5_3" or model == "ode_8_3":
                fixed: dict[str, float] = parameters.ode_5_2()

        self.lineargs.update(kwargs)
        ax.vlines([fixed[var]], 0, ymax, **self.lineargs)

    def plot_hist_fixed(
        self,
        ax: axe.Axes,
        model: str,
        var: str,
        ymax: float,
        parameters: dict | GillespieFixed | ODEFixed,
        **kwargs,
    ):
        if model in ["5_2", "5_3"]:
            if type(parameters) is dict:
                parameters = GillespieFixed(**parameters)
            assert type(parameters) is GillespieFixed, "You fucked the parameters somehow"

        else:
            if type(parameters) is dict:
                parameters = ODEFixed(**parameters)
            assert type(parameters) is ODEFixed, "You fucked the parameters somehow"

        self._plot_vline(ax, model, var, ymax, parameters, **kwargs)

    def plot_walk_fixed(
        self,
        ax: axe.Axes,
        model: str,
        var: str,
        xmax: float,
        parameters: dict | GillespieFixed | ODEFixed,
        **kwargs,
    ):
        if model in ["5_2", "5_3"]:
            if type(parameters) is dict:
                parameters = GillespieFixed(**parameters)
            assert type(parameters) is GillespieFixed, "You fucked the parameters somehow"

        else:
            if type(parameters) is dict:
                parameters = ODEFixed(**parameters)
            assert type(parameters) is ODEFixed, "You fucked the parameters somehow"

        self._plot_hline(ax, model, var, xmax, parameters, **kwargs)
