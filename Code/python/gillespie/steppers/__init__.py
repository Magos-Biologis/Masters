from importlib import import_module
import numpy as np

import gillespie.propensities as dgp


class ssa_stepper:
    step_module: str = "gillespie.steppers"
    steps: int = 10_000

    vj: dict[str, list[np.ndarray[tuple[int], np.dtype[np.int_]]]] = dgp.transitions

    def __init__(
        self,
        variation: str,
        x0: np.ndarray[tuple[int], np.dtype[np.int_]],
        k: np.ndarray[tuple[int, int], np.dtype[np.float64]],
        # v: list[np.ndarray[tuple[int], np.dtype[np.int_]]],
    ) -> None:
        self.sub_module = self.step_module + ".step_function_" + variation
        self.state_changes = self.vj["vj_" + variation]

        self.x0 = x0
        self.k = k

    def step_function(
        self,
        steps: int | np.int_ = steps,
    ) -> tuple[
        np.ndarray[tuple[int], np.dtype[np.float64]],
        np.ndarray[tuple[int, int], np.dtype[np.int_]],
    ]:
        modu = import_module(self.sub_module)
        step_function = getattr(modu, "main")

        return step_function(steps, self.x0, self.state_changes, self.k)
