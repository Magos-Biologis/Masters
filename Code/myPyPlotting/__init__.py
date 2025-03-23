from myPyPlotting import *

import os


from sys import exit
from dataclasses import dataclass


# import pprint


@dataclass
class parameters:
    pass


@dataclass
class base:
    def __init__(
        self,
        figure_name: str,
        subdirectory: str | None = None,
        figure_enviroment: str | None = None
        if type(os.getenv("THESIS_FIGURE_PATH")) is None
        else os.getenv("THESIS_FIGURE_PATH"),
    ) -> None:
        self.figure_name = figure_name
        self.subdirectory = subdirectory
        self.figure_env = figure_enviroment
