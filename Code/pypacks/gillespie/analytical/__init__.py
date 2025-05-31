import os
from dataclasses import dataclass

import numpy as np
from numba import njit

from .analytical_parameter_class import *
from .overcomplicated_two_chemical_system import *
from .simple_two_chemical_system import *
