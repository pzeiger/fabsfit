# -*- coding: utf-8 -*-

import inspect
from numba import njit
import numpy as np


def get_current_function_name():
    return inspect.stack()[1].function


