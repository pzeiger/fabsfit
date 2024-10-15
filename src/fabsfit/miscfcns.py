# -*- coding: utf-8 -*-

import inspect


def get_current_function_name():
    return inspect.stack()[1].function


def get_list_func_args(func):
    sig = inspect.signature(func)
    return list(sig.parameters)


def get_number_params(func):
    sig = inspect.signature(func)
    params = sig.parameters
    return len(params) - 1

