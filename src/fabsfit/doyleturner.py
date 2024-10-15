import numpy as np
from numba import njit
from fabsfit.data import datafiles


default_datafile = datafiles['peng_high']


def elastic_scattering_factor_wrapper(p: np.ndarray):
    
    @njit
    def inner(s: np.ndarray):
        return elastic_scattering_factor(s, p)
    
    return inner


@njit
def elastic_scattering_factor(s: np.ndarray, p: np.ndarray):
    return p[0, 0] * np.exp(-p[1, 0] * s**2) + \
           p[0, 1] * np.exp(-p[1, 1] * s**2) + \
           p[0, 2] * np.exp(-p[1, 2] * s**2) + \
           p[0, 3] * np.exp(-p[1, 3] * s**2) + \
           p[0, 4] * np.exp(-p[1, 4] * s**2)


@njit
def elastic_scattering_factor_fit(s: np.ndarray,
                                  a1: np.double, a2: np.double,
                                  a3: np.double, a4: np.double,
                                  a5: np.double, b1: np.double,
                                  b2: np.double, b3: np.double,
                                  b4: np.double, b5: np.double):
    p = np.array([[a1, a2, a3, a4, a5], [b1, b2, b3, b4, b5]])
    return elastic_scattering_factor(s, p)
