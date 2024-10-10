import numpy as np
from numba import njit
from importlib import resources as impresources
from . import data
from functools import partial


default_datafile = impresources.files(data) / 'peng_high.json'
asfc = 0.023933754  # asymptotic scattering factor constant

@njit
def elastic_scattering_factor(s: np.ndarray, p: np.ndarray):
    return p[0, 0] * np.exp(-p[1, 0] * s**2) + \
           p[0, 1] * np.exp(-p[1, 1] * s**2) + \
           p[0, 2] * np.exp(-p[1, 2] * s**2) + \
           p[0, 3] * np.exp(-p[1, 3] * s**2) + \
           p[0, 4] * np.exp(-p[1, 4] * s**2)


@njit
def elastic_scattering_factor_fit_wrapper(s: np.ndarray, 
                                          a1: np.double, a2: np.double,
                                          a3: np.double, a4: np.double,
                                          a5: np.double, b1: np.double,
                                          b2: np.double, b3: np.double,
                                          b4: np.double, b5: np.double):
    p = np.array([[a1, a2, a3, a4, a5], [b1, b2, b3, b4, b5]])
    return elastic_scattering_factor(s, p)


@njit
def elastic_scattering_factor_extended(s: np.ndarray, smax: np.double, 
                                       alpha: np.double, Z: np.uint8,
                                       p: np.ndarray):
    fel = np.zeros(s.shape)
    mask_le = s <= smax
    mask_gt = s > smax
    fel[mask_le] = elastic_scattering_factor(s[mask_le], p)
    fel[mask_gt] = _asymptotic_scattering_factor(s[mask_gt], alpha, Z)
    return fel


@njit
def _asymptotic_scattering_factor(s: np.ndarray, alpha: np.double, 
                                  Z:np.uint8):
    return asfc * (Z / (s**2 + alpha**2))


@njit
def absorptive_scattering_factor_integrand(theta: np.double, sprim: np.double, 
                                            s: np.ndarray, Biso: np.double,
                                            smax: np.double, alpha: np.double,
                                            Z: np.uint8, p: np.ndarray):
    argp = np.sqrt(0.25 * s**2 + sprim **2 + s * sprim * np.cos(theta))
    argm = np.sqrt(0.25 * s**2 + sprim **2 - s * sprim * np.cos(theta))
    return elastic_scattering_factor_extended(argp, smax, alpha, Z, p) * \
           elastic_scattering_factor_extended(argm, smax, alpha, Z, p) * \
           (1. - np.exp(-2 * Biso * (sprim**2 - .25 * s**2))) * sprim
   

# def absorptive_scattering_factor_integrand_v2(xarray: np.ndarray, 
#                                               s: np.ndarray, Biso: np.double,
#                                               smax: np.double, alpha: np.double,
#                                               Z: np.uint8, p: np.ndarray):
#     sprim = xarray[:, 0]
#     theta = xarray[:, 1]
#     argp = np.sqrt(0.25 * s**2 + sprim **2 + s * sprim * np.cos(theta))
#     argm = np.sqrt(0.25 * s**2 + sprim **2 - s * sprim * np.cos(theta))
#     return elastic_scattering_factor_extended(argp, smax, alpha, Z, p) * \
#            elastic_scattering_factor_extended(argm, smax, alpha, Z, p) * \
#            (1. - np.exp(-2 * Biso * (sprim**2 - .25 * s**2))) * sprim


# def _absorptive_scattering_factor_integrand_pycuba(
#                                               ndim: np.integer,
#                                               xx: np.ndarray,
#                                               ncomp,
#                                               ff,
#                                               userdata,
#                                               sprim_lim, 
#                                               theta_lim,
#                                               s: np.ndarray, Biso: np.double, 
#                                               smax: np.double, alpha: np.double,
#                                               Z: np.integer, p: np.double,
#                                               ):
    
#     sprim = xx[0] * sprim_lim
#     theta = xx[1] * theta_lim
    
#     argp = np.sqrt(0.25 * s**2 + sprim **2 + s * sprim * np.cos(theta))
#     argm = np.sqrt(0.25 * s**2 + sprim **2 - s * sprim * np.cos(theta))
#     result = elastic_scattering_factor_extended(argp, smax, alpha, Z, p) * \
#              elastic_scattering_factor_extended(argm, smax, alpha, Z, p) * \
#              (1. - np.exp(-2 * Biso * (sprim**2 - .25 * s**2))) * sprim  * \
#              theta_lim * sprim_lim
#     ff[0] = result
#     return 0



# # Python code to illustrate 
# # Decorators with parameters in Python 
# def absorptive_scattering_factor_integrand_pycuba(sprim_lim, 
#                                                   theta_lim,
#                                                   s: np.ndarray,
#                                                   Biso: np.double, smax: np.double,
#                                                   alpha: np.double, Z: np.integer,
#                                                   p: np.double):
        
#     return partial(_absorptive_scattering_factor_integrand_pycuba, sprim_lim=sprim_lim, theta_lim=theta_lim, s=s, Biso=Biso, smax=smax, alpha=alpha, Z=Z, p=p)



