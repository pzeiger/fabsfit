from fabsfit.parametrization import BaseParametrization
import numpy as np
from importlib import resources as impresources
from . import data

class Parametrization(BaseParametrization):
    
    
    def _get_datafile(self):
        return impresources.files(data) / 'peng_high.json'
        
    
    def elastic_scattering_factor(self, s: np.ndarray, p: np.ndarray = None):
        if not p:
            p = self.p
        return p[0, 0] * np.exp(-p[1, 0] * s**2.0) + \
               p[0, 1] * np.exp(-p[1, 1] * s**2.0) + \
               p[0, 2] * np.exp(-p[1, 2] * s**2.0) + \
               p[0, 3] * np.exp(-p[1, 3] * s**2.0) + \
               p[0, 4] * np.exp(-p[1, 4] * s**2.0)

    
 



