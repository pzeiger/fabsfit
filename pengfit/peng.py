from parametrization import Parametrization
import numpy as np

class Parametrization(Parametrization):
    
    
    
    def scattering_factor(self, s, p):
        return p[0, 0] * np.exp(-p[1, 0] * s**2.0) \
               + p[0, 1] * np.exp(-p[1, 1] * s**2.0) \
               + p[0, 2] * np.exp(-p[1, 2] * s**2.0) \
               + p[0, 3] * np.exp(-p[1, 3] * s**2.0) \
               + p[0, 4] * np.exp(-p[1, 4] * s**2.0)


 



