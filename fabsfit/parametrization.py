import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import dblquad
from fabsfit.constants import hbarc, mec2
from periodictable import elements
from .miscfcns import get_current_function_name
import multiprocessing as mp
import importlib


class Parametrization():

    def __init__(self, elparam: str, absparam: str, element: str,
                 Ekin: float, Biso: float, datafile: str = ''):
        
        self.elparam = importlib.import_module('.'+elparam,
                                               package='fabsfit')
        
        self.absparam = absparam
        
        if datafile:
            self.datafile = datafile
        else:
            self.datafile = self.elparam.default_datafile
        
        self.Ekin = Ekin
        
        # Initialize the element specific data
        self.element = element                        
        self.Z = np.uint8(elements.__getattribute__(element).number)   # atomic number
        
        # Set some more values and constants
        self.Biso = np.double(Biso)         # in Å**2
        
        
        # Get all the preparations done
        self.p = self._load_params_elastic_scattering_factor()
        self.smax = np.double(self._smax())
        self.alpha, df = self._alpha()
    
    
    def fit(self):
        s = np.linspace(.0, 6., 61)
        fabs, fabserr = self.absorptive_scattering_factor(s)
        np.savetxt(f'data_fabs_{self.element}.txt', np.array([s, fabs, fabserr]).T)
        popt, pcov = curve_fit(
            self.elparam.elastic_scattering_factor_fit_wrapper,
            s, fabs*np.exp(-.5 * self.Biso * s**2),
            maxfev=40000,
        )
        print(popt)
        print(pcov)
        return popt, pcov
    
    
    def absorptive_scattering_factor(self, s: np.ndarray):
        
        args = [(
             self.elparam.absorptive_scattering_factor_integrand, 
             .0, np.inf, #2.*self._k0(),
             .0, 2.*np.pi,
             (np.array([x,]), self.Biso, self.smax, self.alpha, self.Z, self.p),
            ) for x in s]
        print(args)
        
#        for arg in args:
#            print(arg)
#            dblquad(*arg)
        
        with mp.Pool(8) as pool:
            tmp = pool.starmap(dblquad, args)
        fabs = np.array([x[0] for x in tmp]) * self._conversion_factor()
        err = np.array([x[1] for x in tmp]) * self._conversion_factor()
        
        return fabs, err
    
    
    def _load_params_elastic_scattering_factor(self):
        with open(self.datafile, 'r') as fh:
            data = json.load(fh)
        return np.array(data[self.element])
    
    
    def _smax(self):
        if np.uint8(self.Z) == 1:
            return np.double(1.5) * np.ones((1,))
        else:
            return np.double(6.0) * np.ones((1,))
    
    
    def _alpha(self):
        tmp1 = self.elparam.elastic_scattering_factor(self.smax, self.p)
        tmp1 /= self.elparam.asfc * self.Z
        tmp2 = 1. - tmp1 * self.smax**2
        alpha = np.sqrt(tmp2 / tmp1)
        df = np.abs(self.elparam.elastic_scattering_factor(self.smax, self.p) - \
                    self.elparam._asymptotic_scattering_factor(self.smax, alpha, self.Z))
        print('Found alpha = % 2.8e with a mismatch of df = % 2.8e' \
              % (alpha, df))
        return alpha, df
    
    
    def _get_datafile(self):
        raise NotImplementedError(
                  'Method %s not implemented.' % get_current_function_name())
    
    
    def _v_over_c(self):
        return np.sqrt(1. - (mec2 / (self.Ekin + mec2))**2)
    
    
    def _k0(self):
        return self._v_over_c() * mec2 * self._lorentz_factor() / hbarc / 2. / np.pi
    
    
    def _conversion_factor(self):
        return 4. * np.pi * hbarc / mec2 / self._v_over_c()
    
    
    def _lorentz_factor(self):
        return 1. / np.sqrt(1. - self._v_over_c()**2)
    
    