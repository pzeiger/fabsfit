import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad
from fabsfit.constants import hbar, me, c
from periodictable import elements
from .miscfcns import get_current_function_name
import multiprocessing as mp


class BaseParametrization():

    def __init__(self, element: str, Ekin: float, Biso: float, datafile: str = ''):
        
        self.Ekin = Ekin
        
        # Initialize the element specific data
        self.element = element                        
        self.Z = elements.__getattribute__(element).number   # atomic number
        
        # Set some more values and constants
        self.Biso = Biso         # in Å**2
        self.asfc = 0.023933754  # asymptotic scattering factor constant
        
        if datafile:
            self.datafile = datafile
        else:
            self.datafile = self._get_datafile()
        
        # Get all the preparations done
        self.p = self._load_params_elastic_scattering_factor()
        self.smax = self._determine_smax()
        self.alpha, df = self._determine_alpha()
        
    
    def fit(self):
        s = np.linspace(.0, 10., 251)
        fabs, fabserr = self.absorptive_scattering_factor(s)
        np.savetxt(f'data_fabs_{self.element}.txt', np.array([fabs, fabserr]).T)
        popt, pcov = curve_fit(self.elastic_scattering_factor, s, fabs)
        print(popt)
        print(pcov)
        return popt, pcov
    
    def _load_params_elastic_scattering_factor(self):
        with open(self.datafile, 'r') as fh:
            data = json.load(fh)
        return np.array(data[self.element])
    
    
    def elastic_scattering_factor_extended(self, s: np.ndarray):
        fel = np.zeros(s.shape)
        mask_le = s <= self.smax
        mask_gt = s > self.smax
        fel[mask_le] = self.elastic_scattering_factor(s[mask_le])
        fel[mask_gt] = self.asymptotic_scattering_factor(s[mask_gt])
        return fel
    
    
    def elastic_scattering_factor(self, s: np.ndarray):
        raise NotImplementedError(
                  'Method %s not implemented.' % get_current_function_name())
    
    
    def absorptive_scattering_factor(self, s: np.ndarray):
        
        args = [(self._fabs_integrand, -np.inf, np.inf, (np.array([x,]),)) for x in s]
        print(args)
        with mp.Pool(8) as pool:
            tmp = pool.starmap(quad, args)
            print(tmp)
        fabs = np.ndarray([x[0] for x in tmp])
        fabserr = np.ndarray([x[1] for x in tmp])
        return fabs * self._conversion_factor(), fabserr * self._conversion_factor()
    
    
    def _fabs_integrand(self, sprim, s: np.ndarray):
        return self.elastic_scattering_factor_extended(np.abs(s/2 + sprim)) * \
               self.elastic_scattering_factor_extended(np.abs(s/2 - sprim)) * \
               (1 - np.exp(-2 * self.Biso * (sprim**2 - s**2/4)))
    
    
    def asymptotic_scattering_factor(self, s: np.ndarray, alpha=None):
        if not alpha:
            alpha = self.alpha
        return self.asfc * (self.Z / (s**2 + alpha**2))
    
    
    def _determine_smax(self):
        if self.Z == 1:
            return 1.5 * np.ones((1,))
        else:
            return 6.0 * np.ones((1,))
    
    
    def _determine_alpha(self):
        tmp1 = self.elastic_scattering_factor(self.smax) / self.asfc / self.Z
        tmp2 = 1 - tmp1 * self.smax**2
        alpha = np.sqrt(tmp2 / tmp1)
        df = np.abs(self.elastic_scattering_factor(self.smax) - \
                    self.asymptotic_scattering_factor(self.smax, alpha))
        print('Found alpha = % 2.8e with a mismatch of df = % 2.8e' \
              % (alpha, df))
        return alpha, df
    
    
    def _get_datafile(self):
        raise NotImplementedError(
                  'Method %s not implemented.' % get_current_function_name())
        
        
    def _c_over_v(self):
        return np.sqrt(1 - me / (self.Ekin + me))
        
    
    def _conversion_factor(self):
        return 4 * np.pi * hbar / me / self._c_over_v()
        
    