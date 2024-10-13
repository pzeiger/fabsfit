import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import dblquad
from fabsfit.constants import hbarc, mec2
from periodictable import elements
from fabsfit.miscfcns import get_current_function_name
import multiprocessing as mp
import importlib
from numba import njit
from inspect import signature


# asymptotic scattering factor constant
_asfc = 0.023933754


def _elastic_scattering_factor_extended(elastic_scattering_factor,
                                        asymptotic_scattering_factor,
                                        smax: np.double):
    
    f_esf = njit(elastic_scattering_factor)
    f_asf = njit(asymptotic_scattering_factor)
    
    @njit
    def inner(s: np.ndarray):
        fel = np.zeros(s.shape)
        mask_le = s <= smax
        mask_gt = s > smax
        fel[mask_le] = f_esf(s[mask_le])
        fel[mask_gt] = f_asf(s[mask_gt])
        return fel
        
    return inner


def _absorptive_scattering_factor_integrand(elastic_scattering_factor_extended,
                                            Biso: np.double):
    f = njit(elastic_scattering_factor_extended)
    @njit
    def inner(theta: np.double, sprim: np.double, s: np.ndarray):
        argp = np.sqrt(0.25 * s**2 + sprim **2 + s * sprim * np.cos(theta))
        argm = np.sqrt(0.25 * s**2 + sprim **2 - s * sprim * np.cos(theta))
        return  sprim * f(argp) * f(argm) * \
                (1. - np.exp(-2 * Biso * (sprim**2 - .25 * s**2)))
    return inner


def _asymptotic_scattering_factor(alpha: np.double, Z: np.uint8):
    @njit
    def inner(s: np.ndarray):
        return _asfc * (Z / (s**2 + alpha**2))
    return inner


class Parametrization():

    def __init__(self, elparam: str, absparam: str, element: str,
                 Ekin: float, Biso: float, datafile: str = '', debug=False):
        
        self.debug = debug
        
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
        self._set_elastic_scattering_factor()
        self.smax = np.double(self._smax())
        
        # Set all functions
        self._set_asymptotic_scattering_factor()
        if self.debug:
            self.check_asymptotic_scattering_factor()
        
        self._set_elastic_scattering_factor_fit()
        self._set_elastic_scattering_factor_extended()
        self._set_absorptive_scattering_factor_integrand()
    
    
    def fit(self):
        s = np.linspace(.0, 6., 61)
        fabs, fabserr = self.absorptive_scattering_factor(s)
        
        # 
        fabs_dwf = fabs * np.exp(-.5 * self.Biso * s**2)
        
        np.savetxt(f'data_fabs_{self.element}.txt', 
                   np.array([s, fabs, fabs_dwf, fabserr]).T)
        
        sig = signature(self.elastic_scattering_factor_fit)
        params = sig.parameters
        nparams = len(params) - 1
        
        popt, pcov, infodict, mesg, ier = curve_fit(
            self.elastic_scattering_factor_fit,
            s, fabs_dwf,
            p0=[.5 for i in range(nparams)],
            maxfev=20000,
            full_output=True,
        )
        print(popt)
        print(pcov)
        
        self._set_fitted_absorptive_scattering_factor(popt.reshape((2,-1)))
        
        fabs_model = self.fitted_absorptive_scattering_factor(s)
        RSS = np.sum((fabs_dwf - fabs_model)**2)
        TSS = np.sum((fabs_dwf - fabs_dwf.mean())**2)
        R2 = 1 - RSS/TSS
        R2_adj = 1 - (RSS * (s.size-1)) / (TSS * (s.size-nparams))
        
        return popt, pcov, infodict, mesg, ier, R2, R2_adj
        
    
    
    def absorptive_scattering_factor(self, s: np.ndarray):
        args = [(
             self.absorptive_scattering_factor_integrand, 
             .0, 2.*np.inf, #self._k0(),
             .0, 2.*np.pi,
             (np.array([x,]),),
             ) for x in s]
        
        # for arg in args:
        #     print(arg)
        #     dblquad(*arg)
        
        if self.debug:
            print('******************************************************')
            print('Starting integration of absorptive scattering factor')
            print(args)
        
        # Integrate absorptive scattering factor
        with mp.Pool(8) as pool:
            tmp = pool.starmap(dblquad, args)
            
        if self.debug:
            print('******************************************************')
        
        fabs = np.array([x[0] for x in tmp]) * self._conversion_factor()
        err = np.array([x[1] for x in tmp]) * self._conversion_factor()
        
        return fabs, err
    
    
    def _set_elastic_scattering_factor(self):
        self.elastic_scattering_factor = \
            self.elparam.elastic_scattering_factor(self.p)
    
    
    def _set_asymptotic_scattering_factor(self):
        self.alpha = self._alpha()
        self.asymptotic_scattering_factor = \
            _asymptotic_scattering_factor(self.alpha, self.Z)
    
    
    def _set_elastic_scattering_factor_fit(self):
        self.elastic_scattering_factor_fit = \
            self.elparam.elastic_scattering_factor_fit
    
    
    def _set_elastic_scattering_factor_extended(self):
        self.elastic_scattering_factor_extended = \
            _elastic_scattering_factor_extended(
                self.elastic_scattering_factor.py_func,
                self.asymptotic_scattering_factor.py_func,
                self.smax)
    
    
    def _set_absorptive_scattering_factor_integrand(self):
        self.absorptive_scattering_factor_integrand = \
            _absorptive_scattering_factor_integrand(
                self.elastic_scattering_factor_extended.py_func,
                self.Biso)
    
    
    def _set_fitted_absorptive_scattering_factor(self, p):
        self.fitted_absorptive_scattering_factor = \
            self.elparam.elastic_scattering_factor(p=p)
    
    
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
        tmp1 = self.elastic_scattering_factor(self.smax)
        tmp1 /= _asfc * self.Z
        tmp2 = 1. - tmp1 * self.smax**2
        alpha = np.sqrt(tmp2 / tmp1)
        return alpha
    
    
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
    
    
    def check_asymptotic_scattering_factor(self):
        df = np.abs(self.elastic_scattering_factor(self.smax) - \
                    self.asymptotic_scattering_factor(self.smax))
        print('Found alpha = % 2.8e with a mismatch of df = % 2.8e' \
              % (self.alpha, df))
        return df
    
    
    