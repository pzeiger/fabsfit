import json
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import nquad
from fabsfit.constants import hbarc, mec2
from periodictable import elements
import fabsfit.miscfcns as miscfcns
import multiprocessing as mp
import importlib
from numba import njit
import io
from contextlib import redirect_stdout, redirect_stderr
import sys


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


def _integrate_absorptive_scattering_factor(integrand, s: np.ndarray, k0):
    
    out = io.StringIO()
    err = io.StringIO()
    
    with redirect_stdout(out), redirect_stderr(err):
        print(f'Value of s: {s}')
        ab = [(.0, 2.*np.pi), (.0, np.inf)]
        fabs, fabs_err = nquad(integrand, ab, args=(s,),
                               opts={'limit': 1000})
    
    return fabs, fabs_err, out.getvalue(), err.getvalue()



class Parametrization():

    def __init__(self, elparam: str, absparam: str, element: str,
                 Ekin: float, Biso: float, datafile: str = '', 
                 verbose=False, debug=False):
        
        self.debug = debug
        self.verbose = verbose
        
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
        
        self._set_scattering_factor_absfit()
        self._set_elastic_scattering_factor_extended()
        self._set_absorptive_scattering_factor_integrand()
        
        self.smin_fit = .0
        self.smax_fit = 6.
        self.ns = 121
    
    
    def fit(self):
        s = np.linspace(self.smin_fit, self.smax_fit, self.ns)
        
        fabs, fabs_err, comments = self.absorptive_scattering_factor(s)
        
        # multiply by DWF
        fabs_dwf = fabs * np.exp(-.5 * self.Biso * s**2)
        
        if self.debug:
            print(s.shape)
            print(fabs.shape)
            print(fabs_err.shape)
            print(comments.shape)
        
        self.fitdata = np.array([s, fabs, fabs_err, fabs_dwf, comments]).T
        
        header = 'columns: s    fabs    fabs_err    ' + \
            'fabs*exp(-.5*Biso*s**2)    comment'
        
        np.savetxt(f'data_fabs_{self.save_str()}.txt',
                   self.fitdata,
                   fmt='      %.3f %.10e %.10e %.10e %s',
                   header=header)
        
        if self.debug or self.verbose:
            print('******************************************************')
            print('Starting fit of absorptive scattering factor ...')
        
        bounds = np.zeros((2, self.nparams_absfit))
        bounds[0,:] = np.inf
        bounds[1,:5] = -np.inf
        
        if self.debug:
            print(bounds)
        
        popt, pcov, infodict, mesg, ier = curve_fit(
            self._scattering_factor_absfit,
            s, fabs_dwf,
            p0=[np.random.uniform() for i in range(self.nparams_absfit)],
            maxfev=20000,
            full_output=True,
        )
        
        if self.debug or self.verbose:
            print('Done!')
            print('******************************************************')
        
        if self.debug:
            print(popt)
            print(pcov)
        
        self._set_fitted_absorptive_scattering_factor(popt.reshape((2,-1)))
        
        fabs_model = self.fitted_absorptive_scattering_factor(s)
        RSS = np.sum((fabs_dwf - fabs_model)**2)
        TSS = np.sum((fabs_dwf - fabs_dwf.mean())**2)
        R2 = 1 - RSS/TSS
        R2_adj = 1 - (RSS * (s.size-1)) / (TSS * (s.size-self.nparams_absfit))
        
        return popt, pcov, infodict, mesg, ier, R2, R2_adj
        

    def absorptive_scattering_factor(self, s: np.ndarray):
        args = [(self.absorptive_scattering_factor_integrand, 
                 np.array([x,]), self._k0()) for x in s]
                
        # for arg in args:
        #     print(arg)
        #     dblquad(*arg)
        
        if self.debug or self.verbose:
            print('******************************************************')
            print('Starting integration of absorptive scattering factor ...')
        
        if self.debug:
            print(args)
        
        # Integrate absorptive scattering factor
        with mp.Pool(8) as pool:
            tmp = pool.starmap(_integrate_absorptive_scattering_factor, args)
            
        if self.debug or self.verbose:
            print('Done!')
            print('******************************************************')
        
        fabs = np.array([x[0] for x in tmp]) * self._conversion_factor()
        fabs_err = np.array([x[1] for x in tmp]) * self._conversion_factor()
        stdout = [x[2] for x in tmp]
        stderr = [x[3] for x in tmp]
        
        comments = self._process_integration_warnings(stdout, stderr)
        
        return fabs, fabs_err, comments
    
    
    def check_asymptotic_scattering_factor(self):
        df = np.abs(self.elastic_scattering_factor(self.smax) - \
                    self.asymptotic_scattering_factor(self.smax))
        print('Found alpha = % 2.8e with a mismatch of df = % 2.8e' \
              % (self.alpha, df))
        return df
    
    
    def save_str(self):
        return f'{self.element}_Biso{self.Biso}_Ekin{self.Ekin}'
    
    
    def _process_integration_warnings(self, stdout, stderr):
        comments = []
        
        for o, e in zip(stdout, stderr):
            if o[:11] != 'Value of s:':
                print(o)
            if e == '':
                comments.append('')
            elif 'IntegrationWarning' in e and 'roundoff error' in e:
                comments.append('roundoff error occurred, ' + \
                               'fab_err may be underestimated')
            else:
                print(e)
                sys.exit()
        return np.array(comments, dtype=object)
    
    
    def _set_elastic_scattering_factor(self):
        self.elastic_scattering_factor = \
            self.elparam.elastic_scattering_factor_wrapper(self.p)
    
    
    def _set_asymptotic_scattering_factor(self):
        self.alpha = self._alpha()
        self.asymptotic_scattering_factor = \
            _asymptotic_scattering_factor(self.alpha, self.Z)
    
    
    def _set_scattering_factor_absfit(self):
        self._scattering_factor_absfit = \
            self.elparam.elastic_scattering_factor_fit
        self.nparams_absfit = miscfcns.get_number_params(
            self._scattering_factor_absfit)
        x = miscfcns.get_list_func_args(self._scattering_factor_absfit)
        self.list_params_absfit = x[1:]
    
    
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
            self.elparam.elastic_scattering_factor_wrapper(p=p)
    
    
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
                  f'Method {miscfcns.get_current_function_name()} not implemented.')
    
    
    def _v_over_c(self):
        return np.sqrt(1. - (mec2 / (self.Ekin + mec2))**2)
    
    
    def _k0(self):
        return self._v_over_c() * mec2 * self._lorentz_factor() / hbarc / 2. / np.pi
    
    
    def _conversion_factor(self):
        return 4. * np.pi * hbarc / mec2 / self._v_over_c()
    
    
    def _lorentz_factor(self):
        return 1. / np.sqrt(1. - self._v_over_c()**2)
    