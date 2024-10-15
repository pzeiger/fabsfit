from fabsfit.parametrization import Parametrization
from fabsfit.cui.plot import plot_fit
import numpy as np
import pickle


def fit_absorptive_scattering_factor(elparam, absparam, element, Ekin,
                                     Biso, showplot=False,
                                     verbose=False, debug=False):
    
    
    
    param = Parametrization(elparam, absparam, element, Ekin, Biso, debug=debug)
     
    popt, pcov, infodict, mesg, ier, R2, R2_adj = param.fit()
    
    fitinfo = {
        'popt':      popt,
        'pcov':      pcov,
        'infodict':  infodict,
        'mesg':      mesg,
        'ier':       ier,
        'R2':        R2,
        'R2_adj':    R2_adj,
    }
    
    with open(f'fabsfit_fitinfo_{param.save_str()}.pkl', 'wb') as fh:
        pickle.dump(fitinfo, fh)
    
    header = ' '.join(param.list_params_absfit)
    np.savetxt(f'fabsfit_params_{param.save_str()}.txt',
               popt.reshape(1,param.nparams_absfit), 
               header=header)
    
    
    fig, ax = plot_fit(param.fitdata[:,0], param.fitdata[:,3], 
                       param.fitted_absorptive_scattering_factor)
    
    fig.savefig(f'fabsfit_{param.save_str()}.pdf', dpi=300)