from fabsfit.parametrization import Parametrization

def fit_absorptive_scattering_factor(elparam, absparam, element, Ekin,
                                     Biso, verbose=False, debug=False):
    
    
    
    param = Parametrization(elparam, absparam, element, Ekin, Biso, debug=debug)
     
    popt, pcov, infodict, mesg, ier, R2, R2_adj = param.fit()
     
    fig, ax = plot_fit()
     
    fig.savefig('fabsfit_{args.element}.pdf', dpi=300)