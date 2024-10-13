#!/usr/bin/env python3

import sys
import argparse
import fabsfit
from fabsfit.plot import plot_fit
from fabsfit.parametrization import Parametrization

def parse_args():
    parser = argparse.ArgumentParser(
                        prog='fabsfit',
                        description='Fit absorptive scattering factor',
                        epilog='Text at the bottom of help')
    
    parser.add_argument('-v', '--verbose',
                        action='store_true')  # on/off flag
     
    parser.add_argument('element',
                        nargs=1,
                        help='Element abbreviation.')
    
    parser.add_argument('Biso',
                        nargs=1,
                        help='Isotropic Debye-Waller factor exponent in Å**2.')
    
    parser.add_argument('Ekin',
                        nargs=1,
                        default=100,
                        help='Kinetic energy of beam electrons in kV')
    
    parser.add_argument('-e',
                        '--elparam',
                        default='peng',
                        nargs='1',
                        help='Which parametrisation to use for the absorptive scattering factor')
    
    parser.add_argument('-a',
                        '--absparam',
                        default='peng',
                        nargs='1',
                        help='Which parametrisation to use for the absorptive scattering factor')
    
    parser.add_argument('-d',
                        '--datafile',
                        nargs='1',
                        default='peng_high.json',
                        help='Filename of the elastic scattering factor parameter data.')
    
    parser.add_argument('-m', 
                        '--maxfev',
                        default=10000,
                        help='Maximum number of evaluations of fitting function. Try to increase if fit does not converge.')
    
    return parser.parse_args()


def main():
    """
    """
    
    args = parse_args()
    print(args)
    
    param = Parametrization(args.elparam, args.absparam,
                            args.element, args.Ekin,
                            args.Biso, debug=args.debug)
    
    popt, pcov, infodict, mesg, ier, R2, R2_adj = param.fit()
    
    fig, ax = plot_fit()
    
    fig.savefig('fabsfit_{args.element}.pdf', dpi=300)


if __name__ == '__main__':
    main()
