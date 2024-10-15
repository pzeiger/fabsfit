#!/usr/bin/env python3

from fabsfit.plot import plot_fit
from fabsfit.cui.download import download_abtem
from fabsfit.cui.fit import fit_absorptive_scattering_factor
import argparse
import fabsfit.data
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        prog='fabsfit',
        description='Fitting absorptive scattering factors',)
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true')  # on/off flag
    
    parser.add_argument(
        '-D', '--debug',
        action='store_true')  # run with debug info
    
    subparsers = parser.add_subparsers(
        required=True,
        help='subcommand help',)
    
    
    #*************************************
    # download command arguments
    #*************************************
    parser_download = subparsers.add_parser(
        'download',
        help='fabsfit download -- Download Parametrization data from abTEM')
    
    parser_download.add_argument(
        'path',
        nargs='?',
        default=Path(fabsfit.data.__file__).parent,
        help='Specify the target directory to download the ' + \
            'data files to.')
        
    parser_download.set_defaults(func=download)
    
    
    #*************************************
    # fit command arguments
    #*************************************
    parser_fit = subparsers.add_parser(
        'fit',
        help='fabsfit fit -- Fit absorptive scattering factor')
    
    parser_fit.add_argument(
        'element',
        help='Element symbol.')
    
    parser_fit.add_argument(
        'Biso',
        type=float,
        #nargs='?',
        help='Isotropic Debye-Waller factor exponent in Å**2.')
    
    parser_fit.add_argument(
        '-k',
        '--Ekin',
        type=float,
        default=100,
        help='Kinetic energy of beam electrons in keV')
    
    parser_fit.add_argument(
        '-e',
        '--elparam',
        default='peng',
        help='Which parametrisation to use for the absorptive scattering factor')
    
    parser_fit.add_argument(
        '-a',
        '--absparam',
        default='peng',
        help='Which parametrisation to use for the absorptive scattering factor')
    
    parser_fit.add_argument(
        '-f',
        '--datafile',
        nargs=1,
        default='peng_high.json',
        help='Path to the elastic scattering factor parameter data.')
    
    parser_fit.add_argument(
        '-m',
        '--maxfev',
        default=10000,
        help='Maximum number of evaluations of fitting function. Try to increase if fit does not converge.')
    
    parser_fit.add_argument(
        '-d',
        '--download',
        action='store_true',
        help='Download parameter file to data dir')
    
    parser_fit.set_defaults(func=fit)
    
    return parser.parse_args()



def download(args):
    def inner():
        return download_abtem(args.path, args.verbose, args.debug)
    return inner()



def fit(args):
    def inner():
        return fit_absorptive_scattering_factor(args.elparam, args.absparam,
                                                args.element, args.Ekin,
                                                args.Biso, 
                                                verbose=args.verbose,
                                                debug=args.debug)
    return inner()



def main():
    """
    """
    args = parse_args()
    
    if args.debug:
        print(args)
    
    return args.func(args)
    


if __name__ == '__main__':
    main()
