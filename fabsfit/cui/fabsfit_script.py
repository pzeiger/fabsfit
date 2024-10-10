#!/usr/bin/env python3

import sys
import argparse
import fabsfit

def parse_args():
    parser = argparse.ArgumentParser(
                        prog='fabsfit',
                        description='Fit absorptive scattering factor',
                        epilog='Text at the bottom of help')
    
    parser.add_argument('-v', '--verbose',
                        action='store_true')  # on/off flag
    parser.add_argument('Ekin',
                        nargs='?',
                        default=100,
                        help='Kinetic energy of beam electrons in kV')
    parser.add_argument('element',
                        nargs='?',
                        default='H',
                        help='Element abbreviation.')
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
    parser.add_argument('-B', 
                        '--Biso',
                        default=0.0001,
                        help='Isotropic Debye-Waller factor exponent in Å**2.')
    
    return parser.parse_args()


def main():
    """
    """
    
    args = parse_args()
    
    Parametrization()
    
    print(args)
    
    fabsfit.peng.Parametrization(args.element, args.Ekin)


if __name__ == '__main__':
    main()
