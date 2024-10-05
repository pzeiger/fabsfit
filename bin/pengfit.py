#!/usr/bin/env python3

import pengfit
import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
                        prog='ProgramName',
                        description='What the program does',
                        epilog='Text at the bottom of help')
    parser.add_argument('filename', )           # positional argument
    parser.add_argument('-c', '--count')      # option that takes a value
    parser.add_argument('-v', '--verbose',
                        action='store_true')  # on/off flag
    parser.add_argument('-B', '--Biso')
    
    return parser.parse_args() 



def main(args):
    """
    """
    print(args)
    
    print(pengfit.__file__)



if __name__ == '__main__':
    main(parse_args())
