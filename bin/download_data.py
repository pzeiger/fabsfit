#!/usr/bin/env python3

"""
This script downloads
"""

import requests
import argparse

def download_data():
    files = ('kirkland.json',
             'lobato.json',
             'peng_high.json',
             'peng_ionic.json',
             'peng_low.json',
             'waasmaier_kirfel.json')

    with requests.Session() as s:
        for file in files:
            url = f'https://raw.githubusercontent.com/abTEM/abTEM/main/abtem/parametrizations/data/{file}'
            r = s.get(url)
            print(r.headers.keys())
            print(r.__dict__.keys())
            print(r._content)
            ofile = f'../pengfit/data/{file}'
            with open(ofile,'w') as f:
                f.write(r.text)
                f.write('\n')


if __name__ == __main__:
    download_data()


