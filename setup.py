#!/usr/bin/env python3

import os

from setuptools import setup, find_packages


here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "README.md")) as fd:
    long_description = fd.read()

setup(name='fabsfit',
      version='0.1',
      description='Fitting absorptive scattering factors',
      license='GNU General Public License v3.0',
      author='Paul Zeiger',
      author_email='paul.zeiger@physics.uu.se',
      url='https://github.com/pzeiger/fabsfit',
      scripts=['bin/fabsfit', 'bin/'],
      package_dir={"":"src"},
      packages=find_packages(where='src',),
#      python_requires='>=3.10',
      install_requires=open('requirements.txt').read().splitlines(),
      package_data={'': ['*.json']},
      include_package_data=False,
      exclude_package_data={'': ['*.json']},
      )


