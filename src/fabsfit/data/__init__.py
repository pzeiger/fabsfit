import importlib.resources as resources
from pathlib import Path


"""
Here we make all the .json files found in the data dir available.
"""

files = resources.files(__name__)
p = files.glob('*.json')

datafiles = {}
for file in p:
    fstr = Path(file).name
    datafiles[fstr.rstrip('.json')] = file
