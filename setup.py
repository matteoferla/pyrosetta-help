from warnings import warn
from importlib import util
import sys

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# ---------- Non pip modules  ------------------------------------------------------------------------------------------

if not util.find_spec('rdkit'):
    warn('This 3.6+ script **requires** pyrosetta, which has to be downloaded from ' +
         'the Rosetta software site due to licencing.')

from setuptools import setup, find_packages

setup(
    name='pyrosetta_help',
    version='0.1',
    packages=find_packages(),
    url='https://github.com/matteoferla/pyrosetta_help',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description=open('README.md').read()
)
