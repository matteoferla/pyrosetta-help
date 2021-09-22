from warnings import warn
from importlib import util
import sys

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# ---------- Non pip modules  ------------------------------------------------------------------------------------------

if not util.find_spec('pyrosetta'):
    warn('This 3.6+ script **requires** pyrosetta, which has to be downloaded from ' +
         'the Rosetta software site due to licencing.')

# ---------- Setup  ------------------------------------------------------------------------------------------


from setuptools import setup, find_packages

import os
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    __doc__ = f.read()

description = 'A variety of functions to make working with Pyrosetta easier.'

setup(
    name='pyrosetta_help',
    version='0.3.8',
    packages=find_packages(),
    install_requires=['rdkit-to-params'],
    url='https://github.com/matteoferla/pyrosetta_help',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description=description,
    long_description=__doc__,
    long_description_content_type='text/markdown'
)
