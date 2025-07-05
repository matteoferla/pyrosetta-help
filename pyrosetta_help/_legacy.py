"""
These are a series of legacy monkey patches to ensure that
dependencies work as expected despite changes in the libraries.
"""

import numpy as np

# make sure np.alltrue is available
# this was removed in numpy 1.25, but is still used in PyRosetta
if not hasattr(np, 'alltrue'):
    np.alltrue = np.all

# Older Bio
try:
    from Bio import Align
except ImportError:
    print('Please update your Biopython. Any functions that use Bio.Align will not work until you do so.')
    import Bio
    from types import ModuleType
    Bio.Align = ModuleType('Bio.Align')

import Bio.Align
for name in ('PairwiseAligner', 'PairwiseAlignments', 'Alignment'):
    if not hasattr(Bio.Align, name):
        setattr(Bio.Align, name, None)
        print('Please update your Biopython. Any functions that use Bio.Align will not work until you do so.')