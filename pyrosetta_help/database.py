# database/chemical
__all__ = ['database', 'patches', 'residue_types']
"""
This submodule takes a second to parse contents of the folders, therefore it is not imported into the main module.
It is shorthard way for interacting with the Rosetta Database within PyRosetta, but pre-parses all the folders
in order to make tab completion work, see a ``DBEntry`` instance for more.

.. code-block:: python
    from pyrosetta_help.database import database, residue_types, patches 
"""


import pkg_resources
import os
import pkg_resources
from typing import List


class DBEntry:
    """
    A curious shorthard way for interacting with the Rosetta Database within PyRosetta.

    Namely the ``database`` object (type: ``DBEntry``) is the root and pressing tab in Jupyter
    on this object will suggest a few completions, including the files within it it if a directory
    and:

    :ivar text: the content of the file (if not a directory)
    :ivar filenames: is the filenames with that directory (if a directory)
    :ivar is_dir: is boolean of is a directory
    :ivar parts: is the parts of the path within PyRosetta package (List[str])
    :ivar path: is the above but in path format (str)
    :ivar fullpath: is the absolute filepath (i.e. inclusive of ``$CONDA_PREFIX`` if conda based)

    .. code-block:: python
        import pyrosetta_help as ph
        from rdkit_to_params import Params
        p = Params.load(ph.residue_types.)
    """

    def __init__(self, *parts: str):
        self.parts = parts
        self.path = os.path.join(*self.parts)
        if self.is_dir:
            for fn in self.filenames:
                setattr(self, fn, self.__class__(*self.parts, fn))

    @property
    def is_dir(self):
        return pkg_resources.resource_isdir(pyrosetta.__name__, self.path)

    @property
    def fullpath(self):
        return pkg_resources.resource_filename(pyrosetta.__name__, self.path)

    @property
    def text(self) -> str:
        if self.is_dir:
            raise TypeError(f'{self.parts[-1]} is a directory')
        return pkg_resources.resource_string(pyrosetta.__name__, self.path)

    @property
    def filenames(self) -> List[str]:
        if not self.is_dir:
            raise TypeError(f'{self.parts[-1]} is not a directory')
        return pkg_resources.resource_listdir(pyrosetta.__name__, self.path)


database = DBEntry('database')
patches: DBEntry = database.chemical.residue_type_sets.fa_standard.patches
residue_types: DBEntry = database.chemical.residue_type_sets.fa_standard.residue_types
