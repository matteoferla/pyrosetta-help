# mixin classes are private
from ._init import BlueprinterInit as _Init
from ._common import BlueprinterCommon as _Common
from ._subscripted import BlueprinterSubscripted as _Subscripted
from ._expected import BlueprinterExpected as _Expected
from ._remodel import Remodel as _Remodel
from ._pdb_info import BlueCopier as _Copier
# public
from ._pdb_info import ResInfo


class Blueprinter(_Init, _Subscripted, _Common, _Expected, _Remodel, _Copier):
    """
    Make a blueprint file for Rosetta Remodel and load it into the options (``.set(fn)``).
    The rows property is a list of lists, this is what the operations manipulate.
    Where each row is something like [20, 'G', '.']
    subscript assignment allows the 4th field to be altered. The SS is based on the pose.

    * ``.from_pose(pose)`` initialises it from pose, while the reg. init is from sequence and ss.

    Whereas

    * ``.del_span(10, 29)`` deletes a span
    * ``.wobble_span(100, 105)`` wobbles a span (NATAA)
    * ``.mutate(5, 'W')`` mutates resi i to aa
    * ``.insert(3, 'PIKAA W')`` or ``.insert(3, ['PIKAA W', 'PIKAA A'])`` inserts.

    a more flexible approach is using subscripts and slices:

    >>> blue = Blueprinter.from_pose(pose)
    >>> blue[10:14] = 'NATAA' # preceding loop
    >>> del blue[15:20]
    >>> blue[20:25] = 'NATAA' # following loop
    >>> blue[22] = 'PIKAA W'

    blue.set(bluprint_filename)

    remember that if remodelmover is already initialised to call ``rm.register_options()``

    About the subscript setter, there are three (weird) things to note:

    * A star will be convered into the native amino acid. `PIKAA F*CK` on a glutamine, will result in `PIKAA FECK`
    * The ranges are human style: deletion of resi 10-15 means that 6 residues are missing, including 15. Where normally in python 10:15 is 5, 15 is excluded.
    * The value is a string but there are no safegards or checks, even if there are a limited amount of possibilities and it would be really nice seeing if the non-canonical for `EMPTY NC XAA`
    """
    ResInfo = ResInfo

    # ========= init ===========================
    # from BlueprinterInit
    # adding init and classmethods

    # ========= subscripted ===========================
    # from BlueprinterSubscripted
    # adds subscription methods

    # === common operations =============================================
    # from BlueprinterCommon
    # adds some common tasks

    # === expected  =============================================
    # from BlueprinterExpected
    # adds expected_seq dynamic property

