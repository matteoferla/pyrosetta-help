"""
Adds to the NGLView Widget the method ``.add_selector``

>>> view = nv.view_rosetta(pose)
>>> view.add_selector(pose, selector)
>>> view
"""

import pyrosetta
import nglview as nv

def add_selector(self: nv.widget.NGLWidget,
                  pose: pyrosetta.Pose,
                  selector: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector,
                  representation_name: str = 'hyperball',
                  color: str = 'grey',
                  **other):
    """
    Add a representation of type ``representation_name`` (def. 'hyperball') and center
    based upon the Pyrosetta residue selector

    :param self:
    :param pose:
    :param selector:
    :param representation_name:
    :param color:
    :param other:
    :return:
    """
    pdb_info = pose.pdb_info()
    # should probably `deal with pdb_info.segmentID(1).strip()`
    ResidueVector = pyrosetta.rosetta.core.select.residue_selector.ResidueVector
    selections = [f'{pdb_info.number(r)}:{pdb_info.chain(r)}' for r in ResidueVector(selector.apply(pose))]
    selection = ' or '.join(selections)
    self.add_representation(representation_name,
                            colorValue=color,
                            selection=selection,
                            **other)
    self.center(selection)


nv.widget.NGLWidget.add_selector = add_selector
