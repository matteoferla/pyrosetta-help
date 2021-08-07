"""
Adds to the NGLView Widget the method ``.add_selector``

>>> view = nv.view_rosetta(pose)
>>> view.add_selector(pose, selector)
>>> view
"""

import pyrosetta
import nglview
from io import StringIO

def selector_to_ngl(self: nglview.widget.NGLWidget,
                    pose: pyrosetta.Pose,
                    selector: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector):
    pdb_info = pose.pdb_info()
    # should probably `deal with pdb_info.segmentID(1).strip()`
    ResidueVector = pyrosetta.rosetta.core.select.residue_selector.ResidueVector
    selections = [f'{pdb_info.number(r)}:{pdb_info.chain(r)}' for r in ResidueVector(selector.apply(pose))]
    selection = ' or '.join(selections)
    return selection

def add_selector(self: nglview.widget.NGLWidget,
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
    selection = self.selector_to_ngl(pose, selector)
    self.add_representation(representation_name,
                            colorValue=color,
                            selection=selection,
                            **other)
    self.center(selection)

def add_rosetta(self: nglview.widget.NGLWidget,
                pose: pyrosetta.Pose):
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    fh = StringIO(buffer.str())
    c = self.add_component(fh, ext='pdb')
    # c.update_cartoon(color='bfactor')
    return c


def make_pose_comparison(first_pose: pyrosetta.Pose,
                         second_pose: pyrosetta.Pose,
                         first_color: str = '#F8766D',
                         second_color: str = '#00B4C4'):
    self = nglview.widget.NGLWidget()
    c0 = self.add_rosetta(first_pose)
    c0.update_cartoon(color=first_color, smoothSheet=True)
    c1 = self.add_rosetta(second_pose)
    c1.update_cartoon(color=second_color, smoothSheet=True)

# ======= Monkey patch ==========================================

nglview.widget.NGLWidget.add_selector = selector_to_ngl
nglview.widget.NGLWidget.add_selector = add_selector
nglview.widget.NGLWidget.add_rosetta = add_rosetta
nglview.make_pose_comparison = make_pose_comparison