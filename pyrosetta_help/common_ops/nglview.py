"""
Adds to the NGLView Widget the method ``.add_selector``

>>> view = nv.view_rosetta(pose)
>>> view.add_selector(pose, selector)
>>> view
"""

import pyrosetta
import nglview
from io import StringIO
from typing import *
from .constraints import get_NGL_selection_from_AtomID
from .utils import get_pdbstr


def selector_to_ngl(self: nglview.widget.NGLWidget,
                    pose: pyrosetta.Pose,
                    selector: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector):
    """
    Given a pose and a selector return the selection string for NGL.

    :param self:
    :param pose:
    :param selector:
    :return:
    """
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
                pose: pyrosetta.Pose, color: Optional[str] = None):
    """
    The module method ``show_rosetta`` creates an NGLWidget
    This is a monkeypatched bound method.

    :param self:
    :param pose:
    :param bfactor:
    :return:
    """
    fh = StringIO(get_pdbstr(pose))
    c = self.add_component(fh, ext='pdb')
    if color:
        c.update_cartoon(color='color', smoothSheet=True)
    return c


def make_pose_comparison(self: nglview.widget.NGLWidget,
                         first_pose: pyrosetta.Pose,
                         second_pose: pyrosetta.Pose,
                         first_color: str = '#00B4C4',
                         second_color: str = '#F8766D'):
    """
    Adds two objects, the first colored by default in #00B4C4, which is turquoise,
     while the second #F8766D, which is salmon.
     The poses are assumed aligned.

    :param self:
    :param first_pose:
    :param second_pose:
    :param first_color:
    :param second_color:
    :return:
    """
    c0 = self.add_rosetta(first_pose)
    c0.update_cartoon(color=first_color, smoothSheet=True)
    c1 = self.add_rosetta(second_pose)
    c1.update_cartoon(color=second_color, smoothSheet=True)


def add_constraints(self: nglview.widget.NGLWidget,
                    pose: pyrosetta.Pose,
                    component=0,
                    color="skyblue"):
    atom_pairs = [[get_NGL_selection_from_AtomID(pose, con.atom1(), named=False),
                   get_NGL_selection_from_AtomID(pose, con.atom2(), named=False)
                   ]
                  for con in pose.constraint_set().get_all_constraints()
                  ]
    atom_pairs = [[a1, a2] for a1, a2 in atom_pairs if a1 != a2]
    component = getattr(self, f'component_{component}')
    component.add_representation("distance",
                                 atomPair=atom_pairs,
                                 color=color
                                 )
    return atom_pairs


# ======= Monkey patch ==========================================

nglview.widget.NGLWidget.selector_to_ngl = selector_to_ngl
nglview.widget.NGLWidget.add_selector = add_selector
nglview.widget.NGLWidget.add_rosetta = add_rosetta
nglview.widget.NGLWidget.add_constraints = add_constraints
nglview.widget.NGLWidget.make_pose_comparison = make_pose_comparison
