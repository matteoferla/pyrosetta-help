__all__ = ['RingSelector',
           'AlteredSelector',
           'UnalteredSelector',
           'OrListSelector']

# will be fixed.
__all__.append('get_bfactor_vector')

import pyrosetta


class RingSelector:
    """
    Select all residues in the "ring"
    based upon 12A from origin.
    There is probably a saner way.

    NB> This is not actually a residue selector. The logical selectors will not accept it.
    """

    def __init__(self, radius=12):
        self.radius = radius

    def apply(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.utility.vector1_bool:
        sele = pyrosetta.rosetta.utility.vector1_bool(pose.total_residue())
        for r in range(1, pose.total_residue() + 1):
            if abs(pose.residue(r).xyz(1).x) < 12:
                sele[r] = 1
        return sele


class AlteredSelector:
    """
    Select residues that were altered in the threading.
    NB> This is not actually a residue selector. The logical selectors will not accept it.
    """

    def __init__(self, threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover):
        self.threader = threader

    def apply(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.utility.vector1_bool:
        sele = pyrosetta.rosetta.utility.vector1_bool(pose.total_residue())
        mapping = self.threader.get_qt_mapping(pose).mapping()
        for r in range(1, pose.total_residue() + 1):
            if mapping[r] != 0:
                sele[r] = 1
        return sele


class UnalteredSelector:
    """
    Select residues that were unaltered in the threading.
    NB> This is not actually a residue selector. The logical selectors will not accept it.
    """

    def __init__(self, threader: pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover):
        self.threader = threader

    def apply(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.utility.vector1_bool:
        sele = pyrosetta.rosetta.utility.vector1_bool(pose.total_residue())
        mapping = self.threader.get_qt_mapping(pose).mapping()
        for r in range(1, pose.total_residue() + 1):
            if mapping[r] == 0:
                sele[r] = 1
        return sele


def OrListSelector(*selectors) -> pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector:
    """
    OrResidueSelector but 2+
    (not a class, but returns a Or
    :param selectors:
    :return:
    """
    sele = pyrosetta.rosetta.core.select.residue_selector.FalseResidueSelector()
    for subsele in selectors:
        sele = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector(subsele, sele)
    return sele


# --------------------------------------------------------------------------------------------
# not a a faux selector

def get_bfactor_vector(pose: pyrosetta.Pose, cutoff: float, above=True) -> pyrosetta.rosetta.utility.vector1_bool:
    """
    Return a selection vector based on b-factors.
    above = get all above. So to select bad b-factors above is ``True``,
    but to select AF2 bad ones. above is ``False``
    """
    pdb_info = pose.pdb_info()
    vector = pyrosetta.rosetta.utility.vector1_bool(pose.total_residue())
    for r in range(1, pose.total_residue() + 1):
        try:
            atom_index = pose.residue(r).atom_index('CA')
        except AttributeError:
            atom_index = 1
        bfactor = pdb_info.bfactor(r, atom_index)
        if above and bfactor >= cutoff:
            vector[r] = True
        elif not above and bfactor <= cutoff:
            vector[r] = True
        else:
            pass  # vector[r] = False
    return vector
