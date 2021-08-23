__all__ = ['superimpose_by_pLDDT']

# ============= pLDDT superimposition ======================================
# modified from https://gist.github.com/asford/c2404c8b045700f016fda8893325c807
import difflib
import pyrosetta
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID
from pyrosetta.rosetta.core.id import AtomID


def paired_residue_inds(a: pyrosetta.Pose, b: pyrosetta.Pose):
    """Get paired indicies of common residues in two structures."""
    aseq = a.sequence()
    bseq = b.sequence()
    astart, bstart, align_len = difflib.SequenceMatcher(a=aseq, b=bseq).find_longest_match(0, len(aseq), 0, len(bseq))
    return [(astart + i, bstart + i) for i in range(align_len)]


def superimpose_by_pLDDT(pose: pyrosetta.Pose,
                         original: pyrosetta.Pose,
                         cutoff=70,
                         pose_range=None) -> map_core_id_AtomID_core_id_AtomID:
    """
    Superimpose two poses, based on residues with pLDDT above a given threshold.

    :param pose:
    :param original:
    :param cutoff: %
    :param pose_range: optional argument to subset (start:int, end:int)
    :return:
    """
    ca_map = map_core_id_AtomID_core_id_AtomID()
    for w, a in paired_residue_inds(pose, original):
        if original.pdb_info().bfactor(a + 1, 1) <= cutoff:
            continue
        if pose_range is not None and (w < pose_range[0] or w > pose_range[1]):
            continue
        ca_map[AtomID(pose.residue(w + 1).atom_index("CA"), w + 1)] = AtomID(
            original.residue(a + 1).atom_index("CA"), a + 1
        )
    assert len(ca_map), 'No atoms greater than cutoff'
    pyrosetta.rosetta.core.scoring.superimpose_pose(pose, original, ca_map)
    return ca_map
