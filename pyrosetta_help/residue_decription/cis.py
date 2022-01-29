from typing import (List)

import pyrosetta

from ..common_ops import pose_range


def get_cis_residues(pose: pyrosetta.Pose) -> List[int]:
    """
    Returns the pose indices of residues in cis (omega of zero)
    """
    # the cis residue is actually the one after (also a filter is not a generator)
    cis_gen = filter(lambda ri: abs(pose.omega(ri)) < 90., pose_range(pose))
    # citri is not really the plural of cis as cis is a preposition, citer is its adjective.
    # but best not discuss Latin grammar _ulterior_-ly.
    citri: List[int] = []
    for i in cis_gen:  #: int
        if i + 1 > pose.total_residue():
            # terminal residue
            continue
        if not pose.residue(i).connected_residue_at_upper():
            # gap
            continue
        citri.append(i + 1)
    return citri
