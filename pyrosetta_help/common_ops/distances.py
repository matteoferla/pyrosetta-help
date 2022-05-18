__all__ = ['measure_distance_matrix',
           'measure_ligand_distances',
           'measure_inter_residue_distance']

from typing import (Optional, List)
from types import ModuleType

import numpy as np
import pyrosetta
import itertools

residue_selector = pyrosetta.rosetta.core.select.residue_selector  # ModuleType
utility = pyrosetta.rosetta.utility  # noqa: F821
chemical = pyrosetta.rosetta.core.chemical  # ModuleType


def measure_distance_matrix(pose) -> np.ndarray:
    """
    Note the distance matrix is zero indexed as it would be confusing using numpy with one indexed data.

    :param pose:
    :return:
    """
    distances = np.zeros((pose.total_residue(), pose.total_residue())) * np.nan
    residue_iter = map(lambda r: pose.residue(r + 1), range(pose.total_residue()))
    ca_xyzs = [res.xyz('CA') if res.is_protein() else None for res in residue_iter]
    for i in range(len(ca_xyzs)):
        if i is None:
            continue
        for j in range(i):
            d = ca_xyzs[i].distance(ca_xyzs[j])
            distances[i, j] = d
            distances[j, i] = d
    return distances


def measure_ligand_distances(pose: pyrosetta.Pose, target_residue_idx: int) -> List[dict]:
    """
    Get the distances to the ligand from the target residue —closest atom, not centroid.
    Returns a list of dictionaries, one for each ligand, like:
    ``{'ligand_name': ' CA', 'ligand_idx': 430, 'distance': 23.06283024797271}``
    example:

    .. code-block:: python
        import operator
        import pyrosetta_help as help
        distances:List[dict] = ph.get_ligand_distances(pose, 20)
        closest = sorted(distances, key=operator.itemgetter('distance'))[0]
    """
    lig_vector: utility.vector1_bool = residue_selector.ResiduePropertySelector(chemical.ResidueProperty.LIGAND) \
        .apply(pose)
    lig_resis: utility.vector1_unsigned_long = residue_selector.ResidueVector(lig_vector)
    target_residue = pose.residue(target_residue_idx)
    distances = []
    for ligand in map(pose.residue, lig_resis):  #: pyrosetta.Residue
        distances.append(dict(ligand_name=ligand.name3(),
                              ligand_idx=ligand.seqpos(),
                              distance=measure_inter_residue_distance(pose, target_residue_idx, ligand.seqpos())
                              )
                         )
    return distances


def measure_inter_residue_distance(pose: pyrosetta.Pose, query_residue_idx: int, target_residue_idx: int) -> float:
    """
    Get the distances between two residues —closest atom, not centroid.
    Virtual residues may cause problems.
    """
    fortran_range = lambda max_: range(1, max_ + 1)  # noqa: E731 is stupid
    get_xyzs = lambda residue: [residue.xyz(i) for i in fortran_range(residue.natoms())]  # noqa: E731 is stupid
    target_residue = pose.residue(target_residue_idx)
    query_residue = pose.residue(query_residue_idx)
    xyz_gen = itertools.product(get_xyzs(target_residue), get_xyzs(query_residue))
    return min([l_xyz.distance(t_xyz) for l_xyz, t_xyz in xyz_gen])
