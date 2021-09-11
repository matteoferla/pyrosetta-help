__all__ = ['add_pae_constraints', 'add_stretch_constraint', 'get_distance_matrix', 'make_pae_constraint']

import pyrosetta
import numpy as np
from typing import *


def add_pae_constraints(pose: pyrosetta.Pose,
                        errors: np.ndarray,
                        cutoff: float = 12,
                        tolerance: Optional[float] = None,
                        weight: float = 1,
                        adjecency_threshold=5,
                        blank: bool = True):
    """

    Add constrains to the pose based on the errors matrix.
    NB. this matrix is a reshaped version of what AF2 returns.

    A harmonic function is added to CA atoms that are in residues with the error under a specified cutoff.
    The mu is the current distance and the standard deviation of the harmonic is the error times ``weight``.

    ``blank`` as False keeps the current constaints.

    To find out how many were added:

    >>> len(pose.constraint_set().get_all_constraints())

    :param pose:
    :param errors:
    :param cutoff:
    :param tolerance: if None Harmonic, if value, tollerance of FlatHarmonic
    :param weight: this is added to the SD part so squared inverse.
    :param adjecency_threshold: min residue separation of sequence neighbours
    :param blank:
    :return:
    """
    cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
    if not blank:
        previous = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet().get_all_constraints()
        for con in previous:
            cs.add_constraint(con)
    for r1_idx, r2_idx in np.argwhere(errors < cutoff):
        if abs(r1_idx - r2_idx) < adjecency_threshold:
            continue # skip neighbours
        elif r1_idx <= r2_idx:
            continue # add once.
        d_error = errors[r1_idx, r2_idx]
        apc = make_pae_constraint(pose=pose,
                                  residue1_pose_idx=r1_idx + 1,
                                  residue2_pose_idx=r2_idx + 1,
                                  error=d_error,
                                  weight=weight,
                                  tolerance=tolerance)
        cs.add_constraint(apc)
    setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    setup.constraint_set(cs)
    setup.apply(pose)
    return cs

def make_pae_constraint(pose,
                        residue1_pose_idx:int, # one indexed.
                        residue2_pose_idx:int, # one indexed.
                        error:float,
                        tolerance: Optional[float] = None,
                        weight:float=1):
    get_ca = lambda r, i: pyrosetta.AtomID(atomno_in=r.atom_index('CA'), rsd_in=i)
    FlatHarmonicFunc = pyrosetta.rosetta.core.scoring.func.FlatHarmonicFunc
    HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc

    AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
    residue1 = pose.residue(residue1_pose_idx)
    ca1_atom = get_ca(residue1, residue1_pose_idx)
    residue2 = pose.residue(residue2_pose_idx)
    ca2_atom = get_ca(residue2, residue2_pose_idx)
    ca1_xyz = residue1.xyz(ca1_atom.atomno())
    ca2_xyz = residue2.xyz(ca2_atom.atomno())
    d = (ca1_xyz - ca2_xyz).norm()
    if not tolerance:
        fun = HarmonicFunc(x0_in=d, sd_in=error * weight)
    else:
        fun = FlatHarmonicFunc(x0_in=d, sd_in=error * weight, tol_in=tolerance)
    return AtomPairConstraint(ca1_atom, ca2_atom, fun)


def add_stretch_constraint(pose: pyrosetta.Pose,
                           weight: float = 5,
                           slope_in: float = -0.05,
                           residue_index_A: int = 1,
                           residue_index_B: int = -1,
                           distance: Optional[
                               float] = None) -> pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint:
    """
    Add a constraint to "stretch out" the model, because ``slope_in`` is negative.

    :param pose: Pose to add constraint to
    :param weight: how strength of constraint (max of 0.5 for ``SigmoidFunc``)
    :param slope_in: negative number to stretch
    :param residue_index_A: first residue?
    :param residue_index_B: last residue is "-1"
    :param distance: if omitted, the midpoint of Sigmoid will be the current distance
    :return:
    """
    # get current length
    first_ca = pyrosetta.AtomID(atomno_in=pose.residue(residue_index_A).atom_index('CA'),
                                rsd_in=residue_index_A)
    if residue_index_B == -1:
        residue_index_B = pose.total_residue()
    last_ca = pyrosetta.AtomID(atomno_in=pose.residue(residue_index_B).atom_index('CA'),
                               rsd_in=residue_index_B)
    first_ca_xyz = pose.residue(1).xyz(first_ca.atomno())
    last_ca_xyz = pose.residue(pose.total_residue()).xyz(last_ca.atomno())
    if distance is None:
        distance = (first_ca_xyz - last_ca_xyz).norm()
    # make & add con
    sf = pyrosetta.rosetta.core.scoring.func
    AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
    fun = sf.ScalarWeightedFunc(weight, sf.SigmoidFunc(x0_in=distance, slope_in=slope_in))
    con = AtomPairConstraint(first_ca, last_ca, fun)
    pose.constraint_set().add_constraint(con)
    return con

def get_distance_matrix(pose):
    distances = np.zeros((pose.total_residue(), pose.total_residue()))
    ca_xyzs = [pose.residue(r).xyz('CA') for r in range(1, pose.total_residue() + 1)]
    for i in range(len(ca_xyzs)):
        for j in range(i):
            d = ca_xyzs[i].distance(ca_xyzs[j])
            distances[i, j] = d
            distances[j, i] = d
    return distances
