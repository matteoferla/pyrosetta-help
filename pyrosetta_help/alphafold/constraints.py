__all__ = ['add_pae_constraints',
           'add_interchain_pae_constraints',
           'add_stretch_constraint',
           'make_pae_constraint',]

from typing import (Optional, List)

import numpy as np
import pyrosetta
from ..common_ops.distances import measure_distance_matrix


def add_pae_constraints(pose: pyrosetta.Pose,
                        errors: np.ndarray,
                        cutoff: float = 12,
                        tolerance: Optional[float] = None,
                        weight: float = 1,
                        adjecency_threshold=5) -> None:
    """

    Add constrains to the pose based on the errors matrix.
    NB. this matrix is a reshaped version of what AF2 returns.

    A harmonic function is added to CA atoms that are in residues with the error under a specified cutoff.
    The mu is the current distance and the standard deviation of the harmonic is the error times ``weight``.

    To find out how many were added:

    >>> len(pose.constraint_set().get_all_constraints())

    :param pose:
    :param errors:
    :param cutoff:
    :param tolerance: if None Harmonic, if value, tollerance of FlatHarmonic
    :param weight: this is added to the SD part so squared inverse.
    :param adjecency_threshold: min residue separation of sequence neighbours
    :return:
    """
    for r1_idx, r2_idx in np.argwhere(errors < cutoff):
        if abs(r1_idx - r2_idx) < adjecency_threshold:
            continue  # skip neighbours
        elif r1_idx <= r2_idx:
            continue  # add once.
        d_error = errors[r1_idx, r2_idx]
        apc = make_pae_constraint(pose=pose,
                                  residue1_pose_idx=r1_idx + 1,
                                  residue2_pose_idx=r2_idx + 1,
                                  error=d_error,
                                  weight=weight,
                                  tolerance=tolerance)
        pose.add_constraint(apc)


def make_pae_constraint(pose,
                        residue1_pose_idx: int,  # one indexed.
                        residue2_pose_idx: int,  # one indexed.
                        error: float,
                        tolerance: Optional[float] = None,
                        weight: float = 1):
    """
    Add a constraint between two residues based on the PAE error
    from AlphaFold2 (the colourful heatmap in EBI-AF2).

    :param pose:
    :param residue1_pose_idx:
    :param residue2_pose_idx:
    :param error:
    :param tolerance:
    :param weight:
    :return:
    """
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


def add_interchain_pae_constraints(pose, errors, cutoff=15):
    """
    Add constraints between residues that are interacting according to the PAE error matrix
    but are in different 'chains' (sensu PyRosetta FoldTree).

    :param pose:
    :param errors:
    :param cutoff:
    :return:
    """
    xdistances = measure_distance_matrix(pose)
    for c in (1, 2):
        xdistances[pose.chain_begin(c) - 1: pose.chain_end(c),
        pose.chain_begin(c) - 1: pose.chain_end(c)] = np.nan
    with np.errstate(invalid='ignore'):
        mask = xdistances < cutoff
    for r1_idx, r2_idx in np.argwhere(mask):
        apc = make_pae_constraint(pose=pose,
                                  residue1_pose_idx=r1_idx + 1,
                                  residue2_pose_idx=r2_idx + 1,
                                  error=errors[r1_idx, r2_idx],
                                  tolerance=3)  # flatharmonic
        pose.add_constraint(apc)


def add_stretch_constraint(pose: pyrosetta.Pose,
                           weight: float = 5,
                           slope_in: float = -0.05,
                           residue_index_A: int = 1,
                           residue_index_B: int = -1,
                           distance: Optional[float] = None,
                           sigmoid: bool = True
                           ) -> pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint:
    """
    Add a constraint to "stretch out" the model, because ``slope_in`` is negative.
    The weight needs to be negative for sigmoid=False or it will attractive

    :param pose: Pose to add constraint to
    :param weight: how strength of constraint (max of 0.5 for ``SigmoidFunc``)
    :param slope_in: negative number to stretch
    :param residue_index_A: first residue?
    :param residue_index_B: last residue is "-1"
    :param distance: if omitted, the midpoint of Sigmoid will be the current distance
    :param sigmoid: use sigmoid or identity/linear (bad idea)
    :return:
    """
    # get current length
    if residue_index_B == -1:
        residue_index_B = pose.total_residue()
    assert pose.residue(residue_index_A).is_protein, f'residue idx {residue_index_A} is not an AA'
    assert pose.residue(residue_index_B).is_protein, f'residue idx {residue_index_B} is not an AA'
    first_ca = pyrosetta.AtomID(atomno_in=pose.residue(residue_index_A).atom_index('CA'),
                                rsd_in=residue_index_A)
    last_ca = pyrosetta.AtomID(atomno_in=pose.residue(residue_index_B).atom_index('CA'),
                               rsd_in=residue_index_B)
    first_ca_xyz = pose.residue(residue_index_A).xyz(first_ca.atomno())
    last_ca_xyz = pose.residue(residue_index_B).xyz(last_ca.atomno())
    if distance is None:
        distance = (first_ca_xyz - last_ca_xyz).norm()
    # make & add con
    sf = pyrosetta.rosetta.core.scoring.func
    AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint  # noqa
    if sigmoid:
        fun = sf.ScalarWeightedFunc(weight, sf.SigmoidFunc(x0_in=distance, slope_in=slope_in))
    else:
        fun = sf.ScalarWeightedFunc(weight / distance, sf.IdentityFunc())
    con = AtomPairConstraint(first_ca, last_ca, fun)
    pose.add_constraint(con)
    return con

