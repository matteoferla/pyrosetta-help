__all__ = ['constrain_distances']

import pyrosetta
import numpy as np

def constrain_distances(pose: pyrosetta.Pose,
                        errors: np.ndarray,
                        cutoff: float = 5,
                        weight: float = 1,
                        blank: bool = True):
    """
    Add constrains to the pose based on the errors matrix.
    NB. this matrix is a reshaped version of what AF2 returns.

    A harmonic function is added to CA atoms that are in residues with the error under a specified cutoff.
    The mu is the current distance and the standard deviation of the harmonic is the error times ``weight``.

    ``blank`` as False keeps the current constaints.

    To find out how many were added:

    >>> len(pose.constraint_set().get_all_constraints())
    """
    get_ca = lambda r, i: pyrosetta.AtomID(atomno_in=r.atom_index('CA'), rsd_in=i)
    HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
    AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
    cs = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet()
    if not blank:
        previous = pyrosetta.rosetta.core.scoring.constraints.ConstraintSet().get_all_constraints()
        for con in previous:
            cs.add_constraint(con)
    for r1_idx, r2_idx in np.argwhere(errors < cutoff):
        d_error = errors[r1_idx, r2_idx]
        residue1 = pose.residue(r1_idx + 1)
        ca1_atom = get_ca(residue1, r1_idx + 1)
        residue2 = pose.residue(r2_idx + 1)
        ca2_atom = get_ca(residue2, r2_idx + 1)
        ca1_xyz = residue1.xyz(ca1_atom.atomno())
        ca2_xyz = residue2.xyz(ca2_atom.atomno())
        d = (ca1_xyz - ca2_xyz).norm()
        apc = AtomPairConstraint(ca1_atom, ca2_atom, HarmonicFunc(x0_in=d, sd_in=d_error * weight))
        cs.add_constraint(apc)
    setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
    setup.constraint_set(cs)
    setup.apply(pose)
    return cs
