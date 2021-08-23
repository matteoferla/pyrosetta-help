"""
There is likely a mover that does this already and better...
"""


def make_ss(pose, begin: int = 1, end: int = -1, phi: float = 180, psi: float = 180):
    """
    Given a pose, likely one made via ``pyrosetta.pose_from_sequence``,
    make the fortran-indexed range [begin, end] (inclusive) into that SS.
    """
    if end == -1:
        end = pose.total_residue()
    for r in range(begin, end + 1):
        pose.set_phi(r, phi)
        pose.set_psi(r, psi)


def make_alpha_helical(pose, begin: int = 1, end: int = -1):
    make_ss(pose, begin, end, phi=-57.8, psi=-47.0)


def make_310_helical(pose, begin: int = 1, end: int = -1):
    make_ss(pose, begin, end, phi=-74.0, psi=-4.0)


def make_pi_helical(pose, begin: int = 1, end: int = -1):
    make_ss(pose, begin, end, phi=-57.1, psi=-69.7)


def make_sheet(pose, begin: int = 1, end: int = -1):
    make_ss(pose, begin, end, phi=-139, psi=+135)

# PyCharm is licking a mandrake again. Ignore the warning of bytes + str
make_alpha_helical.__doc__ = make_ss.__doc__ + '\nIn this case alpha-helix with phi=-57.8, psi=-47.0'
make_310_helical.__doc__ = make_ss.__doc__ + '\nIn this case 3.10-helix with phi=-74.0, psi=-4.0'
make_pi_helical.__doc__ = make_ss.__doc__ + '\nIn this case pi-helix with phi=-57.1, psi=-69.7'
make_sheet.__doc__ = make_ss.__doc__ + '\nIn this case sheet with phi=-139, psi=+135'
