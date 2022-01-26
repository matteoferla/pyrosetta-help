import pyrosetta

def get_ss(pose: pyrosetta.Pose) -> str:
    """
    returns a string of the SS of types H, S, L
    This is not the ALBEGO classifier
    """
    DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()  # noqa
    DSSP.apply(pose)
    return pose.secstruct()