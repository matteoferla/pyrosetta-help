import pyrosetta

Residue = pyrosetta.rosetta.core.conformation.Residue  # typehinting only


def is_xlinked(residue: Residue) -> bool:
    return get_xlink_idx(residue, raise_on_bond=False) != 0


def get_xlink_idx(residue: Residue, raise_on_bond: bool = True) -> int:
    cxi: int = sum((residue.has_lower_connect(), residue.has_upper_connect())) + 1
    if residue.n_current_residue_connections() >= cxi:
        return cxi
    elif raise_on_bond:
        raise ValueError('Not a crosslinked residue')
    else:
        return 0


def get_xlink_details(residue_i: int, pose: pyrosetta.Pose) -> dict:
    """
    cross-link (SSBOND, LINK etc.) => conn3
    By conn3 I mean a connection that is not
    lower (own N) or upper (own C).
    So if a residue has no connection on one or both of these,
    technically it will be conn2 or conn1, but that's being pedantic.

    NB. Does not check for incomplete connections (which are a product of fancy meddling)

    return
    :param residue_i:
    :param pose:
    :return: dict of keys other_idx, other_name3, other_is_protein own_atom_name other_atom_name
    """
    residue: Residue = pose.residue(residue_i)
    cxi: int = get_xlink_idx(residue)
    other_idx: int = residue.connected_residue_at_resconn(cxi)
    other: Residue = pose.residue(other_idx)
    atom_idx: int = residue.residue_connect_atom_index(cxi)
    # there **will** be a conn: (unless residue.has_incomplete_connection(3) == True)
    other_cxi: int = other.connections_to_residue(residue_i)[1]
    other_atom_idx: int = residue.residue_connect_atom_index(cxi)
    return dict(other_idx=other_idx,
                other_name3=other.name3(),
                other_is_protein=other.is_protein(),
                own_atom_name=residue.atom_name(atom_idx).strip(),
                other_atom_name=other.atom_name(other_atom_idx).strip(),
                )
