import re
from typing import (Optional)

import pyrosetta

__all__ = ['get_NGL_selection_from_AtomID',
           'print_constraint_score',
           'print_constraint_scores',
           'get_AtomID',
           'get_AtomID_by_NGL_sele',
           'get_AtomID_from_pymol_line',
           'get_constraint_score_data',
           'make_constraint_from_pymol_line',
           'print_bad_constraint_scores',
           'constraints2pandas']

def _get_cons(pose):
    """
    KofNConstraint and others are nested. This splits them up
    """
    cons = []
    cs = pose.constraint_set()
    for con in cs.get_all_constraints():
        if hasattr(con, 'member_constraints'):
            cons.extend(con.member_constraints())
        else:
            cons.append(con)
    return con

def get_NGL_selection_from_AtomID(pose: pyrosetta.Pose, atom_id: pyrosetta.AtomID, named: bool = False):
    """
    Given a pyrosetta AtomID give an NGL selection.
    NB ``named`` gives the residue name (``'[SER]3:A.CA'``) but is not a valid selection.

    :param pose:
    :param atom_id:
    :param named:
    :return:
    """
    pose_resi = atom_id.rsd()
    residue = pose.residue(pose_resi)
    atom_name = residue.atom_name(atom_id.atomno()).strip()
    pdb_resi, chain = pose.pdb_info().pose2pdb(pose_resi).strip().split()
    if named:
        return f'[{residue.name3().strip()}]{pdb_resi}:{chain}.{atom_name}'
    else:
        return f'{pdb_resi}:{chain}.{atom_name}'


def print_constraint_score(pose: pyrosetta.Pose, con):
    """
    Print the constraint details. atoms and score

    :param pose:
    :param con:
    :return:
    """
    if hasattr(con, 'member_constraints'):
        print(f'Nested: {con.__class__.__name__}')
        cons = con.member_constraints()
    else:
        cons = [con]
    for con in cons:
        data = get_constraint_score_data(pose, con)
        print('{constraint_name}: {atom_A} – {atom_B} {function_name}={funpart} --> {score:.2}'.format(**data))


def get_constraint_score_data(pose: pyrosetta.Pose, con) -> dict:
    """
    Print the constraint details. atoms and score

    :param pose:
    :param con:
    :return:
    """
    atoms = [get_NGL_selection_from_AtomID(pose, con.atom(i)) for i in range(1, con.natoms())]
    try:
        fun = con.get_func()
        function_name = fun.__class__.__name__
        if hasattr(fun, 'x0') and hasattr(fun, 'sd'):
            funpart = f'{fun.x0():.2f}/{fun.sd():.2f}'
        elif hasattr(fun, 'x0'):
            funpart = f'{fun.x0():.2f}'
        else:
            funpart = 'NA'
    except RuntimeError:
        function_name= 'unknown'
        funpart = 'NA'
    score = con.score(pose)
    info = dict(constraint_name=con.__class__.__name__,
                function_name=function_name,
                funpart=funpart,
                score=score
                )
    for i, atom in enumerate(atoms):
        info[f'atom_{i+1}'] = atom
    return info


def print_constraint_scores(pose: pyrosetta.Pose):
    """
    Prints the scores for each constraint in the pose

    :param pose:
    :return:
    """
    cs = pose.constraint_set()
    for con in cs.get_all_constraints():
        print_constraint_score(pose, con)


def constraints2pandas(pose):
    import pandas as pd
    rows = []
    for con in _get_cons(pose):
        rows.append(get_constraint_score_data(pose, con))
    return pd.DataFrame(rows).fillna('')


def print_bad_constraint_scores(pose, cutoff=0.5):
    for con in _get_cons(pose):
        data = get_constraint_score_data(pose, con)
        if data['score'] > cutoff:
            print('{constraint_name}: {atom_A} – {atom_B} {function_name}={funpart} --> {score:.2}'.format(**data))


def get_AtomID(pose: pyrosetta.Pose, chain: str, resi: int, atomname: str) -> pyrosetta.AtomID:
    r = pose.pdb_info().pdb2pose(res=resi, chain=chain)
    assert r != 0, f'{resi}:{chain} is absent'
    residue = pose.residue(r)
    return pyrosetta.AtomID(atomno_in=residue.atom_index(atomname), rsd_in=r)


def get_AtomID_by_NGL_sele(pose: pyrosetta.Pose, selection: str) -> pyrosetta.AtomID:
    """
    23:A.CA
    """
    if ' ' in selection.strip():
        raise ValueError('single atom selection')
    # chain
    if ':' in selection:
        chain = re.search(r':(\w)', selection).group(1)
    else:
        chain = 'A'
    # atom name
    if ':' in selection:
        name = re.search(r'\.(\w+)', selection).group(1)
    else:
        name = ''
    # residue name
    if '[' in selection:
        resn = re.search(r'\[(\w+)]', selection).group(1)
    else:
        resn = ''
    # residue index
    if re.search(r'\d', selection):
        # sub pattern was incorrect, '' added, but I am not sure what it does.
        resi = int(re.match(r'(\d+)', re.sub(r'\[.*]', '', selection)).group(1))
    else:
        resi = float('nan')
    # assert
    if str(resi) == 'nan' or chain == '' or name == '':
        raise ValueError(f'selection {selection} is not like `23:A.CA`.')
    # return
    return get_AtomID(pose, chain=chain, resi=resi, atomname=name)


def get_AtomID_from_pymol_line(pose: pyrosetta.Pose, line: Optional[str] = None) -> pyrosetta.rosetta.core.id.AtomID:
    """
    Given a copypaste from the console in pymol following an atom selection in edit mode:
    (``You clicked /1amq/B/A/PMP`413/N1 -> (pk2)``) returns that atom in PyRosetta.

    If the line argument is blank the clipboard is read.

    :param pose:
    :param line:
    :return:
    """
    # You clicked /1amq/B/A/PMP`413/N1 -> (pk2)
    if line is None:
        import xerox
        line = xerox.paste()
    rex = re.search(r'\/[\w_]+/\w?/(\w)/\w{3}`(\d+)/(\w+)', line)
    chain, res, atomname = rex.groups()
    return get_AtomID(pose, chain, int(res), atomname)


def make_constraint_from_pymol_line(pose: pyrosetta.Pose,
                                    lines: str) \
        -> pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint:
    """
    two atoms clicked...
    You clicked /1amq/A/A/ASP`222/OD2 -> (pk1)
    You clicked /1amq/B/A/PMP`413/N1 -> (pk2)
    distance measured in PyRosetta hence the pose.

    :param lines:
    :param pose:
    :return:
    """
    fore, aft = lines.strip().split('\n')
    HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
    AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
    fore_atom = get_AtomID_from_pymol_line(pose, fore)
    aft_atom = get_AtomID_from_pymol_line(pose, aft)
    fore_xyz = pose.residue(fore_atom.rsd()).xyz(fore_atom.atomno())
    aft_xyz = pose.residue(aft_atom.rsd()).xyz(aft_atom.atomno())
    d = (fore_xyz - aft_xyz).norm()
    return AtomPairConstraint(fore_atom, aft_atom, HarmonicFunc(x0_in=d, sd_in=0.2))
