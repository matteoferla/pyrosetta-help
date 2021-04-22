import pyrosetta
import re


def get_NGL_selection_from_AtomID(pose:pyrosetta.Pose, atom_id: pyrosetta.AtomID):
    pose_resi = atom_id.rsd()
    residue = pose.residue(pose_resi)
    atom_name = residue.atom_name(atom_id.atomno()).strip()
    pdb_resi, chain = pose.pdb_info().pose2pdb(pose_resi).strip().split()
    return f'[{residue.name3().strip()}]{pdb_resi}:{chain}.{atom_name}'


def print_constraint_score(pose:pyrosetta.Pose, con):
    a = get_NGL_selection_from_AtomID(pose, con.atom1())
    b = get_NGL_selection_from_AtomID(pose, con.atom2())
    fun = con.get_func()
    funpart = f'{fun.__class__.__name__}={fun.x0()}±{fun.sd()}'
    score = con.score(pose)
    print(f'{a} – {b} {funpart} --> {score:.2}')


def print_constraint_scores(pose:pyrosetta.Pose):
    cs = pose.constraint_set()
    for con in cs.get_all_constraints():
        print_constraint_score(pose, con)

def get_AtomID(pose, chain:str, resi:int, atomname:str) -> pyrosetta.rosetta.core.id.AtomID:
    r = pose.pdb_info().pdb2pose(res=resi, chain=chain)
    assert r != 0, f'{resi}:{chain} is absent'
    residue = pose.residue(r)
    return pyrosetta.rosetta.core.id.AtomID(atomno_in=residue.atom_index(atomname), rsd_in=r)

def get_AtomID_by_NGL_sele(pose, selection:str) -> pyrosetta.AtomID:
    """
    23:A.CA
    """
    if ' ' in selection.strip():
        raise ValueError('single atom selection')
    # chain
    if ':' in selection:
        chain = re.match(':(\w)', selection).group(1)
    else:
        chain = 'A'
    # atom name
    if ':' in selection:
        name = re.match('\.(\w+)', selection).group(1)
    else:
        name = ''
    # residue name
    if '[' in selection:
        resn = re.match('\[(\w+)\]', selection).group(1)
    else:
        resn = ''
    # residue index
    if re.search('\d', selection):
        resi = int(re.match('(\d+)', re.sub('\[.*\]', selection)).group(1))
    else:
        resi = float('nan')
    # assert
    if str(resi) == 'nan' or chain == '' or name == '':
        raise ValueError(f'selection {selection} is not like `23:A.CA`.')
    # return
    return get_AtomID(pose, chain=chain, resi=resi, atomname=name)