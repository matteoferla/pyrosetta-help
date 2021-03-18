from typing import *
import pyrosetta
import pandas as pd
import numpy as np

def pose_from_file(pdb_filename: str,
                   params_filenames: Optional[Union[pyrosetta.rosetta.utility.vector1_string, List[str]]] = None) \
        -> pyrosetta.Pose:
    """
    Return a pose like pose_from_file but with params.

    :param pdb_filename:
    :param params_filenames:
    :return:
    """
    pose = pyrosetta.Pose()
    if params_filenames and isinstance(params_filenames, pyrosetta.rosetta.utility.vector1_string):
        pyrosetta.generate_nonstandard_residue_set(pose, params_filenames)
    if params_filenames and isinstance(params_filenames, list):
        params_filenames2 = pyrosetta.rosetta.utility.vector1_string()
        params_filenames2.extend(params_filenames)
        pyrosetta.generate_nonstandard_residue_set(pose, params_filenames2)
    else:
        pass
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdb_filename)
    return pose

def get_local_scorefxn() -> pyrosetta.ScoreFunction:
    ## local scorefxn w/ ED
    # ---------- Local weights ------------------
    # these are mostly the same except for the line:
    # <Set scale_sc_dens_byres="E:0.56,D:0.56,R:0.76,K:0.76,M:0.76,C:0.81,Q:0.81,H:0.81,N:0.81,T:0.81,S:0.81,
    # Y:0.88,W:0.88,A:0.88,F:0.88,P:0.88,I:0.88,L:0.88,V:0.88"/>
    # which is an utter mystery.
    weights = {"cart_bonded_length": 0.5,
               "cart_bonded_torsion": 0.5,
               "cart_bonded_angle": 1.0,
               "pro_close": 0.0,
               "fa_sol": 0.0,
               "elec_dens_fast": 30,  # <-- ED
               "rama": 0.0,
               "rama_prepro": 1.0}
    scorefxn_local = pyrosetta.create_score_function('ref2015_cart')
    stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
    for name, value in weights.items():
        scorefxn_local.set_weight(stm.score_type_from_name(name), value)
    return scorefxn_local

def prep_ED(pose: pyrosetta.Pose, map_filename: str) -> pyrosetta.rosetta.core.scoring.electron_density.ElectronDensity:
    # rmsd & ED fit
    rmsd = pyrosetta.rosetta.core.simple_metrics.metrics.RMSDMetric(pose)
    ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(map_filename)
    initial_fit = ED.matchPose(pose)
    # This is redundant with ED.
    map_mover = pyrosetta.rosetta.protocols.cryst.LoadDensityMapMover(map_filename)
    map_mover.apply(pose)
    # This is redundant with map
    sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
    sdsm.apply(pose)
    return ED


def get_local_relax(scorefxn: Optional[pyrosetta.ScoreFunction] = None,
                    ncyc: int = 3,
                    nexp: int = 2) -> pyrosetta.rosetta.protocols.relax.LocalRelax:
    if scorefxn is None:
        scorefxn = get_local_scorefxn()
    # <LocalRelax name="local_rlx" scorefxn="dens" max_iter="100" ncyc="1" ramp_cart="0" K="16" nexp="2"/>
    relax = pyrosetta.rosetta.protocols.relax.LocalRelax()
    relax.set_sfxn(scorefxn)
    relax.set_K(16)
    relax.set_max_iter(100)
    relax.set_ncyc(ncyc)
    relax.set_nexp(nexp)
    return relax

def do_local_relax(pose: pyrosetta.Pose, scorefxn: Optional[pyrosetta.ScoreFunction] = None) -> None:
    get_local_relax(scorefxn).apply(pose)


def do_chainwise_relax(pose: pyrosetta.Pose,
                       scorefxn: Optional[pyrosetta.ScoreFunction] = None,
                       cycles: int = 5) -> None:
    if scorefxn is None:
        scorefxn = pyrosetta.get_fa_scorefxn()
    for chain_i in range(1, pose.num_chains() + 1):
        chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_i)
        chain_vector = chain_sele.apply(pose)
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(allow_bb=chain_vector)
        movemap.set_chi(allow_chi=chain_vector)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap(movemap)
        relax.apply(pose)


def pose2pandas(pose: pyrosetta.Pose, scorefxn: pyrosetta.ScoreFunction) -> pd.DataFrame:
    """
    Return a pandas dataframe from the scores of the pose
    
    :param pose:
    :return:
    """
    pose.energies().clear_energies()
    scorefxn(pose)
    scores = pd.DataFrame(pose.energies().residue_total_energies_array())
    pi = pose.pdb_info()
    scores['residue'] = scores.index.to_series()\
                              .apply(lambda r: pose.residue(r+1)\
                                                   .name1() + pi.pose2pdb(r+1)
                                    )
    return scores

def add_bfactor_from_score(pose: pyrosetta.Pose):
    """
    Adds the bfactors from total_score
    ``replace_res_remap_bfactors`` may have been a cleaner strategy. This was quicker to write.

    Check this is not the segfaulting version first!!
    """
    # scores
    energies = pose.energies()
    total_score = pyrosetta.rosetta.core.scoring.ScoreType.total_score
    # add to pose
    pdb_info = pose.pdb_info()
    pdb2pose = pdb_info.pdb2pose
    total_scores = np.array([float('nan')] + [energies.residue_total_energies(res)[total_score] for res in
                                              range(1, pose.total_residue() + 1)])
    total_scores -= np.nanmin(total_scores)
    total_scores *= 100 / np.nanmax(total_scores)
    total_scores = np.nan_to_num(total_scores, nan=100)
    for res in range(1, pose.total_residue() + 1):
        for i in range(pose.residue(res).natoms()):
            pdb_info.bfactor(res, i, total_scores[res])

def get_last_res_in_chain(pose, chain='A') -> int:
    """
    Returns last residue index in a chain. THere is probably a mover that does this.

    :param pose:
    :param chain: letter or number
    :return:
    """
    cv = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain).apply(pose)
    rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(cv)
    return max(rv)

def clarify_selector(selector:  pyrosetta.rosetta.core.select.residue_selector.ResidueSelector,
                   pose: pyrosetta.Pose) -> List['str']:
    """
    Given a selector and pose return a list of residues in NGL selection format
    Example, [CMP]787:H

    :param selector:
    :param pose:
    :return: list of residues in NGL selection format
    """
    pose2pdb = pose.pdb_info().pose2pdb
    vector = selector.apply(pose)
    rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(vector)
    return [f'[{pose.residue(r).name3()}]{pose2pdb(r).strip().replace(" ",":")}' for r in rv]


def download_map(code: str):
    import shutil
    import urllib.request as request
    from contextlib import closing
    ftp_path = f'ftp://ftp.ebi.ac.uk/pub/databases/emdb/structures/{code}/map/emd_{code.replace("EMD-", "")}.map.gz'
    file_path = f'{code}.map.gz'
    with closing(request.urlopen(ftp_path)) as r, open(file_path, 'wb') as f:
        shutil.copyfileobj(r, f)



def download_cif(code: str):
    import requests
    import shutil
    """
    Download CIF. Pyrosetta has issues importing Cifs (hetatms).
    This is just a note to self as PyMOL fetch can do it.
    """
    http_path = f'http://www.ebi.ac.uk/pdbe/static/entry/download/{code.lower()}-assembly-1.cif.gz'
    file_path = f'{code}.cif.gz'
    r = requests.get(http_path, stream=True)
    if r.status_code == 200:
        with open(file_path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
        return True
    else:
        return False

def download_opm(code: str):
    """
    Download OPM PDB
    """
    import shutil, requests
    http_path = f'https://opm-assets.storage.googleapis.com/pdb/{code.lower()}.pdb'
    file_path = f'{code}_OMP.pdb'
    r = requests.get(http_path, stream=True)
    if r.status_code == 200:
        with open(file_path, 'wb') as f:
            r.raw.decode_content = True
            shutil.copyfileobj(r.raw, f)
        return True
    else:
        return False