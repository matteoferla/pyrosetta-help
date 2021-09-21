__all__ = ['pose_from_file',
           'pose2pandas',
           'make_blank_pose',
           'add_bfactor_from_score',
           'get_last_res_in_chain',
           'clarify_selector',
           'count_ligands',
           'correct_numbering',
           'get_pdbstr']

from typing import *
import pyrosetta
import pandas as pd
import numpy as np
from warnings import warn
from collections import Counter


def pose_from_file(pdb_filename: str,
                   params_filenames: Optional[Union[pyrosetta.rosetta.utility.vector1_string, List[str]]] = None) \
        -> pyrosetta.Pose:
    """
    Return a pose like pose_from_file but with params.

    :param pdb_filename:
    :param params_filenames:
    :return:
    """
    pose = make_blank_pose(params_filenames)
    pyrosetta.rosetta.core.import_pose.pose_from_file(pose, pdb_filename)
    return pose


vector1_string = pyrosetta.rosetta.utility.vector1_string


def make_blank_pose(params_filenames: Optional[Union[vector1_string, List[str]]] = None) \
        -> pyrosetta.Pose:
    """
    Returns an empty pose, but with params from files.

    :param params_filenames:
    :return:
    """
    pose = pyrosetta.Pose()
    if params_filenames and isinstance(params_filenames, vector1_string):
        pyrosetta.generate_nonstandard_residue_set(pose, params_filenames)
    if params_filenames and isinstance(params_filenames, list):
        params_filenames2 = vector1_string()
        params_filenames2.extend([f for f in params_filenames if f])
        pyrosetta.generate_nonstandard_residue_set(pose, params_filenames2)
    elif params_filenames:
        raise TypeError(f'Unexpected type {type(params_filenames)}')
    else:
        pass
    return pose


def pose2pandas(pose: pyrosetta.Pose, scorefxn: pyrosetta.ScoreFunction) -> pd.DataFrame:
    """
    Return a pandas dataframe from the scores of the pose

    :param pose:
    :return:
    """
    pose.energies().clear_energies()
    scorefxn.weights() # neccessary?
    emopts = pyrosetta.rosetta.core.scoring.methods.EnergyMethodOptions(scorefxn.energy_method_options())
    emopts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    scorefxn.set_energy_method_options(emopts)
    scorefxn(pose)
    scores = pd.DataFrame(pose.energies().residue_total_energies_array())
    pi = pose.pdb_info()
    scores['residue'] = scores.index.to_series() \
        .apply(lambda r: pose.residue(r + 1) \
               .name1() + pi.pose2pdb(r + 1)
               )
    return scores


def add_bfactor_from_score(pose: pyrosetta.Pose):
    """
    Adds the bfactors from total_score.
    Snippet for testing in Jupyter

    >>> import nglview as nv
    >>> view = nv.show_rosetta(pose)
    >>> # view = nv.show_file('test.cif')
    >>> view.clear_representations()
    >>> view.add_tube(radiusType="bfactor", color="bfactor", radiusScale=0.10, colorScale="RdYlBu")
    >>> view

    ``replace_res_remap_bfactors`` may have been a cleaner strategy. This was quicker to write.

    If this fails, it may be because the pose was not scored first.
    """
    if pose.pdb_info().obsolete():
        raise ValueError('Pose pdb_info is flagged as obsolete (change `pose.pdb_info().obsolete(False)`)')
    # scores
    energies = pose.energies()

    def get_res_score(res):
        total_score = pyrosetta.rosetta.core.scoring.ScoreType.total_score
        # if pose.residue(res).is_polymer()
        try:
            return energies.residue_total_energies(res)[total_score]
        except:
            return float('nan')

    # the array goes from zero (nan) to n_residues
    total_scores = np.array([float('nan')] + [get_res_score(res) for res in range(1, pose.total_residue() + 1)])
    mask = np.isnan(total_scores)
    total_scores -= np.nanmin(total_scores)
    total_scores *= 100 / np.nanmax(total_scores)
    total_scores = np.nan_to_num(total_scores, nan=100)
    total_scores[mask] = 0.
    # add to pose
    pdb_info = pose.pdb_info()
    for res in range(pose.total_residue()):
        for i in range(pose.residue(res + 1).natoms()):
            pdb_info.bfactor(res + 1, i + 1, total_scores[res + 1])


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


def clarify_selector(selector: pyrosetta.rosetta.core.select.residue_selector.ResidueSelector,
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
    return [f'[{pose.residue(r).name3()}]{pose2pdb(r).strip().replace(" ", ":")}' for r in rv]

def count_ligands(pose: pyrosetta.Pose) -> List[Tuple[str, int]]:
    LIGAND = pyrosetta.rosetta.core.chemical.ResidueProperty.LIGAND
    pr_rs = pyrosetta.rosetta.core.select.residue_selector
    lig_sele = pr_rs.ResiduePropertySelector(LIGAND)
    count = Counter([pose.residue(i).name3() for i in pr_rs.ResidueVector(lig_sele.apply(pose))])
    return count.most_common()


def correct_numbering(pose):
    """
    A fresh PDBInfo has the PDB residue number the same as the pose one
    as opposed to restarting per chain.

    :param pose:
    :return:
    """
    pdb_info = pose.pdb_info()
    if pdb_info is None:
        pdb_mover = pyrosetta.rosetta.protocols.simple_moves.AddPDBInfoMover()
        pdb_mover.apply(pose)
        pdb_info = pose.pdb_info()
    current_chain = 'A'
    r = 1
    for i in range(1, pose.total_residue() +1):
        if pdb_info.chain(i) != current_chain:
            current_chain = pdb_info.chain(i)
            r = 1
            pdb_info.number(i, 1)
        else:
            pdb_info.number(i, r)
            r += 1

def get_pdbstr(pose):
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    return buffer.str()

