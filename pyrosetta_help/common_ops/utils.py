__all__ = ['pose_from_file',
           'pose2pandas',
           'make_blank_pose',
           'add_bfactor_from_score',
           'get_last_res_in_chain',
           'clarify_selector',
           'count_ligands',
           'correct_numbering',
           'get_pdbstr',
           'pose_range',
           'fix_offset',
           'assign_chain_letter',
           'what_is_chain']

from collections import Counter
from typing import (Optional, Tuple, Union, Iterable, Counter, List)
from Bio import pairwise2
from IPython.display import display, HTML
import string
import numpy as np
import pandas as pd
import pyrosetta
import functools


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

def correct_numbering(pose, chain_references:Optional[Union[None, str]]=None, allow_insertions:bool=False):
    """
    A fresh PDBInfo has the PDB residue number the same as the pose one
    as opposed to restarting per chain this fixes it.

    Additionaly, if the list ``chain_references`` is provided, any element as None is ignored,
    but any element with a sequence (str) where the first char is pos 1, will be used
    for numbering by calling ``fix_offset``.

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
    if not chain_references:
        return
    for i, ref_sequence in enumerate(chain_references):
        chain_i = i + 1
        if ref_sequence is None:
            continue
        assert isinstance(ref_sequence, str)
        fix_offset(pose, chain_i, ref_sequence, allow_insertions)


def fix_offset(pose, chain_i: int, ref_sequence: str, allow_insertions:bool=False) -> HTML:
    """
    Fix the pose numbering on ``chain_i`` based on ``ref_sequence``
    No gaps accepted in ref_sequence if ``allow_insertions`` is False, so pad with Xs appropriately if only termini.
    If ``allow_insertions`` is True, insertion codes are used.
    """
    pdb_info = pose.pdb_info()
    assert not pdb_info.obsolete(), "pdb info is marked as obsolete"
    # ## Align
    alignment = pairwise2.align.globalms(ref_sequence,
                                         pose.chain_sequence(chain_i),
                                         1, -2,  # penalise mismatches by x2
                                         -1, -0.5  # penalise extension by 1/2
                                         )[0]
    # this is relative to alignment and not unaligned seqA, so may have gaps if insertions are allowed:
    resi_map = [i + 1 for i, l in enumerate(alignment.seqB) if l != '-']
    if alignment.seqA.find('-') == -1:  # the ref seq aligns w/o gaps
        pass
    elif not allow_insertions:
        raise ValueError('There is an insertion in the pose, which is not allowed by attribute `allow_insertions` == False')
    else:
        gaps = [i + 1 for i, l in enumerate(alignment.seqA) if l == '-']
        for start, end in ph.rangify(gaps):
            resi_map = resi_map[:start] + [(start - 1, letter) for letter in string.ascii_uppercase[:end + 1 - start]] + resi_map[start:]

    pdb_info = pose.pdb_info()
    chain_resi_offset = -1
    for resi in ph.pose_range(pose):
        residue = pose.residue(resi)
        if residue.chain() != chain_i:
            continue
        if chain_resi_offset == -1:
            chain_resi_offset = resi - 1
        neoresi = resi_map[resi - chain_resi_offset - 1]
        if isinstance(neoresi, int):
            pdb_info.number(resi, neoresi)
        else: # isinstance(neoresi, tuple)
            neoresi, icode = neoresi
            pdb_info.number(resi, neoresi)
            pdb_info.icode(resi, icode)
    # ## Output
    formatted = pairwise2.format_alignment(*alignment)
    a, gap, b, score = formatted.strip().split('\n')
    gap = ''.join(['.' if c == '|' else '*' for c in gap])
    return HTML(f'<div style="font-family:monospace; display: inline-block; white-space: nowrap;">' +
                f'{a}<br/>{gap.replace(" ", "*")}<br/>{b}<br/>{score}</div>')

def assign_chain_letter(pose, chain_i: int, new_letter: str):
    """
    Assign a new letter to chain number ``chain_i``
    """
    pdb_info = pose.pdb_info()
    assert not pdb_info.obsolete(), "pdb info is marked as obsolete"
    for resi in ph.pose_range(pose):
        residue = pose.residue(resi)
        if residue.chain() != chain_i:
            continue
        pdb_info.chain(resi, new_letter)

def get_pdbstr(pose):
    buffer = pyrosetta.rosetta.std.stringbuf()
    pose.dump_pdb(pyrosetta.rosetta.std.ostream(buffer))
    return buffer.str()


def pose_range(pose: pyrosetta.Pose, protein_only=True) -> Iterable:
    """
    range but for the pose...
    For now naive of chains, non-amino acid residues etc.
    """
    all_iter = range(1, pose.total_residue()+1)
    if protein_only:
        return filter(lambda i: pose.residue(i).is_protein(), all_iter)
    return all_iter


def what_is_chain(pose: pyrosetta.Pose, chain: Union[str, int]):
    """
    Given a pose and a chain integer or character,
    return the character or integer of that chain

    :param pose:
    :param chain: str --> pdb_info based. int --> pyrosetta.Pose based
    :return:
    """
    pdb_info = pose.pdb_info()
    if isinstance(chain, int):
        for resi in ph.pose_range(pose):
            residue = pose.residue(resi)
            if residue.chain() != chain:
                continue
            return pdb_info.chain(resi)
    elif isinstance(chain, str):
        for resi in ph.pose_range(pose):
            if pdb_info.chain(resi) != chain:
                continue
            residue = pose.residue(resi)
            return residue.chain()
    else:
        raise TypeError(f'Chain can be int or str not {type(chain)}')
    raise ValueError(f'No chain {chain}')

