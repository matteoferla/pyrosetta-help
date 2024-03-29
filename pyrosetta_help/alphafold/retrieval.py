__all__ = ['pose_from_alphafold2', 'get_alphafold2_error', 'reshape_errors']

from typing import (Union, Dict, List)

import numpy as np
import pyrosetta
import requests


def pose_from_alphafold2(uniprot: str) -> pyrosetta.Pose:
    """
    Returns a pose from the alphafold2 server.

    :param uniprot: uniprot id ("accession"), not gene name or uniprot name.
    :return:
    """
    reply = requests.get(f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v2.pdb')
    reply.raise_for_status()
    pdbblock = reply.text
    pose = pyrosetta.Pose()
    pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdbblock)
    return pose


# ==== errors based methods ====================================================
def get_alphafold2_error(uniprot: str, reshaped=True) -> Union[np.ndarray, list]:
    """
    Returns the distances errors either as numpy matrix (``reshaped=True``)
    or as the weird format from AF2-EBI —see ``help(pyrosetta_help.alphafold.retrieval.reshape_errors)`` for more.

    Remember that the matrix is zero indexed and that these values are in Ångström
    and are not pLDDT, which are stored as b-factors.
    """
    # https://alphafold.ebi.ac.uk/files/AF-Q00341-F1-predicted_aligned_error_v1.json
    reply = requests.get(f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-predicted_aligned_error_v2.json')
    reply.raise_for_status()
    errors = reply.json()
    if reshaped:
        return reshape_errors(errors)
    else:
        return errors


def reshape_errors(errors: List[Dict[str, list]]) -> np.array:
    """
    The JSON from AF2 has a single element list.
    the sole element is a dictionary with keys 'residue1', 'residue2' and 'distance'.
    This method returns a matrix of distances reshaped based on the stated residue indices.
    This is rather unlikely to differ from a regular reshape...
    but idiotically I am not taking changes assuming it is always sorted.
    """
    n_residues = int(np.sqrt(len(errors[0]['distance'])))
    error_matrix = np.zeros((n_residues, n_residues)) * np.nan
    for i, d in enumerate(errors[0]['distance']):
        error_matrix[errors[0]['residue1'][i] - 1, errors[0]['residue2'][i] - 1] = d
    return error_matrix
