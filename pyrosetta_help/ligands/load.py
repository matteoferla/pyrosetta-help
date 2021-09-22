import pyrosetta, warnings, requests
from ..common_ops.utils import make_blank_pose
from rdkit_to_params import Params, neutralise
from typing import *

__all__ = ['parameterised_pose_from_file', 'parameterised_pose_from_pdbblock', 'get_smiles']


def parameterised_pose_from_file(pdb_filename,
                                 wanted_ligands: Union[List[str], Dict[str, Union[str, None]]] = (),
                                 force_parameterisation: bool = False,
                                 neutralise_params: bool = True,
                                 save_params: bool = True,
                                 overriding_params=()) -> pyrosetta.Pose:
    """
    pose loading, the circutous way to not loose ligand or use PDB_component.
    Assumes all mystery components are PDB residues.
    Works best with ignore_unrecognized_res False

    :param pdb_filename:
    :param wanted_ligands: a list of three letter codes or a dictionary of three letter codes to None or smiles
    :param force_parameterisation:
    :param neutralise_params: protonated for pH 7
    :param save_params:
    :param overriding_params: list of params filenames
    :return:
    """
    with open(pdb_filename, 'r') as fh:
        pdbblock = fh.read()
    return parameterised_pose_from_pdbblock(pdbblock,
                                            wanted_ligands=wanted_ligands,
                                            force_parameterisation=force_parameterisation,
                                            neutralise_params=neutralise_params,
                                            save_params=save_params,
                                            overriding_params=overriding_params)


def parameterised_pose_from_pdbblock(pdbblock: str,
                                     wanted_ligands: Union[List[str], Dict[str, Union[str, None]]] = (),
                                     force_parameterisation: bool = False,
                                     neutralise_params: bool = True,
                                     save_params: bool = True,
                                     overriding_params=()) -> pyrosetta.Pose:
    """
    pose loading, the circutous way to not loose ligand or use PDB_component.
    Assumes all mystery components are PDB residues.
    Works best with ignore_unrecognized_res False

    :param pdb_filename:
    :param wanted_ligands: a list of three letter codes or a dictionary of three letter codes to None or smiles
    :param force_parameterisation:
    :param neutralise_params: protonated for pH 7
    :param save_params:
    :param overriding_params: list of params filenames
    :return:
    """

    if pyrosetta.rosetta.basic.options.get_boolean_option('in:file:load_PDB_components'):
        raise ValueError('load_PDB_components is True. Run ``pyrosetta.pose_from_filename`` then.')
    # check if ignore_unrecognized_res ?
    if isinstance(wanted_ligands, dict):
        needed_ligands = dict(wanted_ligands)
    else:
        needed_ligands = {lig: None for lig in wanted_ligands}
    pose = _prep_pose(pdbblock=pdbblock,
                      ligand_dex=needed_ligands,
                      force_parameterisation=force_parameterisation,
                      neutralise_params=neutralise_params,
                      save_params=save_params,
                      overriding_params=overriding_params)
    assert pose.sequence()
    return pose


def _prep_pose(pdbblock: str,
               ligand_dex: Dict,
               force_parameterisation: bool = False,
               neutralise_params: bool = True,
               save_params: bool = True,
               overriding_params: List = ()):
    """
    Parameterise if needed the needed_ligands and add to the list if there's an additional missing PDB residue.
    """
    pose = make_blank_pose(overriding_params)
    rts = pose.residue_type_set_for_pose()
    try:
        for target_ligand in ligand_dex:
            if not rts.has_name3(target_ligand) or force_parameterisation:
                smiles = ligand_dex[target_ligand]
                if not smiles:
                    smiles = get_smiles(target_ligand)
                params = parameterise(pdb_block=pdbblock,
                                      target_ligand=target_ligand,
                                      smiles=smiles,
                                      neutral=neutralise_params,
                                      save=save_params)
                rts = params.add_residuetype(pose)
                pose.conformation().reset_residue_type_set_for_conf(rts)
                rts = pose.residue_type_set_for_pose()  # playing it ubersafe
                assert rts.has_name3(target_ligand)
        pyrosetta.rosetta.core.import_pose.pose_from_pdbstring(pose, pdbblock)
        # pyrosetta.pose_from_file(pose, pdb_filename)
        return pose
    except RuntimeError as err:
        missing = err.args[0].strip()[-3:]
        warnings.warn(f'adding {missing}')
        ligand_dex[missing] = None
        return _prep_pose(pdbblock, ligand_dex, force_parameterisation, neutralise_params, save_params,
                          overriding_params)


def parameterise(pdb_block: str, target_ligand: str, smiles: str, neutral: bool = True, save: bool = True) -> Params:
    params = Params.from_smiles_w_pdbblock(pdb_block=pdb_block,
                                           smiles=smiles,
                                           proximityBonding=False,
                                           name=target_ligand)
    if neutral:
        mol = neutralise(params.mol)
        params = Params.from_mol(mol)
    if save:
        params.dump(f'{target_ligand}.params')  # safekeeping. Not used.
    return params


def get_smiles(ligand_code: str) -> str:
    """
    Get the smiles of a ligand.
    Remember that PDBe smiles need to charged to pH 7.
    """
    ligand_code = ligand_code.upper()
    ligand_data = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{ligand_code}').json()
    return ligand_data[ligand_code][0]['smiles'][0]['name']
