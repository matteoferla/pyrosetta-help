import pyrosetta
import sys

if sys.version_info < (3, 8):
    from typing_extensions import TypedDict
else:
    from typing import TypedDict
from typing import Dict, List
from ..common_ops import pose_range

BondDataType = TypedDict('BondDataType',
                         {'distance':     float,
                          'energy':       float,
                          'acc_resi':     int,
                          'acc_resn':     str,
                          'acc_atm_name': str,
                          'don_resi':     int,
                          'don_resn':     str,
                          'don_atm_name': str,
                          })

def get_hbond_dicts(pose: pyrosetta.Pose) -> Dict[int, List[BondDataType]]:
    hbondset = pose.get_hbonds()
    return {i: [hbond2dict(pose, bond) for bond in hbondset.residue_hbonds(i)] for i in pose_range(pose)}


def hbond2dict(pose: pyrosetta.Pose, bond: pyrosetta.rosetta.core.scoring.hbonds.HBond) -> BondDataType:
    return BondDataType(distance=bond.get_HAdist(pose),
                        energy=bond.energy(),
                        **_get_acceptor(pose, bond),
                        **_get_donor(pose, bond))


def _parse_part(pose: pyrosetta.Pose, resi: int, atom: int, prefix: str = '') -> dict:
    residue: pyrosetta.rosetta.core.conformation.Residue = pose.residue(resi)
    resn: str = residue.name3()
    atm_name: str = residue.atom_name(atom).strip()
    data = dict(resi=resi, resn=resn, atm_name=atm_name)
    return {f'{prefix}_{k}': v for k, v in data.items()} if prefix else data


def _get_acceptor(pose: pyrosetta.Pose, bond: pyrosetta.rosetta.core.scoring.hbonds.HBond) -> dict:
    return _parse_part(pose, resi=bond.acc_res(), atom=bond.acc_atm(), prefix='acc')


def _get_donor(pose: pyrosetta.Pose, bond: pyrosetta.rosetta.core.scoring.hbonds.HBond) -> dict:
    return _parse_part(pose, resi=bond.don_res(), atom=bond.don_hatm(), prefix='don')
