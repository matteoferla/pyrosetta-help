from functools import partial
from typing import (Callable, List, Dict)
import pyrosetta

from .ss import get_ss
from ..common_ops import pose_range

class _BetaTurnPresent:
    """
    Was formerly

        beta_turn_present: Callable[[int], bool] = partial(btd.beta_turn_present, pose)

    But some ligands were causing issues.
    :return:
    """
    def __init__(self, pose):
        self.btd = pyrosetta.rosetta.protocols.features.BetaTurnDetection()
        self.pose=pose

    def __call__(self, resi) -> bool:
        try:
            return self.btd.beta_turn_present(self.pose, resi)
        except RuntimeError:
            return False

def get_betaturns(pose: pyrosetta.Pose) -> Dict[int, str]:
    """
    There is a wee problem in that alpha-helices get classified as beta-turns
    And I am not sure all are beta-turns...
    """

    get_ss(pose)  # fill secstruct
    loop_iter = filter(lambda i: pose.secstruct(i) == 'L', pose_range(pose, protein_only=True))
    maxed = list(loop_iter)[:-3]
    beta_iter = filter(_BetaTurnPresent(pose), maxed)
    btd = pyrosetta.rosetta.protocols.features.BetaTurnDetection()
    return {i: btd.beta_turn_type(pose, i) for i in beta_iter}