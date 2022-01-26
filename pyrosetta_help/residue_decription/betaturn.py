from functools import partial
from typing import (Callable, List, Dict)
import pyrosetta

from .ss import get_ss
from ..common_ops import pose_range

def get_betaturns(pose: pyrosetta.Pose) -> Dict[int, str]:
    """
    There is a wee problem in that alpha-helices get classified as beta-turns
    And I am not sure all are beta-turns...
    """

    get_ss(pose)  # fill secstruct
    btd = pyrosetta.rosetta.protocols.features.BetaTurnDetection()
    beta_turn_present: Callable[[int], bool] = partial(btd.beta_turn_present, pose)
    loop_iter = filter(lambda i: pose.secstruct(i) == 'L', pose_range(pose, protein_only=True))
    maxed = list(loop_iter)[:-3]
    beta_iter = filter(beta_turn_present, maxed)
    return {i: btd.beta_turn_type(pose, i) for i in beta_iter}