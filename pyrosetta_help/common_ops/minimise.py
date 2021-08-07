from typing import *
import pyrosetta
import os

__all__ = ['get_local_scorefxn',
           'prep_ED',
           'get_local_relax',
           'do_local_relax',
           'do_chainwise_relax']

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
    if os.path.splitext(map_filename) == '.gz':
        raise Exception('The file is gzipped')
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