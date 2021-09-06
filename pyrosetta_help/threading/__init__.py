from typing import *
import pyrosetta

__all__ = ['get_alignment', 'write_grishin', 'thread', 'rangify', 'get_nonprotein_pose', 'oligomer_thread']


def get_alignment(target: str, template: str) -> Dict[str, str]:
    """
    Returns alignments.
    """
    from Bio import pairwise2
    alignments = pairwise2.align.globalxs(target,
                                          template,
                                          -1,  # open
                                          -0.1  # extend
                                          )
    return dict(zip(['target', 'template', 'score', 'begin', 'end'], alignments[0]))


def write_grishin(target_name, target_sequence, template_name, template_sequence, outfile):
    with open(outfile, 'w') as w:
        w.write(f'## {target_name} {template_name}\n')
        w.write(f'#\n')
        w.write('scores_from_program: 0\n')
        w.write(f'0 {target_sequence}\n')
        w.write(f'0 {template_sequence}\n')
        w.write('--\n')


def thread(target_sequence: str,
           template_pose: pyrosetta.Pose,
           target_pose: Optional[pyrosetta.Pose] = None,
           target_name: str = 'target',
           template_name: str = 'template') -> Tuple[pyrosetta.Pose,
                                                     pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover,
                                                     pyrosetta.rosetta.utility.vector1_bool]:
    """
    Actual operation.

    :param target_sequence:
    :param template_pose:
    :param target_pose:
    :param target_name:
    :param template_name:
    :return:
    """
    alignment = get_alignment(target_sequence, template_pose.sequence())
    aln_file = f'{target_name}.aln'
    write_grishin(target_name,
                  alignment['target'],
                  template_name,
                  alignment['template'],
                  aln_file)
    align = pyrosetta.rosetta.core.sequence.read_aln(format='grishin', filename=aln_file)
    # pyrosetta.rosetta.protocols.comparative_modeling.PartialThreadingMover
    threader = pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover(align=align[1],
                                                                               template_pose=template_pose)
    if target_pose is None:
        target_pose = pyrosetta.Pose()
        pyrosetta.rosetta.core.pose.make_pose_from_sequence(target_pose,
                                                            target_sequence,
                                                            'fa_standard')
    threader.apply(target_pose)
    qt = threader.get_qt_mapping(target_pose)
    steal = pyrosetta.rosetta.protocols.comparative_modeling.StealSideChainsMover(template_pose, qt)
    steal.apply(target_pose)
    # v = AlteredSelector(threader).apply(itpr3_pose)
    vector = pyrosetta.rosetta.utility.vector1_bool(target_pose.total_residue())
    mapping = qt.mapping()
    for r in range(1, target_pose.total_residue() + 1):
        if mapping[r] != 0:
            vector[r] = 1
    return target_pose, threader, vector


def rangify(values):
    """
    Given a list of integers, returns a list of tuples of ranges (interger pairs).

    :param values:
    :return:
    """
    previous = None
    start = None
    ranges = []
    for r in values:
        if previous is None:
            previous = r
            start = r
        elif r == previous + 1:
            pass
        else:  # r != previous + 1
            ranges.append((start, previous))
            start = r
        previous = r
    ranges.append((start, previous))
    return ranges


def get_nonprotein_pose(pose):
    """
    pyrosetta.rosetta.protocols.comparative_modeling.StealLigandMover requires some weird things.
    This does the same, makes a ligand only pose.

    :param pose:
    :return:
    """
    protein_prop = pyrosetta.rosetta.core.chemical.ResidueProperty.PROTEIN
    prot_sele = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(protein_prop)
    np_sele = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(prot_sele)
    rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(np_sele.apply(pose))
    if len(rv) == 0:
        return None
    poses = []
    for begin, end in rangify(rv):
        poses.append(pyrosetta.rosetta.protocols.grafting.return_region(pose, begin, end))
    neo = poses[0]
    for pose in poses[1:]:
        neo.append_pose_by_jump(pose, neo.num_jump() + 1)
    return neo


def oligomer_thread(pose, sequence):
    faux = pose.clone()
    pyrosetta.rosetta.core.pose.remove_nonprotein_residues(faux)
    threaded = None
    threaders = []
    master_vector = pyrosetta.rosetta.utility.vector1_bool()
    for chain_pose in faux.split_by_chain():
        target_pose, threader, vector = thread(target_sequence=sequence, template_pose=chain_pose)
        threaders.append(threader)
        if threaded is None:
            threaded = target_pose
        else:
            threaded.append_pose_by_jump(target_pose, threaded.num_jump() + 1)
            master_vector.extend(vector)
    threaded.append_pose_by_jump(get_nonprotein_pose(pose), threaded.num_jump() + 1)
    return threaded, threaders, master_vector
