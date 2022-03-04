from typing import (Optional, Tuple, Iterable, Dict)

import pyrosetta

__all__ = ['get_alignment', 'write_grishin', 'thread', 'rangify', 'steal_ligands',
           'get_nonprotein_pose', 'oligomer_thread', 'make_fragment_sets']


def get_alignment(target: str, template: str) -> Dict[str, str]:
    """
    Returns alignments using ``pairwise2.align.globalxs``
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
           template_name: str = 'template',
           fragment_sets: Optional[pyrosetta.rosetta.utility.vector1_std_shared_ptr_core_fragment_FragSet_t] = None,
           align: Optional[pyrosetta.rosetta.core.sequence.SequenceAlignment] = None
           ) -> Tuple[pyrosetta.Pose,
                      pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover,
                      pyrosetta.rosetta.utility.vector1_bool]:
    """
    Given the target sequence and the optional blank target pose (in case there's some reason for it),
    thread the sequence against the template pose — which is assumed to be a single chain.
    Optionally using fragments from fragment_sets.
    The three outputs are the target_pose, the threader instance and a vector of the residues threaded.

    >>> print(threader.frag_libs()[1].nr_frames())
    
    >>> qt = threader.get_qt_mapping(threaded)
    >>> print(f'Template residue {21} is residue {q[21]} in target')

    :param target_sequence:
    :param template_pose:
    :param target_pose:
    :param target_name:
    :param template_name:
    :param fragment_sets:
    :param align:
    :return:
    """
    ## Make target_pose
    if target_pose is None:
        target_pose = pyrosetta.Pose()
        pyrosetta.rosetta.core.pose.make_pose_from_sequence(target_pose,
                                                            target_sequence,
                                                            'fa_standard')
    ## Unterminate and remove ligands
    original_template_pose = template_pose
    clean_template_pose = template_pose.clone()
    pyrosetta.rosetta.core.pose.remove_nonprotein_residues(clean_template_pose)
    ### find
    lowers = pyrosetta.rosetta.utility.vector1_std_pair_unsigned_long_protocols_sic_dock_Vec3_t()
    uppers = pyrosetta.rosetta.utility.vector1_std_pair_unsigned_long_protocols_sic_dock_Vec3_t()
    pyrosetta.rosetta.protocols.sic_dock.get_termini_from_pose(clean_template_pose, lowers, uppers)
    ### remove
    rm_upper = pyrosetta.rosetta.core.conformation.remove_upper_terminus_type_from_conformation_residue
    rm_lower = pyrosetta.rosetta.core.conformation.remove_lower_terminus_type_from_conformation_residue
    for upper, _ in uppers:
        rm_upper(clean_template_pose.conformation(), upper)
    for lower, _ in lowers:
        rm_lower(clean_template_pose.conformation(), lower)
    ## Align
    if align is None:
        alignment = get_alignment(target_sequence, clean_template_pose.sequence())
        aln_file = f'{target_name}.aln'
        write_grishin(target_name,
                      alignment['target'],
                      template_name,
                      alignment['template'],
                      aln_file)
        align = pyrosetta.rosetta.core.sequence.read_aln(format='grishin', filename=aln_file)[1]
    ## Thread
    # not pyrosetta.rosetta.protocols.comparative_modeling.PartialThreadingMover
    threader = pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover(align=align,
                                                                               template_pose=clean_template_pose)
    ### Input frags
    if fragment_sets is not None:
        threader.build_loops(True)
        threader.randomize_loop_coords(True)  # default
        threader.frag_libs(fragment_sets)
    ### Apply
    threader.apply(target_pose)
    ## Steal sidechains
    qt = threader.get_qt_mapping(target_pose)
    steal = pyrosetta.rosetta.protocols.comparative_modeling.StealSideChainsMover(clean_template_pose, qt)
    steal.apply(target_pose)
    ## Make vector of what was changed
    vector = pyrosetta.rosetta.utility.vector1_bool(target_pose.total_residue())
    mapping = qt.mapping()
    for r in range(1, target_pose.total_residue() + 1):
        if mapping[r] != 0:
            vector[r] = 1
    return target_pose, threader, vector


vector1_FragSet_t = pyrosetta.rosetta.utility.vector1_std_shared_ptr_core_fragment_FragSet_t


def make_fragment_sets(*poses: Iterable[pyrosetta.Pose], lengths: Iterable[int] = (3,)) -> vector1_FragSet_t:
    fragsets = vector1_FragSet_t(len(lengths))
    for i, l in enumerate(lengths):
        fragsets[i+1] = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(l)
        for pose in poses:
            pyrosetta.rosetta.core.fragment.steal_constant_length_frag_set_from_pose(pose, fragsets[i+1])
    return fragsets

def steal_ligands(donor_pose, acceptor_pose) -> None:
    """
    Steals non-Protein residues from donor_pose and adds them to acceptor_pose

    Do not use with nucleic acid polymers.

    :param donor_pose:
    :param acceptor_pose:
    :return:
    """
    PROTEIN = pyrosetta.rosetta.core.chemical.ResidueProperty.PROTEIN
    prot_sele = pyrosetta.rosetta.core.select.residue_selector.ResiduePropertySelector(PROTEIN)
    not_sele = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(prot_sele)
    rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(not_sele.apply(donor_pose))
    # if it were DNA...
    # for from_res, to_res in rangify(rv):
    #     pyrosetta.rosetta.core.pose.append_subpose_to_pose(acceptor_pose, donor_pose, from_res, to_res, True)
    for res in rv:
        pyrosetta.rosetta.core.pose.append_subpose_to_pose(acceptor_pose, donor_pose, res, res, True)


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
