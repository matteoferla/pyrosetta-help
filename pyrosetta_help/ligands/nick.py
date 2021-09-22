import pyrosetta

try:
    from rdkit_to_params import Params, neutralise
except ModuleNotFoundError:
    pass
from typing import *
from ..common_ops.downloads import download_pdb
from ..common_ops.utils import make_blank_pose
from .load import parameterised_pose_from_file

import warnings, requests

pr_rs = pyrosetta.rosetta.core.select.residue_selector


def chain_letter_to_number(letter, pose):
    warnings.warn('this is a shitty interim way of going from letter to number...' + \
                  ' I am sure there"s a method somewhere that does this.')
    pdb_info = pose.pdb_info()
    for i in range(1, pose.num_chains() + 1):
        b = pose.chain_begin(i)
        if pdb_info.chain(b) == letter:
            return i
    else:
        raise ValueError


class LigandNicker:
    """
    Given a pdb_file (regular initialisation) or code (``.from_pdbcode`` classmethod)
    and a list of 3-letter codes of wanted residues (``wanted_ligands`` argument),
    it loads it as Pyrosetta Pose (``donor_pose``), ready for ``migrate`` to loads the acceptor_pose and nick
    the residues that are wanted.

    If ``force_parameterisation`` is on it or the residue is novel it parameterises it.
    """

    def __init__(self,
                 pdb_filename: str,
                 chain: str,
                 wanted_ligands: List[str],
                 force_parameterisation: bool = False,
                 neutralise_params: bool = True,
                 save_params: bool = True,
                 overriding_params=()):
        """
        Initialisation loads the donor. ``migrate`` loads the acceptor.

        :param pdb_filename:
        :param chain:
        :param wanted_ligands:
        :param force_parameterisation:
        :param neutralise_params: pH 7 protonation
        :param save_params:
        :param overriding_params: overide paramaterisation and use provide params
        """
        self.pdb_filename = pdb_filename
        self.donor_chain = chain
        self.wanted_ligands = wanted_ligands
        self.extra_params = []
        self.donor_pose = parameterised_pose_from_file(pdb_filename=pdb_filename,
                                                       force_parameterisation=force_parameterisation,
                                                       neutralise_params=neutralise_params,
                                                       save_params=save_params,
                                                       overriding_params=overriding_params)
        self.acceptor_pose = None

    @classmethod
    def from_pdbcode(cls, pdb_code: str, chain: str, *args, **kvargs):
        pdb_filename = download_pdb(pdb_code)
        return cls(pdb_filename=pdb_filename, chain=chain, *args, **kvargs)

    def migrate(self, acceptor_pose: pyrosetta.Pose, acceptor_chain: str = 'A',
                constrained: bool = True, relaxed: bool = True,
                relax_radius: int = 20, relax_cycles: int = 3):
        """
        The acceptor pose is the non-empty pose.

        This method aligns the sequences of the acceptor and donor pose.
        It finds the mapping of the neighbourhood of the wanted residues of the donor_pose
        It superimposes the poses by those residues.
        It adds the residues that need nicking.
        It adds constraints (optionally) based on the hydrogen bonding of the residues around the wanted residues.
        onto the ``.acceptor_pose``. To check:

        >>> print( len(self.acceptor_pose.constraint_set().get_all_constraints())  )

        It then optionally relaxes the neighbourhood.

        :param acceptor_pose:
        :param acceptor_chain:
        :param constrained:
        :param relaxed:
        :param relax_radius:
        :param relax_cycles:
        :return:
        """
        # store
        self.acceptor_pose = acceptor_pose
        self.acceptor_chain = acceptor_chain
        # get wanted residues
        wanted_selector = self.get_wanted_selector()
        donor_neighbors = self.get_surrounding_residue(self.donor_pose,
                                                       self.donor_chain,
                                                       wanted_selector)
        # get the neighbouring residues & superpose w/ them
        dex = self.get_mapping_between_poses(donor_neighbors)
        mapping = self.make_atomID_map(dex, self.donor_pose, self.acceptor_pose)
        pyrosetta.rosetta.core.scoring.superimpose_pose(self.donor_pose,
                                                        self.acceptor_pose,
                                                        mapping)
        wanted_vector = wanted_selector.apply(self.donor_pose)
        # append the ligands. n nucleotides have issues with the grafting command.
        # hence this way.
        added_indices = []
        for r in pr_rs.ResidueVector(wanted_vector):
            r2 = self.acceptor_pose.total_residue() + 1
            dex[r] = r2
            added_indices.append(r2)
            # append_subpose_to_pose is cool with novel residues.
            pyrosetta.rosetta.core.pose.append_subpose_to_pose(self.acceptor_pose,
                                                               self.donor_pose,
                                                               r, r, True)
        # make selector of newly added residues
        self.added_selector = pr_rs.ResidueIndexSelector()
        for r in added_indices:
            self.added_selector.append_index(r)
        # constrain
        if constrained:
            self.constrain_migrated(wanted_vector, dex)
        # relax
        if relaxed:
            self.relax_migrated(distance=relax_radius, cycles=relax_cycles)

    # ---- migrate dependent methods

    def get_wanted_selector(self):
        # select the ligands that are wanted.
        # set_residue_names does not like custom residue types.
        wanted_ligands = {c.rjust(3) for c in self.wanted_ligands}
        # ResidueNameSelector seems to not like custom residues if they arent in the pose
        try:
            resn_sele = pr_rs.ResidueNameSelector()
            resn_sele.set_residue_name3(','.join(wanted_ligands))
            resn_sele.apply(self.donor_pose)
        except RuntimeError as err:
            raise ValueError(f'Residue not in the pose. {err}')
        # end of shitty hack.
        # filter for those within 3A of chain of interest
        # define extended neighbours
        chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(self.donor_chain)
        cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
        cc_sele.central_residue_group_selector(chain_sele)
        cc_sele.threshold(3)
        # intersection:
        and_sele = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(cc_sele, resn_sele)
        return and_sele

    def get_surrounding_residue(self, pose, chain_filter, wanted_selector):
        cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
        cc_sele.central_residue_group_selector(wanted_selector)
        cc_sele.threshold(3)
        prop_sele = pr_rs.ResiduePropertySelector()
        prop_sele.add_property(pyrosetta.rosetta.core.chemical.ResidueProperty.ALPHA_AA)
        neigh_sele = pr_rs.AndResidueSelector(prop_sele, cc_sele)
        neigh_vector = neigh_sele.apply(pose)
        pdb_info = pose.pdb_info()
        rs = pr_rs.ResidueVector(neigh_vector)
        neigh_resis = [r for r in rs if pdb_info.chain(r) == chain_filter]
        return neigh_resis  # these are the pose indices

    def _make_map(self, al, pose_offset=0):
        # pose index (fortran-style) to MSA index (C++-style)
        gap_map = [i for i, r in enumerate(al) if r != '-']
        return {pose_offset + ungap_i: gap_i for ungap_i, gap_i in enumerate(gap_map)}

    def get_mapping_between_poses(self, index_list: List[int]) -> Dict[int, int]:
        """
        Given a list of indices for one pose, return their aligned equivalents in the second.
        Modded to be donor_pose --> acceptor_pose
        """
        from Bio import pairwise2
        acc_ch_idx = chain_letter_to_number(self.acceptor_chain, self.acceptor_pose)
        don_ch_idx = chain_letter_to_number(self.donor_chain, self.donor_pose)
        acc_seq = self.acceptor_pose.chain_sequence(acc_ch_idx)
        don_seq = self.donor_pose.chain_sequence(don_ch_idx)
        alignments = pairwise2.align.globalxs(acc_seq,
                                              don_seq,
                                              -1,  # open
                                              -0.1  # extend
                                              )
        alignments = dict(zip(['target', 'template', 'score', 'begin', 'end'], alignments[0]))

        donor_pose_to_msa_mapping = self._make_map(alignments['template'],
                                                   self.donor_pose.chain_begin(don_ch_idx))
        acceptor_pose_to_msa_mapping = self._make_map(alignments['target'],
                                                      self.acceptor_pose.chain_begin(acc_ch_idx))
        msa_to_acceptor_pose_mapping = dict(zip(acceptor_pose_to_msa_mapping.values(),
                                                acceptor_pose_to_msa_mapping.keys()))
        mapping = {}
        for r in index_list:
            msa_i = donor_pose_to_msa_mapping[r]
            trans = msa_to_acceptor_pose_mapping[msa_i]
            #             print(r, self.donor_pose.residue(r).name1(),
            #                   msa_i, alignments['template'][msa_i],
            #                   trans, self.acceptor_pose.residue(trans).name1() )
            mapping[r] = trans
        return mapping

    def make_atomID_map(self, dex, query_pose, target_pose):
        mapping = pyrosetta.rosetta.std.map_core_id_AtomID_core_id_AtomID()
        get_id = lambda pose, r, atomname: pyrosetta.AtomID(atomno_in=pose.residue(r).atom_index(atomname), rsd_in=r)
        for q, t in dex.items():
            for atomname in ('CA',):
                mapping[get_id(query_pose, q, atomname)] = get_id(target_pose, t, atomname)
        return mapping

    def constrain_migrated(self, wanted_vector, dex):
        hbond_set = self.donor_pose.get_hbonds()
        wanted_indices = pr_rs.ResidueVector(wanted_vector)
        for hbond in hbond_set.hbonds():
            if hbond.don_res() in wanted_indices or hbond.acc_res() in wanted_indices:
                con = self.make_constraint_foreign_hbond(hbond, dex)
                if con:
                    self.acceptor_pose.add_constraint(con)

    def make_constraint_foreign_hbond(self, hbond, dex: Dict[int, int]) \
            -> Union[None, pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint]:
        # acc = hydrogen acceptor residue w/Acceptor atom
        # don = donor residue w/ Bounded and Hydrogen atoms.
        # ------------- sort out donor -------------
        don_res = hbond.don_res()
        don = self.donor_pose.residue(don_res)
        if don_res not in dex:
            warnings.warn(
                f'The H-bond–donor residue in the donor pose has residue {don_res} with no equivalent - why was it absent in CC selector?')
            return
        trans_don_res = dex[don_res]  # acceptor_pose
        trans_don = self.acceptor_pose.residue(trans_don_res)
        don_name3 = don.name3()
        trans_don_name3 = trans_don.name3()
        if not hbond.don_hatm_is_backbone() and don_name3 != trans_don_name3:
            return  # skip as it is a sidechain Hbond of a changed residue.
        elif not hbond.don_hatm_is_backbone() and don_name3 == 'HIS':
            return  # avoiding HID vs HIE pain
        don_atomname = don.atom_name(hbond.don_hatm())
        donor_atom = pyrosetta.AtomID(atomno_in=trans_don.atom_index(don_atomname),
                                      rsd_in=trans_don_res)
        # ------------- sort out acceptor -------------
        acc_res = hbond.acc_res()
        acc = self.donor_pose.residue(acc_res)
        trans_acc_res = dex[acc_res]  # acceptor_pose
        acc_name3 = acc.name3()
        trans_acc = self.acceptor_pose.residue(trans_acc_res)
        trans_acc_name3 = trans_acc.name3()
        if not hbond.acc_atm_is_backbone() and acc_name3 != trans_acc_name3:
            return  # skip as it is a sidechain Hbond of a changed residue.
        elif not hbond.acc_atm_is_backbone() and acc_name3 == 'HIS':
            return  # avoiding HID vs HIE pain
        acc_atomname = acc.atom_name(hbond.acc_atm())
        acceptor_atom = pyrosetta.AtomID(atomno_in=trans_acc.atom_index(acc_atomname),
                                         rsd_in=trans_acc_res)
        # ------------- make constraint -------------
        HarmonicFunc = pyrosetta.rosetta.core.scoring.func.HarmonicFunc
        AtomPairConstraint = pyrosetta.rosetta.core.scoring.constraints.AtomPairConstraint
        d = hbond.get_HAdist(self.donor_pose)
        return AtomPairConstraint(donor_atom,
                                  acceptor_atom,
                                  HarmonicFunc(x0_in=d, sd_in=0.2))

    def relax_migrated(self, distance: int = 20, cycles: int = 3, atom_pair_weight: int = 5):
        scorefxn = pyrosetta.get_fa_scorefxn()
        atom_pair = pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint
        scorefxn.set_weight(atom_pair, atom_pair_weight)
        neigh_sele = pr_rs.NeighborhoodResidueSelector(self.added_selector, distance=distance,
                                                       include_focus_in_subset=True)
        movemap = pyrosetta.MoveMap()
        n = neigh_sele.apply(self.acceptor_pose)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(allow_chi=n)
        ft = self.acceptor_pose.fold_tree()
        for r in pr_rs.ResidueVector(self.added_selector.apply(self.acceptor_pose)):
            j = ft.get_jump_that_builds_residue(r)
            movemap.set_jump(j, True)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.apply(self.acceptor_pose)
