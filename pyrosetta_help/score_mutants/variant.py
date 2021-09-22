from __future__ import annotations

import pyrosetta
import re, os, csv, json
from typing import Optional, List, Dict, Union, Any, Callable, Tuple
from .mutation import Mutation


class MutantScorer:
    """
    Copy pasted from PI4KA <- GNB2 <- SnoopCatcher
    """

    # ===== init ===========================================================================

    def __init__(self,
                 pose: pyrosetta.Pose,
                 modelname: str,
                 scorefxn: Optional[pyrosetta.ScoreFunction] = None,
                 strict_about_starting_residue: bool = True,
                 verbose: bool = False):
        self.pose = pose
        self.modelname = modelname
        self.strict_about_starting_residue = bool(strict_about_starting_residue)
        if scorefxn is None:
            self.scorefxn = pyrosetta.get_fa_scorefxn()
        else:
            self.scorefxn = scorefxn
        self.verbose = verbose
        self.output_folder = 'variants'

    @classmethod
    def from_file(cls, filename: str, params_filenames: Optional[List[str]] = None, **kwargs):
        return cls(pose=cls._load_pose_from_file(filename, params_filenames), **kwargs)

    @classmethod
    def _load_pose_from_file(cls, filename: str, params_filenames: Optional[List[str]] = None) -> pyrosetta.Pose:
        """
        Loads a pose from filename with the params in the params_folder
        :param filename:
        :param params_filenames:
        :return:
        """
        pose = pyrosetta.Pose()
        if params_filenames:
            params_paths = pyrosetta.rosetta.utility.vector1_string()
            params_paths.extend(params_filenames)
            pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)
        return pose

    # ===== score =========================================================================
    def score_mutation(self,
                       mutation_name: str,
                       chains: str,
                       distance: int,
                       cycles: int,
                       interfaces,
                       ref_interface_dG: Dict[str, float],
                       final_func: Optional[Callable] = None,
                       preminimise: bool = False,
                       movement: bool = False) -> Tuple[Dict[str, float], pyrosetta.Pose, pyrosetta.Pose]:
        """
        Scores the mutation ``mutation_name`` (str or Mutation instance)
        returning three objects: a dict of scores, the wt (may differ from pose if preminimise=True) and mutant pose

        :param mutation_name:
        :param chains:
        :param distance:
        :param cycles:
        :param interfaces:
        :param ref_interface_dG: premade if no preminimise.
        :param final_func:
        :param preminimise:
        :return:
        """
        if preminimise:
            premutant = self.pose.clone()
        else:
            premutant = self.pose
        assert len(chains) > 0, 'Specify at least one chain e.g. "A"'
        for chain in chains:
            mutation = self.parse_mutation(mutation_name, chain)
            if self.verbose:
                print(mutation)
            if not self.does_contain(mutation):
                raise ValueError('Absent')
            if preminimise:
                self.relax_around_mover(premutant,
                                        mutation=mutation,
                                        distance=distance,
                                        cycles=cycles)
                if self.verbose:
                    print('preminimisation complete')
        variant = premutant.clone()
        for chain in chains:
            mutation = self.parse_mutation(mutation_name, chain)
            variant = self.make_mutant(variant,
                                       mutation=mutation,
                                       distance=distance,
                                       cycles=cycles,
                                       inplace=True)
        if self.verbose:
            print('mutant made')
        variant.dump_scored_pdb(f'{self.output_folder}/{self.modelname}.{mutation}.pdb', self.scorefxn)
        # score proper
        data = self.score_only(variant=variant,
                               reference=premutant,
                               mutation=mutation,
                               chains=chains,
                               distance=distance,
                               interfaces=interfaces,
                               ref_interface_dG=ref_interface_dG if not preminimise else dict(),
                               final_func=final_func,
                               movement=movement)
        # if ref_interface_dG (above) is empty score_only calcutates it, so it makes no diff why its empty.
        # ie. the ref_interface_dG variable here is empty or an empty argument is passed.
        return data, premutant, variant

    def score_only(self,
                   variant: pyrosetta.Pose,
                   reference: pyrosetta.Pose,
                   mutation: Mutation,
                   chains: str,
                   distance: int,
                   interfaces: List[Tuple[str, str]],
                   ref_interface_dG: Dict,
                   final_func: Optional[Callable] = None,
                   movement: bool=False) -> dict:
        unweighted_scorefxn = self.get_unweighted_scorefxn()
        wt_score = unweighted_scorefxn(reference)
        mut_score = unweighted_scorefxn(variant)
        neigh_vector = self.get_neighbour_vector(pose=reference,
                                                 resi=mutation.pose_resi,
                                                 chain=None,  # if chain is present, it assumes resi is pdb number.
                                                 distance=distance)
        neigh_wt_score = unweighted_scorefxn.get_sub_score(reference, neigh_vector)
        neigh_mut_score = unweighted_scorefxn.get_sub_score(variant, neigh_vector)
        data = {'model': self.modelname,
                'mutation': str(mutation),
                'complex_ddG': mut_score - wt_score,
                'complex_native_dG': wt_score,
                'complex_mutant_dG': mut_score,
                'complex_ddG_neigh_only': neigh_mut_score - neigh_wt_score,
                'complex_native_dG_only': neigh_wt_score,
                'complex_mutant_dG_only': neigh_mut_score,
                'FA_RMSD': self.FA_RMSD(self.pose,
                                        variant,
                                        resi=mutation.pose_resi,
                                        chain=None,   # None becase pose_resi is provided.
                                        distance=distance),
                'CA_RMSD': self.CA_RMSD(self.pose,
                                        variant,
                                        resi=mutation.pose_resi,
                                        chain=None,  # None because pose_resi is provided.
                                        distance=distance)
                }
        # interfaces
        for interface_name, interface_scheme in interfaces:
            if self.has_interface(variant, interface_scheme):
                if interface_name not in ref_interface_dG:
                    ref_interface_dG[interface_name] = self.score_interface(reference,
                                                                            interface_scheme)['interface_dG']
                if self.verbose:
                    print(f'{interface_name} ({interface_scheme}) applicable to {self.modelname}')
                i = self.score_interface(variant, interface_scheme)['interface_dG']
            else:
                print(f'{interface_name} ({interface_scheme}) not applicable to {self.modelname}')
                i = float('nan')
            data[f'{interface_name}_interface_native_dG'] = ref_interface_dG[interface_name]
            data[f'{interface_name}_interface_mutant_dG'] = i
            data[f'{interface_name}_interface_ddG'] = i - ref_interface_dG[interface_name]
        if self.verbose:
            print('interface scored')
        # raw
        wt_scoredex = self.get_wscoredict(reference)
        mut_scoredex = self.get_wscoredict(variant)
        delta_scoredex = self.delta_scoredict(mut_scoredex, wt_scoredex)
        if self.verbose:
            print('scores stored')
        # movement
        if movement:
            data['wt_rmsd'] = self.movement(original=reference, resi=mutation.pdb_resi, chain=chains[0], distance=distance)
            data['mut_rmsd'] = self.movement(original=reference, resi=mutation.pdb_resi, chain=chains[0], distance=distance)
            data['ratio_rmsd'] = data['mut_rmsd'] / data['wt_rmsd']
            if self.verbose:
                print('movement assessed')
        data = {**data,
                **self.prefix_dict(wt_scoredex, 'wt'),
                **self.prefix_dict(mut_scoredex, 'mut'),
                **self.prefix_dict(delta_scoredex, 'delta')}
        if final_func is not None:  # callable
            final_func(data, reference, variant)
            if self.verbose:
                print('extra step done.')
        return data

    def make_output_folder(self):
        if not os.path.exists(self.output_folder):
            os.mkdir(self.output_folder)

    def get_unweighted_scorefxn(self):
        unweighted_scorefxn = self.scorefxn.clone()
        ST = pyrosetta.rosetta.core.scoring.ScoreType
        unweighted_scorefxn.set_weight(ST.atom_pair_constraint, 0)
        unweighted_scorefxn.set_weight(ST.angle_constraint, 0)
        unweighted_scorefxn.set_weight(ST.coordinate_constraint, 0)
        unweighted_scorefxn.set_weight(ST.dihedral_constrain, 0)
        return unweighted_scorefxn


    def score_mutations(self,
                        mutations,
                        chains='A',
                        interfaces=(),  # list of two: name, scheme
                        preminimise=False,
                        distance=10,
                        cycles=5,
                        final_func: Optional[Callable] = None) -> List[Dict[str, Union[float, str]]]:
        self.make_output_folder()
        data = []
        ## wt
        ref_interface_dG = {}
        scores = {}  # not written to csv file.
        if not preminimise:
            # touch the energies:
            self.get_unweighted_scorefxn()(self.pose)
            for interface_name, interface_scheme in interfaces:
                ref_interface_dG[interface_name] = self.score_interface(self.pose, interface_scheme)['interface_dG']
        else:
            pass  # calculated for each.
        ## muts
        for mutation_name in mutations:
            try:
                datum, premutant, mutant = self.score_mutation(mutation_name=mutation_name,
                                                               chains=chains,
                                                               distance=distance,
                                                               cycles=cycles,
                                                               interfaces=interfaces,
                                                               ref_interface_dG=ref_interface_dG,
                                                               final_func=final_func,
                                                               preminimise=preminimise)
            except Exception as error:
                msg = f"{error.__class__.__name__}: {error}"
                print(msg)
                datum = {'model': self.modelname,
                         'mutation': str(mutation_name),
                         'complex_ddG': msg
                         }
            data.append(datum)
        return data

    # ===== mutation ======================================================================

    def parse_mutation(self, mutation: Union[str, Mutation], chain, pose: pyrosetta.Pose = None):
        if pose is None:
            pose = self.pose
        if isinstance(mutation, str):
            mutant = Mutation(mutation, chain, pose)
        elif isinstance(mutation, Mutation):
            mutant = mutation
        else:
            raise TypeError(f'Does not accept mutation of type {mutation.__class__.__name__}')
        if mutant.pose_resi == 0:
            raise ValueError('Not in pose')
        if self.strict_about_starting_residue:
            mutant.assert_valid()
        return mutant

    def make_mutant(self,
                    pose: pyrosetta.Pose,
                    mutation: Union[str, Mutation],
                    chain='A',
                    distance: int = 10,
                    cycles: int = 5,
                    inplace: bool = False,
                    ) -> pyrosetta.Pose:
        """
        Make a point mutant (``A23D``).

        :param pose: pose
        :param mutation:
        :param chain:
        :param inplace:
        :return: a copy if inplace is false
        """
        if pose is None and inplace:
            raise ValueError('inplace parameter = True works only if pose is provide')
        elif pose is None:
            mutant = self.pose.clone()
        elif not inplace:
            mutant = pose.clone()
        else:
            mutant = pose
        if isinstance(mutation, str):
            mutation = Mutation(mutation, chain, mutant)
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=mutation.pose_resi, new_res=mutation.to_resn3).apply(mutant)
        self.relax_around_mover(mutant,
                                mutation=mutation,
                                distance=distance,
                                cycles=cycles,
                                own_chain_only=False)
        return mutant

    def does_contain(self, mutation: Union[Mutation, str], chain: Optional[str] = None) -> bool:
        assert isinstance(mutation, Mutation) or chain is not None, 'mutation as str requires chain.'
        if isinstance(mutation, str):
            mutation = Mutation(mutation, chain, self.pose)
        if mutation.pose_resi == 0:
            return False
        if mutation.pose_resn1 not in mutation._name3.keys():
            return False
        else:
            return True

    def relax_around_mover(self,
                           pose: pyrosetta.Pose,
                           mutation: Optional[Mutation] = None,
                           resi: int = None, chain: str = None, cycles=5, distance=5,
                           own_chain_only=False) -> None:
        """
        Relaxes pose ``distance`` around resi:chain or mutation

        :param resi: PDB residue number.
        :param chain:
        :param pose:
        :param cycles: of relax (3 quick, 15 thorough)
        :param distance:
        :return:
        """
        if mutation is None and resi is None:
            raise ValueError('mutation or resi+chain required')
        elif mutation is not None:
            resi = mutation.pose_resi
            chain = None  # chain None is for pose resi
        else:
            pass
        if pose is None:
            pose = self.pose
        movemap = pyrosetta.MoveMap()
        ####
        n = self.get_neighbour_vector(pose=pose, resi=resi, chain=chain, distance=distance,
                                      own_chain_only=own_chain_only)
        # print(pyrosetta.rosetta.core.select.residue_selector.ResidueVector(n))
        movemap.set_bb(False)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(False)
        movemap.set_chi(allow_chi=n)
        movemap.set_jump(False)

        relax = pyrosetta.rosetta.protocols.relax.FastRelax(self.scorefxn, cycles)
        relax.set_movemap(movemap)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        if self.scorefxn.get_weight(pyrosetta.rosetta.core.scoring.ScoreType.cart_bonded) > 0:
            # it's cartesian!
            relax.cartesian(True)
            relax.minimize_bond_angles(True)
            relax.minimize_bond_lengths(True)
        else:
            relax.cartesian(False)
        relax.apply(pose)

    def get_neighbour_vector(self, pose: pyrosetta.Pose, resi: int, chain: str, distance: int,
                             include_focus_in_subset: bool = True,
                             own_chain_only: bool = False) -> pyrosetta.rosetta.utility.vector1_bool:
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        if chain is None:  # pose numbering.
            resi_sele.set_index(resi)
        else:
            resi_sele.set_index(pose.pdb_info().pdb2pose(chain=chain, res=resi))
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        neigh_sele = NeighborhoodResidueSelector(resi_sele, distance=distance,
                                                 include_focus_in_subset=include_focus_in_subset)
        if own_chain_only and chain is not None:
            chain_sele = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain)
            and_sele = pyrosetta.rosetta.core.select.residue_selector.AndResidueSelector(neigh_sele, chain_sele)
            return and_sele.apply(pose)
        else:
            return neigh_sele.apply(pose)

    # ===== interface ======================================================================

    def score_interface(self, pose: pyrosetta.Pose, interface: str) -> Dict[str, float]:
        if pose is None:
            pose = self.pose
        assert self.has_interface(pose, interface), f'There is no {interface}'
        ia = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(interface)
        ia.apply(pose)
        return {'complex_energy': ia.get_complex_energy(),
                'separated_interface_energy': ia.get_separated_interface_energy(),
                'complexed_sasa': ia.get_complexed_sasa(),
                'crossterm_interface_energy': ia.get_crossterm_interface_energy(),
                'interface_dG': ia.get_interface_dG(),
                'interface_delta_sasa': ia.get_interface_delta_sasa()}

    def has_interface(self, pose: pyrosetta.Pose, interface: str) -> bool:
        have_chains = self.get_present_chains(pose)
        want_chains = set(interface.replace('_', ''))
        return have_chains == want_chains

    def get_present_chains(self, pose: Optional[pyrosetta.Pose] = None):
        if pose is None:
            pose = self.pose
        # pose2pdb = pose.pdb_info().pose2pdb
        # return {pose2pdb(r).split()[1] for r in range(1, pose.total_residue() + 1)}
        pdb_info = pose.pdb_info()
        return {pdb_info.chain(res + 1) for res in range(pose.total_residue())}

    def has_residue(self, pose: pyrosetta.Pose, resi: int, chain: str) -> bool:
        if pose is None:
            pose = self.pose
        pdb2pose = pose.pdb_info().pdb2pose
        r = pdb2pose(res=resi, chain=chain)
        return r != 0

    # ========= RMSD relative to native ===========================================

    def vector2list(self, vector: pyrosetta.rosetta.utility.vector1_bool) -> pyrosetta.rosetta.std.list_unsigned_long_t:
        rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(vector)
        x = pyrosetta.rosetta.std.list_unsigned_long_t()
        assert len(rv) > 0, 'Vector is empty!'
        for w in rv:
            x.append(w)
        return x

    def CA_RMSD(self, poseA: pyrosetta.Pose, poseB: pyrosetta.Pose, resi: int, chain: Union[str, None], distance: int) -> float:
        n = self.get_neighbour_vector(pose=poseA, resi=resi, chain=chain, distance=distance)
        residues = self.vector2list(n)
        return pyrosetta.rosetta.core.scoring.CA_rmsd(poseA, poseB, residues)

    def FA_RMSD(self, poseA: pyrosetta.Pose, poseB: pyrosetta.Pose, resi: int, chain: Union[str, None], distance: int) -> float:
        n = self.get_neighbour_vector(pose=poseA, resi=resi, chain=chain, distance=distance,
                                      include_focus_in_subset=False)
        residues = self.vector2list(n)
        # pyrosetta.rosetta.core.scoring.automorphic_rmsd(residueA, residueB, False)
        return pyrosetta.rosetta.core.scoring.all_atom_rmsd(poseA, poseB, residues)

    # ======= score dictionary operations =================================================
    def get_scoredict(self, pose: pyrosetta.Pose) -> Dict[str, float]:
        """
        Given a pose get the global scores.
        """
        a = pose.energies().total_energies_array()
        return dict(zip(a.dtype.fields.keys(), a.tolist()[0]))

    def get_wscoredict(self, pose: pyrosetta.Pose) -> Dict[str, float]:
        scoredex = self.get_scoredict(pose)
        stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        return {k: scoredex[k] * self.scorefxn.get_weight(stm.score_type_from_name(k)) for k in scoredex}

    def delta_scoredict(self, minuend: Dict[str, float], subtrahend: Dict[str, float]) -> Dict[str, float]:
        """
        minuend - subtrahend = difference
        given two dict return the difference -without using pandas.
        """
        assert all([isinstance(v, float) for v in [*minuend.values(), *subtrahend.values()]]), 'Float only, please...'
        minuend_keys = set(minuend.keys())
        subtrahend_keys = set(subtrahend.keys())
        common_keys = minuend_keys & subtrahend_keys
        minuend_unique_keys = minuend_keys - subtrahend_keys
        subtrahend_unique_keys = minuend_keys - subtrahend_keys
        return {**{k: minuend[k] - subtrahend[k] for k in common_keys},  # common
                **{k: minuend[k] - 0 for k in minuend_unique_keys},  # unique
                **{k: 0 - subtrahend[k] for k in subtrahend_unique_keys}  # unique
                }

    def prefix_dict(self, dex: Dict[str, Any], prefix: str) -> Dict[str, Any]:
        return {f'{prefix}_{k}': v for k, v in dex.items()}

    # ====== movement ==============================================================================

    def movement(self, original: pyrosetta.Pose, resi: int, chain: str, distance: int,
                 trials: int = 50, temperature: int = 1.0, replicate_number: int = 10):
        """
        This method adapted from a notebook of mine, but not from an official source, is not well written.
        It should be a filter and score combo.

        It returns the largest bb_rmsd of the pdb residue resi following backrub.
        """
        # this code is experimental

        n = self.get_neighbour_vector(pose=original, resi=resi, chain=chain, distance=distance,
                                      own_chain_only=False)
        # resi
        if chain is None:  # pose numbering.
            target_res = resi
        else:
            target_res = original.pdb_info().pdb2pose(chain=chain, res=resi)
        # prep
        rv = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(n)
        backrub = pyrosetta.rosetta.protocols.backrub.BackrubMover()
        backrub.set_pivot_residues(rv)
        # https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/GenericMonteCarloMover
        monégasque = pyrosetta.rosetta.protocols.monte_carlo.GenericMonteCarloMover(maxtrials=trials,
                                                                                    max_accepted_trials=trials,
                                                                                    # gen.max_accepted_trials() = 0
                                                                                    task_scaling=5,
                                                                                    # gen.task_scaling()
                                                                                    mover=backrub,
                                                                                    temperature=temperature,
                                                                                    sample_type='low',
                                                                                    drift=True)
        monégasque.set_scorefxn(self.scorefxn)
        # monégasque.add_filter(filters , False , 0.005 , 'low'  , True )
        # define the first 4 atoms (N C CA O)
        am = pyrosetta.rosetta.utility.vector1_unsigned_long(4)
        for i in range(1, 5):
            am[i] = i
        # find most deviant
        best_r = 0
        for i in range(replicate_number):
            variant = original.clone()
            monégasque.apply(variant)
            if monégasque.accept_counter() > 0:
                variant = monégasque.last_accepted_pose()  # pretty sure redundant
                # bb_rmsd is all residues: pyrosetta.rosetta.core.scoring.bb_rmsd(pose, ori)
                r = pyrosetta.rosetta.core.scoring.residue_rmsd_nosuper(variant.residue(target_res),
                                                                        original.residue(target_res), am)
                if r > best_r:
                    best_r = r
        return best_r
