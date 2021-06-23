from typing import *
import pyrosetta


class NeighbourInteractions:
    """
    Gets the per atom energies for the interactions.

    >>> ni = NeighbourInteractions(pose, 1)  # pose residue index 1
    >>> print(ni.describe_best())
    >>> ni.interactions[(' N  ', 2, ' N  ')]

    Unfortunately, bonding is not taken into account therefore the total is more favourable
    as these are ignored.

    >>> ni.total, ni.expected_total
    """

    score_types = ['fa_atr', 'fa_rep', 'fa_sol', 'fa_elec']
    term_relevance_cutoff = 0.5

    def __init__(self,
                 pose: pyrosetta.Pose,
                 target_idx: int,
                 threshold: int = 3,
                 scorefxn: Optional[pyrosetta.ScoreFunction] = None,
                 weighted: bool = True,
                 halved: bool = False):
        # --- input
        self.pose = pose
        self.target_idx = target_idx
        self.target_residue = self.pose.residue(self.target_idx)
        assert self.target_residue, f'{self.target_idx} not in pose'
        self.threshold = threshold
        if scorefxn is None:
            self.scorefxn = pyrosetta.get_fa_scorefxn()
        else:
            self.scorefxn = scorefxn
        self.weighted = weighted
        stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        get_weights = lambda name: self.scorefxn.get_weight(stm.score_type_from_name(name))
        if not weighted:
            self.weights = {name: 1. for name in self.score_types}
        else:
            self.weights = {name: get_weights(st_name) for name, st_name in
                            zip(self.score_types, ['fa_atr', 'fa_rep', 'fa_sol', 'fa_elec'])}
        self.halved = halved
        # --- derived
        self.neighbours = self._get_neighbours()  # ResidueVector
        self.interactions = self._get_interactions()  # Dict[Tuple(str, int, str), Dict[str, float]]

    def _get_interactions(self) -> Dict[Tuple[str, int, str], Dict[str, float]]:
        interactions = {}
        ratio = 2 if self.halved else 1
        # Iterate per target residue's atom per all other residues' atoms
        for i in range(1, self.target_residue.natoms() + 1):
            if self.target_residue.atom_type(i).element() == 'X':
                continue
            iname = self.target_residue.atom_name(i)
            interactions[iname] = {}
            for r in self.neighbours:
                if r == self.target_idx:
                    continue
                other = self.pose.residue(r)
                interactions[iname][r] = {}
                for o in range(1, other.natoms() + 1):
                    oname = other.atom_name(o)
                    if r == self.target_idx and o == i:
                        continue  # self to self should be zero...!
                    elif other.atom_type(o).element() == 'X':
                        continue
                    score = pyrosetta.toolbox.atom_pair_energy.etable_atom_pair_energies(self.target_residue,
                                                                                         i,
                                                                                         other,
                                                                                         o,
                                                                                         self.scorefxn)
                    # interactions[iname][r][oname] = dict(zip(score_types, score))
                    interactions[iname][r][oname] = {st: s * self.weights[st] / ratio
                                                     for st, s in zip(self.score_types, score)}
        # correct fa_rep for bonded
        for atomname, other_residue_index, other_atomname in self._get_connections():
            # interactions[atomname][other_residue_index][other_atomname]['fa_rep'] = float('nan')
            del interactions[atomname][other_residue_index][other_atomname]
        # reshape
        reshaped = {(target_atomname, other_resi, atomname): interactions[target_atomname][other_resi][atomname] for
                    target_atomname in interactions for other_resi in interactions[target_atomname] for atomname in
                    interactions[target_atomname][other_resi]}
        return reshaped

    def _get_neighbours(self) -> pyrosetta.rosetta.core.select.residue_selector.ResidueVector:
        cc_sele = self.get_cc_selector()
        neighs = pyrosetta.rosetta.core.select.residue_selector.ResidueVector(cc_sele.apply(self.pose))
        return neighs

    def get_cc_selector(self):
        resi_sele = self.get_target_selector()
        cc_sele = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
        cc_sele.central_residue_group_selector(resi_sele)
        cc_sele.threshold(self.threshold)
        return cc_sele

    def get_target_selector(self):
        return pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(self.target_idx)

    def _get_connections(self):
        residue = self.pose.residue(self.target_idx)
        for conn_id in range(1, residue.n_current_residue_connections() + 1):
            atomno = residue.residue_connect_atom_index(conn_id)
            atomname = residue.atom_name(atomno)
            adjecent_atomnos = residue.bonded_neighbor(atomno)
            other_residue_index = residue.residue_connection_partner(conn_id)
            other_residue = self.pose.residue(other_residue_index)
            other_atomname, other_atomno = self._get_other_connecting_atom(conn_id)
            other_adjecent_atomnos = other_residue.bonded_neighbor(other_atomno)
            return [(atomname, other_residue_index, other_atomname), ] + \
                   [(residue.atom_name(no), other_residue_index, other_atomname) for no in adjecent_atomnos] + \
                   [(atomname, other_residue_index, other_residue.atom_name(no)) for no in other_adjecent_atomnos]

    def _get_other_connecting_atom(self, conn_id) -> Tuple[str, int]:
        other_residue_index = self.target_residue.residue_connection_partner(conn_id)
        other_residue = self.pose.residue(other_residue_index)
        for other_conn_id in range(1, other_residue.n_current_residue_connections() + 1):
            if other_residue.residue_connection_partner(other_conn_id) == self.target_idx:
                other_atom_index = other_residue.residue_connect_atom_index(other_conn_id)
                return other_residue.atom_name(other_atom_index), other_atom_index
        else:
            raise ValueError('No connection found?!')


    @property
    def best_interactions(self):
        # dict is ordered...
        get_max = lambda k: max(map(abs, self.interactions[k].values()))
        return {k: self.interactions[k] for k in sorted(self.interactions, key=get_max, reverse=True)
                if get_max(k) >= self.term_relevance_cutoff}

    def describe_atom(self, residue, atomname):
        atomno = residue.atom_index(atomname)
        atomtype = residue.atom_type(atomno)
        verdict = {'heavyatom': atomtype.is_heavyatom(),
                   'polar hydrogen': atomtype.is_polar_hydrogen(),
                   'H-acceptor': atomtype.is_acceptor(),
                   'H-donor': atomtype.is_donor(),
                   'aromatic': atomtype.is_aromatic()}
        return ', '.join([atomtype.atom_type_name(), atomtype.element()] +
                         [k for k, v in verdict.items() if v])

    def describe_interaction(self, target_atomname, other_resi, other_atomname):
        other_residue = self.pose.residue(other_resi)
        target_text = self.describe_atom(self.target_residue, target_atomname)
        other_text = self.describe_atom(other_residue, other_atomname)
        # [('fa_sol', 0.7810698032540209), ('fa_elec', 0.37053877404394236)]
        all_terms = self.interactions[target_atomname, other_resi, other_atomname]
        # sort & filter | > 0.25 |
        terms = sorted(filter(lambda t: abs(t[1]) >= self.term_relevance_cutoff,
                              all_terms.items(),
                              ),
                       key=lambda t: -abs(t[1]))
        term_text = ' + '.join([f'{score_type} ({value:.1f} kcal/mol)' for score_type, value in terms])
        return f'{self.target_idx}.{target_atomname.strip()} ({target_text}) - ' + \
               f'{other_resi}.{other_atomname.strip()} ({other_text}) ' + \
               term_text

    def describe_best(self):
        return '\n'.join([self.describe_interaction(*k) for k in self.best_interactions])

    @property
    def expected_total(self):
        self.scorefxn(self.pose)
        return self.scorefxn.get_sub_score(self.pose, self.get_target_selector().apply(self.pose))

    @property
    def total(self):
        return sum([sum(d.values()) for d in self.interactions.values()])



