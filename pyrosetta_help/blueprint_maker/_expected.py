from typing import List
import pyrosetta
from Bio.Align import PairwiseAligner, PairwiseAlignments, Alignment
from IPython.display import display, HTML

class BlueprinterExpected:
    def expected_seq(self):
        """
        Returns the sequence as expected by the blueprint.
        Calls iteratively ``get_expected_aa_from_row``
        """
        return ''.join([self.get_expected_aa_from_row(row) for row in self])

    def get_expected_aa_from_row(self, row: List[str]) -> str:
        """
        Give a row of the blueprint, return the expected amino acid.
        Called by ``expected_seq``.
        """
        if len(row) == 0: # impossible
            return ''
        elif len(row) == 3 and row[2] == '.': # "30 K ."
            return row[1]
        elif len(row) == 3:  # "30 K H" ???
            return '*'  # H/L/E/D in lowercase would be ambiguous
        elif len(row) == 4 and row[3] in ('NATAA', 'NATRO'):  # "30 K H NATAA"
            return row[1]
        elif len(row) == 4 and 'PIKAA' in row[3] and len(row[3]) == len('PIKAA ')+1:  # "30 K H PIKAA R"
            return row[3][-1]
        else: # ALLAA etc..
            return '*'

    def show_aligned(self, pose) -> None:
        """
        notebook output alignment of pose seq and blueprint expectation
        """
        self.show_seq_aligned(pose.sequence(), self.expected_seq())

    def _show_alignment(self, alignment: Alignment) -> None:
        """
        Display the formatted alignment.
        """
        # Extract aligned sequences
        a = ''.join(alignment.target[i] if i != -1 else '-' for i in alignment.aligned[0])
        b = ''.join(alignment.query[i] if i != -1 else '-' for i in alignment.aligned[1])
        # Create a gap line
        gap_line = ''.join('|' if a[i] == b[i] else ' ' for i in range(len(a)))
        formatted_gap = ''.join(['.' if c == '|' else '*' for c in gap_line])
        # Format the score
        score = f"Score={alignment.score}"
        # Display the formatted alignment
        display(HTML(f'<div style="font-family:monospace; display: inline-block; white-space: nowrap;">' +
                     f'{a}<br/>{formatted_gap}<br/>{b}<br/>{score}</div>'))

    def show_poses_aligned(self, pose_A, pose_B):
        """
        Show the sequence alignment of two poses.
        """
        self.show_seq_aligned(pose_A.sequence(), pose_B.sequence())

    def show_seq_aligned(self, seq_A: str, seq_B:str) -> List[Alignment]:
        """
        Show the alignment of two sequences.
        """
        # Align using PairwiseAligner
        aligner = PairwiseAligner()
        alignments: PairwiseAlignments = aligner.align(seq_A, seq_B)
        # Assuming you want the best alignment
        best_alignment: Alignment = alignments[0]

        # Call _show_alignment to display the alignment
        self._show_alignment(best_alignment)
        return [best_alignment]

    def correct(self, pose: pyrosetta.Pose) -> pyrosetta.rosetta.core.select.residue_selector.ResidueSelector:
        """
        This is not a great thing to do. So it is best to relax the neighbours of the vector afterwards.
        If there are valines instead of the intended sequence that is bad.
        Altering the blueprint is required.

        :param pose:
        :return:
        """
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        one2three = pyrosetta.rosetta.protocols.motifs.name3_from_oneletter
        ex_seq = self.expected_seq()
        altered = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        for i, (expected, current) in enumerate(zip(ex_seq, pose.chain_sequence(1))):
            if expected != current:
                altered.append_index(i + 1)
                MutateResidue(target=i + 1, new_res=one2three(expected)).apply(pose)
        return altered

    def correct_and_relax(self, pose: pyrosetta.Pose):
        """
        Runs ``correct`` and then relaxes the pose.
        The former is the last resort fix of a bad blueprint result.

        :param pose:
        :return:
        """
        altered = self.correct(pose)
        movemap = pyrosetta.MoveMap()
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        neigh_sele = NeighborhoodResidueSelector(altered, distance=15, include_focus_in_subset=True)
        n = neigh_sele.apply(pose)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(allow_chi=n)
        scorefxn = pyrosetta.get_fa_scorefxn()
        # cartesian
        scorefxn_cart = pyrosetta.create_score_function('ref2015_cart')
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_cart, 5)
        relax.cartesian(True)
        relax.minimize_bond_angles(True)
        relax.minimize_bond_lengths(True)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.apply(pose)
        # dual
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.apply(pose)





