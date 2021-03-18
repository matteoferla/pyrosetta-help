from typing import List
import pyrosetta

class BlueprinterExpected:
    def expected_seq(self):
        """

        :return:
        """
        return ''.join([self.get_expected_aa_from_row(row) for row in self])

    def get_expected_aa_from_row(self, row: List[str]) -> str:
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

    def show_aligned(self, pose):
        """
        notebook output alignment of pose seq and blueprint expectation
        :param pose:
        :return:
        """
        self.show_seq_aligned(pose.sequence(), self.expected_seq())

    def _show_alignment(self, alignment: list):
        from Bio import pairwise2
        from IPython.display import display, HTML
        formatted = pairwise2.format_alignment(*alignment)
        a, gap, b, score = formatted.strip().split('\n')
        gap = ''.join(['.' if c == '|' else '*' for c in gap])
        display(HTML(f'<div style="font-family:monospace; display: inline-block; white-space: nowrap;">'+
                     '{a}<br/>{gap}<br/>{b}<br/>{score}</div>'))  #ignore pycharm. typehint for display is wrong.

    def show_poses_aligned(self, pose_A, pose_B):
        self.show_seq_aligned(pose_A.sequence(), pose_B.sequence())

    def show_seq_aligned(self, seq_A: str, seq_B:str):
        from Bio import pairwise2
        alignment = pairwise2.align.globalxx(seq_A, seq_B)[0]
        self._show_alignment(alignment)
        return alignment

    def correct(self, pose: pyrosetta.Pose):
        """
        If there are valines instead of the intended sequence that is bad.
        Altering the blueprint is required.

        :param pose:
        :return:
        """
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        one2three = pyrosetta.rosetta.protocols.motifs.name3_from_oneletter
        ex_seq = self.expected_seq()
        altered = []
        for i, (expected, curr) in enumerate(zip(ex_seq, pose.chain_sequence(1))):
            if expected != curr:
                altered.append(i + 1)
                MutateResidue(target=i + 1, new_res=one2three(expected)).apply(pose)
        return altered





