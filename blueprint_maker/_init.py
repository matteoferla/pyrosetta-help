import pyrosetta

class BlueprinterInit:
    def __init__(self, seq: str, ss: str = None):
        self.seq = seq
        self.ss = ss
        self.rows = self._seq2rows(seq)

    def _seq2rows(self, seq):
        rows = []
        for i, aa in enumerate(seq):
            rows.append([i + 1, aa, '.'])
        return rows

    @classmethod
    def from_pose(cls, pose: pyrosetta.Pose):
        seq = pose.chain_sequence(1)
        DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
        DSSP.apply(pose)
        ss = pose.secstruct()  # for this case, no slicing.
        return cls(seq, ss)

    def get_ss(self, idx: int):
        if self.ss is not None and len(self.ss) > idx:
            return self.ss[idx - 1]
        else:
            return 'D'