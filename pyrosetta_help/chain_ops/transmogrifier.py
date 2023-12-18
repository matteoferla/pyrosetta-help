from typing import Tuple, List, Union
import re
from .chain_ops import ChainOps
from Bio.Align import PairwiseAligner


class Transmogrifier(ChainOps):
    """
    This is to convert a mutation from one species (``owned_label``) to another (``wanted_label``)
    based on provided sequences in the chains list of dict

    """

    def __init__(self, chains: List[dict], wanted_label: str, owned_label: str):
        super().__init__(chains)  # self.chains
        self.wanted_label = wanted_label
        self.owned_label = owned_label

    @classmethod
    def from_chain_ops(cls, chain_ops: ChainOps, wanted_label: str, owned_label: str):
        return cls(chain_ops.chains, wanted_label, owned_label)

    def align_seqs(self, chain_selection) -> Tuple[str, str]:
        """
        Align the human seq to the mouse and store it the chain dict
        """
        chain = self.chains[chain_selection]
        wanted_sequence = chain[f'{self.wanted_label}_sequence']
        owned_sequence = chain[f'{self.owned_label}_sequence']
        # Align using PairwiseAligner
        aligner = PairwiseAligner()
        aligner.open_gap_score = -1  # Gap open penalty
        aligner.extend_gap_score = -0.1  # Gap extension penalty
        # Perform the alignment
        alignment = aligner.align(wanted_sequence, owned_sequence)[0]  # Assuming the best alignment
        # Get the aligned sequences
        aligned_wanted = ''.join(wanted_sequence[i] if i != -1 else '-' for i in alignment.aligned[0])
        aligned_owned = ''.join(owned_sequence[i] if i != -1 else '-' for i in alignment.aligned[1])
        return aligned_wanted, aligned_owned

    def covert_A2B(self, seqA: str, seqB: str, resiA: int) -> int:
        """
        Given seqA and seqB as two gap aligned sequences,
        and an off-by-one residue number of seqA without counting gaps,
        return an off-by-one residue number of seqB without counting gaps.
        """
        assert resiA > 0, 'Negative number is no.'
        assert isinstance(resiA, int), 'Float is no'
        # first step
        get_resi2aln_map = lambda seq: [i for i, r in enumerate(seq) if r != '-']
        aln_pos = get_resi2aln_map(seqA)[resiA - 1]
        # second
        return get_resi2aln_map(seqB).index(aln_pos) + 1

    def transmogrify(self, mutation: str, chain_selection: Union[str, int, dict]) -> str:
        chain = self[chain_selection]
        human_resi = int(re.search(r'(\d+)', mutation).group(1))
        mouse_resi = self.covert_A2B(chain[f'{self.owned_label}_aln_sequence'],
                                     chain[f'{self.wanted_label}_aln_sequence'],
                                     human_resi)
        return re.sub(str(human_resi), str(mouse_resi), mutation)


class Murinizer(Transmogrifier):
    """
    Human --> Mouse
    """

    def __init__(self, chains: List[dict]):
        super().__init__(chains, 'mouse', 'human')

# from functools import partial
# med27_murinize = partial(murinize, entry=med27_chain)
