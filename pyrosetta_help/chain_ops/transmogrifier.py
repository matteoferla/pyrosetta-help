from Bio import pairwise2
from typing import Tuple, List, Union
import re
from .chain_ops import ChainOps


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
        alignments = pairwise2.align.globalxs(chain[f'{self.wanted_label}_sequence'],
                                              chain[f'{self.owned_label}_sequence'],
                                              -1,  # open
                                              -0.1  # extend
                                              )
        al = alignments[0]
        chain[f'{self.wanted_label}_aln_sequence'] = al[0]
        chain[f'{self.owned_label}_aln_sequence'] = al[1]
        return al[0], al[1]

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
