from io import StringIO
import requests
from typing import *
from Bio.Blast.NCBIWWW import qblast
from Bio import SearchIO
import pandas as pd
from collections import Counter


class LigandHunter:
    """
    Given a sequence find homologues (Blast) and find if they have ligands (PDBe API).
    The definition of what ligand is a cofactor comes from PDBe and does not count ions or
    triphospho-nucleotides as cofactors, but does count NADH adducts.

    * ``.data`` is a list of dictionaries.
    * ``.to_dataframe()`` converts into a pandas dataframe for further analysis.
    * ``.candidate_ligands`` list of ligand residue 3-letter codes
    """
    cofactor_reference = requests.get('https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors').json()
    _grouped_cofactor_codes = {name: value[0]['cofactors'] for name, value in cofactor_reference.items()}
    cofactor_codes = [code for codes in _grouped_cofactor_codes.values() for code in codes]

    def __init__(self, sequence: str):
        """

        :param sequence: is assumed clean protein sequence
        """
        self.sequence = sequence
        self._blast_result = self._retrieve_homologues()  #:StringIO
        self.data = self._parse_blast_result()
        self._retrieve_ligands()  # inplace.
        self._ligand_data = None

    def to_dataframe(self):
        """
        Converts ``.data`` to a pandas dataframe
        """
        return pd.DataFrame(self.data).transpose()

    @property
    def candidate_ligands(self) -> List[str]:
        return list(set([code for datum in self.data.values() for code in datum['ligand_codes']]))

    @property  # cached the old way
    def ligand_data(self):
        """
        Data for the ligands.

        Note that LigandNicker has a get smiles method, that works like this, but is unrelated.
        """
        if self._ligand_data is None:
            self._ligand_data = requests.post('https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/',
                                              data=','.join(self.candidate_ligands)).json()
        return self._ligand_data

    def get_most_common_ligands(self) -> List[Tuple[str, int]]:
        """
        Uses collections.counter

        :return:
        """
        c = Counter([code for datum in self.data.values() for code in datum['ligand_codes']])
        return c.most_common()

    def get_pdb_entry_by_ligand(self, ligand_code: str) -> dict:
        """
        get pdb **entry** by ligand.
        Returns the first, which should be lowest e-value
        """
        for datum in self.data.values():
            if ligand_code in datum['ligand_codes']:
                return datum
        else:
            raise ValueError(f'{ligand_code} not found in any of the {len(self.data)} hits')

    # ------------ initialisation methods ------------------------

    def _retrieve_homologues(self) -> StringIO:
        return qblast('blastp', 'pdb', self.sequence)

    def _parse_blast_result(self) -> Dict[str, dict]:
        results = {}
        self._blast_result.seek(0)
        for query in SearchIO.parse(self._blast_result, 'blast-xml'):
            for hit in query:
                datum = dict(accession=hit.accession,
                             description=hit.description,
                             evalue=hit.hsps[0].evalue)
                datum['pdb_code'], datum['chain'] = hit.accession.split('_')
                results[datum['pdb_code'].upper()] = datum
        return results

    def _retrieve_ligands(self) -> None:
        query = ','.join(self.data.keys())
        reply = requests.post(url='https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/',
                              data=query)
        for code, datum in reply.json().items():
            entry = self.data[code.upper()]
            entry['ligands'] = datum
            entry['ligand_codes'] = [inner['chem_comp_id'] for inner in datum]
            entry['cofactor_codes'] = [code for code in entry['ligand_codes'] if code in self.cofactor_codes]
            entry['has_cofactor'] = bool(entry['cofactor_codes'])


LigandHunter.cofactor_codes += ['ATP', 'GTP', 'CA', 'MG', 'W']
