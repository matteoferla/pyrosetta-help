import pyrosetta
import re


class Mutation:
    """
    A mutation is an object that has all the details of the mutation.
    A variant, as interpreted in ``Model.score_mutations`` is a pose with a mutation.
    """
    _name3 = {'A': 'ALA',
              'C': 'CYS',
              'D': 'ASP',
              'E': 'GLU',
              'F': 'PHE',
              'G': 'GLY',
              'H': 'HIS',
              'I': 'ILE',
              'L': 'LEU',
              'K': 'LYS',
              'M': 'MET',
              'N': 'ASN',
              'P': 'PRO',
              'Q': 'GLN',
              'R': 'ARG',
              'S': 'SER',
              'T': 'THR',
              'V': 'VAL',
              'W': 'TRP',
              'Y': 'TYR'}

    def __init__(self, mutation_name: str, chain: str, pose: pyrosetta.Pose):
        self.mutation = self.parse_mutation(mutation_name)
        rex = re.match('(\w)(\d+)(\w)', self.mutation)
        self.pdb_resi = int(rex.group(2))
        self.chain = chain
        self.from_resn1 = rex.group(1)
        self.from_resn3 = self._name3[rex.group(1)]
        self.to_resn1 = rex.group(3)
        self.to_resn3 = self._name3[rex.group(3)]
        pose2pdb = pose.pdb_info().pdb2pose
        self.pose_resi = pose2pdb(res=self.pdb_resi, chain=self.chain)
        if self.pose_resi != 0:
            self.pose_residue = pose.residue(self.pose_resi)
            self.pose_resn1 = self.pose_residue.name1()
            self.pose_resn3 = self.pose_residue.name3()
        else:
            self.pose_residue = None
            self.pose_resn1 = None
            self.pose_resn3 = None

    def parse_mutation(self, mutation: str):
        if mutation[:2] == 'p.':
            mutation = mutation.replace('p.', '')
        if mutation[1].isdigit():
            return mutation
        else:
            value2key = lambda value: list(self._name3.keys())[list(self._name3.values()).index(value.upper())]
            return value2key(mutation[:3]) + mutation[3:-3] + value2key(mutation[-3:])

    def is_valid(self):
        return self.pose_resn1 == self.from_resn1

    def assert_valid(self):
        assert self.is_valid(), f'residue {self.pose_resi}(pose)/{self.pdb_resi}(pdb) ' + \
                                f'is a {self.pose_resn3}, not a {self.from_resn3}'

    def __str__(self):
        return self.mutation