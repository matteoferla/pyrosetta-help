from typing import *
import re, json

class ChainOps:
    """
    Works on the list of dict (metadata.json) with among others keys

    * number (pose number)
    * chain (chain letter)
    * gene_name (gene name)

    """

    def __init__(self, chains: List[dict]):
        self.chains = chains

    def load(self, json_filename: str = 'metadata.json'):
        with open(json_filename, 'r') as fh:
            self.chains = json.load(fh)

    def dump(self, json_filename: str = 'metadata.json'):
        with open(json_filename, 'w') as fh:
            json.dump(self.chains, fh)

    def get_entry_of_key(self, value, key):
        return [chain for chain in self.chains if chain[key] == value][0]

    def get_entry(self, value):
        """
        guess key...

        :param value:
        :return:
        """
        if isinstance(value, int):
            return self.get_entry_of_key(value, 'number')
        elif isinstance(value, str) and len(value) == 1:
            return self.get_entry_of_key(value, 'chain')
        elif isinstance(value, str):
            return self.get_entry_of_key(value, 'gene_name')
        elif isinstance(value, dict) and 'gene_name' in value:  # an entry was passed
            return value

    def get_pose_of_chain(self, pose, value):  # pose is a pyrosetta.Pose
        return pose.split_by_chain()[self.get_entry(value)['number']]

    def __getitem__(self, value):
        return self.get_entry(value)

    def __iter__(self):
        return self.chains
