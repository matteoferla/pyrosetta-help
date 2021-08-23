__all__ = ['AF2NotebookAnalyser']

import numpy as np
import os, re, json
from typing import *
import pandas as pd
import pyrosetta
from .retrieval import reshape_errors
from .constraints import add_pae_constraints


class AF2NotebookAnalyser:

    def __init__(self, folder: str, load_poses: bool = True):
        self.folder = folder
        self.scores = self.make_AF2_dataframe()
        self._add_settings()
        self.poses = None
        if load_poses:
            self.poses = self.get_poses()
        self.errors = self.get_errors()


    def make_AF2_dataframe(self) -> pd.DataFrame:
        """
        Given a folder form ColabsFold return a dictionary
        with key rank index
        and value a dictionary of details

        This is convoluted, but it may have been altered by a human.
        """
        filenames = os.listdir(self.folder)
        # group files
        ranked_filenames = {}
        for filename in filenames:
            if not re.match('rank_\d+.*\.pdb', filename):
                continue
            rex = re.match('rank_(\d+)_model_(\d+).*seed_(\d+)(.*)\.pdb', filename)
            rank = int(rex.group(1))
            model = int(rex.group(2))
            seed = int(rex.group(3))
            other = rex.group(4)
            if rank in ranked_filenames and ranked_filenames[rank]['relaxed'] == True:
                continue
            data = dict(name=filename,
                        path=os.path.join(self.folder, filename),
                        rank=rank,
                        model=model,
                        seed=seed,
                        relaxed='_relaxed' in other
                        )
            ranked_filenames[rank] = data
        # make dataframe
        return pd.DataFrame(list(ranked_filenames.values()))

    def _add_settings(self):
        # add data from settings.
        pLDDTs = {}
        pTMscores = {}
        with open(os.path.join(self.folder, 'settings.txt'), 'r') as fh:
            for match in re.findall('rank_(\d+).*pLDDT\:(\d+\.\d+)\ pTMscore:(\d+\.\d+)', fh.read()):
                pLDDTs[int(match[0])] = float(match[1])
                pTMscores[int(match[0])] = float(match[2])
        self.scores['pLDDT'] = self.scores['rank'].map(pLDDTs)
        self.scores['pTMscore'] = self.scores['rank'].map(pTMscores)

    def get_poses(self) -> pyrosetta.rosetta.utility.vector1_core_pose_Pose:
        poses = pyrosetta.rosetta.utility.vector1_core_pose_Pose(len(self.scores))
        for i, row in self.scores.iterrows():
            poses[row['rank']] = pyrosetta.pose_from_file(row['path'])
        return poses

    def get_errors(self) -> dict:
        errors = dict()
        for i, row in self.scores.iterrows():
            filename = row['path'].replace('.pdb', '.json') \
                .replace('_unrelaxed', '_pae') \
                .replace('_relaxed', '_pae')
            with open(filename, 'r') as fh:
                errors[row['rank']] = reshape_errors(json.load(fh))
        return errors

    def constrain_and_relax(self):
        scorefxn = pyrosetta.get_fa_scorefxn()
        ap_st = pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint
        scorefxn.set_weight(ap_st, 1)
        if self.poses is None:
            raise ValueError('Load poses first.')
        for pose, error in zip(self.poses, self.errors):
            add_pae_constraints(pose, error)
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 3)
            relax.apply(pose)
        # Add dG
        self.scores['dG'] = self.scores['rank'].apply(lambda rank: scorefxn(self.poses[rank]))

    def calculate_interface(self, interface = 'A_B'):
        self.get_median_interface_bfactors()
        # interface score.
        if 'dG' not in self.scores:
            return
        newdata = []
        for rank in self.scores['rank']:
            ia = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(interface)
            ia.apply(self.poses[rank])
            newdata.append({'complex_energy': ia.get_complex_energy(),
                    'separated_interface_energy': ia.get_separated_interface_energy(),
                    'complexed_sasa': ia.get_complexed_sasa(),
                    'crossterm_interface_energy': ia.get_crossterm_interface_energy(),
                    'interface_dG': ia.get_interface_dG(),
                    'interface_delta_sasa': ia.get_interface_delta_sasa()})
        # adding multiple columns, hence why not apply route.
        # the order is not chanced, so all good
        newdata = pd.DataFrame(newdata)
        for column in newdata.columns:
            self.scores[column] = newdata[column]

    def _get_median_interface_bfactor(self, pose, residues):
        pbd_info = pose.pdb_info()
        bfactors = [pbd_info.bfactor(r, pose.residue(r).atom_index('CA')) for r in residues]
        return np.array(bfactors).median()

    def get_median_interface_bfactors(self):
        for i, row in self.details['interface_pLDDT'].iterrows():
            residues = list(row.interchain_residues_1) + list(row.interchain_residues_2)
            pose = self.poses[row['rank']] # ``row.rank`` is a function, just like ``row.name``.
            return self._get_median_interface_bfactor(pose, residues)


