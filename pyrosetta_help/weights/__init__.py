import pyrosetta
import os
import pandas as pd
from typing import *
from .terms import term_meanings


class WeightWatcher:
    """
    A class to easily look at the weights for the terms of different scorefunctions.
    It has various functionalities.

    List available scorefunctions:

    .. code-block:: python
      ww = WeightWatcher()
      print( ww.possible_scorefxn_names ) # dynamic attribute

    Find a scorefunction that mentions a word:

    .. code-block:: python
      ww.find_metion('foo')

    Get the comment block of a scorefunction:

    .. code-block:: python
      print(ww.get_scorefxn_comments('ref2015'))

    Get a scorefunction by name by calling the appropriate options first.

    .. code-block:: python
      scorefxn = ww.get_scorefxn('beta_nov16')

    Get weights (including ref) for a scorefunction or a name of one

    .. code-block:: python
      weights = ww.get_weights('beta_nov16')

    Get a pandas table of the weights

    .. code-block:: python
      weight_table = ww.compare(['ref2015', 'beta_nov16'], different_only=True)

    """
    folder = os.path.join(os.path.split(pyrosetta.__file__)[0], 'database', 'scoring', 'weights')
    term_meanings = term_meanings

    # def __init__(self):
    #     pyrosetta.distributed.maybe_init(extra_options='-mute all')

    def get_weight(self,
                   scorefxn: pyrosetta.ScoreFunction,
                   score_type: pyrosetta.rosetta.core.scoring.ScoreType
                   ) -> float:
        return scorefxn.get_weight(getattr(pyrosetta.rosetta.core.scoring.ScoreType, score_type.name))

    def get_weights(self,
                    scorefxn: Union[str, pyrosetta.ScoreFunction]
                    ) -> Dict[str, float]:
        """
        Gets weights for scorefxn
        """
        if isinstance(scorefxn, str):
            # scorefxn = pyrosetta.create_score_function(scorefxn)
            # mod:
            scorefxn = self.get_scorefxn(scorefxn)
        return {'name': scorefxn.get_name(),
                **{score_type.name: self.get_weight(scorefxn, score_type) for score_type in
                   scorefxn.get_nonzero_weighted_scoretypes()},
                **self.get_ref_values_badly(scorefxn, prefix=True)}

    def get_ref_values_badly(self,
                             scorefxn: pyrosetta.ScoreFunction,
                             prefix: Optional[bool] = True
                             ) -> Dict[str, float]:
        """
        I could not find a direct way to get_method_weights.
        Here this horrid way.
        """
        # similarity order 'IVLFCMAGTSWYPHNDEQKR'
        # name1 alphabetical order 'ACDEFGHIKLMNPQRSTVWY'
        # name3 alphabetical order 'ARNDCQEGHILKMFPSTWYV'
        pose = pyrosetta.pose_from_sequence('ARNDCQEGHILKMFPSTWYV')
        scorefxn(pose)
        aa_ref = {}
        for res in range(1, pose.total_residue() + 1):
            name = pose.residue(res).name3()
            value = pose.energies().residue_total_energies(res)[pyrosetta.rosetta.core.scoring.ScoreType.ref]
            if prefix:
                aa_ref[f'ref_{name}'] = value
            else:
                aa_ref[name] = value
        return aa_ref

    def get_scorefxn(self, scorefxn_name: str) -> pyrosetta.ScoreFunction:
        """
        Gets the scorefxn with appropriate corrections.
        """
        if isinstance(scorefxn_name, pyrosetta.ScoreFunction):
            return scorefxn_name  # it's already a scorefxn (not a name)
        corrections = {'beta_july15': False,
                       'beta_nov16': False,
                       'gen_potential': False,
                       'restore_talaris_behavior': False
                       }
        if 'beta_july15' in scorefxn_name or 'beta_nov15' in scorefxn_name:
            # beta_july15 is ref2015
            corrections['beta_july15'] = True
        elif 'beta_nov16' in scorefxn_name:
            corrections['beta_nov16'] = True
        elif 'genpot' in scorefxn_name:
            corrections['gen_potential'] = True
            pyrosetta.rosetta.basic.options.set_boolean_option('corrections:beta_july15', True)
        elif 'talaris' in scorefxn_name:  # 2013 and 2014
            corrections['restore_talaris_behavior'] = True
        else:
            pass
        for corr, value in corrections.items():
            pyrosetta.rosetta.basic.options.set_boolean_option(f'corrections:{corr}', value)
        return pyrosetta.create_score_function(scorefxn_name)

    @property
    def possible_scorefxn_names(self) -> List[str]:
        """
        Returns the scorefxn names
        """
        return sorted([fn.replace('.wts', '') for fn in os.listdir(self.folder) if '.wts' == os.path.splitext(fn)[1]])

    def get_scorefxn_block(self, scorefxn_name: str) -> str:
        """
        Read the file given a scorefxn name
        """
        if '/' in scorefxn_name:
            filename = scorefxn_name
        elif '.wts' in scorefxn_name:
            filename = os.path.join(self.folder, f'{scorefxn_name}')
        else:
            filename = os.path.join(self.folder, f'{scorefxn_name}.wts')
        with open(filename, 'r') as fh:
            return fh.read()

    def get_scorefxn_comments(self, scorefxn_name: str) -> str:
        """
        Get only the comments
        """
        lines = self.get_scorefxn_block(scorefxn_name).split('\n')
        return '\n'.join([line for line in lines if '#' in line])

    def find_metion(self, word: str) -> List[str]:
        """
        Find all the scorefxn names whose files contain the string word (case insensitive).

        >>> find_metion('spades')

        returns ``['hydrate_score12']``
        """
        mentionants = []
        for scorefxn_name in self.get_possible_scorefxn_names:
            block = self.get_scorefxn_block(scorefxn_name)
            if word.lower() in block.lower():
                mentionants.append(scorefxn_name)
        return mentionants

    def compare(self, scorefxn_names: List[str], different_only: bool = False) -> pd.DataFrame:
        terms = pd.DataFrame(map(self.get_weights, scorefxn_names)).fillna(0).set_index('name')
        if different_only:
            unique_cols = [col for col in terms.columns if len(terms[col].unique()) > 1]
            return terms[unique_cols]
        else:
            return terms