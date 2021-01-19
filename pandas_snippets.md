## Pandas

### MutantScorer specific

Given a table of scored variants:

    from score_mutants import MutantScorer
    import pandas as pd
    
    model = Variant(pose, modelname='holo')   
    mutations = 'p.A10W p.D20W'.split()
    data = model.score_mutations(mutations,
                                chain='A',
                                interfaces=(('isolated', 'A_B'),),
                                preminimise=False,
                                distance=12,
                                cycles=5)
    scores = pd.DataFrame(data)

The following methods are handy as the table contains weighted contributions of each term, 
and these find the largest contribution
either the lower or the highest.

    from typing import *

    def get_lowest_contributor(row: pd.Series) -> Tuple[str, float]:
        """
        Not smallest/infinitesimal in abs contribution, but lowest number (i.e. most negative)
        """
        delta_names = [col for col in row.index if 'delta' in col]
        srow = row[delta_names].sort_values(ascending=True)
        return srow.index[0], srow[0]
    
    def get_highest_contributor(row: pd.Series) -> Tuple[str, float]:
        """
        Not largest in abs contribution, but highest number (i.e. most positive)
        """
        delta_names = [col for col in row.index if 'delta' in col]
        srow = row[delta_names].sort_values(ascending=False)
        return srow.index[0], srow[0]
    
    def get_largest_contributor(row: pd.Series) -> Tuple[str, float]:
        """
        Largest in abs amount as per the confusing fact that a very low negative number is large.
        """
        delta_names = [col for col in row.index if 'delta' in col]
        srow = row[delta_names].abs().sort_values(ascending=False)
        return srow.index[0], srow[0]
    
    
Applied to a dataframe:

    scores['highest_contributor'] = scores.apply(lambda row: get_highest_contributor(row)[0].replace('delta_', ''), 1)
    scores['highest_contributor_value'] = scores.apply(lambda row: row['delta_'+row.highest_contributor], 1)
    scores['highest_contributor_wordy'] = scores.apply(lambda row: MutantScorer.term_meanings[row.highest_contributor], 1)
    
    scores['lowest_contributor'] = scores.apply(lambda row: get_lowest_contributor(row)[0].replace('delta_', ''), 1)
    scores['lowest_contributor_value'] = scores.apply(lambda row: row['delta_'+row.lowest_contributor], 1)
    scores['lowest_contributor_wordy'] = scores.apply(lambda row: MutantScorer.term_meanings[row.lowest_contributor], 1)
    
The latter, uses `MutantScorer.term_meanings` dictionary. This could be handy by itself.
As it converts "fa_atr" -> "Lennard-Jones attractive between atoms in different residues (r^6 term, London dispersion forces)." etc.