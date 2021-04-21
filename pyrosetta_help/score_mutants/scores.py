from typing import *
from ..weights import term_meanings
import pandas as pd

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


def extend_scores(scores: pd.DataFrame):
    """
    Adds the following fields:

    * highest/lowest_contributor
    * highest/lowest_contributor_value
    * highest/lowest_contributor_wordy

    :param scores: pd.DataFrame(output_of_variants)
    :return:
    """
    scores['highest_contributor'] = scores.apply(lambda row: get_highest_contributor(row)[0].replace('delta_', ''), 1)
    scores['highest_contributor_value'] = scores.apply(lambda row: row['delta_' + row.highest_contributor], 1)
    scores['highest_contributor_wordy'] = scores.apply(lambda row: term_meanings[row.highest_contributor], 1)

    scores['lowest_contributor'] = scores.apply(lambda row: get_lowest_contributor(row)[0].replace('delta_', ''), 1)
    scores['lowest_contributor_value'] = scores.apply(lambda row: row['delta_' + row.lowest_contributor], 1)
    scores['lowest_contributor_wordy'] = scores.apply(lambda row: term_meanings[row.lowest_contributor], 1)
