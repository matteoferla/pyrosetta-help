from typing import *
import pandas as pd

# https://www.rosettacommons.org/docs/latest/rosetta_basics/scoring/score-types
term_meanings = {
    "fa_atr": "Lennard-Jones attractive between atoms in different residues (r^6 term, London dispersion forces).",
    "fa_rep": "Lennard-Jones repulsive between atoms in different residues (r^12 term, Pauli repulsion forces).",
    "fa_sol": "Lazaridis-Karplus solvation energy.",
    "fa_intra_rep": "Lennard-Jones repulsive between atoms in the same residue.",
    "fa_elec": "Coulombic electrostatic potential with a distance-dependent dielectric.",
    "pro_close": "Proline ring closure energy and energy of psi angle of preceding residue.",
    "hbond_sr_bb": "Backbone-backbone hbonds close in primary sequence.",
    "hbond_lr_bb": "Backbone-backbone hbonds distant in primary sequence.",
    "hbond_bb_sc": "Sidechain-backbone hydrogen bond energy.",
    "hbond_sc": "Sidechain-sidechain hydrogen bond energy.",
    "dslf_fa13": "Disulfide geometry potential.",
    "rama": "Ramachandran preferences.",
    "omega": "Omega dihedral in the backbone. A Harmonic constraint on planarity with standard deviation of ~6 deg.",
    "fa_dun": "Internal energy of sidechain rotamers as derived from Dunbrack's statistics (2010 Rotamer Library used in Talaris2013).",
    "p_aa_pp": "Probability of amino acid at Φ/Ψ.",
    "ref": "Reference energy for each amino acid. Balances internal energy of amino acid terms.  Plays role in design.",
    "METHOD_WEIGHTS": "Not an energy term itself, but the parameters for each amino acid used by the ref energy term.",
    "lk_ball": "Anisotropic contribution to the solvation.",
    "lk_ball_iso": "Same as fa_sol; see below.",
    "lk_ball_wtd": "weighted sum of lk_ball & lk_ball_iso (w1*lk_ball + w2*lk_ball_iso); w2 is negative so that anisotropic contribution(lk_ball) replaces some portion of isotropic contribution (fa_sol=lk_ball_iso).",
    "lk_ball_bridge": "Bonus to solvation coming from bridging waters, measured by overlap of the 'balls' from two interacting polar atoms.",
    "lk_ball_bridge_uncpl": "Same as lk_ball_bridge, but the value is uncoupled with dGfree (i.e. constant bonus, whereas lk_ball_bridge is proportional to dGfree values).",
    "fa_intra_atr_xover4": "Intra-residue LJ attraction, counted for the atom-pairs beyond torsion-relationship.",
    "fa_intra_rep_xover4": "Intra-residue LJ repulsion, counted for the atom-pairs beyond torsion-relationship.",
    "fa_intra_sol_xover4": "Intra-residue LK solvation, counted for the atom-pairs beyond torsion-relationship.",
    "fa_intra_elec": "Intra-residue Coulombic interaction, counted for the atom-pairs beyond torsion-relationship.",
    "rama_prepro": "Backbone torsion preference term that takes into account of whether preceding amono acid is Proline or not.",
    "hxl_tors": "Sidechain hydroxyl group torsion preference for Ser/Thr/Tyr, supersedes yhh_planarity (that covers L- and D-Tyr only).",
    "yhh_planarity": "Sidechain hydroxyl group torsion preference for Tyr, superseded by hxl_tors"
}

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
