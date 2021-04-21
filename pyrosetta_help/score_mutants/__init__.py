# This has several variants...

# this is a mutation, which is validated against a pose (G10W will fail if resi 10 is not glycine).
from .mutation import Mutation

# this generates and scores variant poses
# called Variant in some older repos
from .variant import MutantScorer
from .scores import extend_scores
