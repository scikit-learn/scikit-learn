"""The :mod:`sklearn.inspection` module includes tools for model inspection."""

# TODO: remove me in 0.24 (as well as the noqa markers) and
# import the partial_dependence func directly from the
# ._partial_dependence module instead.
# Pre-cache the import of the deprecated module so that import
# sklearn.inspection.partial_dependence returns the function as in
# 0.21, instead of the module
# https://github.com/scikit-learn/scikit-learn/issues/15842
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore", category=FutureWarning)
    from .partial_dependence import partial_dependence

from ._partial_dependence import individual_conditional_expectation
from ._partial_dependence import plot_individual_conditional_expectation
from ._partial_dependence import IndividualConditionalExpectationDisplay
from ._partial_dependence import plot_partial_dependence
from ._partial_dependence import PartialDependenceDisplay
from ._permutation_importance import permutation_importance


__all__ = [
    'individual_conditional_expectation',
    'partial_dependence',
    'permutation_importance',
    'plot_individual_conditional_expectation',
    'plot_partial_dependence',
    'IndividualConditionalExpectationDisplay',
    'PartialDependenceDisplay'
]
