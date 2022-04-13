"""
This module implements multioutput regression and classification.

The estimators provided in this module are meta-estimators: they require
a base estimator to be provided in their constructor. The meta-estimator
extends single output estimators to multioutput estimators.
"""

# Author: Tim Head <betatim@gmail.com>
# Author: Hugo Bowne-Anderson <hugobowne@gmail.com>
# Author: Chris Rivera <chris.richard.rivera@gmail.com>
# Author: Michael Williamson
# Author: James Ashton Nichols <james.ashton.nichols@gmail.com>
#
# License: BSD 3 clause

from ._multioutput import (
    MultiOutputRegressor,
    MultiOutputClassifier,
    RegressorChain,
    ClassifierChain
)


__all__ = [
    "MultiOutputRegressor",
    "MultiOutputClassifier",
    "RegressorChain",
    "ClassifierChain",
]
