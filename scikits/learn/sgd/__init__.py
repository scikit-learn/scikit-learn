"""
Module that implements Stochastic Gradient Descent related algorithms.

See http://scikit-learn.sourceforge.net/modules/sgd.html for complete
documentation.
"""

from . import sparse
from .sgd import SGD
from .base import Log, ModifiedHuber, Hinge
