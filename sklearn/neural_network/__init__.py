"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

# Licence: BSD 3 clause

from .rbm import BernoulliRBM

from .extreme_learning_machines import ELMClassifier
from .extreme_learning_machines import ELMRegressor

__all__ = ["BernoulliRBM",
           "ELMClassifier",
           "ELMRegressor"]
