"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

# Licence: BSD 3 clause

from .rbm import BernoulliRBM

from .multilayer_perceptron import MLPClassifier
from .multilayer_perceptron import MLPRegressor

from .random_basis_function import RandomBasisFunction

__all__ = ["BernoulliRBM",
           "MLPClassifier",
           "MLPRegressor",
           "RandomBasisFunction"]
