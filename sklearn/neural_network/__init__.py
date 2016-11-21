"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

# License: BSD 3 clause

from .rbm import BernoulliRBM, GaussianBernoulliRBM

from .multilayer_perceptron import MLPClassifier
from .multilayer_perceptron import MLPRegressor

__all__ = ["BernoulliRBM",
           "GaussianBernoulliRBM",
           "MLPClassifier",
           "MLPRegressor"]
