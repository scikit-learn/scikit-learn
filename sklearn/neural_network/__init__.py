"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

# Licence: BSD 3 clause

from .rbm import BernoulliRBM

__all__ = ['BernoulliRBM']

from .mlp import MLPClassifier
from .mlp import MultilayerPerceptronClassifier

from .multilayer_perceptron import MultilayerPerceptronClassifier
from .multilayer_perceptron import MultilayerPerceptronRegressor

__all__ = ["BernoulliRBM",
           "MultilayerPerceptronClassifier",
           "MultilayerPerceptronRegressor"]
