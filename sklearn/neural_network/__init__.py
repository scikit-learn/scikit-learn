HEAD
"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

from .rbm import BernoulliRBM

__all__ = ['BernoulliRBM']

from .mlp import MLPClassifier
(WIP) Added Multi-layer perceptron (MLP)
