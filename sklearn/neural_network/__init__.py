<<<<<<< HEAD
HEAD
=======
>>>>>>> squashed last 40 commits
"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.
"""

# Licence: BSD 3 clause

from .rbm import BernoulliRBM
<<<<<<< HEAD
<<<<<<< HEAD

__all__ = ['BernoulliRBM']

from .mlp import MLPClassifier
(WIP) Added Multi-layer perceptron (MLP)
=======
from .mlp import MultilayerPerceptronClassifier
>>>>>>> rebased
=======
from .multilayer_perceptron import MultilayerPerceptronClassifier
from .multilayer_perceptron import MultilayerPerceptronRegressor

__all__ = ["BernoulliRBM",
           "MultilayerPerceptronClassifier",
           "MultilayerPerceptronRegressor"]
>>>>>>> squashed last 40 commits
