"""
The :mod:`sklearn.neural_network` module includes models based on neural
networks.

**User guide.** See the :ref:`neural_networks_supervised` and
:ref:`neural_networks_unsupervised` sections for further details.
"""

# License: BSD 3 clause

from ._multilayer_perceptron import MLPClassifier, MLPRegressor
from ._rbm import BernoulliRBM

__all__ = ["BernoulliRBM", "MLPClassifier", "MLPRegressor"]
