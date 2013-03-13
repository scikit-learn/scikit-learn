"""
The :mod:`sklearn.neural_networks` module implements neural network models.
It includes the Extreme Learning Machine regressor and classifier, and
the Random Hidden Layer transformers.
"""

from .elm import (ELMRegressor, SimpleELMRegressor,
                  ELMClassifier, SimpleELMClassifier)
from .random_hidden_layer import (RBFRandomHiddenLayer,
                                  SimpleRandomHiddenLayer)

__all__ = ['ELMRegressor',
           'ELMClassifier',
           'SimpleELMRegressor',
           'SimpleELMClassifier',
           'RBFRandomHiddenLayer',
           'SimpleRandomHiddenLayer']
