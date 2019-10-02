from ._multilayer_perceptron import *  # noqa
from ..utils.deprecation import _raise_dep_warning_if_not_pytest


# TODO: remove entire file in 0.24
deprecated_path = 'sklearn.neural_network.multilayer_perceptron'
correct_path = 'sklearn.neural_network'

_raise_dep_warning_if_not_pytest(deprecated_path, correct_path)
