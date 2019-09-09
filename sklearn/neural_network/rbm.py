import warnings

from ._rbm import *  # noqa
from ..utils.deprecation import _get_deprecation_message


msg = _get_deprecation_message(
       deprecated_path='sklearn.neural_network.rbm',
       correct_path='sklearn.neural_network'
)

# TODO: remove entire file in 0.24
warnings.warn(msg, DeprecationWarning)
