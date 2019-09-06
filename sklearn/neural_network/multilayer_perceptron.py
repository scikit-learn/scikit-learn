import warnings

from ._multilayer_perceptron import *  # noqa


msg = ("Any import from sklearn.neural_nework.multilayer_perceptron is "
       "deprecated in version "
       "0.22 and will not work anymore in version 0.24. You should only "
       "import from sklearn.neural_network directly.")

# TODO: remove entire file in 0.24
warnings.warn(msg, DeprecationWarning)
