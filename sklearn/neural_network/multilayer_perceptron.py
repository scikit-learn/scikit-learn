import warnings

from ._multilayer_perceptron import *  # noqa


msg = ("The `sklearn.neural_nework.multilayer_perceptron` module is "
       "deprecated in version 0.22 and will be removed in version 0.24. "
       "The corresponding classes / functions "
       "should instead be imported from sklearn.neural_nework. "
       "Anything that cannot be imported from sklearn.neural_network is now "
       "part of the private API.")

# TODO: remove entire file in 0.24
warnings.warn(msg, DeprecationWarning)
