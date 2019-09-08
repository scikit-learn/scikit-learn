import warnings

from ._rbm import *  # noqa


msg = ("The `sklearn.neural_nework.rbm` module is  deprecated in version "
       "0.22 and will be removed in version 0.24. "
       "The corresponding classes / functions "
       "should instead be imported from sklearn.neural_nework")

# TODO: remove entire file in 0.24
warnings.warn(msg, DeprecationWarning)
