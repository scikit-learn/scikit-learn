"""
Base class for all estimators.

"""
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>

# License: BSD Style

import numpy as np

# XXX: This is work in progress and will probably need refactoring

# TODO:
#       * It might be worth the while to use the inspect module to
#         check at instanciationg the __init__ signature and check that
#         we are not creating a __init__ that is not consistent with the
#         _params

################################################################################
class BaseEstimator(object):
    """ Base class for all estimators in the scikit learn

        _params: parameters of the model, but not parameters of the
        estimation algorithm (such as maxiter)
    """

    def __init__(self, **params):
        assert hasattr(self, '_params'), \
                    'Estimator class without parameter definition'
        self._set_params(**params)


    def _set_params(self, **params):
        for key, value in params.iteritems():
            assert key in self._params, ('Specified parameter, %s, unknown'
                            % key)
            # Right now we are punting on type checking: this gives too
            # much work for no good.
            #assert isinstance(value, self.params[key])
            setattr(self, key, value)


    def _get_params(self):
        out = dict()
        for key in self._params:
            out[key] = getattr(self, key)
        return out


    def __repr__(self):
        options = np.get_printoptions()
        np.set_printoptions(precision=5, threshold=64, edgeitems=2)
        class_name = self.__class__.__name__
        params_str = (',\n' + (1+len(class_name))*' ').join(
                                  '%s=%s' % (k, v) 
                                  for k, v in self._get_params().iteritems())
        np.set_printoptions(**options)
        return '%s(%s)' % (
                class_name,
                params_str
            )

