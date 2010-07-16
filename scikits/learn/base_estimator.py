"""
Base class for all estimators.

"""
# XXX: This is work in progress and will probably need refactoring

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
            assert key in self._params, 'Specified parameter, %s, unknown'
            #assert isinstance(value, self.params[key])
            setattr(self, key, value)


    def _get_params(self):
        out = dict()
        for key in self._params:
            out[key] = getattr(self, key)
        return out


    def __repr__(self):
        params_str = '\n   '.join('%s=%s' % (k, v) 
                                  for k, v in self._get_params().iteritems())
        return '%s(%s)' % (
                self.__class__.__name__,
                params_str
            )

