"""
Base class for all estimators.

"""
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>

# License: BSD Style
import inspect

import numpy as np

from .metrics import explained_variance

################################################################################
class BaseEstimator(object):
    """ Base class for all estimators in the scikit learn

        Note
        =====

        All estimators should specify all the parameters that can be set
        at the class level in their __init__ as explicit keyword
        arguments (no *args, **kwargs).

    """

    @classmethod
    def _get_param_names(cls):
        try:
            args, varargs, kw, default = inspect.getargspec(cls.__init__)
            assert varargs is None, (
                'scikit learn estimators should always specify their '
                'parameters in the signature of their init (no varargs).'
                )
            # Remove 'self'
            # XXX: This is going to fail if the init is a staticmethod, but
            # who would do this?
            args.pop(0)
        except TypeError:
            # No explicit __init__
            args = []
        return args


    def _reinit(self):
        """ Constructs a new estimator with the same parameters
        than self. It's a kind of deep copy without actually copying
        attached data.
        """
        klass = self.__class__
        new_object_params = self._get_params()
        new_object = klass(**new_object_params)
        
        return new_object

    def _get_params(self, deep=False):
        """ Get parameters for current estimator

            Parameters
            ==========
            deep: boolean, optional
                If True, will return the parameters for this estimator and
                contained subobjects that are estimators.
        """
        out = dict()
        for key in self._get_param_names():
            value = getattr (self, key)
            if deep and hasattr (value, '_get_params'):
                deep_items = value._get_params().items()
                out.update((key + '__' + k, val) for k, val in deep_items)
            else:
                out[key] = value
        return out

    def _set_params(self, **params):
        """ Set the parameters of the estimator.
        """
        valid_params = self._get_param_names()
        for key, value in params.iteritems():
            split = key.split('__', 1)
            if len(split) > 1:
                getattr (self, split[0])._set_params(**{split[1]:value})
            else:
                assert key in valid_params, ('Invalid parameter %s '
                                             'for estimator %s' %
                                             (key, self.__class__.__name__))
                setattr(self, key, value)


    def __repr__(self):
        options = np.get_printoptions()
        np.set_printoptions(precision=5, threshold=64, edgeitems=2)
        class_name = self.__class__.__name__

        # Do a multi-line justified repr:
        params_list = list()
        this_line_length = len(class_name)
        line_sep = ',\n' + (1+len(class_name)/2)*' '
        for i, (k, v) in enumerate(self._get_params().iteritems()):
            if type(v) is float:
                # use str for representing floating point numbers
                # this way we get consistent representation across
                # architectures and versions.
                this_repr  = '%s=%s' % (k, str(v))
            else:
                # use repr of the rest
                this_repr  = '%s=%s' % (k, repr(v))
            if i > 0: 
                if (this_line_length + len(this_repr) >= 75
                                            or '\n' in this_repr):
                    params_list.append(line_sep)
                    this_line_length += len(line_sep)
                else:
                    params_list.append(', ')
                    this_line_length += 2
            params_list.append(this_repr)
            this_line_length += len(this_repr)

        params_str = ''.join(params_list)
        np.set_printoptions(**options)
        return '%s(%s)' % (
                class_name,
                params_str
            )


################################################################################
class ClassifierMixin(object):
    """ Mixin class for all classifiers in the scikit learn
    """

    def score(self, X, y):
        """ Returns the mean error rate on the given test data and labels.

            Parameters
            ----------
            X : array-like, shape = [n_samples, n_features]
                Training set.

            y : array-like, shape = [n_samples]
                Labels for X.

            Returns
            -------
            z : float
        """
        return np.mean(self.predict(X) == y)


################################################################################
class RegressorMixin(object):
    """ Mixin class for all regression estimators in the scikit learn
    """

    def score(self, X, y):
        """ Returns the explained variance of the prediction

            Parameters
            ----------
            X : array-like, shape = [n_samples, n_features]
                Training set.

            y : array-like, shape = [n_samples]

            Returns
            -------
            z : float
        """
        return explained_variance(y, self.predict(X))
