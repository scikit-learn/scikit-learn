"""Base class for all estimators."""
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD Style

import copy
import inspect
import numpy as np
from .metrics import r2_score


###############################################################################
def clone(estimator, safe=True):
    """ Constructs a new estimator with the same parameters.

    Clone does a deep copy of the model in an estimator
    without actually copying attached data. It yields a new estimator
    with the same parameters that has not been fit on any data.

    Parameters
    ----------
    estimator: estimator object, or list, tuple or set of objects
        The estimator or group of estimators to be cloned

    safe: boolean, optional
        If safe is false, clone will fall back to a deepcopy on objects
        that are not estimators.

    """
    estimator_type = type(estimator)
    # XXX: not handling dictionnaries
    if estimator_type in (list, tuple, set, frozenset):
        return estimator_type([clone(e, safe=safe) for e in estimator])
    elif not hasattr(estimator, '_get_params'):
        if not safe:
            return copy.deepcopy(estimator)
        else:
            raise ValueError("Cannot clone object '%s' (type %s): "
                    "it does not seem to be a scikit-learn estimator as "
                    "it does not implement a '_get_params' methods."
                    % (repr(estimator), type(estimator)))
    klass = estimator.__class__
    new_object_params = estimator._get_params(deep=False)
    for name, param in new_object_params.iteritems():
        new_object_params[name] = clone(param, safe=False)
    new_object = klass(**new_object_params)
    assert new_object._get_params(deep=False) == new_object_params, (
            'Cannot clone object %s, as the constructor does not '
            'seem to set parameters' % estimator
        )

    return new_object


###############################################################################
def _pprint(params, offset=0, printer=repr):
    """ Pretty print the dictionnary 'params'

        Parameters
        ----------
        params: dict
            The dictionnary to pretty print

        offset: int
            The offset in characters to add at the begin of each line.

        printer:
            The function to convert entries to strings, typically
            the builtin str or repr
    """
    # Do a multi-line justified repr:
    options = np.get_printoptions()
    np.set_printoptions(precision=5, threshold=64, edgeitems=2)
    params_list = list()
    this_line_length = offset
    line_sep = ',\n' + (1 + offset / 2) * ' '
    for i, (k, v) in enumerate(params.iteritems()):
        if type(v) is float:
            # use str for representing floating point numbers
            # this way we get consistent representation across
            # architectures and versions.
            this_repr = '%s=%s' % (k, str(v))
        else:
            # use repr of the rest
            this_repr = '%s=%s' % (k, printer(v))
        if len(this_repr) > 500:
            this_repr = this_repr[:300] + '...' + this_repr[-100:]
        if i > 0:
            if (this_line_length + len(this_repr) >= 75
                                        or '\n' in this_repr):
                params_list.append(line_sep)
                this_line_length = len(line_sep)
            else:
                params_list.append(', ')
                this_line_length += 2
        params_list.append(this_repr)
        this_line_length += len(this_repr)

    np.set_printoptions(**options)
    lines = ''.join(params_list)
    # Strip trailing space to avoid nightmare in doctests
    lines = '\n'.join(l.rstrip(' ') for l in lines.split('\n'))
    return lines


###############################################################################
class BaseEstimator(object):
    """ Base class for all estimators in the scikit learn

        Notes
        -----
        All estimators should specify all the parameters that can be set
        at the class level in their __init__ as explicit keyword
        arguments (no *args, **kwargs).

    """

    @classmethod
    def _get_param_names(cls):
        """ Get parameter names for the estimator
        """
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

    def _get_params(self, deep=True):
        """ Get parameters for the estimator

            Parameters
            ----------
            deep: boolean, optional
                If True, will return the parameters for this estimator and
                contained subobjects that are estimators.
        """
        out = dict()
        for key in self._get_param_names():
            value = getattr(self, key)
            if deep and hasattr(value, '_get_params'):
                deep_items = value._get_params().items()
                out.update((key + '__' + k, val) for k, val in deep_items)
            out[key] = value
        return out

    def _set_params(self, **params):
        """ Set the parameters of the estimator.

        The method works on simple estimators as well as on nested
        objects (such as pipelines). The former have parameters of the
        form <component>__<parameter> so that the its possible to
        update each component of the nested object.
        """
        if not params:
            # Simple optimisation to gain speed (inspect is slow)
            return
        valid_params = self._get_params(deep=True)
        for key, value in params.iteritems():
            split = key.split('__', 1)
            if len(split) > 1:
                # nested objects case
                name, sub_name = split
                assert name in valid_params, ('Invalid parameter %s '
                                              'for estimator %s' %
                                             (name, self))
                sub_object = valid_params[name]
                assert hasattr(sub_object, '_get_params'), (
                    'Parameter %s of %s is not an estimator, cannot set '
                    'sub parameter %s' %
                        (sub_name, self.__class__.__name__, sub_name)
                    )
                sub_object._set_params(**{sub_name: value})
            else:
                # simple objects case
                assert key in valid_params, ('Invalid parameter %s '
                                              'for estimator %s' %
                                             (key, self.__class__.__name__))
                setattr(self, key, value)
        return self

    def __repr__(self):
        class_name = self.__class__.__name__
        return '%s(%s)' % (
                class_name,
                _pprint(self._get_params(deep=False),
                        offset=len(class_name),
                ),
            )

    def __str__(self):
        class_name = self.__class__.__name__
        return '%s(%s)' % (
                class_name,
                _pprint(self._get_params(deep=True),
                        offset=len(class_name),
                        printer=str,
                ),
            )


###############################################################################
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


###############################################################################
class RegressorMixin(object):
    """ Mixin class for all regression estimators in the scikit learn
    """

    def score(self, X, y):
        """ Returns the coefficient of determination of the prediction

            Parameters
            ----------
            X : array-like, shape = [n_samples, n_features]
                Training set.

            y : array-like, shape = [n_samples]

            Returns
            -------
            z : float
        """
        return r2_score(y, self.predict(X))


###############################################################################
class TransformerMixin(object):
    """ Mixin class for all transformers in the scikit learn
    """

    def fit_transform(self, X, y=None, **fit_params):
        """Fit model to data and subsequently transform the data

        Sometimes, fit and transform can be implemented more efficiently
        jointly than separately. In those cases, the estimator will typically
        override the method.

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]
            Training set.

        y : numpy array of shape [n_samples]
            Target values.

        Returns
        -------
        self : returns an instance of self.
        """
        if y is None:
            # fit method of arity 1 (unsupervised transformation)
            return self.fit(X, **fit_params).transform(X)
        else:
            # fit method of arity 2 (supervised transformation)
            return self.fit(X, y, **fit_params).transform(X)


###############################################################################
# XXX: Temporary solution to figure out if an estimator is a classifier

def _get_sub_estimator(estimator):
    """ Returns the final estimator if there is any.
    """
    if hasattr(estimator, 'estimator'):
        # GridSearchCV and other CV-tuned estimators
        return _get_sub_estimator(estimator.estimator)
    if hasattr(estimator, 'steps'):
        # Pipeline
        return _get_sub_estimator(estimator.steps[-1][1])
    return estimator


def is_classifier(estimator):
    """ Returns True if the given estimator is (probably) a classifier.
    """
    estimator = _get_sub_estimator(estimator)
    return isinstance(estimator, ClassifierMixin)
