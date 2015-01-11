"""Base classes for all estimators."""
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD 3 clause

import copy
import inspect
import warnings

import numpy as np
from scipy import sparse
from .externals import six


###############################################################################
def clone(estimator, safe=True):
    """Constructs a new estimator with the same parameters.

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
    # XXX: not handling dictionaries
    if estimator_type in (list, tuple, set, frozenset):
        return estimator_type([clone(e, safe=safe) for e in estimator])
    elif not hasattr(estimator, 'get_params'):
        if not safe:
            return copy.deepcopy(estimator)
        else:
            raise TypeError("Cannot clone object '%s' (type %s): "
                            "it does not seem to be a scikit-learn estimator "
                            "it does not implement a 'get_params' methods."
                            % (repr(estimator), type(estimator)))
    klass = estimator.__class__
    new_object_params = estimator.get_params(deep=False)
    for name, param in six.iteritems(new_object_params):
        new_object_params[name] = clone(param, safe=False)
    new_object = klass(**new_object_params)
    params_set = new_object.get_params(deep=False)

    # quick sanity check of the parameters of the clone
    for name in new_object_params:
        param1 = new_object_params[name]
        param2 = params_set[name]
        if isinstance(param1, np.ndarray):
            # For most ndarrays, we do not test for complete equality
            if not isinstance(param2, type(param1)):
                equality_test = False
            elif (param1.ndim > 0
                    and param1.shape[0] > 0
                    and isinstance(param2, np.ndarray)
                    and param2.ndim > 0
                    and param2.shape[0] > 0):
                equality_test = (
                    param1.shape == param2.shape
                    and param1.dtype == param2.dtype
                    # We have to use '.flat' for 2D arrays
                    and param1.flat[0] == param2.flat[0]
                    and param1.flat[-1] == param2.flat[-1]
                )
            else:
                equality_test = np.all(param1 == param2)
        elif sparse.issparse(param1):
            # For sparse matrices equality doesn't work
            if not sparse.issparse(param2):
                equality_test = False
            elif param1.size == 0 or param2.size == 0:
                equality_test = (
                    param1.__class__ == param2.__class__
                    and param1.size == 0
                    and param2.size == 0
                )
            else:
                equality_test = (
                    param1.__class__ == param2.__class__
                    and param1.data[0] == param2.data[0]
                    and param1.data[-1] == param2.data[-1]
                    and param1.nnz == param2.nnz
                    and param1.shape == param2.shape
                )
        else:
            equality_test = new_object_params[name] == params_set[name]
        if not equality_test:
            raise RuntimeError('Cannot clone object %s, as the constructor '
                               'does not seem to set parameter %s' %
                               (estimator, name))

    return new_object


###############################################################################
def _pprint(params, offset=0, printer=repr):
    """Pretty print the dictionary 'params'

    Parameters
    ----------
    params: dict
        The dictionary to pretty print

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
    line_sep = ',\n' + (1 + offset // 2) * ' '
    for i, (k, v) in enumerate(sorted(six.iteritems(params))):
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
            if (this_line_length + len(this_repr) >= 75 or '\n' in this_repr):
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
    """Base class for all estimators in scikit-learn

    Notes
    -----
    All estimators should specify all the parameters that can be set
    at the class level in their ``__init__`` as explicit keyword
    arguments (no ``*args`` or ``**kwargs``).
    """

    @classmethod
    def _get_param_names(cls):
        """Get parameter names for the estimator"""
        # fetch the constructor or the original constructor before
        # deprecation wrapping if any
        init = getattr(cls.__init__, 'deprecated_original', cls.__init__)
        if init is object.__init__:
            # No explicit constructor to introspect
            return []

        # introspect the constructor arguments to find the model parameters
        # to represent
        args, varargs, kw, default = inspect.getargspec(init)
        if varargs is not None:
            raise RuntimeError("scikit-learn estimators should always "
                               "specify their parameters in the signature"
                               " of their __init__ (no varargs)."
                               " %s doesn't follow this convention."
                               % (cls, ))
        # Remove 'self'
        # XXX: This is going to fail if the init is a staticmethod, but
        # who would do this?
        args.pop(0)
        args.sort()
        return args

    def get_params(self, deep=True):
        """Get parameters for this estimator.

        Parameters
        ----------
        deep: boolean, optional
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        out = dict()
        for key in self._get_param_names():
            # We need deprecation warnings to always be on in order to
            # catch deprecated param values.
            # This is set in utils/__init__.py but it gets overwritten
            # when running under python3 somehow.
            warnings.simplefilter("always", DeprecationWarning)
            try:
                with warnings.catch_warnings(record=True) as w:
                    value = getattr(self, key, None)
                if len(w) and w[0].category == DeprecationWarning:
                    # if the parameter is deprecated, don't show it
                    continue
            finally:
                warnings.filters.pop(0)

            # XXX: should we rather test if instance of estimator?
            if deep and hasattr(value, 'get_params'):
                deep_items = value.get_params().items()
                out.update((key + '__' + k, val) for k, val in deep_items)
            out[key] = value
        return out

    def set_params(self, **params):
        """Set the parameters of this estimator.

        The method works on simple estimators as well as on nested objects
        (such as pipelines). The former have parameters of the form
        ``<component>__<parameter>`` so that it's possible to update each
        component of a nested object.

        Returns
        -------
        self
        """
        if not params:
            # Simple optimisation to gain speed (inspect is slow)
            return self
        valid_params = self.get_params(deep=True)
        for key, value in six.iteritems(params):
            split = key.split('__', 1)
            if len(split) > 1:
                # nested objects case
                name, sub_name = split
                if not name in valid_params:
                    raise ValueError('Invalid parameter %s for estimator %s' %
                                     (name, self))
                sub_object = valid_params[name]
                sub_object.set_params(**{sub_name: value})
            else:
                # simple objects case
                if not key in valid_params:
                    raise ValueError('Invalid parameter %s ' 'for estimator %s'
                                     % (key, self.__class__.__name__))
                setattr(self, key, value)
        return self

    def __repr__(self):
        class_name = self.__class__.__name__
        return '%s(%s)' % (class_name, _pprint(self.get_params(deep=False),
                                               offset=len(class_name),),)


###############################################################################
class ClassifierMixin(object):
    """Mixin class for all classifiers in scikit-learn."""

    def score(self, X, y, sample_weight=None):
        """Returns the mean accuracy on the given test data and labels.

        In multi-label classification, this is the subset accuracy
        which is a harsh metric since you require for each sample that
        each label set be correctly predicted.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test samples.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            True labels for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            Mean accuracy of self.predict(X) wrt. y.

        """
        from .metrics import accuracy_score
        return accuracy_score(y, self.predict(X), sample_weight=sample_weight)


###############################################################################
class RegressorMixin(object):
    """Mixin class for all regression estimators in scikit-learn."""

    def score(self, X, y, sample_weight=None):
        """Returns the coefficient of determination R^2 of the prediction.

        The coefficient R^2 is defined as (1 - u/v), where u is the regression
        sum of squares ((y_true - y_pred) ** 2).sum() and v is the residual
        sum of squares ((y_true - y_true.mean()) ** 2).sum().
        Best possible score is 1.0, lower values are worse.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test samples.

        y : array-like, shape = (n_samples) or (n_samples, n_outputs)
            True values for X.

        sample_weight : array-like, shape = [n_samples], optional
            Sample weights.

        Returns
        -------
        score : float
            R^2 of self.predict(X) wrt. y.
        """

        from .metrics import r2_score
        return r2_score(y, self.predict(X), sample_weight=sample_weight)


###############################################################################
class ClusterMixin(object):
    """Mixin class for all cluster estimators in scikit-learn."""
    def fit_predict(self, X, y=None):
        """Performs clustering on X and returns cluster labels.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            Input data.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            cluster labels
        """
        # non-optimized default implementation; override when a better
        # method is possible for a given clustering algorithm
        self.fit(X)
        return self.labels_


class BiclusterMixin(object):
    """Mixin class for all bicluster estimators in scikit-learn"""

    @property
    def biclusters_(self):
        """Convenient way to get row and column indicators together.

        Returns the ``rows_`` and ``columns_`` members.
        """
        return self.rows_, self.columns_

    def get_indices(self, i):
        """Row and column indices of the i'th bicluster.

        Only works if ``rows_`` and ``columns_`` attributes exist.

        Returns
        -------
        row_ind : np.array, dtype=np.intp
            Indices of rows in the dataset that belong to the bicluster.
        col_ind : np.array, dtype=np.intp
            Indices of columns in the dataset that belong to the bicluster.

        """
        rows = self.rows_[i]
        columns = self.columns_[i]
        return np.nonzero(rows)[0], np.nonzero(columns)[0]

    def get_shape(self, i):
        """Shape of the i'th bicluster.

        Returns
        -------
        shape : (int, int)
            Number of rows and columns (resp.) in the bicluster.
        """
        indices = self.get_indices(i)
        return tuple(len(i) for i in indices)

    def get_submatrix(self, i, data):
        """Returns the submatrix corresponding to bicluster `i`.

        Works with sparse matrices. Only works if ``rows_`` and
        ``columns_`` attributes exist.

        """
        from .utils.validation import check_array
        data = check_array(data, accept_sparse='csr')
        row_ind, col_ind = self.get_indices(i)
        return data[row_ind[:, np.newaxis], col_ind]


###############################################################################
class TransformerMixin(object):
    """Mixin class for all transformers in scikit-learn."""

    def fit_transform(self, X, y=None, **fit_params):
        """Fit to data, then transform it.

        Fits transformer to X and y with optional parameters fit_params
        and returns a transformed version of X.

        Parameters
        ----------
        X : numpy array of shape [n_samples, n_features]
            Training set.

        y : numpy array of shape [n_samples]
            Target values.

        Returns
        -------
        X_new : numpy array of shape [n_samples, n_features_new]
            Transformed array.

        """
        # non-optimized default implementation; override when a better
        # method is possible for a given clustering algorithm
        if y is None:
            # fit method of arity 1 (unsupervised transformation)
            return self.fit(X, **fit_params).transform(X)
        else:
            # fit method of arity 2 (supervised transformation)
            return self.fit(X, y, **fit_params).transform(X)


###############################################################################
class MetaEstimatorMixin(object):
    """Mixin class for all meta estimators in scikit-learn."""
    # this is just a tag for the moment


###############################################################################
# XXX: Temporary solution to figure out if an estimator is a classifier

def _get_sub_estimator(estimator):
    """Returns the final estimator if there is any."""
    if hasattr(estimator, 'estimator'):
        # GridSearchCV and other CV-tuned estimators
        return _get_sub_estimator(estimator.estimator)
    if hasattr(estimator, 'steps'):
        # Pipeline
        return _get_sub_estimator(estimator.steps[-1][1])
    return estimator


def is_classifier(estimator):
    """Returns True if the given estimator is (probably) a classifier."""
    estimator = _get_sub_estimator(estimator)
    return isinstance(estimator, ClassifierMixin)
