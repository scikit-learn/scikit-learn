"""
The :mod:`sklearn.experimental._column_transformer` module implements utilities
to work with heterogeneous data and to apply different transformers to
different columns.
"""
# Author: Andreas Mueller
#         Joris Van den Bossche
# License: BSD


import numpy as np
from scipy import sparse

from ..base import clone, TransformerMixin
from ..externals.joblib import Parallel, delayed
from ..externals import six
from ..pipeline import (
    _fit_one_transformer, _fit_transform_one, _transform_one, _name_estimators)
from ..utils import Bunch
from ..utils.metaestimators import _BaseComposition
from ..utils.validation import check_is_fitted


__all__ = ['ColumnTransformer', 'make_column_transformer']


_ERR_MSG_1DCOLUMN = ("1D data passed to a transformer that expects 2D data. "
                     "Try to specify the column selection as a list of one "
                     "item instead of a scalar.")


class ColumnTransformer(_BaseComposition, TransformerMixin):
    """Applies transformers to columns of an array or pandas DataFrame.

    EXPERIMENTAL

    This estimator allows different columns or column subsets of the input
    to be transformed separately and the results combined into a single
    feature space.
    This is useful for heterogeneous or columnar data, to combine several
    feature extraction mechanisms into a single transformer.

    Read more in the :ref:`User Guide <column_transformer>`.

    Parameters
    ----------
    transformers : list of tuples
        List of (name, transformer, column) tuples specifying the transformer
        objects to be applied to subsets of the data.

        name : string
            Like in FeatureUnion and Pipeline, this allows the transformer and
            its parameters to be set using ``set_params`` and searched in grid
            search.
        transformer : estimator or None
            If None, nothing is output for this triple. Estimator must support
            `fit` and `transform`.
        column : string or int, array-like of string or int, slice or boolean \
 mask array
            Indexes the data on its second axis. Integers are interpreted as
            positional columns, while strings can reference DataFrame columns
            by name.  A scalar string or int should be used where
            ``transformer`` expects X to be a 1d array-like (vector),
            otherwise a 2d array will be passed to the transformer.

    passthrough : string or int, array-like of string or int, slice or \
boolean mask array
        Specifies the columns that should be passed through untransformed.
        This selection of columns is concatenated with the output of the
        transformers. See explanation of `column` above for details on the
        format.

    n_jobs : int, optional
        Number of jobs to run in parallel (default 1).

    transformer_weights : dict, optional
        Multiplicative weights for features per transformer. The output of the
        transformer is multiplied by these weights. Keys are transformer names,
        values the weights.

    Attributes
    ----------
    transformers_ : list
        The collection of fitted transformers as tuples of
        (name, fitted_transformer, column).

    named_transformers_ : Bunch object, a dictionary with attribute access
        Read-only attribute to access any transformer by given name.
        Keys are transformer names and values are the fitted transformer
        objects.

    Notes
    -----
    The order of the columns in the transformed feature matrix follows the
    order of how the columns are specified in the `transformers` list.
    Columns of the original feature matrix that are not specified are
    dropped from the resulting transformed feature matrix.

    Examples
    --------
    >>> from sklearn.experimental import ColumnTransformer
    >>> from sklearn.preprocessing import Normalizer
    >>> union = ColumnTransformer(
    ...     [("norm1", Normalizer(norm='l1'), [0, 1]),
    ...      ("norm2", Normalizer(norm='l1'), slice(2, 4))])
    >>> X = np.array([[0., 1., 2., 2.],
    ...               [1., 1., 0., 1.]])
    >>> # Normalizer scales each row of X to unit norm. Therefore, a separate
    >>> # scaling is applied for the two first and two last elements of each
    >>> # row independently.
    >>> union.fit_transform(X)    # doctest: +NORMALIZE_WHITESPACE
    array([[ 0. ,  1. ,  0.5,  0.5],
           [ 0.5,  0.5,  0. ,  1. ]])

    """

    def __init__(self, transformers, passthrough=None, n_jobs=1,
                 transformer_weights=None):
        self.transformers = transformers
        self.passthrough = passthrough
        self.n_jobs = n_jobs
        self.transformer_weights = transformer_weights

    @property
    def _transformers(self):
        """
        Internal list of transformer only containing the name and
        transformers, dropping the columns. This is for the implementation
        of get_params via BaseComposition._get_params which expects lists
        of tuples of len 2.
        """
        return [(name, trans) for name, trans, _ in self.transformers]

    @_transformers.setter
    def _transformers(self, value):
        self.transformers = [
            (name, trans, col) for ((name, trans), (_, _, col))
            in zip(value, self.transformers)]

    def get_params(self, deep=True):
        """Get parameters for this estimator.

        Parameters
        ----------
        deep : boolean, optional
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        return self._get_params('_transformers', deep=deep)

    def set_params(self, **kwargs):
        """Set the parameters of this estimator.

        Valid parameter keys can be listed with ``get_params()``.

        Returns
        -------
        self
        """
        self._set_params('_transformers', **kwargs)
        return self

    def _iter(self, X=None, skip_none=True, fitted=False):
        """Generate (name, trans, column, weight) tuples
        """
        if fitted:
            transformers = self.transformers_
        else:
            transformers = self.transformers
        get_weight = (self.transformer_weights or {}).get
        return ((name, trans, _get_column(X, column) if X is not None else X,
                 get_weight(name))
                for name, trans, column in transformers
                if not skip_none or trans is not None)

    def _validate_transformers(self):
        names, transformers, _, _ = zip(*self._iter(skip_none=False))

        # validate names
        self._validate_names(names)

        # validate estimators
        for t in transformers:
            if t is None:
                continue
            if (not (hasattr(t, "fit") or hasattr(t, "fit_transform")) or not
                    hasattr(t, "transform")):
                raise TypeError("All estimators should implement fit and "
                                "transform. '%s' (type %s) doesn't" %
                                (t, type(t)))

    @property
    def named_transformers_(self):
        """Access the fitted transformer by name.

        Read-only attribute to access any transformer by given name.
        Keys are transformer names and values are the fitted transformer
        objects.

        """
        # Use Bunch object to improve autocomplete
        return Bunch(**dict([(name, trans) for name, trans, _
                             in self.transformers_]))

    def get_feature_names(self):
        """Get feature names from all transformers.

        Returns
        -------
        feature_names : list of strings
            Names of the features produced by transform.
        """
        check_is_fitted(self, 'transformers_')
        feature_names = []
        for name, trans, _, _ in self._iter(fitted=True):
            if not hasattr(trans, 'get_feature_names'):
                raise AttributeError("Transformer %s (type %s) does not "
                                     "provide get_feature_names."
                                     % (str(name), type(trans).__name__))
            feature_names.extend([name + "__" + f for f in
                                  trans.get_feature_names()])
        return feature_names

    def _update_fitted_transformers(self, transformers):
        transformers = iter(transformers)
        self.transformers_ = [
            (name, None if old is None else next(transformers), column)
            for name, old, column in self.transformers
        ]

    def _fit_transform(self, X, y, func, fitted=False):
        """
        Private function to fit and/or transform on demand.

        Return value (transformers and/or transformed X data) depends
        on the passed function.
        ``fitted=True`` ensures the fitted transformers are used.
        """
        try:
            return Parallel(n_jobs=self.n_jobs)(
                delayed(func)(clone(trans) if not fitted else trans,
                              X_sel, y, weight)
                for name, trans, X_sel, weight in self._iter(X=X,
                                                             fitted=fitted))
        except ValueError as e:
            if "Expected 2D array, got 1D array instead" in str(e):
                raise ValueError(_ERR_MSG_1DCOLUMN)
            else:
                raise

    def fit(self, X, y=None):
        """Fit all transformers using X.

        Parameters
        ----------
        X : array-like or DataFrame of shape [n_samples, n_features]
            Input data, of which specified subsets are used to fit the
            transformers.

        y : array-like, shape (n_samples, ...), optional
            Targets for supervised learning.

        Returns
        -------
        self : ColumnTransformer
            This estimator

        """
        self._validate_transformers()

        transformers = self._fit_transform(X, y, _fit_one_transformer)

        self._update_fitted_transformers(transformers)
        return self

    def fit_transform(self, X, y=None):
        """Fit all transformers, transform the data and concatenate results.

        Parameters
        ----------
        X : array-like or DataFrame of shape [n_samples, n_features]
            Input data, of which specified subsets are used to fit the
            transformers.

        y : array-like, shape (n_samples, ...), optional
            Targets for supervised learning.

        Returns
        -------
        X_t : array-like or sparse matrix, shape (n_samples, sum_n_components)
            hstack of results of transformers. sum_n_components is the
            sum of n_components (output dimension) over transformers. If
            any result is a sparse matrix, everything will be converted to
            sparse matrices.

        """
        result = self._fit_transform(X, y, _fit_transform_one)

        if not result:
            # All transformers are None
            return np.zeros((X.shape[0], 0))
        Xs, transformers = zip(*result)

        self._update_fitted_transformers(transformers)

        if self.passthrough:
            Xs = list(Xs) + [_get_column(X, self.passthrough)]

        if any(sparse.issparse(f) for f in Xs):
            Xs = sparse.hstack(Xs).tocsr()
        else:
            Xs = np.hstack(Xs)
        return Xs

    def transform(self, X):
        """Transform X separately by each transformer, concatenate results.

        Parameters
        ----------
        X : array-like or DataFrame of shape [n_samples, n_features]
            Input data, of which specified subsets are used to fit the
            transformers.

        Returns
        -------
        X_t : array-like or sparse matrix, shape (n_samples, sum_n_components)
            hstack of results of transformers. sum_n_components is the
            sum of n_components (output dimension) over transformers. If
            any result is a sparse matrix, everything will be converted to
            sparse matrices.

        """
        check_is_fitted(self, 'transformers_')

        Xs = self._fit_transform(X, None, _transform_one, fitted=True)

        if not Xs:
            # All transformers are None
            if self.passthrough is None:
                return np.zeros((X.shape[0], 0))
            else:
                return _get_column(X, self.passthrough)

        if self.passthrough:
            Xs = list(Xs) + [_get_column(X, self.passthrough)]

        if any(sparse.issparse(f) for f in Xs):
            Xs = sparse.hstack(Xs).tocsr()
        else:
            Xs = np.hstack(Xs)
        return Xs


def _check_key_type(key, superclass):
    """
    Check that scalar, list or slice is of certain type.

    """
    if isinstance(key, superclass):
        return True
    if isinstance(key, slice):
        return (isinstance(key.start, (superclass, type(None))) and
                isinstance(key.stop, (superclass, type(None))))
    if isinstance(key, list):
        return all(isinstance(x, superclass) for x in key)
    if hasattr(key, 'dtype'):
        return key.dtype.kind == 'i'
    return False


def _get_column(X, key):
    """
    Get feature column(s) from input data X.

    Supported input types (X): numpy arrays, sparse arrays and DataFrames

    Supported key types (key):
    - scalar: output is 1D
    - lists, slices, boolean masks: output is 2D

    Supported key data types:

    - integer or boolean mask (positional):
        - supported for arrays, sparse matrices and dataframes
    - string (key-based):
        - only supported for dataframes
        - So no keys other than strings are allowed (while in principle you
          can use any hashable object as key).

    """
    # check whether we have string column names or integers
    if _check_key_type(key, int):
        column_names = False
    elif _check_key_type(key, six.string_types):
        column_names = True
    elif hasattr(key, 'dtype') and np.issubdtype(key.dtype, np.bool):
        # boolean mask
        column_names = False
        if hasattr(X, 'loc'):
            # pandas boolean masks don't work with iloc, so take loc path
            column_names = True
    else:
        raise ValueError("No valid specification of the columns. Only a "
                         "scalar, list or slice of all integers or all "
                         "strings, or boolean mask is allowed")

    if column_names:
        if hasattr(X, 'loc'):
            # pandas dataframes
            return X.loc[:, key]
        else:
            raise ValueError("Specifying the columns using strings is only "
                             "supported for pandas DataFrames")
    else:
        if hasattr(X, 'iloc'):
            # pandas dataframes
            return X.iloc[:, key]
        else:
            # numpy arrays, sparse arrays
            return X[:, key]


def _get_transformer_list(estimators):
    """
    Construct (name, trans, column) tuples from dictionary

    """
    transformers = estimators.keys()
    names = _name_estimators(transformers)

    transformer_list = [(name, trans, estimators[trans])
                        for (name, trans) in names]
    return sorted(transformer_list, key=lambda x: x[0])


def make_column_transformer(transformers, **kwargs):
    """Construct a ColumnTransformer from the given transformers.

    This is a shorthand for the ColumnTransformer constructor; it does not
    require, and does not permit, naming the transformers. Instead, they will
    be given names automatically based on their types. It also does not allow
    weighting.

    Parameters
    ----------
    transformers : dict of estimators
        Dictionary of transformer to column selection.

    n_jobs : int, optional
        Number of jobs to run in parallel (default 1).

    Returns
    -------
    f : ColumnTransformer

    Examples
    --------
    >>> from sklearn.preprocessing import StandardScaler, OneHotEncoder
    >>> from sklearn.experimental import make_column_transformer
    >>> make_column_transformer(
    ...     {StandardScaler(): ['numerical_column'],
    ...      OneHotEncoder(): ['categorical_column']})
    ...     # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    ColumnTransformer(n_jobs=1, passthrough=None, transformer_weights=None,
         transformers=[('onehotencoder',
                        OneHotEncoder(...),
                        ['categorical_column']),
                       ('standardscaler',
                        StandardScaler(...),
                        ['numerical_column'])])

    """
    n_jobs = kwargs.pop('n_jobs', 1)
    if kwargs:
        raise TypeError('Unknown keyword arguments: "{}"'
                        .format(list(kwargs.keys())[0]))
    transformer_list = _get_transformer_list(transformers)
    return ColumnTransformer(transformer_list, n_jobs=n_jobs)
