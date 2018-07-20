"""
The :mod:`sklearn.compose._column_transformer` module implements utilities
to work with heterogeneous data and to apply different transformers to
different columns.
"""
# Author: Andreas Mueller
#         Joris Van den Bossche
# License: BSD
from itertools import chain

import numpy as np
from scipy import sparse

from ..base import clone, TransformerMixin
from ..utils import Parallel, delayed
from ..externals import six
from ..pipeline import (
    _fit_one_transformer, _fit_transform_one, _transform_one, _name_estimators)
from ..preprocessing import FunctionTransformer
from ..utils import Bunch
from ..utils.metaestimators import _BaseComposition
from ..utils.validation import check_is_fitted


__all__ = ['ColumnTransformer', 'make_column_transformer']


_ERR_MSG_1DCOLUMN = ("1D data passed to a transformer that expects 2D data. "
                     "Try to specify the column selection as a list of one "
                     "item instead of a scalar.")


class ColumnTransformer(_BaseComposition, TransformerMixin):
    """Applies transformers to columns of an array or pandas DataFrame.

    EXPERIMENTAL: some behaviors may change between releases without
    deprecation.

    This estimator allows different columns or column subsets of the input
    to be transformed separately and the results combined into a single
    feature space.
    This is useful for heterogeneous or columnar data, to combine several
    feature extraction mechanisms or transformations into a single transformer.

    Read more in the :ref:`User Guide <column_transformer>`.

    .. versionadded:: 0.20

    Parameters
    ----------
    transformers : list of tuples
        List of (name, transformer, column(s)) tuples specifying the
        transformer objects to be applied to subsets of the data.

        name : string
            Like in Pipeline and FeatureUnion, this allows the transformer and
            its parameters to be set using ``set_params`` and searched in grid
            search.
        transformer : estimator or {'passthrough', 'drop'}
            Estimator must support `fit` and `transform`. Special-cased
            strings 'drop' and 'passthrough' are accepted as well, to
            indicate to drop the columns or to pass them through untransformed,
            respectively.
        column(s) : string or int, array-like of string or int, slice, \
boolean mask array or callable
            Indexes the data on its second axis. Integers are interpreted as
            positional columns, while strings can reference DataFrame columns
            by name.  A scalar string or int should be used where
            ``transformer`` expects X to be a 1d array-like (vector),
            otherwise a 2d array will be passed to the transformer.
            A callable is passed the input data `X` and can return any of the
            above.

    remainder : {'passthrough', 'drop'} or estimator, default 'drop'
        By default, only the specified columns in `transformers` are
        transformed and combined in the output, and the non-specified
        columns are dropped. (default of ``'drop'``).
        By specifying ``remainder='passthrough'``, all remaining columns that
        were not specified in `transformers` will be automatically passed
        through. This subset of columns is concatenated with the output of
        the transformers.
        By setting ``remainder`` to be an estimator, the remaining
        non-specified columns will use the ``remainder`` estimator. The
        estimator must support `fit` and `transform`.

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
        (name, fitted_transformer, column). `fitted_transformer` can be an
        estimator, 'drop', or 'passthrough'. If there are remaining columns,
        the final element is a tuple of the form:
        ('remainder', transformer, remaining_columns) corresponding to the
        ``remainder`` parameter. If there are remaining columns, then
        ``len(transformers_)==len(transformers)+1``, otherwise
        ``len(transformers_)==len(transformers)``.

    named_transformers_ : Bunch object, a dictionary with attribute access
        Read-only attribute to access any transformer by given name.
        Keys are transformer names and values are the fitted transformer
        objects.

    Notes
    -----
    The order of the columns in the transformed feature matrix follows the
    order of how the columns are specified in the `transformers` list.
    Columns of the original feature matrix that are not specified are
    dropped from the resulting transformed feature matrix, unless specified
    in the `passthrough` keyword. Those columns specified with `passthrough`
    are added at the right to the output of the transformers.

    See also
    --------
    sklearn.compose.make_column_transformer : convenience function for
        combining the outputs of multiple transformer objects applied to
        column subsets of the original feature space.

    Examples
    --------
    >>> from sklearn.compose import ColumnTransformer
    >>> from sklearn.preprocessing import Normalizer
    >>> ct = ColumnTransformer(
    ...     [("norm1", Normalizer(norm='l1'), [0, 1]),
    ...      ("norm2", Normalizer(norm='l1'), slice(2, 4))])
    >>> X = np.array([[0., 1., 2., 2.],
    ...               [1., 1., 0., 1.]])
    >>> # Normalizer scales each row of X to unit norm. A separate scaling
    >>> # is applied for the two first and two last elements of each
    >>> # row independently.
    >>> ct.fit_transform(X)    # doctest: +NORMALIZE_WHITESPACE
    array([[0. , 1. , 0.5, 0.5],
           [0.5, 0.5, 0. , 1. ]])

    """

    def __init__(self, transformers, remainder='drop', n_jobs=1,
                 transformer_weights=None):
        self.transformers = transformers
        self.remainder = remainder
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

    def _iter(self, X=None, fitted=False, replace_strings=False):
        """Generate (name, trans, column, weight) tuples
        """
        if fitted:
            transformers = self.transformers_
        else:
            transformers = self.transformers
            if self._remainder[2] is not None:
                transformers = chain(transformers, [self._remainder])
        get_weight = (self.transformer_weights or {}).get

        for name, trans, column in transformers:
            sub = None if X is None else _get_column(X, column)

            if replace_strings:
                # replace 'passthrough' with identity transformer and
                # skip in case of 'drop'
                if trans == 'passthrough':
                    trans = FunctionTransformer(
                        validate=False, accept_sparse=True,
                        check_inverse=False)
                elif trans == 'drop':
                    continue

            yield (name, trans, sub, get_weight(name))

    def _validate_transformers(self):
        if not self.transformers:
            return

        names, transformers, _ = zip(*self.transformers)

        # validate names
        self._validate_names(names)

        # validate estimators
        for t in transformers:
            if t in ('drop', 'passthrough'):
                continue
            if (not (hasattr(t, "fit") or hasattr(t, "fit_transform")) or not
                    hasattr(t, "transform")):
                raise TypeError("All estimators should implement fit and "
                                "transform, or can be 'drop' or 'passthrough' "
                                "specifiers. '%s' (type %s) doesn't." %
                                (t, type(t)))

    def _validate_remainder(self, X):
        """
        Validates ``remainder`` and defines ``_remainder`` targeting
        the remaining columns.
        """
        is_transformer = ((hasattr(self.remainder, "fit")
                           or hasattr(self.remainder, "fit_transform"))
                          and hasattr(self.remainder, "transform"))
        if (self.remainder not in ('drop', 'passthrough')
                and not is_transformer):
            raise ValueError(
                "The remainder keyword needs to be one of 'drop', "
                "'passthrough', or estimator. '%s' was passed instead" %
                self.remainder)

        n_columns = X.shape[1]
        cols = []
        for _, _, columns in self.transformers:
            cols.extend(_get_column_indices(X, columns))
        remaining_idx = sorted(list(set(range(n_columns)) - set(cols))) or None

        self._remainder = ('remainder', self.remainder, remaining_idx)

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
            if trans == 'drop':
                continue
            elif trans == 'passthrough':
                raise NotImplementedError(
                    "get_feature_names is not yet supported when using "
                    "a 'passthrough' transformer.")
            elif not hasattr(trans, 'get_feature_names'):
                raise AttributeError("Transformer %s (type %s) does not "
                                     "provide get_feature_names."
                                     % (str(name), type(trans).__name__))
            feature_names.extend([name + "__" + f for f in
                                  trans.get_feature_names()])
        return feature_names

    def _update_fitted_transformers(self, transformers):
        # transformers are fitted; excludes 'drop' cases
        transformers = iter(transformers)
        transformers_ = []

        transformer_iter = self.transformers
        if self._remainder[2] is not None:
            transformer_iter = chain(transformer_iter, [self._remainder])

        for name, old, column in transformer_iter:
            if old == 'drop':
                trans = 'drop'
            elif old == 'passthrough':
                # FunctionTransformer is present in list of transformers,
                # so get next transformer, but save original string
                next(transformers)
                trans = 'passthrough'
            else:
                trans = next(transformers)
            transformers_.append((name, trans, column))

        # sanity check that transformers is exhausted
        assert not list(transformers)
        self.transformers_ = transformers_

    def _validate_output(self, result):
        """
        Ensure that the output of each transformer is 2D. Otherwise
        hstack can raise an error or produce incorrect results.
        """
        names = [name for name, _, _, _ in self._iter(replace_strings=True)]
        for Xs, name in zip(result, names):
            if not getattr(Xs, 'ndim', 0) == 2:
                raise ValueError(
                    "The output of the '{0}' transformer should be 2D (scipy "
                    "matrix, array, or pandas DataFrame).".format(name))

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
                for _, trans, X_sel, weight in self._iter(
                    X=X, fitted=fitted, replace_strings=True))
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
        self._validate_remainder(X)
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
        self._validate_remainder(X)
        self._validate_transformers()

        result = self._fit_transform(X, y, _fit_transform_one)

        if not result:
            # All transformers are None
            return np.zeros((X.shape[0], 0))

        Xs, transformers = zip(*result)

        self._update_fitted_transformers(transformers)
        self._validate_output(Xs)

        return _hstack(list(Xs))

    def transform(self, X):
        """Transform X separately by each transformer, concatenate results.

        Parameters
        ----------
        X : array-like or DataFrame of shape [n_samples, n_features]
            The data to be transformed by subset.

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
        self._validate_output(Xs)

        if not Xs:
            # All transformers are None
            return np.zeros((X.shape[0], 0))

        return _hstack(list(Xs))


def _check_key_type(key, superclass):
    """
    Check that scalar, list or slice is of a certain type.

    This is only used in _get_column and _get_column_indices to check
    if the `key` (column specification) is fully integer or fully string-like.

    Parameters
    ----------
    key : scalar, list, slice, array-like
        The column specification to check
    superclass : int or six.string_types
        The type for which to check the `key`

    """
    if isinstance(key, superclass):
        return True
    if isinstance(key, slice):
        return (isinstance(key.start, (superclass, type(None))) and
                isinstance(key.stop, (superclass, type(None))))
    if isinstance(key, list):
        return all(isinstance(x, superclass) for x in key)
    if hasattr(key, 'dtype'):
        if superclass is int:
            return key.dtype.kind == 'i'
        else:
            # superclass = six.string_types
            return key.dtype.kind in ('O', 'U', 'S')
    return False


def _hstack(X):
    """
    Stacks X horizontally.

    Supports input types (X): list of
        numpy arrays, sparse arrays and DataFrames
    """
    if any(sparse.issparse(f) for f in X):
        return sparse.hstack(X).tocsr()
    else:
        return np.hstack(X)


def _get_column(X, key):
    """
    Get feature column(s) from input data X.

    Supported input types (X): numpy arrays, sparse arrays and DataFrames

    Supported key types (key):
    - scalar: output is 1D
    - lists, slices, boolean masks: output is 2D
    - callable that returns any of the above

    Supported key data types:

    - integer or boolean mask (positional):
        - supported for arrays, sparse matrices and dataframes
    - string (key-based):
        - only supported for dataframes
        - So no keys other than strings are allowed (while in principle you
          can use any hashable object as key).

    """
    if callable(key):
        key = key(X)

    # check whether we have string column names or integers
    if _check_key_type(key, int):
        column_names = False
    elif _check_key_type(key, six.string_types):
        column_names = True
    elif hasattr(key, 'dtype') and np.issubdtype(key.dtype, np.bool_):
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


def _get_column_indices(X, key):
    """
    Get feature column indices for input data X and key.

    For accepted values of `key`, see the docstring of _get_column

    """
    n_columns = X.shape[1]

    if callable(key):
        key = key(X)

    if _check_key_type(key, int):
        if isinstance(key, int):
            return [key]
        elif isinstance(key, slice):
            return list(range(n_columns)[key])
        else:
            return list(key)

    elif _check_key_type(key, six.string_types):
        try:
            all_columns = list(X.columns)
        except AttributeError:
            raise ValueError("Specifying the columns using strings is only "
                             "supported for pandas DataFrames")
        if isinstance(key, six.string_types):
            columns = [key]
        elif isinstance(key, slice):
            start, stop = key.start, key.stop
            if start is not None:
                start = all_columns.index(start)
            if stop is not None:
                # pandas indexing with strings is endpoint included
                stop = all_columns.index(stop) + 1
            else:
                stop = n_columns + 1
            return list(range(n_columns)[slice(start, stop)])
        else:
            columns = list(key)

        return [all_columns.index(col) for col in columns]

    elif hasattr(key, 'dtype') and np.issubdtype(key.dtype, np.bool_):
        # boolean mask
        return list(np.arange(n_columns)[key])
    else:
        raise ValueError("No valid specification of the columns. Only a "
                         "scalar, list or slice of all integers or all "
                         "strings, or boolean mask is allowed")


def _get_transformer_list(estimators):
    """
    Construct (name, trans, column) tuples from list

    """
    transformers = [trans[1] for trans in estimators]
    columns = [trans[0] for trans in estimators]
    names = [trans[0] for trans in _name_estimators(transformers)]

    transformer_list = list(zip(names, transformers, columns))
    return transformer_list


def make_column_transformer(*transformers, **kwargs):
    """Construct a ColumnTransformer from the given transformers.

    This is a shorthand for the ColumnTransformer constructor; it does not
    require, and does not permit, naming the transformers. Instead, they will
    be given names automatically based on their types. It also does not allow
    weighting.

    Parameters
    ----------
    *transformers : tuples of column selections and transformers

    remainder : {'passthrough', 'drop'} or estimator, default 'passthrough'
        By default, all remaining columns that were not specified in
        `transformers` will be automatically passed through (default of
        ``'passthrough'``). This subset of columns is concatenated with the
        output of the transformers.
        By using ``remainder='drop'``, only the specified columns in
        `transformers` are transformed and combined in the output, and the
        non-specified columns are dropped.
        By setting ``remainder`` to be an estimator, the remaining
        non-specified columns will use the ``remainder`` estimator. The
        estimator must support `fit` and `transform`.

    n_jobs : int, optional
        Number of jobs to run in parallel (default 1).

    Returns
    -------
    ct : ColumnTransformer

    See also
    --------
    sklearn.compose.ColumnTransformer : Class that allows combining the
        outputs of multiple transformer objects used on column subsets
        of the data into a single feature space.

    Examples
    --------
    >>> from sklearn.preprocessing import StandardScaler, OneHotEncoder
    >>> from sklearn.compose import make_column_transformer
    >>> make_column_transformer(
    ...     (['numerical_column'], StandardScaler()),
    ...     (['categorical_column'], OneHotEncoder()))
    ...     # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    ColumnTransformer(n_jobs=1, remainder='drop', transformer_weights=None,
             transformers=[('standardscaler',
                            StandardScaler(...),
                            ['numerical_column']),
                           ('onehotencoder',
                            OneHotEncoder(...),
                            ['categorical_column'])])

    """
    n_jobs = kwargs.pop('n_jobs', 1)
    remainder = kwargs.pop('remainder', 'drop')
    if kwargs:
        raise TypeError('Unknown keyword arguments: "{}"'
                        .format(list(kwargs.keys())[0]))
    transformer_list = _get_transformer_list(transformers)
    return ColumnTransformer(transformer_list, n_jobs=n_jobs,
                             remainder=remainder)
