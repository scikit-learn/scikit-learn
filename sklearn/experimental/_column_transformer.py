"""
The :mod:`sklearn.experimental._column_transformer` module implements utilities
to work with heterogeneous data and to apply different transformers to
different columns.
"""
# Author: Andreas Mueller
#         Joris Van den Bossche
# License: BSD


import numpy as np

from ..externals import six
from ..pipeline import FeatureUnion, _name_estimators

__all__ = ['ColumnTransformer', 'make_column_transformer']


_ERR_MSG_1DCOLUMN = ("1D data passed to a transformer that expects 2D data. "
                     "Try to specify the column selection as a list of one "
                     "item instead of a scalar.")


class ColumnTransformer(FeatureUnion):
    """Applies transformers to columns of an array or pandas DataFrame
    (EXPERIMENTAL).

    This estimator applies transformer objects to columns of the input,
    then concatenates the results. This is useful for heterogeneous or
    columnar data, to combine several feature extraction mechanisms into a
    single transformer.

    Read more in the :ref:`User Guide <column_transformer>`.

    Parameters
    ----------
    transformer_list : list of tuples
        List of (name, transformer, column) tuples specifying the transformer
        objects to be applied to subsets of the data. The columns can be
        specified as a scalar or slice/list (for multiple columns) of integer
        or string values, or a boolean mask. Integers are interpreted as the
        positional columns, strings as the column labels of `X`.
        When passing a single column to a transformer that expects 2D input
        data, the column should be specified a list of one element.
        Subselections of a pandas DataFrames are passed through to the
        the transformers as pandas objects.

    n_jobs : int, optional
        Number of jobs to run in parallel (default 1).

    transformer_weights : dict, optional
        Multiplicative weights for features per transformer.
        Keys are transformer names, values the weights.

    Examples
    --------
    >>> from sklearn.experimental import ColumnTransformer
    >>> from sklearn.preprocessing import Normalizer
    >>> union = ColumnTransformer(
    ...     [("norm1", Normalizer(norm='l1'), [0, 1]),
    ...      ("norm2", Normalizer(norm='l1'), slice(2, 4))])
    >>> X = np.array([[0., 1., 2., 2.],
    ...               [1., 1., 0., 1.]])
    >>> union.fit_transform(X)    # doctest: +NORMALIZE_WHITESPACE
    array([[ 0. ,  1. ,  0.5,  0.5],
           [ 0.5,  0.5,  0. ,  1. ]])

    """

    @property
    def _transformers(self):
        return [(name, trans) for name, trans, col in self.transformer_list]

    @_transformers.setter
    def _transformers(self, value):
        self.transformer_list = [
            (name, trans, col) for ((name, trans), (_, _, col))
            in zip(value, self.transformer_list)]

    def _iter(self, X=None, skip_none=True):
        """Generate (name, trans, column, weight) tuples
        """
        get_weight = (self.transformer_weights or {}).get
        return ((name, trans, _get_column(X, column), get_weight(name))
                for name, trans, column in self.transformer_list
                if not skip_none or trans is not None)

    def _update_transformer_list(self, transformers):
        transformers = iter(transformers)
        self.transformer_list[:] = [
            (name, None if old is None else next(transformers), column)
            for name, old, column in self.transformer_list
        ]

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
        try:
            return super(ColumnTransformer, self).fit(X, y=y)
        except ValueError as e:
            if "Got X with X.ndim=1. Reshape your data" in str(e):
                raise ValueError(_ERR_MSG_1DCOLUMN)
            else:
                raise

    def fit_transform(self, X, y=None, **fit_params):
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
            one result is a sparse matrix, everything will be converted to
            sparse matrices.
        """
        try:
            return super(ColumnTransformer, self).fit_transform(X, y=y,
                                                                **fit_params)
        except ValueError as e:
            if "Got X with X.ndim=1. Reshape your data" in str(e):
                raise ValueError(_ERR_MSG_1DCOLUMN)
            else:
                raise

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
            one result is a sparse matrix, everything will be converted to
            sparse matrices.
        """
        try:
            return super(ColumnTransformer, self).transform(X)
        except ValueError as e:
            if "Got X with X.ndim=1. Reshape your data" in str(e):
                raise ValueError(_ERR_MSG_1DCOLUMN)
            else:
                raise


def _get_column(X, key):
    """
    Get feature column(s) from input data X.

    Supported input types (X): numpy arrays, sparse arrays and dataframes

    Supported key types (key):
    - scalar: output is 1D
    - lists, slices, boolean masks: output is 2D

    Supported key data types:

    - integer or boolean mask (positional):
        - supported for (sparse) arrays or dataframes
    - string (key-based):
        - only supported for dataframes
        - So no keys other than strings are allowed (while in principle you
          can use any hashable object as key).

    """
    if X is None:
        return X

    # check whether we have string column names or integers
    if (isinstance(key, int)
            or (isinstance(key, list)
                and all(isinstance(col, int) for col in key))
            or (isinstance(key, slice)
                and isinstance(key.start, (int, type(None)))
                and isinstance(key.stop, (int, type(None))))):
        column_names = False
    elif (isinstance(key, six.string_types)
            or (isinstance(key, list)
                and all(isinstance(col, six.string_types) for col in key))
            or (isinstance(key, slice)
                and isinstance(key.start, (six.string_types, type(None)))
                and isinstance(key.stop, (six.string_types, type(None))))):
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
    ColumnTransformer(n_jobs=1,
         transformer_list=[('onehotencoder',
                            OneHotEncoder(...),
                            ['categorical_column']),
                           ('standardscaler',
                            StandardScaler(...),
                            ['numerical_column'])],
         transformer_weights=None)

    """
    n_jobs = kwargs.pop('n_jobs', 1)
    if kwargs:
        raise TypeError('Unknown keyword arguments: "{}"'
                        .format(list(kwargs.keys())[0]))
    transformer_list = _get_transformer_list(transformers)
    return ColumnTransformer(transformer_list, n_jobs=n_jobs)
