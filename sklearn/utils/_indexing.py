import numbers
import sys
import warnings
from collections import UserList
from itertools import compress, islice

import numpy as np
from scipy.sparse import issparse

from ._array_api import _is_numpy_namespace, get_namespace
from ._param_validation import Interval, validate_params
from .extmath import _approximate_mode
from .validation import (
    _is_arraylike_not_scalar,
    _is_pandas_df,
    _is_polars_df_or_series,
    _use_interchange_protocol,
    check_array,
    check_consistent_length,
    check_random_state,
)


def _array_indexing(array, key, key_dtype, axis):
    """Index an array or scipy.sparse consistently across NumPy version."""
    xp, is_array_api = get_namespace(array)
    if is_array_api:
        return xp.take(array, key, axis=axis)
    if issparse(array) and key_dtype == "bool":
        key = np.asarray(key)
    if isinstance(key, tuple):
        key = list(key)
    return array[key, ...] if axis == 0 else array[:, key]


def _pandas_indexing(X, key, key_dtype, axis):
    """Index a pandas dataframe or a series."""
    if _is_arraylike_not_scalar(key):
        key = np.asarray(key)

    if key_dtype == "int" and not (isinstance(key, slice) or np.isscalar(key)):
        # using take() instead of iloc[] ensures the return value is a "proper"
        # copy that will not raise SettingWithCopyWarning
        return X.take(key, axis=axis)
    else:
        # check whether we should index with loc or iloc
        indexer = X.iloc if key_dtype == "int" else X.loc
        return indexer[:, key] if axis else indexer[key]


def _list_indexing(X, key, key_dtype):
    """Index a Python list."""
    if np.isscalar(key) or isinstance(key, slice):
        # key is a slice or a scalar
        return X[key]
    if key_dtype == "bool":
        # key is a boolean array-like
        return list(compress(X, key))
    # key is a integer array-like of key
    return [X[idx] for idx in key]


def _polars_indexing(X, key, key_dtype, axis):
    """Indexing X with polars interchange protocol."""
    # Polars behavior is more consistent with lists
    if isinstance(key, np.ndarray):
        # Convert each element of the array to a Python scalar
        key = key.tolist()
    elif not (np.isscalar(key) or isinstance(key, slice)):
        key = list(key)

    if axis == 1:
        # Here we are certain to have a polars DataFrame; which can be indexed with
        # integer and string scalar, and list of integer, string and boolean
        return X[:, key]

    if key_dtype == "bool":
        # Boolean mask can be indexed in the same way for Series and DataFrame (axis=0)
        return X.filter(key)

    # Integer scalar and list of integer can be indexed in the same way for Series and
    # DataFrame (axis=0)
    X_indexed = X[key]
    if np.isscalar(key) and len(X.shape) == 2:
        # `X_indexed` is a DataFrame with a single row; we return a Series to be
        # consistent with pandas
        pl = sys.modules["polars"]
        return pl.Series(X_indexed.row(0))
    return X_indexed


def _determine_key_type(key, accept_slice=True):
    """Determine the data type of key.

    Parameters
    ----------
    key : scalar, slice or array-like
        The key from which we want to infer the data type.

    accept_slice : bool, default=True
        Whether or not to raise an error if the key is a slice.

    Returns
    -------
    dtype : {'int', 'str', 'bool', None}
        Returns the data type of key.
    """
    err_msg = (
        "No valid specification of the columns. Only a scalar, list or "
        "slice of all integers or all strings, or boolean mask is "
        "allowed"
    )

    dtype_to_str = {int: "int", str: "str", bool: "bool", np.bool_: "bool"}
    array_dtype_to_str = {
        "i": "int",
        "u": "int",
        "b": "bool",
        "O": "str",
        "U": "str",
        "S": "str",
    }

    if key is None:
        return None
    if isinstance(key, tuple(dtype_to_str.keys())):
        try:
            return dtype_to_str[type(key)]
        except KeyError:
            raise ValueError(err_msg)
    if isinstance(key, slice):
        if not accept_slice:
            raise TypeError(
                "Only array-like or scalar are supported. A Python slice was given."
            )
        if key.start is None and key.stop is None:
            return None
        key_start_type = _determine_key_type(key.start)
        key_stop_type = _determine_key_type(key.stop)
        if key_start_type is not None and key_stop_type is not None:
            if key_start_type != key_stop_type:
                raise ValueError(err_msg)
        if key_start_type is not None:
            return key_start_type
        return key_stop_type
    # TODO(1.9) remove UserList when the force_int_remainder_cols param
    # of ColumnTransformer is removed
    if isinstance(key, (list, tuple, UserList)):
        unique_key = set(key)
        key_type = {_determine_key_type(elt) for elt in unique_key}
        if not key_type:
            return None
        if len(key_type) != 1:
            raise ValueError(err_msg)
        return key_type.pop()
    if hasattr(key, "dtype"):
        xp, is_array_api = get_namespace(key)
        # NumPy arrays are special-cased in their own branch because the Array API
        # cannot handle object/string-based dtypes that are often used to index
        # columns of dataframes by names.
        if is_array_api and not _is_numpy_namespace(xp):
            if xp.isdtype(key.dtype, "bool"):
                return "bool"
            elif xp.isdtype(key.dtype, "integral"):
                return "int"
            else:
                raise ValueError(err_msg)
        else:
            try:
                return array_dtype_to_str[key.dtype.kind]
            except KeyError:
                raise ValueError(err_msg)
    raise ValueError(err_msg)


def _safe_indexing(X, indices, *, axis=0):
    """Return rows, items or columns of X using indices.

    .. warning::

        This utility is documented, but **private**. This means that
        backward compatibility might be broken without any deprecation
        cycle.

    Parameters
    ----------
    X : array-like, sparse-matrix, list, pandas.DataFrame, pandas.Series
        Data from which to sample rows, items or columns. `list` are only
        supported when `axis=0`.
    indices : bool, int, str, slice, array-like
        - If `axis=0`, boolean and integer array-like, integer slice,
          and scalar integer are supported.
        - If `axis=1`:
            - to select a single column, `indices` can be of `int` type for
              all `X` types and `str` only for dataframe. The selected subset
              will be 1D, unless `X` is a sparse matrix in which case it will
              be 2D.
            - to select multiples columns, `indices` can be one of the
              following: `list`, `array`, `slice`. The type used in
              these containers can be one of the following: `int`, 'bool' and
              `str`. However, `str` is only supported when `X` is a dataframe.
              The selected subset will be 2D.
    axis : int, default=0
        The axis along which `X` will be subsampled. `axis=0` will select
        rows while `axis=1` will select columns.

    Returns
    -------
    subset
        Subset of X on axis 0 or 1.

    Notes
    -----
    CSR, CSC, and LIL sparse matrices are supported. COO sparse matrices are
    not supported.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils import _safe_indexing
    >>> data = np.array([[1, 2], [3, 4], [5, 6]])
    >>> _safe_indexing(data, 0, axis=0)  # select the first row
    array([1, 2])
    >>> _safe_indexing(data, 0, axis=1)  # select the first column
    array([1, 3, 5])
    """
    if indices is None:
        return X

    if axis not in (0, 1):
        raise ValueError(
            "'axis' should be either 0 (to index rows) or 1 (to index "
            " column). Got {} instead.".format(axis)
        )

    indices_dtype = _determine_key_type(indices)

    if axis == 0 and indices_dtype == "str":
        raise ValueError("String indexing is not supported with 'axis=0'")

    if axis == 1 and isinstance(X, list):
        raise ValueError("axis=1 is not supported for lists")

    if axis == 1 and hasattr(X, "shape") and len(X.shape) != 2:
        raise ValueError(
            "'X' should be a 2D NumPy array, 2D sparse matrix or "
            "dataframe when indexing the columns (i.e. 'axis=1'). "
            "Got {} instead with {} dimension(s).".format(type(X), len(X.shape))
        )

    if (
        axis == 1
        and indices_dtype == "str"
        and not (_is_pandas_df(X) or _use_interchange_protocol(X))
    ):
        raise ValueError(
            "Specifying the columns using strings is only supported for dataframes."
        )

    if hasattr(X, "iloc"):
        # TODO: we should probably use _is_pandas_df_or_series(X) instead but this
        # would require updating some tests such as test_train_test_split_mock_pandas.
        return _pandas_indexing(X, indices, indices_dtype, axis=axis)
    elif _is_polars_df_or_series(X):
        return _polars_indexing(X, indices, indices_dtype, axis=axis)
    elif hasattr(X, "shape"):
        return _array_indexing(X, indices, indices_dtype, axis=axis)
    else:
        return _list_indexing(X, indices, indices_dtype)


def _safe_assign(X, values, *, row_indexer=None, column_indexer=None):
    """Safe assignment to a numpy array, sparse matrix, or pandas dataframe.

    Parameters
    ----------
    X : {ndarray, sparse-matrix, dataframe}
        Array to be modified. It is expected to be 2-dimensional.

    values : ndarray
        The values to be assigned to `X`.

    row_indexer : array-like, dtype={int, bool}, default=None
        A 1-dimensional array to select the rows of interest. If `None`, all
        rows are selected.

    column_indexer : array-like, dtype={int, bool}, default=None
        A 1-dimensional array to select the columns of interest. If `None`, all
        columns are selected.
    """
    row_indexer = slice(None, None, None) if row_indexer is None else row_indexer
    column_indexer = (
        slice(None, None, None) if column_indexer is None else column_indexer
    )

    if hasattr(X, "iloc"):  # pandas dataframe
        with warnings.catch_warnings():
            # pandas >= 1.5 raises a warning when using iloc to set values in a column
            # that does not have the same type as the column being set. It happens
            # for instance when setting a categorical column with a string.
            # In the future the behavior won't change and the warning should disappear.
            # TODO(1.3): check if the warning is still raised or remove the filter.
            warnings.simplefilter("ignore", FutureWarning)
            X.iloc[row_indexer, column_indexer] = values
    else:  # numpy array or sparse matrix
        X[row_indexer, column_indexer] = values


def _get_column_indices_for_bool_or_int(key, n_columns):
    # Convert key into list of positive integer indexes
    try:
        idx = _safe_indexing(np.arange(n_columns), key)
    except IndexError as e:
        raise ValueError(
            f"all features must be in [0, {n_columns - 1}] or [-{n_columns}, 0]"
        ) from e
    return np.atleast_1d(idx).tolist()


def _get_column_indices(X, key):
    """Get feature column indices for input data X and key.

    For accepted values of `key`, see the docstring of
    :func:`_safe_indexing`.
    """
    key_dtype = _determine_key_type(key)
    if _use_interchange_protocol(X):
        return _get_column_indices_interchange(X.__dataframe__(), key, key_dtype)

    n_columns = X.shape[1]
    if isinstance(key, (list, tuple)) and not key:
        # we get an empty list
        return []
    elif key_dtype in ("bool", "int"):
        return _get_column_indices_for_bool_or_int(key, n_columns)
    else:
        try:
            all_columns = X.columns
        except AttributeError:
            raise ValueError(
                "Specifying the columns using strings is only supported for dataframes."
            )
        if isinstance(key, str):
            columns = [key]
        elif isinstance(key, slice):
            start, stop = key.start, key.stop
            if start is not None:
                start = all_columns.get_loc(start)
            if stop is not None:
                # pandas indexing with strings is endpoint included
                stop = all_columns.get_loc(stop) + 1
            else:
                stop = n_columns + 1
            return list(islice(range(n_columns), start, stop))
        else:
            columns = list(key)

        try:
            column_indices = []
            for col in columns:
                col_idx = all_columns.get_loc(col)
                if not isinstance(col_idx, numbers.Integral):
                    raise ValueError(
                        f"Selected columns, {columns}, are not unique in dataframe"
                    )
                column_indices.append(col_idx)

        except KeyError as e:
            raise ValueError("A given column is not a column of the dataframe") from e

        return column_indices


def _get_column_indices_interchange(X_interchange, key, key_dtype):
    """Same as _get_column_indices but for X with __dataframe__ protocol."""

    n_columns = X_interchange.num_columns()

    if isinstance(key, (list, tuple)) and not key:
        # we get an empty list
        return []
    elif key_dtype in ("bool", "int"):
        return _get_column_indices_for_bool_or_int(key, n_columns)
    else:
        column_names = list(X_interchange.column_names())

        if isinstance(key, slice):
            if key.step not in [1, None]:
                raise NotImplementedError("key.step must be 1 or None")
            start, stop = key.start, key.stop
            if start is not None:
                start = column_names.index(start)

            if stop is not None:
                stop = column_names.index(stop) + 1
            else:
                stop = n_columns + 1
            return list(islice(range(n_columns), start, stop))

        selected_columns = [key] if np.isscalar(key) else key

        try:
            return [column_names.index(col) for col in selected_columns]
        except ValueError as e:
            raise ValueError("A given column is not a column of the dataframe") from e


@validate_params(
    {
        "replace": ["boolean"],
        "n_samples": [Interval(numbers.Integral, 1, None, closed="left"), None],
        "random_state": ["random_state"],
        "stratify": ["array-like", "sparse matrix", None],
    },
    prefer_skip_nested_validation=True,
)
def resample(*arrays, replace=True, n_samples=None, random_state=None, stratify=None):
    """Resample arrays or sparse matrices in a consistent way.

    The default strategy implements one step of the bootstrapping
    procedure.

    Parameters
    ----------
    *arrays : sequence of array-like of shape (n_samples,) or \
            (n_samples, n_outputs)
        Indexable data-structures can be arrays, lists, dataframes or scipy
        sparse matrices with consistent first dimension.

    replace : bool, default=True
        Implements resampling with replacement. If False, this will implement
        (sliced) random permutations.

    n_samples : int, default=None
        Number of samples to generate. If left to None this is
        automatically set to the first dimension of the arrays.
        If replace is False it should not be larger than the length of
        arrays.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for shuffling
        the data.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    stratify : {array-like, sparse matrix} of shape (n_samples,) or \
            (n_samples, n_outputs), default=None
        If not None, data is split in a stratified fashion, using this as
        the class labels.

    Returns
    -------
    resampled_arrays : sequence of array-like of shape (n_samples,) or \
            (n_samples, n_outputs)
        Sequence of resampled copies of the collections. The original arrays
        are not impacted.

    See Also
    --------
    shuffle : Shuffle arrays or sparse matrices in a consistent way.

    Examples
    --------
    It is possible to mix sparse and dense arrays in the same run::

      >>> import numpy as np
      >>> X = np.array([[1., 0.], [2., 1.], [0., 0.]])
      >>> y = np.array([0, 1, 2])

      >>> from scipy.sparse import coo_matrix
      >>> X_sparse = coo_matrix(X)

      >>> from sklearn.utils import resample
      >>> X, X_sparse, y = resample(X, X_sparse, y, random_state=0)
      >>> X
      array([[1., 0.],
             [2., 1.],
             [1., 0.]])

      >>> X_sparse
      <3x2 sparse matrix of type '<... 'numpy.float64'>'
          with 4 stored elements in Compressed Sparse Row format>

      >>> X_sparse.toarray()
      array([[1., 0.],
             [2., 1.],
             [1., 0.]])

      >>> y
      array([0, 1, 0])

      >>> resample(y, n_samples=2, random_state=0)
      array([0, 1])

    Example using stratification::

      >>> y = [0, 0, 1, 1, 1, 1, 1, 1, 1]
      >>> resample(y, n_samples=5, replace=False, stratify=y,
      ...          random_state=0)
      [1, 1, 1, 0, 1]
    """
    max_n_samples = n_samples
    random_state = check_random_state(random_state)

    if len(arrays) == 0:
        return None

    first = arrays[0]
    n_samples = first.shape[0] if hasattr(first, "shape") else len(first)

    if max_n_samples is None:
        max_n_samples = n_samples
    elif (max_n_samples > n_samples) and (not replace):
        raise ValueError(
            "Cannot sample %d out of arrays with dim %d when replace is False"
            % (max_n_samples, n_samples)
        )

    check_consistent_length(*arrays)

    if stratify is None:
        if replace:
            indices = random_state.randint(0, n_samples, size=(max_n_samples,))
        else:
            indices = np.arange(n_samples)
            random_state.shuffle(indices)
            indices = indices[:max_n_samples]
    else:
        # Code adapted from StratifiedShuffleSplit()
        y = check_array(stratify, ensure_2d=False, dtype=None)
        if y.ndim == 2:
            # for multi-label y, map each distinct row to a string repr
            # using join because str(row) uses an ellipsis if len(row) > 1000
            y = np.array([" ".join(row.astype("str")) for row in y])

        classes, y_indices = np.unique(y, return_inverse=True)
        n_classes = classes.shape[0]

        class_counts = np.bincount(y_indices)

        # Find the sorted list of instances for each class:
        # (np.unique above performs a sort, so code is O(n logn) already)
        class_indices = np.split(
            np.argsort(y_indices, kind="mergesort"), np.cumsum(class_counts)[:-1]
        )

        n_i = _approximate_mode(class_counts, max_n_samples, random_state)

        indices = []

        for i in range(n_classes):
            indices_i = random_state.choice(class_indices[i], n_i[i], replace=replace)
            indices.extend(indices_i)

        indices = random_state.permutation(indices)

    # convert sparse matrices to CSR for row-based indexing
    arrays = [a.tocsr() if issparse(a) else a for a in arrays]
    resampled_arrays = [_safe_indexing(a, indices) for a in arrays]
    if len(resampled_arrays) == 1:
        # syntactic sugar for the unit argument case
        return resampled_arrays[0]
    else:
        return resampled_arrays


def shuffle(*arrays, random_state=None, n_samples=None):
    """Shuffle arrays or sparse matrices in a consistent way.

    This is a convenience alias to ``resample(*arrays, replace=False)`` to do
    random permutations of the collections.

    Parameters
    ----------
    *arrays : sequence of indexable data-structures
        Indexable data-structures can be arrays, lists, dataframes or scipy
        sparse matrices with consistent first dimension.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for shuffling
        the data.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    n_samples : int, default=None
        Number of samples to generate. If left to None this is
        automatically set to the first dimension of the arrays.  It should
        not be larger than the length of arrays.

    Returns
    -------
    shuffled_arrays : sequence of indexable data-structures
        Sequence of shuffled copies of the collections. The original arrays
        are not impacted.

    See Also
    --------
    resample : Resample arrays or sparse matrices in a consistent way.

    Examples
    --------
    It is possible to mix sparse and dense arrays in the same run::

      >>> import numpy as np
      >>> X = np.array([[1., 0.], [2., 1.], [0., 0.]])
      >>> y = np.array([0, 1, 2])

      >>> from scipy.sparse import coo_matrix
      >>> X_sparse = coo_matrix(X)

      >>> from sklearn.utils import shuffle
      >>> X, X_sparse, y = shuffle(X, X_sparse, y, random_state=0)
      >>> X
      array([[0., 0.],
             [2., 1.],
             [1., 0.]])

      >>> X_sparse
      <3x2 sparse matrix of type '<... 'numpy.float64'>'
          with 3 stored elements in Compressed Sparse Row format>

      >>> X_sparse.toarray()
      array([[0., 0.],
             [2., 1.],
             [1., 0.]])

      >>> y
      array([2, 1, 0])

      >>> shuffle(y, n_samples=2, random_state=0)
      array([0, 1])
    """
    return resample(
        *arrays, replace=False, n_samples=n_samples, random_state=random_state
    )
