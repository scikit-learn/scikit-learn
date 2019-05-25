# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
#          Joris Van den Bossche <jorisvandenbossche@gmail.com>
# License: BSD 3 clause

import numbers
import warnings

import numpy as np
from scipy import sparse

from .. import get_config as _get_config
from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils import deprecated
from ..utils.fixes import _argmax, _object_dtype_isnan
from ..utils.validation import check_is_fitted

from .base import _transform_selected
from .label import _encode, _encode_check_unknown


__all__ = [
    'OneHotEncoder',
    'OrdinalEncoder'
]


class _BaseEncoder(BaseEstimator, TransformerMixin):
    """
    Base class for encoders that includes the code to categorize and
    transform the input features.

    """

    def _check_X(self, X):
        """
        Perform custom check_array:
        - convert list of strings to object dtype
        - check for missing values for object dtype data (check_array does
          not do that)
        - return list of features (arrays): this list of features is
          constructed feature by feature to preserve the data types
          of pandas DataFrame columns, as otherwise information is lost
          and cannot be used, eg for the `categories_` attribute.

        """
        if not (hasattr(X, 'iloc') and getattr(X, 'ndim', 0) == 2):
            # if not a dataframe, do normal check_array validation
            X_temp = check_array(X, dtype=None)
            if (not hasattr(X, 'dtype')
                    and np.issubdtype(X_temp.dtype, np.str_)):
                X = check_array(X, dtype=np.object)
            else:
                X = X_temp
            needs_validation = False
        else:
            # pandas dataframe, do validation later column by column, in order
            # to keep the dtype information to be used in the encoder.
            needs_validation = True

        n_samples, n_features = X.shape
        X_columns = []

        for i in range(n_features):
            Xi = self._get_feature(X, feature_idx=i)
            Xi = check_array(Xi, ensure_2d=False, dtype=None,
                             force_all_finite=needs_validation)
            X_columns.append(Xi)

        return X_columns, n_samples, n_features

    def _get_feature(self, X, feature_idx):
        if hasattr(X, 'iloc'):
            # pandas dataframes
            return X.iloc[:, feature_idx]
        # numpy arrays, sparse arrays
        return X[:, feature_idx]

    def _fit(self, X, handle_unknown='error'):
        X_list, n_samples, n_features = self._check_X(X)

        if self._categories != 'auto':
            if len(self._categories) != n_features:
                raise ValueError("Shape mismatch: if categories is an array,"
                                 " it has to be of shape (n_features,).")

        self.categories_ = []

        for i in range(n_features):
            Xi = X_list[i]
            if self._categories == 'auto':
                cats = _encode(Xi)
            else:
                cats = np.array(self._categories[i], dtype=Xi.dtype)
                if Xi.dtype != object:
                    if not np.all(np.sort(cats) == cats):
                        raise ValueError("Unsorted categories are not "
                                         "supported for numerical categories")
                if handle_unknown == 'error':
                    diff = _encode_check_unknown(Xi, cats)
                    if diff:
                        msg = ("Found unknown categories {0} in column {1}"
                               " during fit".format(diff, i))
                        raise ValueError(msg)
            self.categories_.append(cats)

    def _transform(self, X, handle_unknown='error'):
        X_list, n_samples, n_features = self._check_X(X)

        X_int = np.zeros((n_samples, n_features), dtype=np.int)
        X_mask = np.ones((n_samples, n_features), dtype=np.bool)

        for i in range(n_features):
            Xi = X_list[i]
            diff, valid_mask = _encode_check_unknown(Xi, self.categories_[i],
                                                     return_mask=True)

            if not np.all(valid_mask):
                if handle_unknown == 'error':
                    msg = ("Found unknown categories {0} in column {1}"
                           " during transform".format(diff, i))
                    raise ValueError(msg)
                else:
                    # Set the problematic rows to an acceptable value and
                    # continue `The rows are marked `X_mask` and will be
                    # removed later.
                    X_mask[:, i] = valid_mask
                    # cast Xi into the largest string type necessary
                    # to handle different lengths of numpy strings
                    if (self.categories_[i].dtype.kind in ('U', 'S')
                            and self.categories_[i].itemsize > Xi.itemsize):
                        Xi = Xi.astype(self.categories_[i].dtype)
                    else:
                        Xi = Xi.copy()

                    Xi[~valid_mask] = self.categories_[i][0]
            # We use check_unknown=False, since _encode_check_unknown was
            # already called above.
            _, encoded = _encode(Xi, self.categories_[i], encode=True,
                                 check_unknown=False)
            X_int[:, i] = encoded

        return X_int, X_mask


class OneHotEncoder(_BaseEncoder):
    """Encode categorical integer features as a one-hot numeric array.

    The input to this transformer should be an array-like of integers or
    strings, denoting the values taken on by categorical (discrete) features.
    The features are encoded using a one-hot (aka 'one-of-K' or 'dummy')
    encoding scheme. This creates a binary column for each category and
    returns a sparse matrix or dense array.

    By default, the encoder derives the categories based on the unique values
    in each feature. Alternatively, you can also specify the `categories`
    manually.
    The OneHotEncoder previously assumed that the input features take on
    values in the range [0, max(values)). This behaviour is deprecated.

    This encoding is needed for feeding categorical data to many scikit-learn
    estimators, notably linear models and SVMs with the standard kernels.

    Note: a one-hot encoding of y labels should use a LabelBinarizer
    instead.

    Read more in the :ref:`User Guide <preprocessing_categorical_features>`.

    Parameters
    ----------
    categories : 'auto' or a list of lists/arrays of values, default='auto'.
        Categories (unique values) per feature:

        - 'auto' : Determine categories automatically from the training data.
        - list : ``categories[i]`` holds the categories expected in the ith
          column. The passed categories should not mix strings and numeric
          values within a single feature, and should be sorted in case of
          numeric values.

        The used categories can be found in the ``categories_`` attribute.

    drop : 'first' or a list/array of shape (n_features,), default=None.
        Specifies a methodology to use to drop one of the categories per
        feature. This is useful in situations where perfectly collinear
        features cause problems, such as when feeding the resulting data
        into a neural network or an unregularized regression.

        - None : retain all features (the default).
        - 'first' : drop the first category in each feature. If only one
          category is present, the feature will be dropped entirely.
        - array : ``drop[i]`` is the category in feature ``X[:, i]`` that
          should be dropped.

    sparse : boolean, default=True
        Will return sparse matrix if set True else will return an array.

    dtype : number type, default=np.float
        Desired dtype of output.

    handle_unknown : 'error' or 'ignore', default='error'.
        Whether to raise an error or ignore if an unknown categorical feature
        is present during transform (default is to raise). When this parameter
        is set to 'ignore' and an unknown category is encountered during
        transform, the resulting one-hot encoded columns for this feature
        will be all zeros. In the inverse transform, an unknown category
        will be denoted as None.

    n_values : 'auto', int or array of ints, default='auto'
        Number of values per feature.

        - 'auto' : determine value range from training data.
        - int : number of categorical values per feature.
                Each feature value should be in ``range(n_values)``
        - array : ``n_values[i]`` is the number of categorical values in
                  ``X[:, i]``. Each feature value should be
                  in ``range(n_values[i])``

        .. deprecated:: 0.20
            The `n_values` keyword was deprecated in version 0.20 and will
            be removed in 0.22. Use `categories` instead.

    categorical_features : 'all' or array of indices or mask, default='all'
        Specify what features are treated as categorical.

        - 'all': All features are treated as categorical.
        - array of indices: Array of categorical feature indices.
        - mask: Array of length n_features and with dtype=bool.

        Non-categorical features are always stacked to the right of the matrix.

        .. deprecated:: 0.20
            The `categorical_features` keyword was deprecated in version
            0.20 and will be removed in 0.22.
            You can use the ``ColumnTransformer`` instead.

    Attributes
    ----------
    categories_ : list of arrays
        The categories of each feature determined during fitting
        (in order of the features in X and corresponding with the output
        of ``transform``). This includes the category specified in ``drop``
        (if any).

    drop_idx_ : array of shape (n_features,)
        ``drop_idx_[i]`` is the index in ``categories_[i]`` of the category to
        be dropped for each feature. None if all the transformed features will
        be retained.

    active_features_ : array
        Indices for active features, meaning values that actually occur
        in the training set. Only available when n_values is ``'auto'``.

        .. deprecated:: 0.20
            The ``active_features_`` attribute was deprecated in version
            0.20 and will be removed in 0.22.

    feature_indices_ : array of shape (n_features,)
        Indices to feature ranges.
        Feature ``i`` in the original data is mapped to features
        from ``feature_indices_[i]`` to ``feature_indices_[i+1]``
        (and then potentially masked by ``active_features_`` afterwards)

        .. deprecated:: 0.20
            The ``feature_indices_`` attribute was deprecated in version
            0.20 and will be removed in 0.22.

    n_values_ : array of shape (n_features,)
        Maximum number of values per feature.

        .. deprecated:: 0.20
            The ``n_values_`` attribute was deprecated in version
            0.20 and will be removed in 0.22.

    Examples
    --------
    Given a dataset with two features, we let the encoder find the unique
    values per feature and transform the data to a binary one-hot encoding.

    >>> from sklearn.preprocessing import OneHotEncoder
    >>> enc = OneHotEncoder(handle_unknown='ignore')
    >>> X = [['Male', 1], ['Female', 3], ['Female', 2]]
    >>> enc.fit(X)
    ... # doctest: +ELLIPSIS
    ... # doctest: +NORMALIZE_WHITESPACE
    OneHotEncoder(categorical_features=None, categories=None, drop=None,
       dtype=<... 'numpy.float64'>, handle_unknown='ignore',
       n_values=None, sparse=True)

    >>> enc.categories_
    [array(['Female', 'Male'], dtype=object), array([1, 2, 3], dtype=object)]
    >>> enc.transform([['Female', 1], ['Male', 4]]).toarray()
    array([[1., 0., 1., 0., 0.],
           [0., 1., 0., 0., 0.]])
    >>> enc.inverse_transform([[0, 1, 1, 0, 0], [0, 0, 0, 1, 0]])
    array([['Male', 1],
           [None, 2]], dtype=object)
    >>> enc.get_feature_names()
    array(['x0_Female', 'x0_Male', 'x1_1', 'x1_2', 'x1_3'], dtype=object)
    >>> drop_enc = OneHotEncoder(drop='first').fit(X)
    >>> drop_enc.categories_
    [array(['Female', 'Male'], dtype=object), array([1, 2, 3], dtype=object)]
    >>> drop_enc.transform([['Female', 1], ['Male', 2]]).toarray()
    array([[0., 0., 0.],
           [1., 1., 0.]])

    See also
    --------
    sklearn.preprocessing.OrdinalEncoder : performs an ordinal (integer)
      encoding of the categorical features.
    sklearn.feature_extraction.DictVectorizer : performs a one-hot encoding of
      dictionary items (also handles string-valued features).
    sklearn.feature_extraction.FeatureHasher : performs an approximate one-hot
      encoding of dictionary items or strings.
    sklearn.preprocessing.LabelBinarizer : binarizes labels in a one-vs-all
      fashion.
    sklearn.preprocessing.MultiLabelBinarizer : transforms between iterable of
      iterables and a multilabel format, e.g. a (samples x classes) binary
      matrix indicating the presence of a class label.
    """

    def __init__(self, n_values=None, categorical_features=None,
                 categories=None, drop=None, sparse=True, dtype=np.float64,
                 handle_unknown='error'):
        self.categories = categories
        self.sparse = sparse
        self.dtype = dtype
        self.handle_unknown = handle_unknown
        self.n_values = n_values
        self.categorical_features = categorical_features
        self.drop = drop

    # Deprecated attributes

    @deprecated("The ``active_features_`` attribute was deprecated in version "
                "0.20 and will be removed 0.22.")
    @property
    def active_features_(self):
        check_is_fitted(self, 'categories_')
        return self._active_features_

    @deprecated("The ``feature_indices_`` attribute was deprecated in version "
                "0.20 and will be removed 0.22.")
    @property
    def feature_indices_(self):
        check_is_fitted(self, 'categories_')
        return self._feature_indices_

    @deprecated("The ``n_values_`` attribute was deprecated in version "
                "0.20 and will be removed 0.22.")
    @property
    def n_values_(self):
        check_is_fitted(self, 'categories_')
        return self._n_values_

    def _handle_deprecations(self, X):
        # internal version of the attributes to handle deprecations
        self._n_values = self.n_values
        self._categories = getattr(self, '_categories', None)
        self._categorical_features = getattr(self, '_categorical_features',
                                             None)

        # user manually set the categories or second fit -> never legacy mode
        if self.categories is not None or self._categories is not None:
            self._legacy_mode = False
            if self.categories is not None:
                self._categories = self.categories

        # categories not set -> infer if we need legacy mode or not
        elif self.n_values is not None and self.n_values != 'auto':
            msg = (
                "Passing 'n_values' is deprecated in version 0.20 and will be "
                "removed in 0.22. You can use the 'categories' keyword "
                "instead. 'n_values=n' corresponds to "
                "'categories=[range(n)] * n_features'."
            )
            warnings.warn(msg, DeprecationWarning)
            self._legacy_mode = True

        else:  # n_values = 'auto'
            # n_values can also be None (default to catch usage), so set
            # _n_values to 'auto' explicitly
            self._n_values = 'auto'
            if self.handle_unknown == 'ignore':
                # no change in behaviour, no need to raise deprecation warning
                self._legacy_mode = False
                self._categories = 'auto'
                if self.n_values == 'auto':
                    # user manually specified this
                    msg = (
                        "Passing 'n_values' is deprecated in version 0.20 and "
                        "will be removed in 0.22. n_values='auto' can be "
                        "replaced with categories='auto'."
                    )
                    warnings.warn(msg, DeprecationWarning)
            else:
                # check if we have integer or categorical input
                try:
                    check_array(X, dtype=np.int)
                except ValueError:
                    self._legacy_mode = False
                    self._categories = 'auto'
                else:
                    if self.drop is None:
                        msg = (
                            "The handling of integer data will change in "
                            "version 0.22. Currently, the categories are "
                            "determined based on the range "
                            "[0, max(values)], while in the future they "
                            "will be determined based on the unique "
                            "values.\nIf you want the future behaviour "
                            "and silence this warning, you can specify "
                            "\"categories='auto'\".\n"
                            "In case you used a LabelEncoder before this "
                            "OneHotEncoder to convert the categories to "
                            "integers, then you can now use the "
                            "OneHotEncoder directly."
                        )
                        warnings.warn(msg, FutureWarning)
                        self._legacy_mode = True
                    else:
                        msg = (
                            "The handling of integer data will change in "
                            "version 0.22. Currently, the categories are "
                            "determined based on the range "
                            "[0, max(values)], while in the future they "
                            "will be determined based on the unique "
                            "values.\n The old behavior is not compatible "
                            "with the `drop` parameter. Instead, you "
                            "must manually specify \"categories='auto'\" "
                            "if you wish to use the `drop` parameter on "
                            "an array of entirely integer data. This will "
                            "enable the future behavior."
                        )
                        raise ValueError(msg)

        # if user specified categorical_features -> always use legacy mode
        if self.categorical_features is not None:
            if (isinstance(self.categorical_features, str)
                    and self.categorical_features == 'all'):
                warnings.warn(
                    "The 'categorical_features' keyword is deprecated in "
                    "version 0.20 and will be removed in 0.22. The passed "
                    "value of 'all' is the default and can simply be removed.",
                    DeprecationWarning)
            else:
                if self.categories is not None:
                    raise ValueError(
                        "The 'categorical_features' keyword is deprecated, "
                        "and cannot be used together with specifying "
                        "'categories'.")
                warnings.warn(
                    "The 'categorical_features' keyword is deprecated in "
                    "version 0.20 and will be removed in 0.22. You can "
                    "use the ColumnTransformer instead.", DeprecationWarning)
                # Set categories_ to empty list if no categorical columns exist
                n_features = X.shape[1]
                sel = np.zeros(n_features, dtype=bool)
                sel[np.asarray(self.categorical_features)] = True
                if sum(sel) == 0:
                    self.categories_ = []
                self._legacy_mode = True
            self._categorical_features = self.categorical_features
        else:
            self._categorical_features = 'all'

        # Prevents new drop functionality from being used in legacy mode
        if self._legacy_mode and self.drop is not None:
            raise ValueError(
                "The `categorical_features` and `n_values` keywords "
                "are deprecated, and cannot be used together "
                "with 'drop'.")

    def fit(self, X, y=None):
        """Fit OneHotEncoder to X.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to determine the categories of each feature.

        Returns
        -------
        self
        """

        self._validate_keywords()

        self._handle_deprecations(X)

        if self._legacy_mode:
            _transform_selected(X, self._legacy_fit_transform, self.dtype,
                                self._categorical_features,
                                copy=True)
            return self
        else:
            self._fit(X, handle_unknown=self.handle_unknown)
            self.drop_idx_ = self._compute_drop_idx()
            return self

    def _compute_drop_idx(self):
        if self.drop is None:
            return None
        elif (isinstance(self.drop, str) and self.drop == 'first'):
            return np.zeros(len(self.categories_), dtype=np.int_)
        elif not isinstance(self.drop, str):
            try:
                self.drop = np.asarray(self.drop, dtype=object)
                droplen = len(self.drop)
            except (ValueError, TypeError):
                msg = ("Wrong input for parameter `drop`. Expected "
                       "'first', None or array of objects, got {}")
                raise ValueError(msg.format(type(self.drop)))
            if droplen != len(self.categories_):
                msg = ("`drop` should have length equal to the number "
                       "of features ({}), got {}")
                raise ValueError(msg.format(len(self.categories_),
                                            len(self.drop)))
            missing_drops = [(i, val) for i, val in enumerate(self.drop)
                             if val not in self.categories_[i]]
            if any(missing_drops):
                msg = ("The following categories were supposed to be "
                       "dropped, but were not found in the training "
                       "data.\n{}".format(
                           "\n".join(
                                ["Category: {}, Feature: {}".format(c, v)
                                    for c, v in missing_drops])))
                raise ValueError(msg)
            return np.array([np.where(cat_list == val)[0][0]
                             for (val, cat_list) in
                             zip(self.drop, self.categories_)], dtype=np.int_)
        else:
            msg = ("Wrong input for parameter `drop`. Expected "
                   "'first', None or array of objects, got {}")
            raise ValueError(msg.format(type(self.drop)))

    def _validate_keywords(self):
        if self.handle_unknown not in ('error', 'ignore'):
            msg = ("handle_unknown should be either 'error' or 'ignore', "
                   "got {0}.".format(self.handle_unknown))
            raise ValueError(msg)
        # If we have both dropped columns and ignored unknown
        # values, there will be ambiguous cells. This creates difficulties
        # in interpreting the model.
        if self.drop is not None and self.handle_unknown != 'error':
            raise ValueError(
                "`handle_unknown` must be 'error' when the drop parameter is "
                "specified, as both would create categories that are all "
                "zero.")

    def _legacy_fit_transform(self, X):
        """Assumes X contains only categorical features."""
        dtype = getattr(X, 'dtype', None)
        X = check_array(X, dtype=np.int)
        if np.any(X < 0):
            raise ValueError("OneHotEncoder in legacy mode cannot handle "
                             "categories encoded as negative integers. "
                             "Please set categories='auto' explicitly to "
                             "be able to use arbitrary integer values as "
                             "category identifiers.")
        n_samples, n_features = X.shape
        if (isinstance(self._n_values, str) and
                self._n_values == 'auto'):
            n_values = np.max(X, axis=0) + 1
        elif isinstance(self._n_values, numbers.Integral):
            if (np.max(X, axis=0) >= self._n_values).any():
                raise ValueError("Feature out of bounds for n_values=%d"
                                 % self._n_values)
            n_values = np.empty(n_features, dtype=np.int)
            n_values.fill(self._n_values)
        else:
            try:
                n_values = np.asarray(self._n_values, dtype=int)
            except (ValueError, TypeError):
                raise TypeError("Wrong type for parameter `n_values`. Expected"
                                " 'auto', int or array of ints, got %r"
                                % type(self._n_values))
            if n_values.ndim < 1 or n_values.shape[0] != X.shape[1]:
                raise ValueError("Shape mismatch: if n_values is an array,"
                                 " it has to be of shape (n_features,).")

        self._n_values_ = n_values
        self.categories_ = [np.arange(n_val - 1, dtype=dtype)
                            for n_val in n_values]
        n_values = np.hstack([[0], n_values])
        indices = np.cumsum(n_values)
        self._feature_indices_ = indices

        column_indices = (X + indices[:-1]).ravel()
        row_indices = np.repeat(np.arange(n_samples, dtype=np.int32),
                                n_features)
        data = np.ones(n_samples * n_features)
        out = sparse.coo_matrix((data, (row_indices, column_indices)),
                                shape=(n_samples, indices[-1]),
                                dtype=self.dtype).tocsr()

        if (isinstance(self._n_values, str) and
                self._n_values == 'auto'):
            mask = np.array(out.sum(axis=0)).ravel() != 0
            active_features = np.where(mask)[0]
            out = out[:, active_features]
            self._active_features_ = active_features

            self.categories_ = [
                np.unique(X[:, i]).astype(dtype) if dtype
                else np.unique(X[:, i]) for i in range(n_features)]

        return out if self.sparse else out.toarray()

    def fit_transform(self, X, y=None):
        """Fit OneHotEncoder to X, then transform X.

        Equivalent to fit(X).transform(X) but more convenient.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to encode.

        Returns
        -------
        X_out : sparse matrix if sparse=True else a 2-d array
            Transformed input.
        """

        self._validate_keywords()

        self._handle_deprecations(X)

        if self._legacy_mode:
            return _transform_selected(
                X, self._legacy_fit_transform, self.dtype,
                self._categorical_features, copy=True)
        else:
            return self.fit(X).transform(X)

    def _legacy_transform(self, X):
        """Assumes X contains only categorical features."""
        X = check_array(X, dtype=np.int)
        if np.any(X < 0):
            raise ValueError("OneHotEncoder in legacy mode cannot handle "
                             "categories encoded as negative integers. "
                             "Please set categories='auto' explicitly to "
                             "be able to use arbitrary integer values as "
                             "category identifiers.")
        n_samples, n_features = X.shape

        indices = self._feature_indices_
        if n_features != indices.shape[0] - 1:
            raise ValueError("X has different shape than during fitting."
                             " Expected %d, got %d."
                             % (indices.shape[0] - 1, n_features))

        # We use only those categorical features of X that are known using fit.
        # i.e lesser than n_values_ using mask.
        # This means, if self.handle_unknown is "ignore", the row_indices and
        # col_indices corresponding to the unknown categorical feature are
        # ignored.
        mask = (X < self._n_values_).ravel()
        if np.any(~mask):
            if self.handle_unknown not in ['error', 'ignore']:
                raise ValueError("handle_unknown should be either error or "
                                 "unknown got %s" % self.handle_unknown)
            if self.handle_unknown == 'error':
                raise ValueError("unknown categorical feature present %s "
                                 "during transform." % X.ravel()[~mask])

        column_indices = (X + indices[:-1]).ravel()[mask]
        row_indices = np.repeat(np.arange(n_samples, dtype=np.int32),
                                n_features)[mask]
        data = np.ones(np.sum(mask))
        out = sparse.coo_matrix((data, (row_indices, column_indices)),
                                shape=(n_samples, indices[-1]),
                                dtype=self.dtype).tocsr()
        if (isinstance(self._n_values, str) and
                self._n_values == 'auto'):
            out = out[:, self._active_features_]

        return out if self.sparse else out.toarray()

    def _transform_new(self, X):
        """New implementation assuming categorical input"""
        # validation of X happens in _check_X called by _transform
        X_int, X_mask = self._transform(X, handle_unknown=self.handle_unknown)

        n_samples, n_features = X_int.shape

        if self.drop is not None:
            to_drop = self.drop_idx_.reshape(1, -1)

            # We remove all the dropped categories from mask, and decrement all
            # categories that occur after them to avoid an empty column.

            keep_cells = X_int != to_drop
            X_mask &= keep_cells
            X_int[X_int > to_drop] -= 1
            n_values = [len(cats) - 1 for cats in self.categories_]
        else:
            n_values = [len(cats) for cats in self.categories_]

        mask = X_mask.ravel()
        n_values = np.array([0] + n_values)
        feature_indices = np.cumsum(n_values)
        indices = (X_int + feature_indices[:-1]).ravel()[mask]
        indptr = X_mask.sum(axis=1).cumsum()
        indptr = np.insert(indptr, 0, 0)
        data = np.ones(n_samples * n_features)[mask]

        out = sparse.csr_matrix((data, indices, indptr),
                                shape=(n_samples, feature_indices[-1]),
                                dtype=self.dtype)
        if not self.sparse:
            return out.toarray()
        else:
            return out

    def transform(self, X):
        """Transform X using one-hot encoding.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to encode.

        Returns
        -------
        X_out : sparse matrix if sparse=True else a 2-d array
            Transformed input.
        """
        check_is_fitted(self, 'categories_')
        if self._legacy_mode:
            return _transform_selected(X, self._legacy_transform, self.dtype,
                                       self._categorical_features,
                                       copy=True)
        else:
            return self._transform_new(X)

    def inverse_transform(self, X):
        """Convert the back data to the original representation.

        In case unknown categories are encountered (all zeros in the
        one-hot encoding), ``None`` is used to represent this category.

        Parameters
        ----------
        X : array-like or sparse matrix, shape [n_samples, n_encoded_features]
            The transformed data.

        Returns
        -------
        X_tr : array-like, shape [n_samples, n_features]
            Inverse transformed array.

        """
        # if self._legacy_mode:
        #     raise ValueError("only supported for categorical features")

        check_is_fitted(self, 'categories_')
        X = check_array(X, accept_sparse='csr')

        n_samples, _ = X.shape
        n_features = len(self.categories_)
        if self.drop is None:
            n_transformed_features = sum(len(cats)
                                         for cats in self.categories_)
        else:
            n_transformed_features = sum(len(cats) - 1
                                         for cats in self.categories_)

        # validate shape of passed X
        msg = ("Shape of the passed X data is not correct. Expected {0} "
               "columns, got {1}.")
        if X.shape[1] != n_transformed_features:
            raise ValueError(msg.format(n_transformed_features, X.shape[1]))

        # create resulting array of appropriate dtype
        dt = np.find_common_type([cat.dtype for cat in self.categories_], [])
        X_tr = np.empty((n_samples, n_features), dtype=dt)

        j = 0
        found_unknown = {}

        for i in range(n_features):
            if self.drop is None:
                cats = self.categories_[i]
            else:
                cats = np.delete(self.categories_[i], self.drop_idx_[i])
            n_categories = len(cats)

            # Only happens if there was a column with a unique
            # category. In this case we just fill the column with this
            # unique category value.
            if n_categories == 0:
                X_tr[:, i] = self.categories_[i][self.drop_idx_[i]]
                j += n_categories
                continue
            sub = X[:, j:j + n_categories]
            # for sparse X argmax returns 2D matrix, ensure 1D array
            labels = np.asarray(_argmax(sub, axis=1)).flatten()
            X_tr[:, i] = cats[labels]
            if self.handle_unknown == 'ignore':
                unknown = np.asarray(sub.sum(axis=1) == 0).flatten()
                # ignored unknown categories: we have a row of all zero
                if unknown.any():
                    found_unknown[i] = unknown
            # drop will either be None or handle_unknown will be error. If
            # self.drop is not None, then we can safely assume that all of
            # the nulls in each column are the dropped value
            elif self.drop is not None:
                dropped = np.asarray(sub.sum(axis=1) == 0).flatten()
                if dropped.any():
                    X_tr[dropped, i] = self.categories_[i][self.drop_idx_[i]]

            j += n_categories

        # if ignored are found: potentially need to upcast result to
        # insert None values
        if found_unknown:
            if X_tr.dtype != object:
                X_tr = X_tr.astype(object)

            for idx, mask in found_unknown.items():
                X_tr[mask, idx] = None

        return X_tr

    def get_feature_names(self, input_features=None):
        """Return feature names for output features.

        Parameters
        ----------
        input_features : list of string, length n_features, optional
            String names for input features if available. By default,
            "x0", "x1", ... "xn_features" is used.

        Returns
        -------
        output_feature_names : array of string, length n_output_features

        """
        check_is_fitted(self, 'categories_')
        cats = self.categories_
        if input_features is None:
            input_features = ['x%d' % i for i in range(len(cats))]
        elif len(input_features) != len(self.categories_):
            raise ValueError(
                "input_features should have length equal to number of "
                "features ({}), got {}".format(len(self.categories_),
                                               len(input_features)))

        feature_names = []
        for i in range(len(cats)):
            names = [
                input_features[i] + '_' + str(t) for t in cats[i]]
            if self.drop is not None:
                names.pop(self.drop_idx_[i])
            feature_names.extend(names)

        return np.array(feature_names, dtype=object)


class OrdinalEncoder(_BaseEncoder):
    """Encode categorical features as an integer array.

    The input to this transformer should be an array-like of integers or
    strings, denoting the values taken on by categorical (discrete) features.
    The features are converted to ordinal integers. This results in
    a single column of integers (0 to n_categories - 1) per feature.

    Read more in the :ref:`User Guide <preprocessing_categorical_features>`.

    Parameters
    ----------
    categories : 'auto' or a list of lists/arrays of values.
        Categories (unique values) per feature:

        - 'auto' : Determine categories automatically from the training data.
        - list : ``categories[i]`` holds the categories expected in the ith
          column. The passed categories should not mix strings and numeric
          values, and should be sorted in case of numeric values.

        The used categories can be found in the ``categories_`` attribute.

    dtype : number type, default np.float64
        Desired dtype of output.

    Attributes
    ----------
    categories_ : list of arrays
        The categories of each feature determined during fitting
        (in order of the features in X and corresponding with the output
        of ``transform``).

    Examples
    --------
    Given a dataset with two features, we let the encoder find the unique
    values per feature and transform the data to an ordinal encoding.

    >>> from sklearn.preprocessing import OrdinalEncoder
    >>> enc = OrdinalEncoder()
    >>> X = [['Male', 1], ['Female', 3], ['Female', 2]]
    >>> enc.fit(X)
    ... # doctest: +ELLIPSIS
    OrdinalEncoder(categories='auto', dtype=<... 'numpy.float64'>)
    >>> enc.categories_
    [array(['Female', 'Male'], dtype=object), array([1, 2, 3], dtype=object)]
    >>> enc.transform([['Female', 3], ['Male', 1]])
    array([[0., 2.],
           [1., 0.]])

    >>> enc.inverse_transform([[1, 0], [0, 1]])
    array([['Male', 1],
           ['Female', 2]], dtype=object)

    See also
    --------
    sklearn.preprocessing.OneHotEncoder : performs a one-hot encoding of
      categorical features.
    sklearn.preprocessing.LabelEncoder : encodes target labels with values
      between 0 and n_classes-1.
    """

    def __init__(self, categories='auto', dtype=np.float64):
        self.categories = categories
        self.dtype = dtype

    def fit(self, X, y=None):
        """Fit the OrdinalEncoder to X.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to determine the categories of each feature.

        Returns
        -------
        self

        """
        # base classes uses _categories to deal with deprecations in
        # OneHoteEncoder: can be removed once deprecations are removed
        self._categories = self.categories
        self._fit(X)

        return self

    def transform(self, X):
        """Transform X to ordinal codes.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to encode.

        Returns
        -------
        X_out : sparse matrix or a 2-d array
            Transformed input.

        """
        X_int, _ = self._transform(X)
        return X_int.astype(self.dtype, copy=False)

    def inverse_transform(self, X):
        """Convert the data back to the original representation.

        Parameters
        ----------
        X : array-like or sparse matrix, shape [n_samples, n_encoded_features]
            The transformed data.

        Returns
        -------
        X_tr : array-like, shape [n_samples, n_features]
            Inverse transformed array.

        """
        check_is_fitted(self, 'categories_')
        X = check_array(X, accept_sparse='csr')

        n_samples, _ = X.shape
        n_features = len(self.categories_)

        # validate shape of passed X
        msg = ("Shape of the passed X data is not correct. Expected {0} "
               "columns, got {1}.")
        if X.shape[1] != n_features:
            raise ValueError(msg.format(n_features, X.shape[1]))

        # create resulting array of appropriate dtype
        dt = np.find_common_type([cat.dtype for cat in self.categories_], [])
        X_tr = np.empty((n_samples, n_features), dtype=dt)

        for i in range(n_features):
            labels = X[:, i].astype('int64', copy=False)
            X_tr[:, i] = self.categories_[i][labels]

        return X_tr

    def _more_tags(self):
        return {'X_types': ['categorical']}
