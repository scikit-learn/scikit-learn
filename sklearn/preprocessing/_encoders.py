# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
#          Joris Van den Bossche <jorisvandenbossche@gmail.com>
# License: BSD 3 clause

import numbers
import warnings

import numpy as np
from scipy import sparse

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array, is_scalar_nan
from ..utils.validation import check_is_fitted
from ..utils.validation import _deprecate_positional_args

from ..utils._encode import _encode, _check_unknown, _unique, _get_counts


__all__ = [
    'OneHotEncoder',
    'OrdinalEncoder'
]


class _BaseEncoder(TransformerMixin, BaseEstimator):
    """
    Base class for encoders that includes the code to categorize and
    transform the input features.

    """

    def _check_X(self, X, force_all_finite=True):
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
            X_temp = check_array(X, dtype=None,
                                 force_all_finite=force_all_finite)
            if (not hasattr(X, 'dtype')
                    and np.issubdtype(X_temp.dtype, np.str_)):
                X = check_array(X, dtype=object,
                                force_all_finite=force_all_finite)
            else:
                X = X_temp
            needs_validation = False
        else:
            # pandas dataframe, do validation later column by column, in order
            # to keep the dtype information to be used in the encoder.
            needs_validation = force_all_finite

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

    def _fit(self, X, handle_unknown='error', force_all_finite=True,
             return_counts=False):
        X_list, n_samples, n_features = self._check_X(
            X, force_all_finite=force_all_finite)

        if self.categories != 'auto':
            if len(self.categories) != n_features:
                raise ValueError("Shape mismatch: if categories is an array,"
                                 " it has to be of shape (n_features,).")

        self.categories_ = []
        category_counts = []

        for i in range(n_features):
            Xi = X_list[i]

            if self.categories == 'auto':
                result = _unique(Xi, return_counts=return_counts)
                if return_counts:
                    cats, counts = result
                    category_counts.append(counts)
                else:
                    cats = result
            else:
                cats = np.array(self.categories[i], dtype=Xi.dtype)
                if Xi.dtype.kind not in 'OU':
                    sorted_cats = np.sort(cats)
                    error_msg = ("Unsorted categories are not "
                                 "supported for numerical categories")
                    # if there are nans, nan should be the last element
                    stop_idx = -1 if np.isnan(sorted_cats[-1]) else None
                    if (np.any(sorted_cats[:stop_idx] != cats[:stop_idx]) or
                        (np.isnan(sorted_cats[-1]) and
                         not np.isnan(sorted_cats[-1]))):
                        raise ValueError(error_msg)

                if handle_unknown == 'error':
                    diff = _check_unknown(Xi, cats)
                    if diff:
                        msg = ("Found unknown categories {0} in column {1}"
                               " during fit".format(diff, i))
                        raise ValueError(msg)
                if return_counts:
                    category_counts.append(_get_counts(Xi, cats))

            self.categories_.append(cats)

        output = {'n_samples': n_samples}
        if return_counts:
            output['category_counts'] = category_counts
        return output

    def _transform(self, X, handle_unknown='error', force_all_finite=True):
        X_list, n_samples, n_features = self._check_X(
            X, force_all_finite=force_all_finite)

        X_int = np.zeros((n_samples, n_features), dtype=int)
        X_mask = np.ones((n_samples, n_features), dtype=bool)

        if n_features != len(self.categories_):
            raise ValueError(
                "The number of features in X is different to the number of "
                "features of the fitted data. The fitted data had {} features "
                "and the X has {} features."
                .format(len(self.categories_,), n_features)
            )

        for i in range(n_features):
            Xi = X_list[i]
            diff, valid_mask = _check_unknown(Xi, self.categories_[i],
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
            # We use check_unknown=False, since _check_unknown was
            # already called above.
            X_int[:, i] = _encode(Xi, uniques=self.categories_[i],
                                  check_unknown=False)

        return X_int, X_mask

    def _more_tags(self):
        return {'X_types': ['categorical']}


class OneHotEncoder(_BaseEncoder):
    """
    Encode categorical features as a one-hot numeric array.

    The input to this transformer should be an array-like of integers or
    strings, denoting the values taken on by categorical (discrete) features.
    The features are encoded using a one-hot (aka 'one-of-K' or 'dummy')
    encoding scheme. This creates a binary column for each category and
    returns a sparse matrix or dense array (depending on the ``sparse``
    parameter)

    By default, the encoder derives the categories based on the unique values
    in each feature. Alternatively, you can also specify the `categories`
    manually.

    This encoding is needed for feeding categorical data to many scikit-learn
    estimators, notably linear models and SVMs with the standard kernels.

    Note: a one-hot encoding of y labels should use a LabelBinarizer
    instead.

    Read more in the :ref:`User Guide <preprocessing_categorical_features>`.

    .. versionchanged:: 0.20

    Parameters
    ----------
    categories : 'auto' or a list of array-like, default='auto'
        Categories (unique values) per feature:

        - 'auto' : Determine categories automatically from the training data.
        - list : ``categories[i]`` holds the categories expected in the ith
          column. The passed categories should not mix strings and numeric
          values within a single feature, and should be sorted in case of
          numeric values.

        The used categories can be found in the ``categories_`` attribute.

        .. versionadded:: 0.20

    drop : {'first', 'if_binary'} or a array-like of shape (n_features,), \
            default=None
        Specifies a methodology to use to drop one of the categories per
        feature. This is useful in situations where perfectly collinear
        features cause problems, such as when feeding the resulting data
        into a neural network or an unregularized regression. Drop is not
        support when `min_frequency` or `max_categories` is set to combine
        infrequent categories.

        However, dropping one category breaks the symmetry of the original
        representation and can therefore induce a bias in downstream models,
        for instance for penalized linear classification or regression models.

        - None : retain all features (the default).
        - 'first' : drop the first category in each feature. If only one
          category is present, the feature will be dropped entirely.
        - 'if_binary' : drop the first category in each feature with two
          categories. Features with 1 or more than 2 categories are
          left intact.
        - array : ``drop[i]`` is the category in feature ``X[:, i]`` that
          should be dropped.

        .. versionchanged:: 0.23
           Added option 'if_binary'.

    sparse : bool, default=True
        Will return sparse matrix if set True else will return an array.

    dtype : number type, default=float
        Desired dtype of output.

    handle_unknown : {'error', 'ignore', 'auto'}, default='error'
        Whether to raise an error or ignore if an unknown categorical feature
        is present during transform (default is to raise). When this parameter
        is set to 'auto' and an unknown category is encountered during
        transform, the resulting one-hot encoded columns for this feature
        will be all zeros. In the inverse transform, an unknown category
        will be denoted as None. Read more in the
        :ref:`User Guide <one_hot_encoder_infrequent_categories>`

        .. versionadded:: 0.24
            `'auto'` was added to automatically handle unknown categories
            and infrequent categories.

        .. deprecated:: 0.24
            `'ignore'` is deprecated in favor of `'auto'`. This option will be
            removed in 0.26.

    min_frequency : int or float, default=1
        Specifies the minimum frequency for a category not to be considered
        infrequent.

            1. If int, categories with a smaller cardinality will be considered
            infrequent.

            2. If float, categories with a smaller cardinality than
            `min_frequency * n_samples`  will be considered infrequent.

        .. versionadded:: 0.24

    max_categories : int, default=None
        Specifies an upper limit to the number of output features for each
        input feature when considering infrequent categories. Note that
        `max_categories` includes the category representing the infrequent
        categories along with the frequent categories. If `None`, there is no
        limit to the number of output features.

        .. versionadded:: 0.24

    Attributes
    ----------
    categories_ : list of arrays
        The categories of each feature determined during fitting
        (in order of the features in X and corresponding with the output
        of ``transform``). This includes the category specified in ``drop``
        (if any).

    drop_idx_ : array of shape (n_features,)
        - ``drop_idx_[i]`` is the index in ``categories_[i]`` of the category
          to be dropped for each feature.
        - ``drop_idx_[i] = None`` if no category is to be dropped from the
          feature with index ``i``, e.g. when `drop='if_binary'` and the
          feature isn't binary.
        - ``drop_idx_ = None`` if all the transformed features will be
          retained.

        .. versionchanged:: 0.23
           Added the possibility to contain `None` values.

    infrequent_indices_ : list of shape (n_features,)
        Defined only when `min_frequency` or `max_categories` is set to a
        non-default value. `infrequent_indices_[i]` is an array of indices
        mapping from `categories_[i]` to the infrequent categories.
        `infrequent_indices_[i]` is None if the ith input feature has no
        infrequent categories.

    See Also
    --------
    OrdinalEncoder : Performs an ordinal (integer)
      encoding of the categorical features.
    sklearn.feature_extraction.DictVectorizer : Performs a one-hot encoding of
      dictionary items (also handles string-valued features).
    sklearn.feature_extraction.FeatureHasher : Performs an approximate one-hot
      encoding of dictionary items or strings.
    LabelBinarizer : Binarizes labels in a one-vs-all
      fashion.
    MultiLabelBinarizer : Transforms between iterable of
      iterables and a multilabel format, e.g. a (samples x classes) binary
      matrix indicating the presence of a class label.

    Examples
    --------
    Given a dataset with two features, we let the encoder find the unique
    values per feature and transform the data to a binary one-hot encoding.

    >>> from sklearn.preprocessing import OneHotEncoder

    One can discard categories not seen during `fit`:

    >>> enc = OneHotEncoder(handle_unknown='auto')
    >>> X = [['Male', 1], ['Female', 3], ['Female', 2]]
    >>> enc.fit(X)
    OneHotEncoder(handle_unknown='auto')
    >>> enc.categories_
    [array(['Female', 'Male'], dtype=object), array([1, 2, 3], dtype=object)]
    >>> enc.transform([['Female', 1], ['Male', 4]]).toarray()
    array([[1., 0., 1., 0., 0.],
           [0., 1., 0., 0., 0.]])
    >>> enc.inverse_transform([[0, 1, 1, 0, 0], [0, 0, 0, 1, 0]])
    array([['Male', 1],
           [None, 2]], dtype=object)
    >>> enc.get_feature_names(['gender', 'group'])
    array(['gender_Female', 'gender_Male', 'group_1', 'group_2', 'group_3'],
      dtype=object)

    One can always drop the first column for each feature:

    >>> drop_enc = OneHotEncoder(drop='first').fit(X)
    >>> drop_enc.categories_
    [array(['Female', 'Male'], dtype=object), array([1, 2, 3], dtype=object)]
    >>> drop_enc.transform([['Female', 1], ['Male', 2]]).toarray()
    array([[0., 0., 0.],
           [1., 1., 0.]])

    Or drop a column for feature only having 2 categories:

    >>> drop_binary_enc = OneHotEncoder(drop='if_binary').fit(X)
    >>> drop_binary_enc.transform([['Female', 1], ['Male', 2]]).toarray()
    array([[0., 1., 0., 0.],
           [1., 0., 1., 0.]])
    """

    @_deprecate_positional_args
    def __init__(self, *, categories='auto', drop=None, sparse=True,
                 dtype=np.float64, handle_unknown='error',
                 min_frequency=1, max_categories=None):
        self.categories = categories
        self.sparse = sparse
        self.dtype = dtype
        self.handle_unknown = handle_unknown
        self.drop = drop
        self.min_frequency = min_frequency
        self.max_categories = max_categories

    def _validate_keywords(self):

        if self.handle_unknown not in ('error', 'ignore', 'auto'):
            msg = (f"handle_unknown should be either 'error', 'ignore', 'auto'"
                   f"got {self.handle_unknown}.")
            raise ValueError(msg)
        # If we have both dropped columns and ignored unknown
        # values, there will be ambiguous cells. This creates difficulties
        # in interpreting the model.
        if self.drop is not None and self.handle_unknown != 'error':
            raise ValueError(
                "`handle_unknown` must be 'error' when the drop parameter is "
                "specified, as both would create categories that are all "
                "zero.")

        if self.max_categories is not None and self.max_categories <= 1:
            raise ValueError("max_categories must be greater than 1")

        if isinstance(self.min_frequency, numbers.Integral):
            if not self.min_frequency >= 1:
                raise ValueError("min_frequency must be an integer at least "
                                 "1 or a float in (0.0, 1.0); got the "
                                 f"integer {self.min_frequency}")
        else:  # float
            if not 0.0 < self.min_frequency < 1.0:
                raise ValueError("min_frequency must be an integer at least "
                                 "1 or a float in (0.0, 1.0); got the "
                                 f"float {self.min_frequency}")

        self._infrequent_enabled = (
            (self.max_categories is not None and self.max_categories > 1)
            or
            (isinstance(self.min_frequency, numbers.Integral) and
             self.min_frequency > 1)
            or
            (isinstance(self.min_frequency, numbers.Real) and
             self.min_frequency < 1.0)
        )

        # validates infrequent category features
        if self.drop is not None and self._infrequent_enabled:
            raise ValueError("infrequent categories are not supported when "
                             "drop is specified")

        # TODO: Remove when handle_unknown='ignore' is deprecated
        if self.handle_unknown == 'ignore':
            warnings.warn("handle_unknown='ignore' is deprecated in favor "
                          "of 'auto' in version 0.24 and will be removed in "
                          "version 0.26", FutureWarning)
            if self._infrequent_enabled:
                raise ValueError("infrequent categories are only supported "
                                 "when handle_unknown is 'error' or 'auto'")

    def _compute_drop_idx(self):
        if self.drop is None:
            return None
        elif isinstance(self.drop, str):
            if self.drop == 'first':
                return np.zeros(len(self.categories_), dtype=object)
            elif self.drop == 'if_binary':
                return np.array([0 if len(cats) == 2 else None
                                for cats in self.categories_], dtype=object)
            else:
                msg = (
                    "Wrong input for parameter `drop`. Expected "
                    "'first', 'if_binary', None or array of objects, got {}"
                    )
                raise ValueError(msg.format(type(self.drop)))

        else:
            try:
                self.drop = np.asarray(self.drop, dtype=object)
                droplen = len(self.drop)
            except (ValueError, TypeError):
                msg = (
                    "Wrong input for parameter `drop`. Expected "
                    "'first', 'if_binary', None or array of objects, got {}"
                    )
                raise ValueError(msg.format(type(self.drop)))
            if droplen != len(self.categories_):
                msg = ("`drop` should have length equal to the number "
                       "of features ({}), got {}")
                raise ValueError(msg.format(len(self.categories_),
                                            len(self.drop)))
            missing_drops = []
            drop_indices = []
            for col_idx, (val, cat_list) in enumerate(zip(self.drop,
                                                          self.categories_)):
                if not is_scalar_nan(val):
                    drop_idx = np.where(cat_list == val)[0]
                    if drop_idx.size:  # found drop idx
                        drop_indices.append(drop_idx[0])
                    else:
                        missing_drops.append((col_idx, val))
                    continue

                # val is nan, find nan in categories manually
                for cat_idx, cat in enumerate(cat_list):
                    if is_scalar_nan(cat):
                        drop_indices.append(cat_idx)
                        break
                else:  # loop did not break thus drop is missing
                    missing_drops.append((col_idx, val))

            if any(missing_drops):
                msg = ("The following categories were supposed to be "
                       "dropped, but were not found in the training "
                       "data.\n{}".format(
                           "\n".join(
                                ["Category: {}, Feature: {}".format(c, v)
                                    for c, v in missing_drops])))
                raise ValueError(msg)
            return np.array(drop_indices, dtype=object)

    def _identify_infrequent(self, category_count, n_samples, col_idx):
        """Compute the infrequent indices

        Parameters
        ----------
        category_count : ndarray of shape (n_cardinality,)
            Category counts.

        n_samples : int
            Number of samples.

        col_idx : int
            Index of the current category. Only used for the error message.

        Returns
        -------
        output : ndarray of shape (n_infrequent_categories,) or None
            If there are infrequent categories, indices of infrequent
            categories. Otherwise None.
        """
        if isinstance(self.min_frequency, numbers.Integral):
            infrequent_mask = category_count < self.min_frequency
        else:  # float
            min_frequency_abs = n_samples * self.min_frequency
            infrequent_mask = category_count < min_frequency_abs

        n_current_features = category_count.size - infrequent_mask.sum() + 1
        if (self.max_categories is not None
                and self.max_categories < n_current_features):
            # stable sort to preserve original count order
            smallest_levels = np.argsort(
                category_count, kind='mergesort'
                )[:-self.max_categories + 1]
            infrequent_mask[smallest_levels] = True

        output = np.flatnonzero(infrequent_mask)
        if output.size == category_count.size:
            raise ValueError(f"All categories in column {col_idx} are "
                             "infrequent, try decreasing min_frequency")
        return output if output.size > 0 else None

    def _fit_infrequent_category_mapping(self, fit_results):
        """Fit infrequent categories.

        Defines the private attribute: `_default_to_infrequent_mappings`. For
        feature `i`, `_default_to_infrequent_mappings[i]` defines the mapping
        from the integer encoding returned by `super().transform()` into
        infrequent categories. If `_default_to_infrequent_mappings[i]` is None,
        there were no infrequent categories in the training set.

        Parameters
        ----------
        fit_results : dict
            return values from `super()._fit()`
        """
        if not self._infrequent_enabled:
            return

        n_samples = fit_results["n_samples"]
        category_counts = fit_results["category_counts"]

        self.infrequent_indices_ = [
            self._identify_infrequent(category_count, n_samples, col_idx)
            for col_idx, category_count in enumerate(category_counts)
        ]

        # compute mapping from default mapping to infrequent mapping
        default_to_infrequent_mappings = []

        for cats, infreq_idx in zip(self.categories_,
                                    self.infrequent_indices_):
            # no infrequent categories
            if infreq_idx is None:
                default_to_infrequent_mappings.append(None)
                continue

            n_cats = len(cats)
            # infrequent indices exist
            mapping = np.empty(n_cats, dtype=int)
            n_infrequent_cats = infreq_idx.size

            # infrequent categories are mapped to the last element.
            n_frequent_cats = n_cats - n_infrequent_cats
            mapping[infreq_idx] = n_frequent_cats

            frequent_indices = np.setdiff1d(np.arange(n_cats), infreq_idx)
            mapping[frequent_indices] = np.arange(n_frequent_cats)

            default_to_infrequent_mappings.append(mapping)

        self._default_to_infrequent_mappings = default_to_infrequent_mappings

    def _map_to_infrequent_categories(self, X_int, X_mask):
        """Map categories to infrequent categories. This modifies X_int
        in-place.

        Parameters
        ----------
        X_int: ndarray of shape (n_samples, n_features)
            Integer encoded categories.
        """
        if not self._infrequent_enabled:
            return

        n_features = X_int.shape[1]
        for col_idx in range(n_features):
            infrequent_idx = self.infrequent_indices_[col_idx]
            if infrequent_idx is None:
                continue

            X_int[~X_mask[:, col_idx], col_idx] = infrequent_idx[0]
            if self.handle_unknown == 'auto':
                # unknown values will be mapped to infrequent in the next for
                # loop
                X_mask[:, col_idx] = True

        for i, mapping in enumerate(self._default_to_infrequent_mappings):
            if mapping is None:
                continue
            X_int[:, i] = np.take(mapping, X_int[:, i])

    def _compute_transformed_categories(self, i):
        """Compute the transformed categories used for column `i`.

        1. Dropped columns are removed.
        2. If there are infrequent categories, the category is named
        'infrequent'. If 'infrequent' is already a category, then then new
        category is called 'infrequent_sklearn'.
        """
        cats = self.categories_[i]

        if self.drop_idx_ is not None:
            if self.drop_idx_[i] is None:
                return cats
            return np.delete(cats, self.drop_idx_[i])

        # drop is None
        if not self._infrequent_enabled:
            return cats

        # infrequent is enabled
        infreq_idx = self.infrequent_indices_[i]
        if infreq_idx is None:
            return cats

        frequent_indices = np.setdiff1d(np.arange(len(cats)), infreq_idx)

        if cats.dtype.kind in 'US' and 'infrequent' in cats:
            infrequent_cat = 'infrequent_sklearn'
        else:
            infrequent_cat = 'infrequent'
        return np.r_[cats[frequent_indices],
                     np.array([infrequent_cat], dtype=object)]

    def _get_transformed_categories(self):
        """Transformed categories."""
        return [self._compute_transformed_categories(i)
                for i in range(len(self.categories_))]

    def _get_n_transformed_features(self):
        """Number of transformed features."""
        if self.drop_idx_ is not None:
            output = []
            for drop_idx, cats in zip(self.drop_idx_, self.categories_):
                if drop_idx is None:
                    output.append(len(cats))
                else:
                    output.append(len(cats) - 1)
            return output

        # drop is None
        output = [len(cats) for cats in self.categories_]

        if not self._infrequent_enabled:
            return output

        # infrequent is enabled
        for i, infreq_idx in enumerate(self.infrequent_indices_):
            if infreq_idx is None:
                continue
            output[i] = output[i] - infreq_idx.size + 1

        return output

    def fit(self, X, y=None):
        """
        Fit OneHotEncoder to X.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to determine the categories of each feature.

        y : None
            Ignored. This parameter exists only for compatibility with
            :class:`~sklearn.pipeline.Pipeline`.

        Returns
        -------
        self
        """
        self._validate_keywords()
        fit_results = self._fit(X, handle_unknown=self.handle_unknown,
                                force_all_finite='allow-nan',
                                return_counts=self._infrequent_enabled)
        self._fit_infrequent_category_mapping(fit_results)
        self.drop_idx_ = self._compute_drop_idx()
        return self

    def fit_transform(self, X, y=None):
        """
        Fit OneHotEncoder to X, then transform X.

        Equivalent to fit(X).transform(X) but more convenient.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to encode.

        y : None
            Ignored. This parameter exists only for compatibility with
            :class:`~sklearn.pipeline.Pipeline`.

        Returns
        -------
        X_out : sparse matrix if sparse=True else a 2-d array
            Transformed input.
        """
        self._validate_keywords()
        return super().fit_transform(X, y)

    def transform(self, X):
        """
        Transform X using one-hot encoding.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to encode.

        Returns
        -------
        X_out : sparse matrix if sparse=True else a 2-d array
            Transformed input.
        """
        check_is_fitted(self)
        # validation of X happens in _check_X called by _transform
        X_int, X_mask = self._transform(X, handle_unknown=self.handle_unknown,
                                        force_all_finite='allow-nan')
        self._map_to_infrequent_categories(X_int, X_mask)

        n_samples, n_features = X_int.shape

        if self.drop_idx_ is not None:
            to_drop = self.drop_idx_.copy()
            # We remove all the dropped categories from mask, and decrement all
            # categories that occur after them to avoid an empty column.
            keep_cells = X_int != to_drop
            for i, cats in enumerate(self.categories_):
                # drop='if_binary' but feature isn't binary
                if to_drop[i] is None:
                    # set to cardinality to not drop from X_int
                    to_drop[i] = len(cats)

            to_drop = to_drop.reshape(1, -1)
            X_int[X_int > to_drop] -= 1
            X_mask &= keep_cells

        n_values = self._get_n_transformed_features()

        mask = X_mask.ravel()
        feature_indices = np.cumsum([0] + n_values)
        indices = (X_int + feature_indices[:-1]).ravel()[mask]

        indptr = np.empty(n_samples + 1, dtype=int)
        indptr[0] = 0
        np.sum(X_mask, axis=1, out=indptr[1:])
        np.cumsum(indptr[1:], out=indptr[1:])
        data = np.ones(indptr[-1])

        out = sparse.csr_matrix((data, indices, indptr),
                                shape=(n_samples, feature_indices[-1]),
                                dtype=self.dtype)
        if not self.sparse:
            return out.toarray()
        else:
            return out

    def inverse_transform(self, X):
        """
        Convert the data back to the original representation.

        In case unknown categories are encountered (all zeros in the
        one-hot encoding), ``None`` is used to represent this category.

        For a given input feature, if there is an infrequent category,
        'infrequent' will be used to represent the infrequent category. If
        'infrequent' is already a category, 'infrequent_sklearn' will be used
        instead.

        Parameters
        ----------
        X : array-like or sparse matrix, shape [n_samples, n_encoded_features]
            The transformed data.

        Returns
        -------
        X_tr : array-like, shape [n_samples, n_features]
            Inverse transformed array.
        """
        check_is_fitted(self)
        X = check_array(X, accept_sparse='csr')

        n_samples, _ = X.shape
        n_features = len(self.categories_)
        transformed_features = self._get_transformed_categories()
        n_features_out = sum(cats.shape[0] for cats in transformed_features)

        # validate shape of passed X
        msg = ("Shape of the passed X data is not correct. Expected {0} "
               "columns, got {1}.")
        if X.shape[1] != n_features_out:
            raise ValueError(msg.format(n_features_out, X.shape[1]))

        # create resulting array of appropriate dtype
        dt = np.find_common_type([cat.dtype
                                  for cat in transformed_features], [])
        X_tr = np.empty((n_samples, n_features), dtype=dt)

        j = 0
        found_unknown = {}

        if self._infrequent_enabled:
            infrequent_indices = self.infrequent_indices_
        else:
            infrequent_indices = [None] * n_features

        for i in range(n_features):
            cats = transformed_features[i]
            n_categories = cats.shape[0]

            # Only happens if there was a column with a unique
            # category. In this case we just fill the column with this
            # unique category value.
            if n_categories == 0:
                X_tr[:, i] = self.categories_[i][self.drop_idx_[i]]
                j += n_categories
                continue
            sub = X[:, j:j + n_categories]
            # for sparse X argmax returns 2D matrix, ensure 1D array
            labels = np.asarray(sub.argmax(axis=1)).flatten()
            X_tr[:, i] = cats[labels]

            if (self.handle_unknown == 'ignore' or
                (self.handle_unknown == 'auto' and
                 infrequent_indices[i] is None)):
                unknown = np.asarray(sub.sum(axis=1) == 0).flatten()
                # ignored unknown categories: we have a row of all zero
                if unknown.any():
                    found_unknown[i] = unknown
            else:
                dropped = np.asarray(sub.sum(axis=1) == 0).flatten()
                if dropped.any():
                    if self.drop_idx_ is None:
                        all_zero_samples = np.flatnonzero(dropped)
                        raise ValueError(
                            f"Samples {all_zero_samples} can not be inverted "
                            "when drop=None and handle_unknown='error' "
                            "because they contain all zeros")
                    # we can safely assume that all of the nulls in each column
                    # are the dropped value
                    X_tr[dropped, i] = self.categories_[i][
                        self.drop_idx_[i]
                    ]

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
        """
        Return feature names for output features.

        For a given input feature, if there is an infrequent category, the most
        'infrequent' will be used as a feature name. If 'infrequent' is already
        a category, 'infrequent_sklearn' will be used instead.

        Parameters
        ----------
        input_features : list of str of shape (n_features,)
            String names for input features if available. By default,
            "x0", "x1", ... "xn_features" is used.

        Returns
        -------
        output_feature_names : ndarray of shape (n_output_features,)
            Array of feature names.
        """
        check_is_fitted(self)
        cats = self._get_transformed_categories()
        if input_features is None:
            input_features = ['x%d' % i for i in range(len(cats))]
        elif len(input_features) != len(cats):
            raise ValueError(
                "input_features should have length equal to number of "
                "features ({}), got {}".format(len(cats), len(input_features)))

        feature_names = []
        for i in range(len(cats)):
            names = [
                input_features[i] + '_' + str(t) for t in cats[i]]
            feature_names.extend(names)

        return np.array(feature_names, dtype=object)


class OrdinalEncoder(_BaseEncoder):
    """
    Encode categorical features as an integer array.

    The input to this transformer should be an array-like of integers or
    strings, denoting the values taken on by categorical (discrete) features.
    The features are converted to ordinal integers. This results in
    a single column of integers (0 to n_categories - 1) per feature.

    Read more in the :ref:`User Guide <preprocessing_categorical_features>`.

    .. versionadded:: 0.20

    Parameters
    ----------
    categories : 'auto' or a list of array-like, default='auto'
        Categories (unique values) per feature:

        - 'auto' : Determine categories automatically from the training data.
        - list : ``categories[i]`` holds the categories expected in the ith
          column. The passed categories should not mix strings and numeric
          values, and should be sorted in case of numeric values.

        The used categories can be found in the ``categories_`` attribute.

    dtype : number type, default np.float64
        Desired dtype of output.

    handle_unknown : {'error', 'use_encoded_value'}, default='error'
        When set to 'error' an error will be raised in case an unknown
        categorical feature is present during transform. When set to
        'use_encoded_value', the encoded value of unknown categories will be
        set to the value given for the parameter `unknown_value`. In
        :meth:`inverse_transform`, an unknown category will be denoted as None.

        .. versionadded:: 0.24

    unknown_value : int or np.nan, default=None
        When the parameter handle_unknown is set to 'use_encoded_value', this
        parameter is required and will set the encoded value of unknown
        categories. It has to be distinct from the values used to encode any of
        the categories in `fit`. If set to np.nan, the `dtype` parameter must
        be a float dtype.

        .. versionadded:: 0.24

    Attributes
    ----------
    categories_ : list of arrays
        The categories of each feature determined during ``fit`` (in order of
        the features in X and corresponding with the output of ``transform``).
        This does not include categories that weren't seen during ``fit``.

    See Also
    --------
    OneHotEncoder : Performs a one-hot encoding of categorical features.
    LabelEncoder : Encodes target labels with values between 0 and
        ``n_classes-1``.

    Examples
    --------
    Given a dataset with two features, we let the encoder find the unique
    values per feature and transform the data to an ordinal encoding.

    >>> from sklearn.preprocessing import OrdinalEncoder
    >>> enc = OrdinalEncoder()
    >>> X = [['Male', 1], ['Female', 3], ['Female', 2]]
    >>> enc.fit(X)
    OrdinalEncoder()
    >>> enc.categories_
    [array(['Female', 'Male'], dtype=object), array([1, 2, 3], dtype=object)]
    >>> enc.transform([['Female', 3], ['Male', 1]])
    array([[0., 2.],
           [1., 0.]])

    >>> enc.inverse_transform([[1, 0], [0, 1]])
    array([['Male', 1],
           ['Female', 2]], dtype=object)
    """

    @_deprecate_positional_args
    def __init__(self, *, categories='auto', dtype=np.float64,
                 handle_unknown='error', unknown_value=None):
        self.categories = categories
        self.dtype = dtype
        self.handle_unknown = handle_unknown
        self.unknown_value = unknown_value

    def fit(self, X, y=None):
        """
        Fit the OrdinalEncoder to X.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to determine the categories of each feature.

        y : None
            Ignored. This parameter exists only for compatibility with
            :class:`~sklearn.pipeline.Pipeline`.

        Returns
        -------
        self
        """
        if self.handle_unknown == 'use_encoded_value':
            if is_scalar_nan(self.unknown_value):
                if np.dtype(self.dtype).kind != 'f':
                    raise ValueError(
                        f"When unknown_value is np.nan, the dtype "
                        "parameter should be "
                        f"a float dtype. Got {self.dtype}."
                    )
            elif not isinstance(self.unknown_value, numbers.Integral):
                raise TypeError(f"unknown_value should be an integer or "
                                f"np.nan when "
                                f"handle_unknown is 'use_encoded_value', "
                                f"got {self.unknown_value}.")
        elif self.unknown_value is not None:
            raise TypeError(f"unknown_value should only be set when "
                            f"handle_unknown is 'use_encoded_value', "
                            f"got {self.unknown_value}.")

        self._fit(X)

        if self.handle_unknown == 'use_encoded_value':
            for feature_cats in self.categories_:
                if 0 <= self.unknown_value < len(feature_cats):
                    raise ValueError(f"The used value for unknown_value "
                                     f"{self.unknown_value} is one of the "
                                     f"values already used for encoding the "
                                     f"seen categories.")

        return self

    def transform(self, X):
        """
        Transform X to ordinal codes.

        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
            The data to encode.

        Returns
        -------
        X_out : sparse matrix or a 2-d array
            Transformed input.
        """
        X_int, X_mask = self._transform(X, handle_unknown=self.handle_unknown)
        X_trans = X_int.astype(self.dtype, copy=False)

        # create separate category for unknown values
        if self.handle_unknown == 'use_encoded_value':
            X_trans[~X_mask] = self.unknown_value
        return X_trans

    def inverse_transform(self, X):
        """
        Convert the data back to the original representation.

        Parameters
        ----------
        X : array-like or sparse matrix, shape [n_samples, n_encoded_features]
            The transformed data.

        Returns
        -------
        X_tr : array-like, shape [n_samples, n_features]
            Inverse transformed array.
        """
        check_is_fitted(self)
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

        found_unknown = {}

        for i in range(n_features):
            labels = X[:, i].astype('int64', copy=False)
            if self.handle_unknown == 'use_encoded_value':
                unknown_labels = labels == self.unknown_value
                X_tr[:, i] = self.categories_[i][np.where(
                    unknown_labels, 0, labels)]
                found_unknown[i] = unknown_labels
            else:
                X_tr[:, i] = self.categories_[i][labels]

        # insert None values for unknown values
        if found_unknown:
            X_tr = X_tr.astype(object, copy=False)

            for idx, mask in found_unknown.items():
                X_tr[mask, idx] = None

        return X_tr
