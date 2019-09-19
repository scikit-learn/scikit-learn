# Authors: Andreas Mueller <amueller@ais.uni-bonn.de>
#          Joris Van den Bossche <jorisvandenbossche@gmail.com>
# License: BSD 3 clause

import numpy as np
from scipy import sparse
from itertools import count

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array
from ..utils import safe_indexing
from ..utils.fixes import _argmax
from ..utils.validation import check_is_fitted

from .label import _encode, _encode_check_unknown


__all__ = [
    'OneHotEncoder',
    'OrdinalEncoder'
]


class _EncoderUnion(TransformerMixin, BaseEstimator):
    """Base class for encoders that includes the code to encode and
    transform the input features one by one.

    The encoders passed to `_fit_list` must define:

    1. `fit(X)` where `X` is a ndarray of shape (n_samples,)
    2. `transform(X)` where `X` is a ndarray of shape (n_samples,)
    3. `inverse_transform(X)` where `X` is an encoded ndarray of sparse array
    """

    @property
    def categories_(self):
        return [encoder.categories_ for encoder in self._single_encoders]

    def _check_X(self, X):
        """Perform custom check_array:
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

        X_columns = []

        for i in range(X.shape[1]):
            Xi = safe_indexing(X, i, axis=1)
            Xi = check_array(Xi, ensure_2d=False, dtype=None,
                             force_all_finite=needs_validation)
            X_columns.append(Xi)

        return X_columns

    def _check_categories(self, n_features):
        """Check categories are consistent with n_features"""
        categories = self.categories

        if (self.categories != 'auto' and
                len(self.categories) != n_features):
            raise ValueError("Shape mismatch: if categories is an array, "
                             "it has to be of shape (n_features,).")

        if isinstance(self.categories, str) and self.categories == 'auto':
            categories = ['auto'] * n_features
        else:
            categories = self.categories

        return categories

    def _fit_list(self, X_list, single_encoders):
        """Fit single_encoders on X_list"""
        for X_col, encoder in zip(X_list, single_encoders):
            encoder.fit(X_col)

        # maps indices from original X to indices in the transformed X
        # this is used in inverse_transform
        orig_idx_to_X_trans_idx = []
        X_trans_idx = 0
        for encoder in single_encoders:
            n_features_out = encoder.n_features_out_
            begin, end = X_trans_idx, X_trans_idx + n_features_out
            orig_idx_to_X_trans_idx.append((begin, end))
            X_trans_idx += n_features_out

        self._orig_idx_to_X_trans_idx = orig_idx_to_X_trans_idx
        self._single_encoders = single_encoders

    def _transform_list(self, X_list):
        """Transform X_list with encoders"""
        n_features = len(X_list)
        if n_features != len(self.categories_):
            raise ValueError(
                "The number of features in X is different to the number of "
                "features of the fitted data. The fitted data had {} features "
                "and the X has {} features."
                .format(len(self.categories_,), n_features)
            )

        X_trans = []
        for encoder, X_col in zip(self._single_encoders, X_list):
            if encoder.n_features_out_ != 0:
                X_trans.append(encoder.transform(X_col))
        return self._hstack(X_trans)

    def _hstack(self, X_trans):
        if any(sparse.issparse(X_tran) for X_tran in X_trans):
            X_trans_stacked = sparse.hstack(X_trans).tocsr()
        else:
            X_trans_stacked = np.hstack(X_trans)
        return X_trans_stacked

    def inverse_transform(self, X):
        """Convert the data back into the original representation.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape \
            (n_samples, n_encoded_features)
            The transformed data.

        Returns
        -------
        X_tr : ndarray of shape (n_samples, n_features)
            Inverse transformed array.
        """
        check_is_fitted(self)
        X = check_array(X, accept_sparse='csr')

        n_features = sum(encoder.n_features_out_ for encoder in
                         self._single_encoders)

        # validate shape of passed X
        msg = ("Shape of the passed X data is not correct. Expected {0} "
               "columns, got {1}.")
        if X.shape[1] != n_features:
            raise ValueError(msg.format(n_features, X.shape[1]))

        X_trs = []
        for encoder, (begin, end) in zip(self._single_encoders,
                                         self._orig_idx_to_X_trans_idx):
            X_slice = safe_indexing(X, slice(begin, end), axis=1)
            X_trs.append(encoder.inverse_transform(X_slice))
        return self._hstack(X_trs)

    def _more_tags(self):
        return {'X_types': ['categorical']}


class OneHotEncoder(_EncoderUnion):
    """Encode categorical features as a one-hot numeric array.

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

    Examples
    --------
    Given a dataset with two features, we let the encoder find the unique
    values per feature and transform the data to a binary one-hot encoding.

    >>> from sklearn.preprocessing import OneHotEncoder
    >>> enc = OneHotEncoder(handle_unknown='ignore')
    >>> X = [['Male', 1], ['Female', 3], ['Female', 2]]
    >>> enc.fit(X)
    OneHotEncoder(handle_unknown='ignore')

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
    def __init__(self, categories='auto', drop=None, sparse=True,
                 dtype=np.float64, handle_unknown='error'):

        self.categories = categories
        self.sparse = sparse
        self.dtype = dtype
        self.handle_unknown = handle_unknown
        self.drop = drop

    @property
    def drop_idx_(self):
        return np.array([encoder.drop_idx_ for encoder in
                         self._single_encoders], dtype=np.int_)

    def _fit(self, X):
        """Validate keywords and fit `X` and return `X_list`."""
        self._validate_keywords()
        X_list = self._check_X(X)
        n_features = len(X_list)

        categories = self._check_categories(n_features)
        drop_kwargs = self._check_drop(n_features)

        encoders = [
            _SingleOneHotEncoder(categories=cat,
                                 dtype=self.dtype,
                                 handle_unknown=self.handle_unknown,
                                 sparse=self.sparse,
                                 feature_idx=idx,
                                 **drop_kwarg)
            for idx, cat, drop_kwarg in zip(count(), categories, drop_kwargs)]

        super()._fit_list(X_list, encoders)

        # validate encoders
        missing_drops = []
        for idx, encoder in enumerate(encoders):
            drop_idx = encoder.drop_idx_
            if isinstance(drop_idx, str) and drop_idx == 'missing':
                missing_drops.append((idx, encoder.drop_category))

        if any(missing_drops):
            msg = ("The following categories were supposed to be "
                   "dropped, but were not found in the training "
                   "data.\n{}".format("\n".join(
                            ["Category: {}, Feature: {}".format(c, v)
                                for c, v in missing_drops])))
            raise ValueError(msg)

        return X_list

    def fit(self, X, y=None):
        """Fit OneHotEncoder to X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categroies of each feature.

        Returns
        -------
        self
        """
        self._fit(X)
        return self

    def transform(self, X):
        """Transform X to encoding

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to encode.

        Returns
        -------
        X_out : {ndarray, sparse matrix} of shape \
            (n_samples, n_encoded_features)
            Transformed array. When `sparse=True`, `X_out` is a sparse array.
        """
        check_is_fitted(self)
        X_list = self._check_X(X)
        return super()._transform_list(X_list)

    def fit_transform(self, X, y=None):
        """Fit encoder to X and transform X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to encode.

        Returns
        -------
        X_out : {ndarray, sparse matrix} of shape \
            (n_samples, n_encoded_features)
            Transformed array. When `sparse=True`, `X_out` is a sparse array.
        """
        X_list = self._fit(X)
        return super()._transform_list(X_list)

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
        check_is_fitted(self)
        cats = self.categories_
        if input_features is None:
            input_features = ['x%d' % i for i in range(len(cats))]
        elif len(input_features) != len(self.categories_):
            raise ValueError(
                "input_features should have length equal to number of "
                "features ({}), got {}".format(len(self.categories_),
                                               len(input_features)))

        feature_names = []
        for input_feat, cat, encoder in zip(input_features, cats,
                                            self._single_encoders):
            names = [input_feat + '_' + str(t) for t in cat]
            if self.drop is not None:
                names.pop(encoder.drop_idx_)
            feature_names.extend(names)

        return np.array(feature_names, dtype=object)

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

    def _check_drop(self, n_features):
        if self.drop is None:
            return [{'drop_category': None}] * n_features
        if isinstance(self.drop, str) and self.drop == 'first':
            return [{'drop_first': True}] * n_features
        if not isinstance(self.drop, str):
            try:
                drops = np.asarray(self.drop, dtype=object)
                drops_len = len(drops)
            except (ValueError, TypeError):
                msg = ("Wrong input for parameter `drop`. Expected "
                       "'first', None or array of objects, got {}")
                raise ValueError(msg.format(type(self.drop)))
            if drops_len != n_features:
                msg = ("`drop` should have length equal to the number "
                       "of features ({}), got {}")
                raise ValueError(msg.format(n_features,
                                            len(self.drop)))
            return [{'drop_category': drop} for drop in drops]
        else:
            msg = ("Wrong input for parameter `drop`. Expected "
                   "'first', None or array of objects, got {}")
            raise ValueError(msg.format(type(self.drop)))


class _SingleOneHotEncoder(TransformerMixin, BaseEstimator):
    """One hot encoder for a single categorical feature.

    When calling `fit`, the attribute `drop_idx_` will be set to 'missing'
    when the drop category is not found in the dataset or given by
    categories. `drop_idx_` should be checked be the caller to make sure
    the encoder is valid."""
    def __init__(self, categories='auto', drop_category=None, drop_first=False,
                 dtype=np.float64, handle_unknown='error',
                 sparse=True, feature_idx=0):
        self.categories = categories
        self.dtype = dtype
        self.handle_unknown = handle_unknown
        self.drop_category = drop_category
        self.drop_first = drop_first
        self.sparse = sparse
        self.feature_idx = feature_idx

    def fit(self, X):
        """Fit one hot encoder for a single categorical feature.

        Parameters
        ----------
        X : ndarray of shape (n_samples,)
            Categorical feature to encode.

        Returns
        ------
        self
        """
        if isinstance(self.categories, str) and self.categories == 'auto':
            cats = _encode(X)
        else:  # categories were given
            cats = np.array(self.categories, dtype=X.dtype)
            if X.dtype != object and not np.all(np.sort(cats) == cats):
                raise ValueError("Unsorted categories are not "
                                 "supported for numerical categories")

            if self.handle_unknown == 'error':
                diff = _encode_check_unknown(X, cats)
                if diff:
                    msg = ("Found unknown categories {0} in column {1} "
                           "during fit".format(diff, self.feature_idx))
                    raise ValueError(msg)

        self.categories_ = cats

        # compute drop idx
        if self.drop_first:
            self.drop_idx_ = 0
        elif self.drop_category is None:
            self.drop_idx_ = None
        elif self.drop_category not in self.categories_:
            # This is an error state. Caller should check this value and
            # handle accordingly.
            self.drop_idx_ = 'missing'
        else:
            self.drop_idx_ = np.where(self.categories_ ==
                                      self.drop_category)[0][0]

        if self.drop_idx_ is not None:
            self.n_features_out_ = len(self.categories_) - 1
        else:
            self.n_features_out_ = len(self.categories_)
        return self

    def transform(self, X):
        """Transform a single categorical feature.

        Parameters
        ----------
        X : ndarray of shape (n_samples,)
            Categorical feature to encode.

        Returns
        -------
        X_tr : {ndarray, sparse matrix} of shape \
            (n_samples, n_encoded_features)
            Encoded feature. If `sparse=True` a sparse matrix is returned.
        """
        diff, X_mask = _encode_check_unknown(X, self.categories_,
                                             return_mask=True)
        if not np.all(X_mask):
            if self.handle_unknown == 'error':
                msg = ("Found unknown categories {0} in column {1} "
                       "during transform".format(diff, self.feature_idx))
                raise ValueError(msg)

            # cast Xi into the largest string type necessary
            # to handle different lengths of numpy strings
            if (self.categories_.dtype.kind in ('U', 'S')
                    and self.categories_.itemsize > X.itemsize):
                X = X.astype(self.categories_.dtype)
            else:
                X = X.copy()
            X[~X_mask] = self.categories_[0]

        _, X_encoded = _encode(X, self.categories_, encode=True,
                               check_unknown=False)

        if self.drop_idx_ is not None:
            keep_cells = X_encoded != self.drop_idx_
            X_mask &= keep_cells

            # adjust encoding to remove the dropped column
            X_encoded[X_encoded > self.drop_idx_] -= 1

        n_samples = X.shape[0]
        X_mask_non_zero = np.flatnonzero(X_mask)

        out = sparse.csr_matrix((n_samples, self.n_features_out_),
                                dtype=self.dtype)
        out[X_mask_non_zero, X_encoded[X_mask]] = 1

        if self.sparse:
            return out
        else:
            return out.toarray()

    def inverse_transform(self, X):
        """Inverse transform to a single categorical feature.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples,)
            The transformed data of a single feature.

        Returns
        -------
        X_tr : ndarray of shape (n_samples, 1)
            Inverse transform.
        """
        # Only happens if there was a column with a unique
        # category. In this case we just fill the column with this
        # unique category value.
        if self.n_features_out_ == 0:
            value = self.categories_[self.drop_idx_]
            return np.full((X.shape[0], 1), value)

        if self.drop_idx_ is None:
            cats = self.categories_
        else:
            cats = np.delete(self.categories_, self.drop_idx_)

        # for sparse X argmax returns 2D matrix, ensure 1D array
        labels = np.asarray(_argmax(X, axis=1)).flatten()
        X_tr = cats[labels]
        if self.handle_unknown == 'ignore':
            unknown = np.asarray(X.sum(axis=1) == 0).flatten()
            # ignored unknown categories: we have a row of all zero
            if unknown.any():
                if X_tr.dtype != object:
                    X_tr = X_tr.astype(object)
                X_tr[unknown] = None
        # drop will either be None or handle_unknown will be error. If
        # self.drop is not None, then we can safely assume that all of
        # the nulls in each column are the dropped value
        elif self.drop_idx_ is not None:
            dropped = np.asarray(X.sum(axis=1) == 0).flatten()
            if dropped.any():
                X_tr[dropped] = self.categories_[self.drop_idx_]

        return X_tr[:, None]


class OrdinalEncoder(_EncoderUnion):
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
    OrdinalEncoder()
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

    def _fit(self, X):
        """Validate keywords and fit `X` and return `X_list`."""
        X_list = self._check_X(X)
        categories = self._check_categories(len(X_list))

        encoders = [_SingleOrdinalEncoder(categories=cat,
                                          dtype=self.dtype,
                                          feature_idx=idx)
                    for idx, cat in enumerate(categories)]

        super()._fit_list(X_list, encoders)
        return X_list

    def fit(self, X, y=None):
        """Fit OneHotEncoder to X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categroies of each feature.

        Returns
        -------
        self
        """
        self._fit(X)
        return self

    def transform(self, X):
        """Transform X to encoding

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to encode.

        Returns
        -------
        X_out : ndarray of shape (n_samples, n_encoded_features)
            Transformed array.
        """
        check_is_fitted(self)
        X_list = self._check_X(X)
        return super()._transform_list(X_list)

    def fit_transform(self, X, y=None):
        """Fit encoder to X and transform X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to encode.

        Returns
        -------
        X_out : ndarray of shape (n_samples, n_encoded_features)
            Transformed array.
        """
        X_list = self._fit(X)
        return super()._transform_list(X_list)


class _SingleOrdinalEncoder(TransformerMixin, BaseEstimator):
    """Ordinal Encoder for a single categorical feature."""
    def __init__(self, categories='auto', dtype=np.float64, feature_idx=0):
        self.categories = categories
        self.dtype = dtype
        self.feature_idx = feature_idx

    def fit(self, X, y=None):
        """Fit ordinal encoder on a single categorical feature.

        Parameters
        ----------
        X : ndarray of shape (n_samples,)
            Categorical feature to encode.

        Returns
        -------
        self
        """
        if isinstance(self.categories, str) and self.categories == 'auto':
            cats = _encode(X)
        else:  # categories were given
            cats = np.array(self.categories, dtype=X.dtype)
            if X.dtype != object and not np.all(np.sort(cats) == cats):
                raise ValueError("Unsorted categories are not "
                                 "supported for numerical categories")
            diff = _encode_check_unknown(X, cats)
            if diff:
                msg = ("Found unknown categories {0} in column {1} "
                       "during fit".format(diff, self.feature_idx))
                raise ValueError(msg)

        self.categories_ = cats
        self.n_features_out_ = 1
        return self

    def transform(self, X):
        """Transform on an single categorical feature.

        Parameters
        ----------
        X : ndarray of shape (n_samples,)
            Categorical feature to encode.

        Returns
        -------
        X_tr : ndarray of shape (n_samples, 1)
            Encoded categorical feature.
        """
        diff, valid_mask = _encode_check_unknown(X, self.categories_,
                                                 return_mask=True)
        if not np.all(valid_mask):
            if self.handle_unknown == 'error':
                msg = ("Found unknown categories {0} in column {1} "
                       "during transform".format(diff, self.feature_idx))
                raise ValueError(msg)

            # cast Xi into the largest string type necessary
            # to handle different lengths of numpy strings
            if (self.categories_.dtype.kind in ('U', 'S')
                    and self.categories_.itemsize > X.itemsize):
                X = X.astype(self.categories_.dtype)
            else:
                X = X.copy()
            X[~valid_mask] = self.categories_[0]

        _, encoded = _encode(X, self.categories_, encode=True,
                             check_unknown=False)
        return encoded[:, None].astype(self.dtype, copy=False)

    def inverse_transform(self, X):
        """Convert the data back into the original representation.

        Parameters
        ----------
        X : ndarray of shape (n_samples,)
            Transformed data.

        Returns
        -------
        X_tr : ndarry of shape (n_samples, 1)
            Inverse transformed array.
        """
        labels = X.astype('int64', copy=False)
        return self.categories_[labels]
