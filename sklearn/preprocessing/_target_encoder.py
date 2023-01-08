import numpy as np

from numbers import Real

from ._encoders import _BaseEncoder
from ..base import OneToOneFeatureMixin
from ._target_encoder_fast import _fit_encoding_fast
from ..utils.validation import _check_y, check_consistent_length
from ..utils._param_validation import Interval, StrOptions


class TargetRegressorEncoder(OneToOneFeatureMixin, _BaseEncoder):
    """Target Encoder for Regression Targets.

    Each category is encoded based on its effect on the target variable. The encoding
    scheme mixes the global target mean with the target mean conditioned on the
    category.

    Read more in the :ref:`User Guide <target_regressor_encoder>`.

    .. note::
        :class:`TargetRegressorEncoder` uses a cross validation scheme in
        :meth:`fit_transform` to prevent leaking the target during training. In
        :meth:`fit_transform`, categorical encodings are obtained from one split and
        used to encoding the other split. Afterwards, a final categorical encoding is
        obtained from all the training data, which is used to encode data
        during :meth:`transform`. This means that `fit().transform()` does not equal
        `fit_transform()`.

    .. versionadded:: 1.3

    Parameters
    ----------
    categories : 'auto' or a list of array-like, default='auto'
        Categories (unique values) per feature:

        - 'auto' : Determine categories automatically from the training data.
        - list : ``categories[i]`` holds the categories expected in the ith column. The
          passed categories should not mix strings and numeric values within a single
          feature, and should be sorted in case of numeric values.

        The used categories can be found in the ``categories_`` attribute.

    smooth : float, default=30.0
        The amount of mixing the categorical encoding with the global target mean. A
        larger `smooth` value will put more weight on the global target mean.

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy used in
        :meth:`fit_transform`. Possible inputs for cv are:

        - None, to use the default 5-fold cross validation,
        - integer, to specify the number of folds in a `(Stratified)KFold`,
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        Refer :ref:`User Guide <cross_validation>` for the various cross-validation
        strategies that can be used here.

    Attributes
    ----------
    encodings_ : list of shape (n_features,) of ndarray
        For feature `i`, `encodings_[i]` is the encoding matching the
        categories listed in `categories_[i]`.

    categories_ : list of shape (n_features,) of ndarray
        The categories of each feature determined during fitting
        (in order of the features in X and corresponding with the output
        of :meth:`transform`).

    encoding_mean_ : float
        The overall mean of the target.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    OrdinalEncoder : Performs an ordinal (integer) encoding of the categorical features.
    OneHotEncoder : Performs a one-hot encoding of categorical features.
    """

    _parameter_constraints: dict = {
        "categories": [StrOptions({"auto"}), list],
        "smooth": [Interval(Real, 1, None, closed="left")],
        "cv": ["cv_object"],
    }

    def __init__(self, categories="auto", smooth=30.0, cv=5):
        self.categories = categories
        self.smooth = smooth
        self.cv = cv

    def fit(self, X, y):
        """Fit the :class:`TargetRegressorEncoder` to X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categories of each feature.

        y : array-like of shape (n_samples,)
            The target data used to encode the categories.

        Returns
        -------
        self : object
            Fitted encoder.
        """
        self._fit_encodings_all(X, y)
        return self

    def fit_transform(self, X, y):
        """Fit :class:`TargetRegressorEncoder` and transform X with the target encoding.

        .. note::
            :class:`TargetRegressorEncoder` uses a cross validation scheme in
            :meth:`fit_transform` to prevent leaking the target during training. In
            :meth:`fit_transform`, categorical encodings are obtained from one split and
            used to encoding the other split. Afterwards, a final categorical encoding
            is obtained from all the training data, which is used to encode data
            during :meth:`transform`. This means that `fit().transform()` does not equal
            `fit_transform()`.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categories of each feature.

        y : array-like of shape (n_samples,)
            The target data used to encode the categories.

        Returns
        -------
        X_trans : ndarray of shape (n_samples, n_features)
            Transformed input.
        """
        from ..model_selection._split import check_cv  # avoid circular input

        X_int, X_mask, y, n_categories = self._fit_encodings_all(X, y)
        cv = check_cv(self.cv)
        X_out = np.empty_like(X_int, dtype=np.float64)
        X_invalid = ~X_mask

        for train_idx, test_idx in cv.split(X, y):
            X_train, y_train = X_int[train_idx, :], y[train_idx]
            y_mean = np.mean(y_train)
            encodings = _fit_encoding_fast(
                X_train, y_train, n_categories, self.smooth, y_mean
            )
            self._transform_X_int(X_out, X_int, X_invalid, test_idx, encodings, y_mean)
        return X_out

    def transform(self, X):
        """Transform X with the target encoding.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categories of each feature.

        Returns
        -------
        X_trans : ndarray of shape (n_samples, n_features)
            Transformed input.
        """
        X_int, X_mask = self._transform(
            X, handle_unknown="ignore", force_all_finite="allow-nan"
        )
        X_out = np.empty_like(X_int, dtype=np.float64)
        self._transform_X_int(
            X_out, X_int, ~X_mask, slice(None), self.encodings_, self.encoding_mean_
        )
        return X_out

    def _fit_encodings_all(self, X, y):
        """Fit a target encoding with all the data."""
        self._validate_params()
        check_consistent_length(X, y)
        self._fit(X, handle_unknown="ignore", force_all_finite="allow-nan")
        y = _check_y(y, y_numeric=True, estimator=self).astype(np.float64, copy=False)
        self.encoding_mean_ = np.mean(y)

        X_int, X_mask = self._transform(
            X, handle_unknown="ignore", force_all_finite="allow-nan"
        )
        n_categories = np.fromiter(
            (len(cat) for cat in self.categories_),
            dtype=np.int64,
            count=len(self.categories_),
        )
        self.encodings_ = _fit_encoding_fast(
            X_int, y, n_categories, self.smooth, self.encoding_mean_
        )
        return X_int, X_mask, y, n_categories

    def _transform_X_int(self, X_out, X_int, X_invalid, indicies, encodings, y_mean):
        """Transform X_int using encodings."""
        for f_idx, encoding in enumerate(encodings):
            X_out[indicies, f_idx] = np.take(encoding, X_int[indicies, f_idx])
            X_out[X_invalid[:, f_idx], f_idx] = y_mean
