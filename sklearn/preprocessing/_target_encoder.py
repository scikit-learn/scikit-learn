import numpy as np

from numbers import Real

from ._encoders import _BaseEncoder
from ..base import OneToOneFeatureMixin
from ._target_encoder_fast import _fit_encoding_fast
from ._target_encoder_fast import _fit_encoding_fast_auto_smooth
from ..utils.validation import _check_y, check_consistent_length
from ..utils.multiclass import type_of_target
from ..utils._param_validation import Interval, StrOptions


class TargetEncoder(OneToOneFeatureMixin, _BaseEncoder):
    """Target Encoder for regression and classification targets.

    Each category is encoded based on its effect on the target. The encoding
    scheme mixes the global target mean with the target mean conditioned on
    the value of the category. [MIC]_

    Read more in the :ref:`User Guide <target_encoder>`.

    .. note::
        :class:`TargetEncoder` uses a cross validation scheme in
        :meth:`fit_transform` to prevent leaking the target during training. In
        :meth:`fit_transform`, the data is split according to the `cv` parameter.
        Categorical encodings are learned from split and used to encode the other split.
        Afterwards, a final categorical encoding is learned from all the data, which
        is then used to encode data during :meth:`transform`. This means that
        `fit(X, y).transform(X)` does not equal `fit_transform(X, y)`.

    .. versionadded:: 1.3

    Parameters
    ----------
    categories : "auto" or a list of array-like, default="auto"
        Categories (unique values) per feature:

        - `"auto"` : Determine categories automatically from the training data.
        - list : `categories[i]` holds the categories expected in the i-th column. The
          passed categories should not mix strings and numeric values within a single
          feature, and should be sorted in case of numeric values.

        The used categories is stored in the `categories_` fitted attribute.

    target_type : {"auto", "continuous", "binary"}, default="auto"
        Type of target

        - `"auto"` : Type of target is inferred with
          :func:`~sklearn.utils.multiclass.type_of_target`
        - `"continuous"` : Continuous target
        - `"binary"` : Binary target

    smooth : "auto" or float, default="auto"
        The amount of mixing of the categorical encoding with the global target mean. A
        larger `smooth` value will put more weight on the global target mean.
        If `"auto"`, then `smooth` is estimated using an empirical bayes estimate.

    cv : int, cross-validation generator or an iterable, default=None
        Determines the cross-validation splitting strategy used in
        :meth:`fit_transform`. Possible inputs for cv are:

        - `None`, to use the default 5-fold cross validation,
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
        (in order of the features in `X` and corresponding with the output
        of :meth:`transform`).

    target_type_ : str
        Type of target.

    encoding_mean_ : float
        The overall mean of the target. This value is only used in :meth:`transform`
        to encode categories.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    OrdinalEncoder : Performs an ordinal (integer) encoding of the categorical features.
    OneHotEncoder : Performs a one-hot encoding of categorical features.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.preprocessing import TargetEncoder
    >>> X = np.array([["dog"] * 20 + ["cat"] * 30 + ["snake"] * 38], dtype=object).T
    >>> y = [10.3] * 5 + [40.1] * 15 + [20.4] * 5 + [11.1] * 25 + [21.2] * 8 + [49] * 30
    >>> enc = TargetEncoder(smooth=5.0)
    >>> X_trans = enc.fit_transform(X, y)

    References
    ----------
    .. [MIC] :doi:`Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
        categorical attributes in classification and prediction problems"
        SIGKDD Explor. Newsl. 3, 1 (July 2001), 27–32. <10.1145/507533.507538>`
    """

    _parameter_constraints: dict = {
        "categories": [StrOptions({"auto"}), list],
        "target_type": [StrOptions({"auto", "continuous", "binary"})],
        "smooth": [StrOptions({"auto"}), Interval(Real, 1, None, closed="left")],
        "cv": ["cv_object"],
    }

    def __init__(self, categories="auto", target_type="auto", smooth="auto", cv=5):
        self.categories = categories
        self.smooth = smooth
        self.target_type = target_type
        self.cv = cv

    def fit(self, X, y):
        """Fit the :class:`TargetEncoder` to X and y.

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
        """Fit :class:`TargetEncoder` and transform X with the target encoding.

        .. note::
            :class:`TargetEncoder` uses a cross validation scheme in
            :meth:`fit_transform` to prevent leaking the target during training. In
            :meth:`fit_transform`, the data is split according to the `cv` parameter.
            Categorical encodings are learned from split and used to encode the other
            split. Afterwards, a final categorical encoding is learned from all the
            data, which is then used to encode data during :meth:`transform`. This
            means that `fit(X, y).transform(X)` does not equal `fit_transform(X, y)`.

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
        from ..model_selection._split import check_cv  # avoid circular import

        X_int, X_mask, y, n_categories = self._fit_encodings_all(X, y)
        cv = check_cv(self.cv)
        X_out = np.empty_like(X_int, dtype=np.float64)
        X_invalid = ~X_mask

        for train_idx, test_idx in cv.split(X, y):
            X_train, y_train = X_int[train_idx, :], y[train_idx]
            y_mean = np.mean(y_train)

            if self.smooth == "auto":
                y_variance = np.var(y_train)
                encodings = _fit_encoding_fast_auto_smooth(
                    X_train, y_train, n_categories, y_mean, y_variance
                )
            else:
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
        from ..preprocessing import LabelEncoder  # avoid circular import

        self._validate_params()
        check_consistent_length(X, y)
        self._fit(X, handle_unknown="ignore", force_all_finite="allow-nan")

        if self.target_type == "auto":
            accepted_target_types = ("binary", "continuous")
            inferred_type_of_target = type_of_target(y, input_name="y")
            if inferred_type_of_target not in accepted_target_types:
                raise ValueError(
                    f"Target type was inferred to be {inferred_type_of_target!r}. Only"
                    f" {accepted_target_types} are supported."
                )
            self.target_type_ = inferred_type_of_target
        else:
            self.target_type_ = self.target_type

        if self.target_type_ == "binary":
            y = LabelEncoder().fit_transform(y)
        else:  # continuous
            y = _check_y(y, y_numeric=True, estimator=self)

        self.encoding_mean_ = np.mean(y)

        X_int, X_mask = self._transform(
            X, handle_unknown="ignore", force_all_finite="allow-nan"
        )
        n_categories = np.fromiter(
            (len(category_for_feature) for category_for_feature in self.categories_),
            dtype=np.int64,
            count=len(self.categories_),
        )
        if self.smooth == "auto":
            y_variance = np.var(y)
            self.encodings_ = _fit_encoding_fast_auto_smooth(
                X_int, y, n_categories, self.encoding_mean_, y_variance
            )
        else:
            self.encodings_ = _fit_encoding_fast(
                X_int, y, n_categories, self.smooth, self.encoding_mean_
            )

        return X_int, X_mask, y, n_categories

    def _transform_X_int(self, X_out, X_int, X_invalid, indices, encodings, y_mean):
        """Transform X_int using encodings."""
        for f_idx, encoding in enumerate(encodings):
            X_out[indices, f_idx] = encoding[X_int[indices, f_idx]]
            X_out[X_invalid[:, f_idx], f_idx] = y_mean
