from numbers import Integral, Real

import numpy as np

from ..base import OneToOneFeatureMixin, _fit_context
from ..utils._param_validation import Interval, StrOptions
from ..utils.multiclass import type_of_target
from ..utils.validation import (
    _check_feature_names_in,
    _check_y,
    check_consistent_length,
    check_is_fitted,
)
from ._encoders import _BaseEncoder
from ._target_encoder_fast import _fit_encoding_fast, _fit_encoding_fast_auto_smooth


class TargetEncoder(OneToOneFeatureMixin, _BaseEncoder):
    """Target Encoder for regression and classification targets.

    Each category is encoded based on a shrunk estimate of the average target
    values for observations belonging to the category. The encoding scheme mixes
    the global target mean with the target mean conditioned on the value of the
    category. [MIC]_

    :class:`TargetEncoder` considers missing values, such as `np.nan` or `None`,
    as another category and encodes them like any other category. Categories
    that are not seen during :meth:`fit` are encoded with the target mean, i.e.
    `target_mean_`.

    Read more in the :ref:`User Guide <target_encoder>`.

    .. note::
        `fit(X, y).transform(X)` does not equal `fit_transform(X, y)` because a
        cross fitting scheme is used in `fit_transform` for encoding. See the
        :ref:`User Guide <target_encoder>` for details.

    .. versionadded:: 1.3

    Parameters
    ----------
    categories : "auto" or list of shape (n_features,) of array-like, default="auto"
        Categories (unique values) per feature:

        - `"auto"` : Determine categories automatically from the training data.
        - list : `categories[i]` holds the categories expected in the i-th column. The
          passed categories should not mix strings and numeric values within a single
          feature, and should be sorted in case of numeric values.

        The used categories are stored in the `categories_` fitted attribute.

    target_type : {"auto", "continuous", "binary", "multiclass"}, default="auto"
        Type of target.

        - `"auto"` : Type of target is inferred with
          :func:`~sklearn.utils.multiclass.type_of_target`.
        - `"continuous"` : Continuous target
        - `"binary"` : Binary target
        - `"multiclass"` : Multiclass target

        .. note::
            The type of target inferred with `"auto"` may not be the desired target
            type used for modeling. For example, if the target consisted of integers
            between 0 and 100, then :func:`~sklearn.utils.multiclass.type_of_target`
            will infer the target as `"multiclass"`. In this case, setting
            `target_type="continuous"` will specify the target as a regression
            problem. The `target_type_` attribute gives the target type used by the
            encoder.

        .. versionadded:: 1.4
           Added the option 'multiclass'.

    smooth : "auto" or float, default="auto"
        The amount of mixing of the target mean conditioned on the value of the
        category with the global target mean. A larger `smooth` value will put
        more weight on the global target mean.
        If `"auto"`, then `smooth` is set to an empirical Bayes estimate.

    cv : int, default=5
        Determines the number of folds in the cross fitting strategy used in
        :meth:`fit_transform`. For classification targets, `StratifiedKFold` is used
        and for continuous targets, `KFold` is used.

    shuffle : bool, default=True
        Whether to shuffle the data in :meth:`fit_transform` before splitting into
        folds. Note that the samples within each split will not be shuffled.

    random_state : int, RandomState instance or None, default=None
        When `shuffle` is True, `random_state` affects the ordering of the
        indices, which controls the randomness of each fold. Otherwise, this
        parameter has no effect.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    encodings_ : list of shape (n_features,) of ndarray
        Encodings learnt on all of `X`.
        For feature `i`, `encodings_[i]` are the encodings matching the
        categories listed in `categories_[i]`.

    categories_ : list of shape (n_features,) of ndarray
        The categories of each feature determined during fitting or specified
        in `categories`
        (in order of the features in `X` and corresponding with the output
        of :meth:`transform`).

    target_type_ : str
        Type of target.

    target_mean_ : float
        The overall mean of the target. This value is only used in :meth:`transform`
        to encode categories.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    classes_ : ndarray or None
        If `target_type_` is 'multiclass', holds the label for each class,
        otherwise `None`.

    See Also
    --------
    OrdinalEncoder : Performs an ordinal (integer) encoding of the categorical features.
        Contrary to TargetEncoder, this encoding is not supervised. Treating the
        resulting encoding as a numerical features therefore lead arbitrarily
        ordered values and therefore typically lead to lower predictive performance
        when used as preprocessing for a classifier or regressor.
    OneHotEncoder : Performs a one-hot encoding of categorical features. This
        unsupervised encoding is better suited for low cardinality categorical
        variables as it generate one new feature per unique category.

    References
    ----------
    .. [MIC] :doi:`Micci-Barreca, Daniele. "A preprocessing scheme for high-cardinality
       categorical attributes in classification and prediction problems"
       SIGKDD Explor. Newsl. 3, 1 (July 2001), 27–32. <10.1145/507533.507538>`

    Examples
    --------
    With `smooth="auto"`, the smoothing parameter is set to an empirical Bayes estimate:

    >>> import numpy as np
    >>> from sklearn.preprocessing import TargetEncoder
    >>> X = np.array([["dog"] * 20 + ["cat"] * 30 + ["snake"] * 38], dtype=object).T
    >>> y = [90.3] * 5 + [80.1] * 15 + [20.4] * 5 + [20.1] * 25 + [21.2] * 8 + [49] * 30
    >>> enc_auto = TargetEncoder(smooth="auto")
    >>> X_trans = enc_auto.fit_transform(X, y)

    >>> # A high `smooth` parameter puts more weight on global mean on the categorical
    >>> # encodings:
    >>> enc_high_smooth = TargetEncoder(smooth=5000.0).fit(X, y)
    >>> enc_high_smooth.target_mean_
    44...
    >>> enc_high_smooth.encodings_
    [array([44..., 44..., 44...])]

    >>> # On the other hand, a low `smooth` parameter puts more weight on target
    >>> # conditioned on the value of the categorical:
    >>> enc_low_smooth = TargetEncoder(smooth=1.0).fit(X, y)
    >>> enc_low_smooth.encodings_
    [array([20..., 80..., 43...])]
    """

    _parameter_constraints: dict = {
        "categories": [StrOptions({"auto"}), list],
        "target_type": [StrOptions({"auto", "continuous", "binary", "multiclass"})],
        "smooth": [StrOptions({"auto"}), Interval(Real, 0, None, closed="left")],
        "cv": [Interval(Integral, 2, None, closed="left")],
        "shuffle": ["boolean"],
        "random_state": ["random_state"],
    }

    def __init__(
        self,
        categories="auto",
        target_type="auto",
        smooth="auto",
        cv=5,
        shuffle=True,
        random_state=None,
    ):
        self.categories = categories
        self.smooth = smooth
        self.target_type = target_type
        self.cv = cv
        self.shuffle = shuffle
        self.random_state = random_state

    @_fit_context(prefer_skip_nested_validation=True)
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

    @_fit_context(prefer_skip_nested_validation=True)
    def fit_transform(self, X, y):
        """Fit :class:`TargetEncoder` and transform X with the target encoding.

        .. note::
            `fit(X, y).transform(X)` does not equal `fit_transform(X, y)` because a
            cross fitting scheme is used in `fit_transform` for encoding. See the
            :ref:`User Guide <target_encoder>`. for details.

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
        from ..model_selection import KFold, StratifiedKFold  # avoid circular import

        X_ordinal, X_known_mask, y_encoded, n_categories = self._fit_encodings_all(X, y)

        # The cv splitter is voluntarily restricted to *KFold to enforce non
        # overlapping validation folds, otherwise the fit_transform output will
        # not be well-specified.
        if self.target_type_ == "continuous":
            cv = KFold(self.cv, shuffle=self.shuffle, random_state=self.random_state)
        else:
            cv = StratifiedKFold(
                self.cv, shuffle=self.shuffle, random_state=self.random_state
            )

        # If 'multiclass' multiply axis=1 by num classes else keep shape the same
        if self.target_type_ == "multiclass":
            X_out = np.empty(
                (X_ordinal.shape[0], X_ordinal.shape[1] * len(self.classes_)),
                dtype=np.float64,
            )
        else:
            X_out = np.empty_like(X_ordinal, dtype=np.float64)

        for train_idx, test_idx in cv.split(X, y):
            X_train, y_train = X_ordinal[train_idx, :], y_encoded[train_idx]
            y_train_mean = np.mean(y_train, axis=0)

            if self.target_type_ == "multiclass":
                encodings = []
                for i in range(len(self.classes_)):
                    y_class = y_train[:, i]
                    encoding = self._learn_encodings(
                        X_train,
                        y_class,
                        n_categories,
                        y_train_mean[i],
                        save_encodings=False,
                    )
                    encodings += encoding
            else:
                encodings = self._learn_encodings(
                    X_train,
                    y_train,
                    n_categories,
                    y_train_mean,
                    save_encodings=False,
                )
            self._transform_X_ordinal(
                X_out,
                X_ordinal,
                ~X_known_mask,
                test_idx,
                encodings,
                y_train_mean,
            )
        return X_out

    def transform(self, X):
        """Transform X with the target encoding.

        .. note::
            `fit(X, y).transform(X)` does not equal `fit_transform(X, y)` because a
            cross fitting scheme is used in `fit_transform` for encoding. See the
            :ref:`User Guide <target_encoder>`. for details.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categories of each feature.

        Returns
        -------
        X_trans : ndarray of shape (n_samples, n_features)
            Transformed input.
        """
        X_ordinal, X_known_mask = self._transform(
            X, handle_unknown="ignore", force_all_finite="allow-nan"
        )

        # If 'multiclass' multiply axis=1 by num of classes else keep shape the same
        if self.target_type_ == "multiclass":
            X_out = np.empty(
                (X_ordinal.shape[0], X_ordinal.shape[1] * len(self.classes_)),
                dtype=np.float64,
            )
        else:
            X_out = np.empty_like(X_ordinal, dtype=np.float64)

        self._transform_X_ordinal(
            X_out,
            X_ordinal,
            ~X_known_mask,
            slice(None),
            self.encodings_,
            self.target_mean_,
        )
        return X_out

    def _fit_encodings_all(self, X, y):
        """Fit a target encoding with all the data."""
        # avoid circular import
        from ..preprocessing import (
            LabelBinarizer,
            LabelEncoder,
        )

        check_consistent_length(X, y)
        self._fit(X, handle_unknown="ignore", force_all_finite="allow-nan")

        if self.target_type == "auto":
            accepted_target_types = ("binary", "multiclass", "continuous")
            inferred_type_of_target = type_of_target(y, input_name="y")
            if inferred_type_of_target not in accepted_target_types:
                raise ValueError(
                    f"Target type was inferred to be {inferred_type_of_target!r}. Only"
                    f" {accepted_target_types} are supported."
                )
            self.target_type_ = inferred_type_of_target
        else:
            self.target_type_ = self.target_type

        self.classes_ = None
        if self.target_type_ == "binary":
            y = LabelEncoder().fit_transform(y)
        elif self.target_type_ == "multiclass":
            label_binarizer = LabelBinarizer()
            y = label_binarizer.fit_transform(y)
            self.classes_ = label_binarizer.classes_
        else:  # continuous
            y = _check_y(y, y_numeric=True, estimator=self)

        self.target_mean_ = np.mean(y, axis=0)

        X_ordinal, X_known_mask = self._transform(
            X, handle_unknown="ignore", force_all_finite="allow-nan"
        )
        n_categories = np.fromiter(
            (len(category_for_feature) for category_for_feature in self.categories_),
            dtype=np.int64,
            count=len(self.categories_),
        )
        if self.target_type_ == "multiclass":
            encodings = []
            for i in range(len(self.classes_)):
                y_class = y[:, i]
                encoding = self._learn_encodings(
                    X_ordinal,
                    y_class,
                    n_categories,
                    self.target_mean_[i],
                    save_encodings=False,
                )
                encodings += encoding
            self.encodings_ = encodings
        else:
            self._learn_encodings(
                X_ordinal,
                y,
                n_categories,
                self.target_mean_,
                save_encodings=True,
            )

        return X_ordinal, X_known_mask, y, n_categories

    def _learn_encodings(self, X_ordinal, y, n_categories, target_mean, save_encodings):
        """Learn target encodings, optionally saving as attribute."""
        if self.smooth == "auto":
            y_variance = np.var(y)
            encodings = _fit_encoding_fast_auto_smooth(
                X_ordinal,
                y,
                n_categories,
                target_mean,
                y_variance,
            )
        else:
            encodings = _fit_encoding_fast(
                X_ordinal,
                y,
                n_categories,
                self.smooth,
                target_mean,
            )

        if save_encodings:
            self.encodings_ = encodings
        return encodings

    def _transform_X_ordinal(
        self,
        X_out,
        X_ordinal,
        X_unknown_mask,
        indices,
        encodings,
        target_mean,
    ):
        """Transform X_ordinal using encodings.

        In the multiclass case, `X_ordinal` and `X_unknown_mask` have column
        (axis=1) size `n_features`, while column of `X_out` and length of
        `encodings` length is of size `n_features * n_classes`. `multi_idx`
        deals with this by repeating feature indicies for each class.
        E.g., for 3 features, 2 classes: 0,1,2,0,1,2

        Additionally, `X_out` needs to order columns such that features (f)
        are grouped together:
        f0_c0, f0_c1, f1_c0, f1_c1, f2_c0, f2_c1
        instead of classes (c) being grouped together, as in `encodings`:
        f0_c0, f1_c0, f2_c0, f0_c1, f1_c1, f2_c1
        This is handled by `out_indx`.
        """

        if self.target_type_ == "multiclass":
            n_classes = len(self.classes_)
            n_features = self.n_features_in_

        for e_idx, encoding in enumerate(encodings):
            if self.target_type_ == "multiclass":
                # `X_out` should have features grouped together
                out_indx = (e_idx // n_classes) + ((e_idx % n_classes) * n_features)
                # Repeat feature indicies for each class
                feat_idx = e_idx - (n_features * (e_idx // n_features))
                mean_idx = e_idx // n_features
                X_out[indices, out_indx] = encoding[X_ordinal[indices, feat_idx]]
                X_out[X_unknown_mask[:, feat_idx], e_idx] = target_mean[mean_idx]
            else:
                X_out[indices, e_idx] = encoding[X_ordinal[indices, e_idx]]
                X_out[X_unknown_mask[:, e_idx], e_idx] = target_mean

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Not used, present here for API consistency by convention.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Transformed feature names. `feature_names_in_` is used unless it is
            not defined, in which case the following input feature names are
            generated: `["x0", "x1", ..., "x(n_features_in_ - 1)"]`.
            When `type_of_target_` is "multiclass" the names are of the format
            '<feature_name>_<class_name>'.
        """
        check_is_fitted(self, "n_features_in_")
        feature_names = _check_feature_names_in(self, input_features)
        if self.target_type_ == "multiclass":
            feature_names = [
                f"{feature_name}_{class_name}"
                for feature_name in feature_names
                for class_name in self.classes_
            ]
            return np.asarray(feature_names, dtype=object)
        else:
            return feature_names

    def _more_tags(self):
        return {"requires_y": True}
