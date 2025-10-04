# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections.abc import Hashable
from numbers import Integral, Real

import numpy as np

from sklearn.base import OneToOneFeatureMixin, _fit_context
from sklearn.preprocessing._encoders import _BaseEncoder
from sklearn.preprocessing._target_encoder_fast import (
    _fit_encoding_fast,
    _fit_encoding_fast_auto_smooth,
)
from sklearn.utils._param_validation import Interval, StrOptions
from sklearn.utils.multiclass import type_of_target
from sklearn.utils.validation import (
    _check_feature_names_in,
    _check_y,
    check_consistent_length,
    check_is_fitted,
    validate_data,
)

# Distinct, hashable sentinels for each NA-like category (do not unify!)
_NONE_SENTINEL = object()  # exactly None
_NAN_SENTINEL = object()  # float NaN (np.nan, math.nan)
_PD_NA_SENTINEL = object()  # pandas.NA
_PD_NAT_SENTINEL = object()  # pandas.NaT
_NAT_SENTINEL = object()  # NumPy datetime/timedelta NaT (np.isnat == True)


def _looks_like_pandas_na(x) -> bool:
    """Return True if x is pandas.NA (without importing pandas)."""
    t = type(x)
    mod = getattr(t, "__module__", "")
    return t.__name__ == "NAType" and mod.startswith("pandas")


def _looks_like_pandas_nat(x) -> bool:
    """Return True if x is pandas.NaT (without importing pandas)."""
    t = type(x)
    mod = getattr(t, "__module__", "")
    # Some pandas versions use "NaTType", others expose "NaT" type name.
    return t.__name__ in {"NaTType", "NaT"} and mod.startswith("pandas")


def _norm_key(x):
    """Return a dict key that preserves category identity for NA-like values.

    We keep None, float NaN, pandas.NA, pandas.NaT, and NumPy NaT as distinct
    keys so the small-batch fast path matches the vectorized reference exactly.
    """
    # 1) exactly None
    if x is None:
        return _NONE_SENTINEL

    # 2) plain Python/NumPy float NaN
    try:
        # isinstance(x, float) cheaply catches np.float64 and Python float
        if isinstance(x, float) and np.isnan(x):
            return _NAN_SENTINEL
    except Exception:
        pass

    # 3) pandas extension sentinels (duck-typed; no pandas import required)
    try:
        if _looks_like_pandas_na(x):
            return _PD_NA_SENTINEL
        if _looks_like_pandas_nat(x):
            return _PD_NAT_SENTINEL
    except Exception:
        # Be defensive—if any weird object makes detection fail, fall through.
        pass

    # 4) NumPy datetime/timedelta NaT (non-float)
    try:
        if np.isnat(x):
            return _NAT_SENTINEL
    except Exception:
        pass

    # 5) everything else as-is
    return x


class TargetEncoder(OneToOneFeatureMixin, _BaseEncoder):
    """Target Encoder for regression and classification targets.

    Each category is encoded based on a shrunk estimate of the average target
    values for observations belonging to the category. The encoding scheme mixes
    the global target mean with the target mean conditioned on the value of the
    category (see [MIC]_).

    When the target type is "multiclass", encodings are based
    on the conditional probability estimate for each class. The target is first
    binarized using the "one-vs-all" scheme via
    :class:`~sklearn.preprocessing.LabelBinarizer`, then the average target
    value for each class and each category is used for encoding, resulting in
    `n_features` * `n_classes` encoded output features.

    :class:`TargetEncoder` considers missing values, such as `np.nan` or `None`,
    as another category and encodes them like any other category. Categories
    that are not seen during :meth:`fit` are encoded with the target mean, i.e.
    `target_mean_`.

    For a demo on the importance of the `TargetEncoder` internal cross-fitting,
    see
    :ref:`sphx_glr_auto_examples_preprocessing_plot_target_encoder_cross_val.py`.
    For a comparison of different encoders, refer to
    :ref:`sphx_glr_auto_examples_preprocessing_plot_target_encoder.py`. Read
    more in the :ref:`User Guide <target_encoder>`.

    .. note::
        `fit(X, y).transform(X)` does not equal `fit_transform(X, y)` because a
        :term:`cross fitting` scheme is used in `fit_transform` for encoding.
        See the :ref:`User Guide <target_encoder>` for details.

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

        .. versionchanged:: 1.4
           Added the option 'multiclass'.

    smooth : "auto" or float, default="auto"
        The amount of mixing of the target mean conditioned on the value of the
        category with the global target mean. A larger `smooth` value will put
        more weight on the global target mean.
        If `"auto"`, then `smooth` is set to an empirical Bayes estimate.

    cv : int, default=5
        Determines the number of folds in the :term:`cross fitting` strategy used in
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
    encodings_ : list of shape (n_features,) or (n_features * n_classes) of \
                    ndarray
        Encodings learnt on all of `X`.
        For feature `i`, `encodings_[i]` are the encodings matching the
        categories listed in `categories_[i]`. When `target_type_` is
        "multiclass", the encoding for feature `i` and class `j` is stored in
        `encodings_[j + (i * len(classes_))]`. E.g., for 2 features (f) and
        3 classes (c), encodings are ordered:
        f0_c0, f0_c1, f0_c2, f1_c0, f1_c1, f1_c2,

    categories_ : list of shape (n_features,) of ndarray
        The categories of each input feature determined during fitting or
        specified in `categories`
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
        If `target_type_` is 'binary' or 'multiclass', holds the label for each class,
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
    np.float64(44.3)
    >>> enc_high_smooth.encodings_
    [array([44.1, 44.4, 44.3])]

    >>> # On the other hand, a low `smooth` parameter puts more weight on target
    >>> # conditioned on the value of the categorical:
    >>> enc_low_smooth = TargetEncoder(smooth=1.0).fit(X, y)
    >>> enc_low_smooth.encodings_
    [array([21, 80.8, 43.2])]
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
        # --- existing public params ---
        self.categories = categories
        self.smooth = smooth
        self.target_type = target_type
        self.cv = cv
        self.shuffle = shuffle
        self.random_state = random_state

        # --- private: small-batch fast-path controls & lazy caches ---
        # Heuristic threshold for tiny batches (internal, not user-facing).
        # We only trigger the dict-lookup fast path when n_samples <= this,
        # or when categories >> rows (see transform()).
        self._small_batch_threshold = 256

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
        # in TargetEncoder.fit(...)
        X_df = X if hasattr(X, "columns") else None

        # validate once; this sets n_features_in_ and catches shape/type issues
        X, y = validate_data(
            self,
            X,
            y,
            ensure_2d=True,
            dtype=None,
            reset=True,
            ensure_all_finite="allow-nan",
        )

        # if X_df exists, make sure it has the same shape as validated X
        if X_df is not None:
            try:
                if getattr(X_df, "shape", None) != X.shape:
                    X_df = None  # shape mismatch; fall back to ndarray
            except Exception:
                X_df = None

        # set feature_names_in_ only if DataFrame *and* all string columns
        if X_df is not None:
            cols = np.asarray(X_df.columns, dtype=object)
            if cols.size == X.shape[1] and all(isinstance(c, str) for c in cols):
                self.feature_names_in_ = cols

        # now call the internal fitter with a DataFrame when available,
        # otherwise use the validated ndarray. This lets _BaseEncoder._fit
        # discover DataFrame metadata when present.
        self._fit_encodings_all(X_df if X_df is not None else X, y)

        # ---- Lazy fast-path initialization ----
        # We only allocate placeholders; actual per-feature maps/views are built
        # on the first small-batch transform that needs them.
        # categories_ is a list-like of length n_features after fitting.
        try:
            n_features = len(self.categories_)
        except Exception as e:  # be defensive if upstream changes
            raise AttributeError(
                "TargetEncoder.fit must set 'categories_' before fast path."
            ) from e

        # ---- initialize lazy fast-path placeholders (underscore attrs) ----
        self._te_index_maps_ = [None] * n_features
        self._te_enc_vecs_ = [None] * n_features
        self._te_enc_blocks_ = [None] * n_features
        self._te_defaults_ = [None] * n_features
        self._te_is_multiclass_ = getattr(self, "target_type_", None) == "multiclass"

        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def fit_transform(self, X, y):
        """Fit :class:`TargetEncoder` and transform X with the target encoding.

        .. note::
            `fit(X, y).transform(X)` does not equal `fit_transform(X, y)` because a
            :term:`cross fitting` scheme is used in `fit_transform` for encoding.
            See the :ref:`User Guide <target_encoder>`. for details.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data to determine the categories of each feature.

        y : array-like of shape (n_samples,)
            The target data used to encode the categories.

        Returns
        -------
        X_trans : ndarray of shape (n_samples, n_features) or \
                    (n_samples, (n_features * n_classes))
            Transformed input.
        """
        from sklearn.model_selection import (  # avoid circular import
            KFold,
            StratifiedKFold,
        )

        X_ordinal, X_known_mask, y_encoded, n_categories = self._fit_encodings_all(X, y)
        # ---- Lazy fast-path initialization (same as fit) ----
        try:
            n_features = len(self.categories_)
        except Exception as e:
            raise AttributeError(
                "TargetEncoder.fit_transform must set 'categories_' before fast path."
            ) from e

        self._te_index_maps_ = [None] * n_features
        self._te_enc_vecs_ = [None] * n_features
        self._te_enc_blocks_ = [None] * n_features
        self._te_defaults_ = [None] * n_features
        self._te_is_multiclass_ = getattr(self, "target_type_", None) == "multiclass"

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
                encodings = self._fit_encoding_multiclass(
                    X_train,
                    y_train,
                    n_categories,
                    y_train_mean,
                )
            else:
                encodings = self._fit_encoding_binary_or_continuous(
                    X_train,
                    y_train,
                    n_categories,
                    y_train_mean,
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

    def _ensure_fastpath_structs_for_feature(self, j: int) -> None:
        """Lazily build fast-path caches for feature j.
        Caches:
        - self._te_index_maps_[j]: dict normalized_category -> index in categories_[j]
        - self._te_enc_vecs_[j]:   1D view (n_cats_j,) for regression/binary
        - self._te_enc_blocks_[j]: 2D view (n_classes, n_cats_j) for multiclass
          ( class order == self.classes_ )
        - self._te_defaults_[j]:   scalar (binary/regression) or 1D vector (n_classes,)
          ( for unseen categories )
        """
        # already built
        maps = getattr(self, "_te_index_maps_", None)
        if maps is not None and maps[j] is not None:
            return

        # sanity (fit must have happened)
        if not hasattr(self, "categories_") or not hasattr(self, "encodings_"):
            raise AttributeError(
                "TargetEncoder must be fitted before using the fast path."
            )

        # ensure placeholders exist (fit created them)
        if self._te_index_maps_ is None:
            n_features = len(self.categories_)
            self._te_index_maps_ = [None] * n_features
            self._te_enc_vecs_ = [None] * n_features
            self._te_enc_blocks_ = [None] * n_features
            self._te_defaults_ = [None] * n_features
            self._te_is_multiclass_ = (
                getattr(self, "target_type_", None) == "multiclass"
            )

        cats_j = self.categories_[j]
        n_cats_j = len(cats_j)
        is_multi = bool(self._te_is_multiclass_)

        # -------- index map: category -> index (distinct NA-like sentinels) --------
        index_map: dict[Hashable, int] = {}
        for i in range(n_cats_j):
            key = _norm_key(cats_j[i])
            # Paranoia: if a collision happens (shouldn’t with distinct sentinels),
            # disable fastpath for this feature by leaving its map as None.
            if key in index_map and index_map[key] != i:
                self._te_index_maps_[j] = None
                self._te_enc_vecs_[j] = None
                self._te_enc_blocks_[j] = None
                self._te_defaults_[j] = None
                return
            index_map[key] = i

        if not is_multi:
            # regression/binary: encodings_[j] is 1D vector (n_cats_j,)
            enc_vec = np.asarray(self.encodings_[j])
            if enc_vec.ndim != 1:
                enc_vec = enc_vec.reshape(-1)
            self._te_enc_vecs_[j] = enc_vec

            # default for unseen: scalar target mean
            default = np.asarray(self.target_mean_, dtype=float)
            if default.ndim != 0:
                # fall back to global mean of enc_vec (should not happen normally)
                default = np.asarray(float(np.mean(enc_vec)))
            self._te_defaults_[j] = default
        else:
            # multiclass: encodings_ is feature-major, class-fast:
            # e_idx for (feature j, class c) == j * n_classes + c
            n_classes = len(self.classes_)
            # stack per-class 1D arrays into a (n_classes, n_cats_j) block
            block = np.empty((n_classes, n_cats_j), dtype=float)
            base = j * n_classes
            for c in range(n_classes):
                block[c, :] = self.encodings_[base + c]
            self._te_enc_blocks_[j] = block

            # default for unseen: 1D vector (n_classes,) in self.classes_ order
            default = np.asarray(self.target_mean_, dtype=float)
            if default.ndim == 0:
                default = np.full((n_classes,), float(default), dtype=float)
            elif default.ndim == 1 and default.shape[0] == n_classes:
                # good as-is
                pass
            else:
                # robust fallback: per-class mean across categories
                default = block.mean(axis=1)
            self._te_defaults_[j] = default

        self._te_index_maps_[j] = index_map

    def transform(self, X):
        """Transform X with the target encoding.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data to encode. Missing values (e.g. ``None`` or ``np.nan``)
            are treated as categories. Categories unseen during :meth:`fit`
            are encoded with ``target_mean_``.

        Returns
        -------
        X_trans : ndarray of shape (n_samples, n_features) or \
                (n_samples, n_features * n_classes)
            Encoded representation of ``X``. For binary and continuous targets,
            one column per input feature is returned. For multiclass targets,
            one column per (feature, class) pair is returned, with classes ordered
            as in ``classes_``.
        """
        # Decide fast path using only cheap metadata;
        # no heavy validation/conversion here.
        check_is_fitted(self)
        X_checked = validate_data(
            self,
            X,
            reset=False,
            ensure_2d=True,
            dtype=None,
            ensure_all_finite="allow-nan",
        )
        n_samples = X_checked.shape[0]
        small_thresh = getattr(self, "_small_batch_threshold", 256)
        use_small_batch = n_samples is not None and n_samples <= small_thresh

        # ---- Small-batch path (tiny inputs only) ---------------------------------
        if (
            use_small_batch
            and getattr(self, "categories_", None) is not None
            and getattr(self, "encodings_", None) is not None
        ):
            X_arr = np.asarray(X_checked, dtype=object)
            n_samples, n_features = X_arr.shape
            is_multi = getattr(self, "target_type_", None) == "multiclass"
            if is_multi:
                n_classes = len(self.classes_)
                X_out = np.empty((n_samples, n_features * n_classes), dtype=np.float64)
            else:
                X_out = np.empty((n_samples, n_features), dtype=np.float64)
                n_classes = 0  # unused below

            norm_key = _norm_key  # local ref for speed

            for j in range(n_features):
                self._ensure_fastpath_structs_for_feature(j)

                idx_map = self._te_index_maps_[j]
                get = idx_map.get
                col = X_arr[:, j]

                if not is_multi:
                    enc_vec = self._te_enc_vecs_[j]
                    default = float(self._te_defaults_[j])
                    out_col = np.empty(n_samples, dtype=float)
                    for i in range(n_samples):
                        idx = get(norm_key(col[i]), -1)
                        out_col[i] = enc_vec[idx] if idx >= 0 else default
                    X_out[:, j] = out_col
                else:
                    enc_block = self._te_enc_blocks_[j]  # (n_classes, n_cats_j)
                    default_v = np.asarray(self._te_defaults_[j], dtype=float)
                    out_block = np.empty((n_samples, n_classes), dtype=float)
                    for i in range(n_samples):
                        idx = get(norm_key(col[i]), -1)
                        out_block[i, :] = enc_block[:, idx] if idx >= 0 else default_v
                    start = j * n_classes
                    X_out[:, start : start + n_classes] = out_block

            return X_out

        # ---- Large/normal batches: baseline vectorized path ----------------------
        X_ordinal, X_known_mask = self._transform(
            X, handle_unknown="ignore", ensure_all_finite="allow-nan"
        )

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
        from sklearn.preprocessing import LabelBinarizer, LabelEncoder

        check_consistent_length(X, y)
        self._fit(X, handle_unknown="ignore", ensure_all_finite="allow-nan")

        if self.target_type == "auto":
            accepted_target_types = ("binary", "multiclass", "continuous")
            inferred_type_of_target = type_of_target(y, input_name="y")
            if inferred_type_of_target not in accepted_target_types:
                raise ValueError(
                    "Unknown label type: Target type was inferred to be "
                    f"{inferred_type_of_target!r}. Only {accepted_target_types} are "
                    "supported."
                )
            self.target_type_ = inferred_type_of_target
        else:
            self.target_type_ = self.target_type

        self.classes_ = None
        if self.target_type_ == "binary":
            label_encoder = LabelEncoder()
            y = label_encoder.fit_transform(y)
            self.classes_ = label_encoder.classes_
        elif self.target_type_ == "multiclass":
            label_binarizer = LabelBinarizer()
            y = label_binarizer.fit_transform(y)
            self.classes_ = label_binarizer.classes_
        else:  # continuous
            y = _check_y(y, y_numeric=True, estimator=self)

        self.target_mean_ = np.mean(y, axis=0)

        X_ordinal, X_known_mask = self._transform(
            X, handle_unknown="ignore", ensure_all_finite="allow-nan"
        )
        n_categories = np.fromiter(
            (len(category_for_feature) for category_for_feature in self.categories_),
            dtype=np.int64,
            count=len(self.categories_),
        )
        if self.target_type_ == "multiclass":
            encodings = self._fit_encoding_multiclass(
                X_ordinal,
                y,
                n_categories,
                self.target_mean_,
            )
        else:
            encodings = self._fit_encoding_binary_or_continuous(
                X_ordinal,
                y,
                n_categories,
                self.target_mean_,
            )
        self.encodings_ = encodings

        return X_ordinal, X_known_mask, y, n_categories

    def _fit_encoding_binary_or_continuous(
        self, X_ordinal, y, n_categories, target_mean
    ):
        """Learn target encodings."""
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
        return encodings

    def _fit_encoding_multiclass(self, X_ordinal, y, n_categories, target_mean):
        """Learn multiclass encodings.

        Learn encodings for each class (c) then reorder encodings such that
        the same features (f) are grouped together. `reorder_index` enables
        converting from:
        f0_c0, f1_c0, f0_c1, f1_c1, f0_c2, f1_c2
        to:
        f0_c0, f0_c1, f0_c2, f1_c0, f1_c1, f1_c2
        """
        n_features = self.n_features_in_
        n_classes = len(self.classes_)

        encodings = []
        for i in range(n_classes):
            y_class = y[:, i]
            encoding = self._fit_encoding_binary_or_continuous(
                X_ordinal,
                y_class,
                n_categories,
                target_mean[i],
            )
            encodings.extend(encoding)

        reorder_index = (
            idx
            for start in range(n_features)
            for idx in range(start, (n_classes * n_features), n_features)
        )
        return [encodings[idx] for idx in reorder_index]

    def _transform_X_ordinal(
        self,
        X_out,
        X_ordinal,
        X_unknown_mask,
        row_indices,
        encodings,
        target_mean,
    ):
        """Transform X_ordinal using encodings.

        In the multiclass case, `X_ordinal` and `X_unknown_mask` have column
        (axis=1) size `n_features`, while `encodings` has length of size
        `n_features * n_classes`. `feat_idx` deals with this by repeating
        feature indices by `n_classes` E.g., for 3 features, 2 classes:
        0,0,1,1,2,2

        Additionally, `target_mean` is of shape (`n_classes`,) so `mean_idx`
        cycles through 0 to `n_classes` - 1, `n_features` times.
        """
        if self.target_type_ == "multiclass":
            n_classes = len(self.classes_)
            for e_idx, encoding in enumerate(encodings):
                # Repeat feature indices by n_classes
                feat_idx = e_idx // n_classes
                # Cycle through each class
                mean_idx = e_idx % n_classes
                X_out[row_indices, e_idx] = encoding[X_ordinal[row_indices, feat_idx]]
                X_out[X_unknown_mask[:, feat_idx], e_idx] = target_mean[mean_idx]
        else:
            for e_idx, encoding in enumerate(encodings):
                X_out[row_indices, e_idx] = encoding[X_ordinal[row_indices, e_idx]]
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

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.target_tags.required = True
        return tags
