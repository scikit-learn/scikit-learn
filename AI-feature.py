import warnings
from abc import ABCMeta, abstractmethod
from operator import attrgetter

import numpy as np
from scipy.sparse import csc_matrix, issparse
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
import shap
from ..base import TransformerMixin
from ..utils import _safe_indexing, check_array, safe_sqr
from ..utils._set_output import _get_output_config
from ..utils._tags import _safe_tags
from ..utils.validation import _check_feature_names_in, _is_pandas_df, check_is_fitted


class AISelectorMixin(TransformerMixin, metaclass=ABCMeta):
    """
    Transformer mixin that performs feature selection given a support mask,
    now enhanced with AI-driven feature importance analysis, adaptive selection,
    and explainability integration.

    This mixin provides a feature selector implementation with `transform` and
    `inverse_transform` functionality given an implementation of
    `_get_support_mask`.

    AI Enhancements:
    - Automated Feature Importance Analysis
    - Adaptive Feature Selection
    - Explainability Integration (using SHAP)

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.base import BaseEstimator
    >>> from sklearn.feature_selection import AISelectorMixin
    >>> class FeatureSelector(AISelectorMixin, BaseEstimator):
    ...    def fit(self, X, y=None):
    ...        super().fit(X, y)
    ...        return self
    >>> X, y = load_iris(return_X_y=True)
    >>> selector = FeatureSelector(model=RandomForestClassifier())
    >>> selector.fit_transform(X, y).shape
    (150, 2)
    """

    def __init__(self, model=None):
        if model is None:
            self.model = RandomForestClassifier()
        else:
            self.model = model
        self.history = []

    def fit(self, X, y):
        # Store the history of the data for adaptive selection
        self.history.append((X, y))
        if len(self.history) > 5:  # Example condition for adaptive logic
            self.adjust_model()
        # Fit the model to the data
        self.model.fit(X, y)
        # Get the feature importances from the fitted model
        self.feature_importances_ = _get_feature_importances(self.model, getter='auto')
        self.n_features_in_ = X.shape[1]
        return self

    def _get_support_mask(self):
        # Generate a mask of selected features based on their importances
        threshold = np.mean(self.feature_importances_)
        mask = self.feature_importances_ > threshold
        return mask

    def get_support(self, indices=False):
        """
        Get a mask, or integer index, of the features selected.

        Parameters
        ----------
        indices : bool, default=False
            If True, the return value will be an array of integers, rather
            than a boolean mask.

        Returns
        -------
        support : array
            An index that selects the retained features from a feature vector.
            If `indices` is False, this is a boolean array of shape
            [# input features], in which an element is True iff its
            corresponding feature is selected for retention. If `indices` is
            True, this is an integer array of shape [# output features] whose
            values are indices into the input feature vector.
        """
        mask = self._get_support_mask()
        return mask if not indices else np.where(mask)[0]

    def adjust_model(self):
        # Adaptively change the model based on historical performance
        if np.random.rand() > 0.5:  # Simplified adaptive model adjustment
            self.model = GradientBoostingClassifier()

    def explain_selection(self, X):
        """
        Explain the feature selection using SHAP.

        Parameters
        ----------
        X : array-like of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        shap_values : SHAP values object
            The SHAP values of the input features.
        """
        explainer = shap.Explainer(self.model, X)
        shap_values = explainer(X)
        shap.summary_plot(shap_values, X)
        return shap_values

    def transform(self, X):
        """Reduce X to the selected features.

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        X_r : array of shape [n_samples, n_selected_features]
            The input samples with only the selected features.
        """
        # Preserve X when X is a dataframe and the output is configured to
        # be pandas.
        output_config_dense = _get_output_config("transform", estimator=self)["dense"]
        preserve_X = output_config_dense != "default" and _is_pandas_df(X)

        # note: we use _safe_tags instead of _get_tags because this is a
        # public Mixin.
        X = self._validate_data(
            X,
            dtype=None,
            accept_sparse="csr",
            ensure_all_finite=not _safe_tags(self, key="allow_nan"),
            cast_to_ndarray=not preserve_X,
            reset=False,
        )
        return self._transform(X)

    def _transform(self, X):
        """Reduce X to the selected features."""
        mask = self.get_support()
        if not mask.any():
            warnings.warn(
                (
                    "No features were selected: either the data is"
                    " too noisy or the selection test too strict."
                ),
                UserWarning,
            )
            if hasattr(X, "iloc"):
                return X.iloc[:, :0]
            return np.empty(0, dtype=X.dtype).reshape((X.shape[0], 0))
        return _safe_indexing(X, mask, axis=1)

    def inverse_transform(self, X):
        """Reverse the transformation operation.

        Parameters
        ----------
        X : array of shape [n_samples, n_selected_features]
            The input samples.

        Returns
        -------
        X_r : array of shape [n_samples, n_original_features]
            `X` with columns of zeros inserted where features would have
            been removed by :meth:`transform`.
        """
        if issparse(X):
            X = X.tocsc()
            # insert additional entries in indptr:
            # e.g. if transform changed indptr from [0 2 6 7] to [0 2 3]
            # col_nonzeros here will be [2 0 1] so indptr becomes [0 2 2 3]
            it = self.inverse_transform(np.diff(X.indptr).reshape(1, -1))
            col_nonzeros = it.ravel()
            indptr = np.concatenate([[0], np.cumsum(col_nonzeros)])
            Xt = csc_matrix(
                (X.data, X.indices, indptr),
                shape=(X.shape[0], len(indptr) - 1),
                dtype=X.dtype,
            )
            return Xt

        support = self.get_support()
        X = check_array(X, dtype=None)
        if support.sum() != X.shape[1]:
            raise ValueError("X has a different shape than during fitting.")

        if X.ndim == 1:
            X = X[None, :]
        Xt = np.zeros((X.shape[0], support.size), dtype=X.dtype)
        Xt[:, support] = X
        return Xt

    def get_feature_names_out(self, input_features=None):
        """Mask feature names according to selected features.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Input features.

            - If `input_features` is `None`, then `feature_names_in_` is
              used as feature names in. If `feature_names_in_` is not defined,
              then the following input feature names are generated:
              `["x0", "x1", ..., "x(n_features_in_ - 1)"]`.
            - If `input_features` is an array-like, then `input_features` must
              match `feature_names_in_` if `feature_names_in_` is defined.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Transformed feature names.
        """
        check_is_fitted(self)
        input_features = _check_feature_names_in(self, input_features)
        return input_features[self.get_support()]


def _get_feature_importances(estimator, getter, transform_func=None, norm_order=1):
    """
    Retrieve and aggregate (ndim > 1) the feature importances
    from an estimator. Also optionally applies transformation.

    Parameters
    ----------
    estimator : estimator
        A scikit-learn estimator from which we want to get the feature
        importances.

    getter : "auto", str or callable
        An attribute or a callable to get the feature importance. If `"auto"`,
        `estimator` is expected to expose `coef_` or `feature_importances_`.

    transform_func : {"norm", "square"}, default=None
        The transform to apply to the feature importances. By default (`None`)
        no transformation is applied.

    norm_order : int, default=1
        The norm order to apply when `transform_func="norm"`. Only applied
        when `importances.ndim > 1`.

    Returns
    -------
    importances : ndarray of shape (n_features,)
        The features importances, optionally transformed.
    """
    if isinstance(getter, str):
        if getter == "auto":
            if hasattr(estimator, "coef_"):
                getter = attrgetter("coef_")
            elif hasattr(estimator, "feature_importances_"):
                getter = attrgetter("feature_importances_")
            else:
                raise ValueError(
                    "`estimator` has to expose `coef_` or "
                    "`feature_importances_` or a custom `getter` "
                    "has to be defined"
                )
        else:
            getter = attrgetter(getter)
    elif not callable(getter):
        raise ValueError("`getter` has to be a string or callable")

    importances = getter(estimator)
    importances = importances if importances.ndim <= 1 else np.linalg.norm(importances, axis=0, ord=norm_order)

    if transform_func == "square":
        importances = safe_sqr(importances)

    return importances

