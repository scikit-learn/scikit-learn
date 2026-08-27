# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

import numpy as np

from sklearn.base import BaseEstimator, RegressorMixin, _fit_context, clone
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import FunctionTransformer
from sklearn.utils import _safe_indexing, check_array
from sklearn.utils._metadata_requests import (
    MetadataRouter,
    MethodMapping,
    _manual_routing,
    _routing_enabled,
    process_routing,
)
from sklearn.utils._param_validation import HasMethods, StrOptions
from sklearn.utils._repr_html.estimator import _VisualBlock
from sklearn.utils._tags import get_tags
from sklearn.utils.validation import check_is_fitted

__all__ = ["TransformedTargetRegressor"]


class TransformedTargetRegressor(RegressorMixin, BaseEstimator):
    """Meta-estimator to regress on a transformed target.

    Useful for applying a non-linear transformation to the target `y` in
    regression problems. This transformation can be given as a Transformer
    such as the :class:`~sklearn.preprocessing.QuantileTransformer` or as a
    function and its inverse such as `np.log` and `np.exp`.

    The computation during :meth:`fit` is::

        regressor.fit(X, func(y))

    or::

        regressor.fit(X, transformer.transform(y))

    The computation during :meth:`predict` is::

        inverse_func(regressor.predict(X))

    or::

        transformer.inverse_transform(regressor.predict(X))

    Read more in the :ref:`User Guide <transformed_target_regressor>`.

    .. versionadded:: 0.20

    Parameters
    ----------
    regressor : object, default=None
        Regressor object such as derived from
        :class:`~sklearn.base.RegressorMixin`. This regressor will
        automatically be cloned each time prior to fitting. If `regressor is
        None`, :class:`~sklearn.linear_model.LinearRegression` is created and used.

    transformer : object, default=None
        Estimator object such as derived from
        :class:`~sklearn.base.TransformerMixin`. Cannot be set at the same time
        as `func` and `inverse_func`. If `transformer is None` as well as
        `func` and `inverse_func`, the transformer will be an identity
        transformer. Note that the transformer will be cloned during fitting.
        Also, the transformer is restricting `y` to be a numpy array.

    func : function, default=None
        Function to apply to `y` before passing to :meth:`fit`. Cannot be set
        at the same time as `transformer`. If `func is None`, the function used will be
        the identity function. If `func` is set, `inverse_func` also needs to be
        provided. The function needs to return a 2-dimensional array.

    inverse_func : function, default=None
        Function to apply to the prediction of the regressor. Cannot be set at
        the same time as `transformer`. The inverse function is used to return
        predictions to the same space of the original training labels. If
        `inverse_func` is set, `func` also needs to be provided. The inverse
        function needs to return a 2-dimensional array.

    check_inverse : bool, default=True
        Whether to check that `transform` followed by `inverse_transform`
        or `func` followed by `inverse_func` leads to the original targets.

    bias_correction : {"additive", "multiplicative", "taylor"} or None, \
            default=None
        Strategy to correct the systematic prediction bias introduced by
        Jensen's inequality when using a nonlinear target transformation.
        Without correction, ``inverse_transform(E[Z|X])`` is a biased
        estimate of ``E[Y|X]``.

        - ``None``: no bias correction (default).
        - ``"additive"``: a global additive constant is computed during
          :meth:`fit` and added to every prediction.
        - ``"multiplicative"``: a global multiplicative factor is computed
          during :meth:`fit` and applied to every prediction.
        - ``"taylor"``: a per-sample correction based on a second-order
          Taylor expansion of the inverse transform. The second derivative
          is computed numerically. This requires a smooth inverse transform
          and will warn if the numerical Hessian is ill-behaved.

        .. versionadded:: 1.10

    Attributes
    ----------
    regressor_ : object
        Fitted regressor.

    transformer_ : object
        Transformer used in :meth:`fit` and :meth:`predict`.

    bias_correction_factor_ : float, ndarray, or None
        The correction factor computed during :meth:`fit`. For
        ``bias_correction="multiplicative"`` this is the multiplicative
        factor; for ``"additive"`` it is the additive constant.
        ``None`` when ``bias_correction`` is ``None`` or ``"taylor"``.

        .. versionadded:: 1.10

    residual_variance_ : float, ndarray, or None
        Variance of the residuals in the transformed target space,
        computed during :meth:`fit`. Used by ``bias_correction="taylor"``
        at prediction time. ``None`` when ``bias_correction`` is not
        ``"taylor"``.

        .. versionadded:: 1.10

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying regressor exposes such an attribute when fit.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    sklearn.preprocessing.FunctionTransformer : Construct a transformer from an
        arbitrary callable.

    Notes
    -----
    Internally, the target `y` is always converted into a 2-dimensional array
    to be used by scikit-learn transformers. At the time of prediction, the
    output will be reshaped to a have the same number of dimensions as `y`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.linear_model import LinearRegression
    >>> from sklearn.compose import TransformedTargetRegressor
    >>> tt = TransformedTargetRegressor(regressor=LinearRegression(),
    ...                                 func=np.log, inverse_func=np.exp)
    >>> X = np.arange(4).reshape(-1, 1)
    >>> y = np.exp(2 * X).ravel()
    >>> tt.fit(X, y)
    TransformedTargetRegressor(...)
    >>> tt.score(X, y)
    1.0
    >>> tt.regressor_.coef_
    array([2.])

    For a more detailed example use case refer to
    :ref:`sphx_glr_auto_examples_compose_plot_transformed_target.py`.
    """

    _parameter_constraints: dict = {
        "regressor": [HasMethods(["fit", "predict"]), None],
        "transformer": [HasMethods("transform"), None],
        "func": [callable, None],
        "inverse_func": [callable, None],
        "check_inverse": ["boolean"],
        "bias_correction": [
            StrOptions({"additive", "multiplicative", "taylor"}),
            None,
        ],
    }

    def __init__(
        self,
        regressor=None,
        *,
        transformer=None,
        func=None,
        inverse_func=None,
        check_inverse=True,
        bias_correction=None,
    ):
        self.regressor = regressor
        self.transformer = transformer
        self.func = func
        self.inverse_func = inverse_func
        self.check_inverse = check_inverse
        self.bias_correction = bias_correction

    def _fit_transformer(self, y):
        """Check transformer and fit transformer.

        Create the default transformer, fit it and make additional inverse
        check on a subset (optional).

        """
        if self.transformer is not None and (
            self.func is not None or self.inverse_func is not None
        ):
            raise ValueError(
                "'transformer' and functions 'func'/'inverse_func' cannot both be set."
            )
        elif self.transformer is not None:
            self.transformer_ = clone(self.transformer)
        else:
            if (self.func is not None and self.inverse_func is None) or (
                self.func is None and self.inverse_func is not None
            ):
                lacking_param, existing_param = (
                    ("func", "inverse_func")
                    if self.func is None
                    else ("inverse_func", "func")
                )
                raise ValueError(
                    f"When '{existing_param}' is provided, '{lacking_param}' must also"
                    f" be provided. If {lacking_param} is supposed to be the default,"
                    " you need to explicitly pass it the identity function."
                )
            self.transformer_ = FunctionTransformer(
                func=self.func,
                inverse_func=self.inverse_func,
                validate=True,
                check_inverse=self.check_inverse,
            )
            # We are transforming the target here and not the features, so we set the
            # output of FunctionTransformer() to be a numpy array (default) and to not
            # depend on the global configuration:
            self.transformer_.set_output(transform="default")
        # XXX: sample_weight is not currently passed to the
        # transformer. However, if transformer starts using sample_weight, the
        # code should be modified accordingly. At the time to consider the
        # sample_prop feature, it is also a good use case to be considered.
        self.transformer_.fit(y)
        if self.check_inverse:
            idx_selected = slice(None, None, max(1, y.shape[0] // 10))
            y_sel = _safe_indexing(y, idx_selected)
            y_sel_t = self.transformer_.transform(y_sel)
            if not np.allclose(y_sel, self.transformer_.inverse_transform(y_sel_t)):
                warnings.warn(
                    (
                        "The provided functions or transformer are"
                        " not strictly inverse of each other. If"
                        " you are sure you want to proceed regardless"
                        ", set 'check_inverse=False'"
                    ),
                    UserWarning,
                )

    @_fit_context(
        # TransformedTargetRegressor.regressor/transformer are not validated yet.
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y, **fit_params):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        **fit_params : dict
            - If `enable_metadata_routing=False` (default): Parameters directly passed
              to the `fit` method of the underlying regressor.

            - If `enable_metadata_routing=True`: Parameters safely routed to the `fit`
              method of the underlying regressor.

            .. versionchanged:: 1.6
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        if y is None:
            raise ValueError(
                f"This {self.__class__.__name__} estimator "
                "requires y to be passed, but the target y is None."
            )
        y = check_array(
            y,
            input_name="y",
            accept_sparse=False,
            ensure_all_finite=True,
            ensure_2d=False,
            dtype="numeric",
            allow_nd=True,
        )

        # store the number of dimension of the target to predict an array of
        # similar shape at predict
        self._training_dim = y.ndim

        # transformers are designed to modify X which is 2d dimensional, we
        # need to modify y accordingly.
        if y.ndim == 1:
            y_2d = y.reshape(-1, 1)
        else:
            y_2d = y
        self._fit_transformer(y_2d)

        # transform y and convert back to 1d array if needed
        y_trans = self.transformer_.transform(y_2d)
        # FIXME: a FunctionTransformer can return a 1D array even when validate
        # is set to True. Therefore, we need to check the number of dimension
        # first.
        if y_trans.ndim == 2 and y_trans.shape[1] == 1 and self._training_dim == 1:
            y_trans = y_trans.squeeze(axis=1)

        self.regressor_ = self._get_regressor(get_clone=True)
        if _routing_enabled():
            routed_params = process_routing(self, "fit", **fit_params)
        else:
            routed_params = _manual_routing({"regressor": {"fit": fit_params}})

        self.regressor_.fit(X, y_trans, **routed_params.regressor.fit)

        self._fit_bias_correction(X, y, y_trans)

        if hasattr(self.regressor_, "feature_names_in_"):
            self.feature_names_in_ = self.regressor_.feature_names_in_

        return self

    def _fit_bias_correction(self, X, y, y_trans):
        """Compute bias correction data needed at prediction time."""
        if self.bias_correction is None:
            self.bias_correction_factor_ = None
            self.residual_variance_ = None
            return

        z_pred_train = self.regressor_.predict(X)
        if z_pred_train.ndim == 1:
            z_pred_train_2d = z_pred_train.reshape(-1, 1)
        else:
            z_pred_train_2d = z_pred_train

        if self.bias_correction in ("multiplicative", "additive"):
            y_pred_train = self.transformer_.inverse_transform(z_pred_train_2d)

            if y.ndim == 1:
                y_for_mean = y.reshape(-1, 1)
            else:
                y_for_mean = y

            y_true_mean = np.mean(y_for_mean, axis=0)
            y_pred_mean = np.mean(y_pred_train, axis=0)

            if self.bias_correction == "multiplicative":
                denom_near_zero = np.abs(y_pred_mean) < 1e-10
                if np.any(denom_near_zero):
                    warnings.warn(
                        "Mean of inverse-transformed training predictions is"
                        " near zero. Multiplicative bias correction factor"
                        " set to 1.0 for affected outputs.",
                        UserWarning,
                    )
                factor = np.where(denom_near_zero, 1.0, y_true_mean / y_pred_mean)
                self.bias_correction_factor_ = (
                    factor.item() if factor.size == 1 else factor
                )
            else:
                factor = y_true_mean - y_pred_mean
                self.bias_correction_factor_ = (
                    factor.item() if factor.size == 1 else factor
                )

            self.residual_variance_ = None

        elif self.bias_correction == "taylor":
            if y_trans.ndim == 1:
                y_trans_for_var = y_trans.reshape(-1, 1)
            else:
                y_trans_for_var = y_trans

            residuals = y_trans_for_var - z_pred_train_2d
            var = np.var(residuals, axis=0)
            self.residual_variance_ = var.item() if var.size == 1 else var
            self.bias_correction_factor_ = None

    def _inverse_transform_array(self, z):
        """Apply inverse_transform on a 2D array and return a 2D result."""
        result = self.transformer_.inverse_transform(z)
        if result.ndim == 1:
            return result.reshape(-1, 1)
        return result

    def _apply_bias_correction(self, z_pred, y_pred):
        """Apply the fitted bias correction to inverse-transformed predictions.

        Parameters
        ----------
        z_pred : ndarray of shape (n_samples, n_outputs)
            Predictions in the transformed space (2D).
        y_pred : ndarray of shape (n_samples, n_outputs)
            Inverse-transformed predictions (2D).

        Returns
        -------
        y_corrected : ndarray of shape (n_samples, n_outputs)
            Bias-corrected predictions.
        """
        if self.bias_correction is None:
            return y_pred

        if self.bias_correction == "multiplicative":
            return y_pred * self.bias_correction_factor_

        if self.bias_correction == "additive":
            return y_pred + self.bias_correction_factor_

        # Taylor expansion: f_inv(z) + (sigma^2 / 2) * f_inv''(z)
        h = np.maximum(np.abs(z_pred) * 1e-4, 1e-4)
        f_plus = self._inverse_transform_array(z_pred + h)
        f_minus = self._inverse_transform_array(z_pred - h)
        hessian = (f_plus - 2 * y_pred + f_minus) / (h**2)

        non_finite_mask = ~np.isfinite(hessian)
        if np.any(non_finite_mask):
            warnings.warn(
                "Numerical second derivative of the inverse transform"
                " produced non-finite values. Bias correction is skipped"
                " for affected samples.",
                UserWarning,
            )
            hessian = np.where(non_finite_mask, 0.0, hessian)

        correction = (self.residual_variance_ / 2) * hessian
        return y_pred + correction

    def predict(self, X, **predict_params):
        """Predict using the base regressor, applying inverse.

        The regressor is used to predict and the `inverse_func` or
        `inverse_transform` is applied before returning the prediction.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Samples.

        **predict_params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters directly passed
              to the `predict` method of the underlying regressor.

            - If `enable_metadata_routing=True`: Parameters safely routed to the
              `predict` method of the underlying regressor.

            .. versionchanged:: 1.6
                See :ref:`Metadata Routing User Guide <metadata_routing>`
                for more details.

        Returns
        -------
        y_hat : ndarray of shape (n_samples,)
            Predicted values.
        """
        check_is_fitted(self)
        if _routing_enabled():
            routed_params = process_routing(self, "predict", **predict_params)
        else:
            routed_params = _manual_routing({"regressor": {"predict": predict_params}})

        pred = self.regressor_.predict(X, **routed_params.regressor.predict)
        if pred.ndim == 1:
            pred_2d = pred.reshape(-1, 1)
        else:
            pred_2d = pred

        pred_trans = self.transformer_.inverse_transform(pred_2d)
        if pred_trans.ndim == 1:
            pred_trans = pred_trans.reshape(-1, 1)

        pred_trans = self._apply_bias_correction(pred_2d, pred_trans)

        if (
            self._training_dim == 1
            and pred_trans.ndim == 2
            and pred_trans.shape[1] == 1
        ):
            pred_trans = pred_trans.squeeze(axis=1)

        return pred_trans

    def __sklearn_tags__(self):
        regressor = self._get_regressor()
        tags = super().__sklearn_tags__()
        tags.regressor_tags.poor_score = True
        tags.input_tags.sparse = get_tags(regressor).input_tags.sparse
        tags.target_tags.multi_output = get_tags(regressor).target_tags.multi_output
        return tags

    @property
    def n_features_in_(self):
        """Number of features seen during :term:`fit`."""
        # For consistency with other estimators we raise an AttributeError so
        # that hasattr() returns False the estimator isn't fitted.
        try:
            check_is_fitted(self)
        except NotFittedError as nfe:
            raise AttributeError(
                f"{self.__class__.__name__} object has no n_features_in_ attribute."
            ) from nfe

        return self.regressor_.n_features_in_

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.6

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = MetadataRouter(owner=self).add(
            regressor=self._get_regressor(),
            method_mapping=MethodMapping()
            .add(caller="fit", callee="fit")
            .add(caller="predict", callee="predict"),
        )
        return router

    def _get_regressor(self, get_clone=False):
        if self.regressor is None:
            if _routing_enabled():
                return LinearRegression().set_fit_request(sample_weight=True)
            else:
                return LinearRegression()

        return clone(self.regressor) if get_clone else self.regressor

    def _sk_visual_block_(self):
        regressor = (
            self.regressor_ if hasattr(self, "regressor_") else self._get_regressor()
        )
        return _VisualBlock(
            "serial",
            [regressor],
            names=[f"regressor: {regressor.__class__.__name__}"],
            name_details=[str(regressor)],
            dash_wrapped=True,
        )
