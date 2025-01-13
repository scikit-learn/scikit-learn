# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

import numpy as np

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin, _fit_context, clone
from ..exceptions import DataConversionWarning, NotFittedError
from ..linear_model import LinearRegression, LogisticRegression
from ..preprocessing import FunctionTransformer
from ..utils import Bunch, _safe_indexing, check_array
from ..utils._metadata_requests import (
    MetadataRouter,
    MethodMapping,
    _routing_enabled,
    process_routing,
)
from ..utils._param_validation import HasMethods, Hidden, StrOptions
from ..utils._tags import get_tags
from ..utils.metaestimators import available_if
from ..utils.validation import check_is_fitted

__all__ = ["TransformedTargetClassifier", "TransformedTargetRegressor"]


def _estimator_has(attr):
    """Check if we can delegate a method to the underlying estimator.

    First, we check the fitted `estimator_` if available, otherwise we check the
    unfitted `estimator`. We raise the original `AttributeError` if `attr` does
    not exist. This function is used together with `available_if`.
    """

    def check(self):
        if hasattr(self, "estimator_"):
            getattr(self.estimator_, attr)
        else:
            getattr(self.estimator, attr)

        return True

    return check


class BaseTransformedTarget(BaseEstimator):
    """Base class for transformed target meta-estimator.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """

    _parameter_constraints: dict = {
        "estimator": [HasMethods(["fit", "predict"]), None],
        "transformer": [HasMethods("transform"), None],
        "func": [callable, None],
        "inverse_func": [callable, None],
        "check_inverse": ["boolean"],
    }

    def __init__(
        self,
        estimator=None,
        *,
        transformer=None,
        func=None,
        inverse_func=None,
        check_inverse=True,
    ):
        self.estimator = estimator
        self.transformer = transformer
        self.func = func
        self.inverse_func = inverse_func
        self.check_inverse = check_inverse

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

            self._check_inverse(y_sel, y_sel_t)

    @_fit_context(
        # BaseTransformedTarget.estimator/transformer are not validated yet.
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y, **params):
        """Fit the model according to the given training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        y : array-like of shape (n_samples,)
            Target values.

        **params : dict
            - If `enable_metadata_routing=False` (default): Parameters directly passed
              to the `fit` method of the underlying estimator.

            - If `enable_metadata_routing=True`: Parameters safely routed to the `fit`
              method of the underlying estimator.

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

        y_2d = self._validate_y(y)

        self._fit_transformer(y_2d)

        # transform y and convert back to 1d array if needed
        y_trans = self.transformer_.transform(y_2d)
        # FIXME: a FunctionTransformer can return a 1D array even when validate
        # is set to True. Therefore, we need to check the number of dimension
        # first.
        if y_trans.ndim == 2 and y_trans.shape[1] == 1:
            y_trans = y_trans.squeeze(axis=1)

        self.estimator_ = self._get_estimator(get_clone=True)
        if _routing_enabled():
            routed_params = process_routing(self, "fit", **params)
        else:
            routed_params = Bunch(estimator=Bunch(fit=params))

        self.estimator_.fit(X, y_trans, **routed_params.estimator.fit)

        if hasattr(self.estimator_, "feature_names_in_"):
            self.feature_names_in_ = self.estimator_.feature_names_in_

        if hasattr(self.estimator_, "classes_"):
            self.classes_ = self.estimator_.classes_
        elif hasattr(self.transformer_, "classes_"):
            self.classes_ = self.transformer_.classes_

        return self

    def predict(self, X, **params):
        """Predict using the base estimator, applying inverse.

        The estimator is used to predict and the `inverse_func` or
        `inverse_transform` is applied before returning the prediction.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Samples.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters directly passed
              to the `predict` method of the underlying estimator.

            - If `enable_metadata_routing=True`: Parameters safely routed to the
              `predict` method of the underlying estimator.

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
            routed_params = process_routing(self, "predict", **params)
        else:
            routed_params = Bunch(estimator=Bunch(predict=params))

        pred = self.estimator_.predict(X, **routed_params.estimator.predict)
        if pred.ndim == 1:
            pred_trans = self.transformer_.inverse_transform(pred.reshape(-1, 1))
        else:
            pred_trans = self.transformer_.inverse_transform(pred)
        if (
            self._training_dim == 1
            and pred_trans.ndim == 2
            and pred_trans.shape[1] == 1
        ):
            pred_trans = pred_trans.squeeze(axis=1)

        return pred_trans

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()

        estimator = self._get_estimator()
        estimator_tags = get_tags(estimator)
        if estimator_tags.classifier_tags is not None:
            tags.classifier_tags.poor_score = True
        if estimator_tags.regressor_tags is not None:
            tags.regressor_tags.poor_score = True

        tags.input_tags.sparse = estimator_tags.input_tags.sparse
        tags.target_tags.multi_output = estimator_tags.target_tags.multi_output
        tags.target_tags.required = True
        return tags

    @property
    def n_features_in_(self):
        """Number of features seen during :term:`fit`."""
        # For consistency with other estimators we raise a AttributeError so
        # that hasattr() returns False the estimator isn't fitted.
        try:
            check_is_fitted(self)
        except NotFittedError as nfe:
            raise AttributeError(
                "{} object has no n_features_in_ attribute.".format(
                    self.__class__.__name__
                )
            ) from nfe

        return self.estimator_.n_features_in_


class TransformedTargetClassifier(ClassifierMixin, BaseTransformedTarget):
    """Meta-estimator to classify based on a transformed target.

    Useful for applying a transformation to the target `y` in
    classification problems. This transformation can be given as a Transformer
    such as the :class:`~sklearn.preprocessing.LabelEncoder` or as a
    function.

    The computation during :meth:`fit` is::

        classifier.fit(X, func(y))

    or::

        classifier.fit(X, transformer.transform(y))

    The computation during :meth:`predict` is::

        inverse_func(classifier.predict(X))

    or::

        transformer.inverse_transform(classifier.predict(X))

    Read more in the :ref:`User Guide <transformed_target_classifier>`.

    Parameters
    ----------
    estimator : object, default=None
        Classifier object such as derived from
        :class:`~sklearn.base.ClassifierMixin`. This classifier will
        automatically be cloned each time prior to fitting. If `estimator is
        None`, :class:`~sklearn.linear_model.LogisticRegression` is created and used.

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
        Function to apply to the prediction of the classifier. Cannot be set at
        the same time as `transformer`. The inverse function is used to return
        predictions to the same space of the original training labels. If
        `inverse_func` is set, `func` also needs to be provided. The inverse
        function needs to return a 2-dimensional array.

    check_inverse : bool, default=True
        Whether to check that `transform` followed by `inverse_transform`
        or `func` followed by `inverse_func` leads to the original targets.

    Attributes
    ----------
    estimator_ : object
        Fitted estimator.

    transformer_ : object
        Transformer used in :meth:`fit` and :meth:`predict`.

    classes_ : ndarray of shape (n_classes,)
        The class labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    TransformedTargetRegressor : Meta-estimator to regress on a
        transformed target.
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
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.compose import TransformedTargetClassifier
    >>> from sklearn.preprocessing import LabelEncoder
    >>> tt = TransformedTargetClassifier(estimator=LogisticRegression(),
    ...                                 transformer=LabelEncoder())
    >>> X = np.arange(4).reshape(-1, 1)
    >>> y = np.array["c_1", "c_1", "c_2", "c_2"]
    >>> tt.fit(X, y)
    TransformedTargetClassifier(...)
    >>> tt.score(X, y)
    1.0
    >>> tt.estimator_.coef_
    array([[0.95826546]])
    """

    def _validate_y(self, y):
        y = check_array(
            y,
            input_name="y",
            accept_sparse=False,
            ensure_all_finite=True,
            ensure_2d=False,
            dtype=None,
            allow_nd=True,
        )

        # store the number of dimension of the target to predict an array of
        # similar shape at predict
        self._training_dim = y.ndim

        # transformers are designed to modify X which is 2d dimensional,
        # but not label transformers as they modify y which is 1d
        # we check the input tags and modify y accordingly

        requires_2d_input = get_tags(self.transformer).input_tags.two_d_array
        if requires_2d_input and y.ndim == 1:
            y_2d = y.reshape(-1, 1)
        else:
            if not get_tags(self._get_estimator()).target_tags.multi_output:
                warnings.warn(
                    (
                        "A column-vector y was passed when a 1d array was"
                        " expected. Please change the shape of y to "
                        "(n_samples,), for example using ravel()."
                    ),
                    DataConversionWarning,
                    stacklevel=2,
                )
            y_2d = y

        return y_2d

    def _check_inverse(self, y, y_t):
        if not np.array_equal(y, self.transformer_.inverse_transform(y_t)):
            warnings.warn(
                (
                    "The provided functions or transformer are"
                    " not strictly inverse of each other. If"
                    " you are sure you want to proceed regardless"
                    ", set 'check_inverse=False'"
                ),
                UserWarning,
            )

    def _get_estimator(self, get_clone=False):
        if self.estimator is None:
            return LogisticRegression()

        return clone(self.estimator) if get_clone else self.estimator

    @available_if(_estimator_has("predict_proba"))
    def predict_proba(self, X, **params):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters directly passed
              to the `predict` method of the underlying estimator.

            - If `enable_metadata_routing=True`: Parameters safely routed to the
              `predict` method of the underlying estimator.

        Returns
        -------
        p : ndarray of shape (n_samples, n_classes)
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        check_is_fitted(self)
        if _routing_enabled():
            routed_params = process_routing(self, "predict_proba", **params)
        else:
            routed_params = Bunch(estimator=Bunch(predict_proba=params))

        return self.estimator_.predict_proba(X, **routed_params.estimator.predict_proba)

    @available_if(_estimator_has("predict_log_proba"))
    def predict_log_proba(self, X, **params):
        """Predict class log-probabilities for X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters directly passed
              to the `predict` method of the underlying estimator.

            - If `enable_metadata_routing=True`: Parameters safely routed to the
              `predict` method of the underlying estimator.

        Returns
        -------
        p : ndarray of shape (n_samples, n_classes)
            The class log-probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        check_is_fitted(self)
        if _routing_enabled():
            routed_params = process_routing(self, "predict_log_proba", **params)
        else:
            routed_params = Bunch(estimator=Bunch(predict_log_proba=params))

        return self.estimator_.predict_log_proba(
            X, **routed_params.estimator.predict_log_proba
        )

    @available_if(_estimator_has("decision_function"))
    def decision_function(self, X, **params):
        """Average of the decision functions of the base classifiers.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters directly passed
              to the `predict` method of the underlying estimator.

            - If `enable_metadata_routing=True`: Parameters safely routed to the
              `predict` method of the underlying estimator.

        Returns
        -------
        score : ndarray of shape (n_samples, k)
            The decision function of the input samples. The columns correspond
            to the classes in sorted order, as they appear in the attribute
            ``classes_``. Regression and binary classification are special
            cases with ``k == 1``, otherwise ``k==n_classes``.
        """
        check_is_fitted(self)
        if _routing_enabled():
            routed_params = process_routing(self, "decision_function", **params)
        else:
            routed_params = Bunch(estimator=Bunch(decision_function=params))

        return self.estimator_.decision_function(
            X, **routed_params.estimator.decision_function
        )

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = MetadataRouter(owner=self.__class__.__name__).add(
            estimator=self._get_estimator(),
            method_mapping=MethodMapping()
            .add(caller="fit", callee="fit")
            .add(caller="predict", callee="predict")
            .add(callee="predict_proba", caller="predict_proba")
            .add(callee="decision_function", caller="decision_function")
            .add(callee="predict_log_proba", caller="predict_log_proba"),
        )
        return router


class TransformedTargetRegressor(RegressorMixin, BaseTransformedTarget):
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
    estimator : object, default=None
        Regressor object such as derived from
        :class:`~sklearn.base.RegressorMixin`. This regressor will
        automatically be cloned each time prior to fitting. If `estimator is
        None`, :class:`~sklearn.linear_model.LinearRegression` is created and used.

        .. versionadded:: 1.6
           `regressor` was renamed to `estimator`.

    regressor : object, default=None
        Regressor object such as derived from
        :class:`~sklearn.base.RegressorMixin`. This regressor will
        automatically be cloned each time prior to fitting. If `estimator is
        None`, :class:`~sklearn.linear_model.LinearRegression` is created and used.

        .. deprecated:: 1.6
            `regressor` was deprecated in 1.6 and will be removed in 1.8.
            Use `estimator` instead.

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

    Attributes
    ----------
    estimator_ : object
        Fitted estimator.

    transformer_ : object
        Transformer used in :meth:`fit` and :meth:`predict`.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    TransformedTargetRegressor : Meta-estimator to classify
        based on a transformed target.
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
    >>> tt = TransformedTargetRegressor(estimator=LinearRegression(),
    ...                                 func=np.log, inverse_func=np.exp)
    >>> X = np.arange(4).reshape(-1, 1)
    >>> y = np.exp(2 * X).ravel()
    >>> tt.fit(X, y)
    TransformedTargetRegressor(...)
    >>> tt.score(X, y)
    1.0
    >>> tt.estimator_.coef_
    array([2.])

    For a more detailed example use case refer to
    :ref:`sphx_glr_auto_examples_compose_plot_transformed_target.py`.
    """

    _parameter_constraints: dict = {
        **BaseTransformedTarget._parameter_constraints,
        "regressor": [
            None,
            HasMethods(["fit", "predict"]),
            Hidden(StrOptions({"deprecated"})),
        ],
    }

    # TODO(1.8) remove
    def __init__(
        self,
        estimator=None,
        regressor="deprecated",
        *,
        transformer=None,
        func=None,
        inverse_func=None,
        check_inverse=True,
    ):
        super().__init__(
            estimator=estimator,
            transformer=transformer,
            func=func,
            inverse_func=inverse_func,
            check_inverse=check_inverse,
        )

        self.regressor = regressor

    def _validate_y(self, y):
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

        return y_2d

    def _check_inverse(self, y, y_t):
        if not np.allclose(y, self.transformer_.inverse_transform(y_t)):
            warnings.warn(
                (
                    "The provided functions or transformer are"
                    " not strictly inverse of each other. If"
                    " you are sure you want to proceed regardless"
                    ", set 'check_inverse=False'"
                ),
                UserWarning,
            )

    def _get_estimator(self, get_clone=False):
        # TODO(1.8): remove
        estimator_ = self.estimator
        if self.estimator is None and self.regressor != "deprecated":
            estimator_ = self.regressor

            warnings.warn(
                (
                    "`regressor` has been deprecated in 1.6 and will be removed"
                    " in 1.8. Please use `estimator` instead."
                ),
                FutureWarning,
            )
        # TODO(1.8) remove
        elif self.estimator is not None and self.regressor != "deprecated":
            raise ValueError(
                "You must pass only one estimator to TransformedTargetRegressor."
                " Use `estimator`."
            )
        # TODO(1.8) remove
        elif self.estimator is None and self.regressor == "deprecated":
            estimator_ = None

        # TODO(1.8) replace estimator_ by self.estimator in remaining code
        if estimator_ is None:
            return LinearRegression()

        return clone(estimator_) if get_clone else estimator_

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
        router = MetadataRouter(owner=self.__class__.__name__).add(
            estimator=self._get_estimator(),
            method_mapping=MethodMapping()
            .add(caller="fit", callee="fit")
            .add(caller="predict", callee="predict"),
        )
        return router
