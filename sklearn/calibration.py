"""Methods for calibrating predicted probabilities."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from functools import partial
from inspect import signature
from math import log
from numbers import Integral, Real

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit

from sklearn._loss import HalfBinomialLoss, HalfMultinomialLoss
from sklearn.base import (
    BaseEstimator,
    ClassifierMixin,
    MetaEstimatorMixin,
    RegressorMixin,
    _fit_context,
    clone,
)
from sklearn.externals import array_api_extra as xpx
from sklearn.frozen import FrozenEstimator
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import LeaveOneOut, check_cv, cross_val_predict
from sklearn.preprocessing import LabelEncoder, label_binarize
from sklearn.svm import LinearSVC
from sklearn.utils import Bunch, _safe_indexing, column_or_1d, get_tags, indexable
from sklearn.utils._array_api import (
    _convert_to_numpy,
    _half_multinomial_loss,
    _is_numpy_namespace,
    get_namespace,
    get_namespace_and_device,
    move_to,
)
from sklearn.utils._param_validation import (
    HasMethods,
    Interval,
    StrOptions,
    validate_params,
)
from sklearn.utils._plotting import (
    _BinaryClassifierCurveDisplayMixin,
    _validate_style_kwargs,
)
from sklearn.utils._response import _get_response_values, _process_predict_proba
from sklearn.utils.extmath import softmax
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _routing_enabled,
    process_routing,
)
from sklearn.utils.multiclass import check_classification_targets
from sklearn.utils.parallel import Parallel, delayed
from sklearn.utils.validation import (
    _check_method_params,
    _check_pos_label_consistency,
    _check_response_method,
    _check_sample_weight,
    _num_samples,
    check_array,
    check_consistent_length,
    check_is_fitted,
)


class CalibratedClassifierCV(ClassifierMixin, MetaEstimatorMixin, BaseEstimator):
    """Calibrate probabilities using isotonic, sigmoid, or temperature scaling.

    This class uses cross-validation to both estimate the parameters of a
    classifier and subsequently calibrate a classifier. With
    `ensemble=True`, for each cv split it
    fits a copy of the base estimator to the training subset, and calibrates it
    using the testing subset. For prediction, predicted probabilities are
    averaged across these individual calibrated classifiers. When
    `ensemble=False`, cross-validation is used to obtain unbiased predictions,
    via :func:`~sklearn.model_selection.cross_val_predict`, which are then
    used for calibration. For prediction, the base estimator, trained using all
    the data, is used. This is the prediction method implemented when
    `probabilities=True` for :class:`~sklearn.svm.SVC` and :class:`~sklearn.svm.NuSVC`
    estimators (see :ref:`User Guide <scores_probabilities>` for details).

    Already fitted classifiers can be calibrated by wrapping the model in a
    :class:`~sklearn.frozen.FrozenEstimator`. In this case all provided
    data is used for calibration. The user has to take care manually that data
    for model fitting and calibration are disjoint.

    The calibration is based on the :term:`decision_function` method of the
    `estimator` if it exists, else on :term:`predict_proba`.

    Read more in the :ref:`User Guide <calibration>`.
    In order to learn more on the CalibratedClassifierCV class, see the
    following calibration examples:
    :ref:`sphx_glr_auto_examples_calibration_plot_calibration.py`,
    :ref:`sphx_glr_auto_examples_calibration_plot_calibration_curve.py`, and
    :ref:`sphx_glr_auto_examples_calibration_plot_calibration_multiclass.py`.

    Parameters
    ----------
    estimator : estimator instance, default=None
        The classifier whose output need to be calibrated to provide more
        accurate `predict_proba` outputs. The default classifier is
        a :class:`~sklearn.svm.LinearSVC`.

        .. versionadded:: 1.2

    method : {'sigmoid', 'isotonic', 'temperature'}, default='sigmoid'
        The method to use for calibration. Can be:

        - 'sigmoid', which corresponds to Platt's method (i.e. a binary logistic
          regression model).
        - 'isotonic', which is a non-parametric approach.
        - 'temperature', temperature scaling.

        Sigmoid and isotonic calibration methods natively support only binary
        classifiers and extend to multi-class classification using a One-vs-Rest (OvR)
        strategy with post-hoc renormalization, i.e., adjusting the probabilities after
        calibration to ensure they sum up to 1.

        In contrast, temperature scaling naturally supports multi-class calibration by
        applying `softmax(classifier_logits/T)` with a value of `T` (temperature)
        that optimizes the log loss.

        For very uncalibrated classifiers on very imbalanced datasets, sigmoid
        calibration might be preferred because it fits an additional intercept
        parameter. This helps shift decision boundaries appropriately when the
        classifier being calibrated is biased towards the majority class.

        Isotonic calibration is not recommended when the number of calibration samples
        is too low ``(≪1000)`` since it then tends to overfit.

        .. versionchanged:: 1.8
           Added option 'temperature'.

    cv : int, cross-validation generator, or iterable, default=None
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`~sklearn.model_selection.StratifiedKFold` is used. If ``y`` is
        neither binary nor multiclass, :class:`~sklearn.model_selection.KFold`
        is used.

        Refer to the :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    n_jobs : int, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors.

        Base estimator clones are fitted in parallel across cross-validation
        iterations.

        See :term:`Glossary <n_jobs>` for more details.

        .. versionadded:: 0.24

    ensemble : bool, or "auto", default="auto"
        Determines how the calibrator is fitted.

        "auto" will use `False` if the `estimator` is a
        :class:`~sklearn.frozen.FrozenEstimator`, and `True` otherwise.

        If `True`, the `estimator` is fitted using training data, and
        calibrated using testing data, for each `cv` fold. The final estimator
        is an ensemble of `n_cv` fitted classifier and calibrator pairs, where
        `n_cv` is the number of cross-validation folds. The output is the
        average predicted probabilities of all pairs.

        If `False`, `cv` is used to compute unbiased predictions, via
        :func:`~sklearn.model_selection.cross_val_predict`, which are then
        used for calibration. At prediction time, the classifier used is the
        `estimator` trained on all the data.
        Note that this method is also internally implemented  in
        :mod:`sklearn.svm` estimators with the `probabilities=True` parameter.

        .. versionadded:: 0.24

        .. versionchanged:: 1.6
            `"auto"` option is added and is the default.

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes,)
        The class labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

        .. versionadded:: 1.0

    calibrated_classifiers_ : list (len() equal to cv or 1 if `ensemble=False`)
        The list of classifier and calibrator pairs.

        - When `ensemble=True`, `n_cv` fitted `estimator` and calibrator pairs.
          `n_cv` is the number of cross-validation folds.
        - When `ensemble=False`, the `estimator`, fitted on all the data, and fitted
          calibrator.

        .. versionchanged:: 0.24
            Single calibrated classifier case when `ensemble=False`.

    See Also
    --------
    calibration_curve : Compute true and predicted probabilities
        for a calibration curve.

    References
    ----------
    .. [1] B. Zadrozny & C. Elkan.
       `Obtaining calibrated probability estimates from decision trees
       and naive Bayesian classifiers
       <https://cseweb.ucsd.edu/~elkan/calibrated.pdf>`_, ICML 2001.

    .. [2] B. Zadrozny & C. Elkan.
       `Transforming Classifier Scores into Accurate Multiclass
       Probability Estimates
       <https://web.archive.org/web/20060720141520id_/http://www.research.ibm.com:80/people/z/zadrozny/kdd2002-Transf.pdf>`_,
       KDD 2002.

    .. [3] J. Platt. `Probabilistic Outputs for Support Vector Machines
       and Comparisons to Regularized Likelihood Methods
       <https://www.researchgate.net/profile/John-Platt-2/publication/2594015_Probabilistic_Outputs_for_Support_Vector_Machines_and_Comparisons_to_Regularized_Likelihood_Methods/links/004635154cff5262d6000000/Probabilistic-Outputs-for-Support-Vector-Machines-and-Comparisons-to-Regularized-Likelihood-Methods.pdf>`_,
       1999.

    .. [4] A. Niculescu-Mizil & R. Caruana.
       `Predicting Good Probabilities with Supervised Learning
       <https://www.cs.cornell.edu/~alexn/papers/calibration.icml05.crc.rev3.pdf>`_,
       ICML 2005.

    .. [5] Chuan Guo, Geoff Pleiss, Yu Sun, Kilian Q. Weinberger.
       :doi:`On Calibration of Modern Neural Networks<10.48550/arXiv.1706.04599>`.
       Proceedings of the 34th International Conference on Machine Learning,
       PMLR 70:1321-1330, 2017.

    Examples
    --------
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.naive_bayes import GaussianNB
    >>> from sklearn.calibration import CalibratedClassifierCV
    >>> X, y = make_classification(n_samples=100, n_features=2,
    ...                            n_redundant=0, random_state=42)
    >>> base_clf = GaussianNB()
    >>> calibrated_clf = CalibratedClassifierCV(base_clf, cv=3)
    >>> calibrated_clf.fit(X, y)
    CalibratedClassifierCV(...)
    >>> len(calibrated_clf.calibrated_classifiers_)
    3
    >>> calibrated_clf.predict_proba(X)[:5, :]
    array([[0.110, 0.889],
           [0.072, 0.927],
           [0.928, 0.072],
           [0.928, 0.072],
           [0.072, 0.928]])
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = make_classification(n_samples=100, n_features=2,
    ...                            n_redundant=0, random_state=42)
    >>> X_train, X_calib, y_train, y_calib = train_test_split(
    ...        X, y, random_state=42
    ... )
    >>> base_clf = GaussianNB()
    >>> base_clf.fit(X_train, y_train)
    GaussianNB()
    >>> from sklearn.frozen import FrozenEstimator
    >>> calibrated_clf = CalibratedClassifierCV(FrozenEstimator(base_clf))
    >>> calibrated_clf.fit(X_calib, y_calib)
    CalibratedClassifierCV(...)
    >>> len(calibrated_clf.calibrated_classifiers_)
    1
    >>> calibrated_clf.predict_proba([[-0.5, 0.5]])
    array([[0.936, 0.063]])
    """

    _parameter_constraints: dict = {
        "estimator": [
            HasMethods(["fit", "predict_proba"]),
            HasMethods(["fit", "decision_function"]),
            None,
        ],
        "method": [StrOptions({"isotonic", "sigmoid", "temperature"})],
        "cv": ["cv_object"],
        "n_jobs": [Integral, None],
        "ensemble": ["boolean", StrOptions({"auto"})],
    }

    def __init__(
        self,
        estimator=None,
        *,
        method="sigmoid",
        cv=None,
        n_jobs=None,
        ensemble="auto",
    ):
        self.estimator = estimator
        self.method = method
        self.cv = cv
        self.n_jobs = n_jobs
        self.ensemble = ensemble

    def _get_estimator(self):
        """Resolve which estimator to return (default is LinearSVC)"""
        if self.estimator is None:
            # we want all classifiers that don't expose a random_state
            # to be deterministic (and we don't want to expose this one).
            estimator = LinearSVC(random_state=0)
            if _routing_enabled():
                estimator.set_fit_request(sample_weight=True)
        else:
            estimator = self.estimator

        return estimator

    @_fit_context(
        # CalibratedClassifierCV.estimator is not validated yet
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y, sample_weight=None, **fit_params):
        """Fit the calibrated model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.

        **fit_params : dict
            Parameters to pass to the `fit` method of the underlying
            classifier.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        check_classification_targets(y)
        X, y = indexable(X, y)
        estimator = self._get_estimator()

        _ensemble = self.ensemble
        if _ensemble == "auto":
            _ensemble = not isinstance(estimator, FrozenEstimator)

        self.calibrated_classifiers_ = []

        # Set `classes_` using all `y`
        label_encoder_ = LabelEncoder().fit(y)
        self.classes_ = label_encoder_.classes_
        if self.method == "temperature" and isinstance(y[0], str):
            # for temperature scaling if `y` contains strings then encode it
            # right here to avoid fitting LabelEncoder again within the
            # `_fit_calibrator` function.
            y = label_encoder_.transform(y=y)

        if _routing_enabled():
            routed_params = process_routing(
                self,
                "fit",
                sample_weight=sample_weight,
                **fit_params,
            )
        else:
            # sample_weight checks
            fit_parameters = signature(estimator.fit).parameters
            supports_sw = "sample_weight" in fit_parameters
            if sample_weight is not None and not supports_sw:
                estimator_name = type(estimator).__name__
                warnings.warn(
                    f"Since {estimator_name} does not appear to accept"
                    " sample_weight, sample weights will only be used for the"
                    " calibration itself. This can be caused by a limitation of"
                    " the current scikit-learn API. See the following issue for"
                    " more details:"
                    " https://github.com/scikit-learn/scikit-learn/issues/21134."
                    " Be warned that the result of the calibration is likely to be"
                    " incorrect."
                )
            routed_params = Bunch()
            routed_params.splitter = Bunch(split={})  # no routing for splitter
            routed_params.estimator = Bunch(fit=fit_params)
            if sample_weight is not None and supports_sw:
                routed_params.estimator.fit["sample_weight"] = sample_weight

        xp, is_array_api, device_ = get_namespace_and_device(X)
        if is_array_api:
            y, sample_weight = move_to(y, sample_weight, xp=xp, device=device_)
        # Check that each cross-validation fold can have at least one
        # example per class
        if isinstance(self.cv, int):
            n_folds = self.cv
        elif hasattr(self.cv, "n_splits"):
            n_folds = self.cv.n_splits
        else:
            n_folds = None
        if n_folds and xp.any(xp.unique_counts(y)[1] < n_folds):
            raise ValueError(
                f"Requesting {n_folds}-fold "
                "cross-validation but provided less than "
                f"{n_folds} examples for at least one class."
            )
        if isinstance(self.cv, LeaveOneOut):
            raise ValueError(
                "LeaveOneOut cross-validation does not allow"
                "all classes to be present in test splits. "
                "Please use a cross-validation generator that allows "
                "all classes to appear in every test and train split."
            )
        cv = check_cv(self.cv, y, classifier=True)

        if _ensemble:
            parallel = Parallel(n_jobs=self.n_jobs)
            self.calibrated_classifiers_ = parallel(
                delayed(_fit_classifier_calibrator_pair)(
                    clone(estimator),
                    X,
                    y,
                    train=train,
                    test=test,
                    method=self.method,
                    classes=self.classes_,
                    xp=xp,
                    sample_weight=sample_weight,
                    fit_params=routed_params.estimator.fit,
                )
                for train, test in cv.split(X, y, **routed_params.splitter.split)
            )
        else:
            this_estimator = clone(estimator)
            method_name = _check_response_method(
                this_estimator,
                ["decision_function", "predict_proba"],
            ).__name__
            predictions = cross_val_predict(
                estimator=this_estimator,
                X=X,
                y=y,
                cv=cv,
                method=method_name,
                n_jobs=self.n_jobs,
                params=routed_params.estimator.fit,
            )
            if self.classes_.shape[0] == 2:
                # Ensure shape (n_samples, 1) in the binary case
                if method_name == "predict_proba":
                    # Select the probability column of the positive class
                    predictions = _process_predict_proba(
                        y_pred=predictions,
                        target_type="binary",
                        classes=self.classes_,
                        pos_label=self.classes_[1],
                    )
                predictions = predictions.reshape(-1, 1)

            if sample_weight is not None:
                # Check that the sample_weight dtype is consistent with the
                # predictions to avoid unintentional upcasts.
                sample_weight = _check_sample_weight(
                    sample_weight, predictions, dtype=predictions.dtype
                )

            this_estimator.fit(X, y, **routed_params.estimator.fit)
            # Note: Here we don't pass on fit_params because the supported
            # calibrators don't support fit_params anyway
            calibrated_classifier = _fit_calibrator(
                this_estimator,
                predictions,
                y,
                self.classes_,
                self.method,
                xp=xp,
                sample_weight=sample_weight,
            )
            self.calibrated_classifiers_.append(calibrated_classifier)

        first_clf = self.calibrated_classifiers_[0].estimator
        if hasattr(first_clf, "n_features_in_"):
            self.n_features_in_ = first_clf.n_features_in_
        if hasattr(first_clf, "feature_names_in_"):
            self.feature_names_in_ = first_clf.feature_names_in_
        return self

    def predict_proba(self, X):
        """Calibrated probabilities of classification.

        This function returns calibrated probabilities of classification
        according to each class on an array of test vectors X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The samples, as accepted by `estimator.predict_proba`.

        Returns
        -------
        C : ndarray of shape (n_samples, n_classes)
            The predicted probas.
        """
        check_is_fitted(self)
        # Compute the arithmetic mean of the predictions of the calibrated
        # classifiers
        xp, _, device_ = get_namespace_and_device(X)
        mean_proba = xp.zeros((_num_samples(X), self.classes_.shape[0]), device=device_)
        for calibrated_classifier in self.calibrated_classifiers_:
            proba = calibrated_classifier.predict_proba(X)
            mean_proba += proba

        mean_proba /= len(self.calibrated_classifiers_)

        return mean_proba

    def predict(self, X):
        """Predict the target of new samples.

        The predicted class is the class that has the highest probability,
        and can thus be different from the prediction of the uncalibrated classifier.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The samples, as accepted by `estimator.predict`.

        Returns
        -------
        C : ndarray of shape (n_samples,)
            The predicted class.
        """
        xp, _ = get_namespace(X)
        check_is_fitted(self)
        class_indices = xp.argmax(self.predict_proba(X), axis=1)
        if isinstance(self.classes_[0], str):
            class_indices = _convert_to_numpy(class_indices, xp=xp)

        return self.classes_[class_indices]

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
        router = (
            MetadataRouter(owner=self)
            .add_self_request(self)
            .add(
                estimator=self._get_estimator(),
                method_mapping=MethodMapping().add(caller="fit", callee="fit"),
            )
            .add(
                splitter=self.cv,
                method_mapping=MethodMapping().add(caller="fit", callee="split"),
            )
        )
        return router

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        estimator_tags = get_tags(self._get_estimator())
        tags.input_tags.sparse = estimator_tags.input_tags.sparse
        tags.array_api_support = (
            estimator_tags.array_api_support and self.method == "temperature"
        )
        return tags


def _fit_classifier_calibrator_pair(
    estimator,
    X,
    y,
    train,
    test,
    method,
    classes,
    xp,
    sample_weight=None,
    fit_params=None,
):
    """Fit a classifier/calibration pair on a given train/test split.

    Fit the classifier on the train set, compute its predictions on the test
    set and use the predictions as input to fit the calibrator along with the
    test labels.

    Parameters
    ----------
    estimator : estimator instance
        Cloned base estimator.

    X : array-like, shape (n_samples, n_features)
        Sample data.

    y : array-like, shape (n_samples,)
        Targets.

    train : ndarray, shape (n_train_indices,)
        Indices of the training subset.

    test : ndarray, shape (n_test_indices,)
        Indices of the testing subset.

    method : {'sigmoid', 'isotonic', 'temperature'}
        Method to use for calibration.

    classes : ndarray, shape (n_classes,)
        The target classes.

    xp : namespace
        Array API namespace.

    sample_weight : array-like, default=None
        Sample weights for `X`.

    fit_params : dict, default=None
        Parameters to pass to the `fit` method of the underlying
        classifier.

    Returns
    -------
    calibrated_classifier : _CalibratedClassifier instance
    """
    fit_params_train = _check_method_params(X, params=fit_params, indices=train)
    X_train, y_train = _safe_indexing(X, train), _safe_indexing(y, train)
    X_test, y_test = _safe_indexing(X, test), _safe_indexing(y, test)

    estimator.fit(X_train, y_train, **fit_params_train)

    predictions, _ = _get_response_values(
        estimator,
        X_test,
        response_method=["decision_function", "predict_proba"],
    )
    if predictions.ndim == 1:
        # Reshape binary output from `(n_samples,)` to `(n_samples, 1)`
        predictions = predictions.reshape(-1, 1)

    if sample_weight is not None:
        # Check that the sample_weight dtype is consistent with the predictions
        # to avoid unintentional upcasts.
        sample_weight = _check_sample_weight(sample_weight, X, dtype=predictions.dtype)
        sw_test = _safe_indexing(sample_weight, test)
    else:
        sw_test = None
    calibrated_classifier = _fit_calibrator(
        estimator,
        predictions,
        y_test,
        classes,
        method,
        xp=xp,
        sample_weight=sw_test,
    )
    return calibrated_classifier


def _fit_calibrator(clf, predictions, y, classes, method, xp, sample_weight=None):
    """Fit calibrator(s) and return a `_CalibratedClassifier`
    instance.

    A separate calibrator is fitted for each of the `n_classes`
    (i.e. `len(clf.classes_)`). However, if `n_classes` is 2 or if
    `method` is 'temperature', only one calibrator is fitted.

    Parameters
    ----------
    clf : estimator instance
        Fitted classifier.

    predictions : array-like, shape (n_samples, n_classes) or (n_samples, 1) \
                    when binary.
        Raw predictions returned by the un-calibrated base classifier.

    y : array-like, shape (n_samples,)
        The targets. For `method="temperature"`, `y` needs to be label encoded.

    classes : ndarray, shape (n_classes,)
        All the prediction classes.

    method : {'sigmoid', 'isotonic', 'temperature'}
        The method to use for calibration.

    xp : namespace
        Array API namespace.

    sample_weight : ndarray, shape (n_samples,), default=None
        Sample weights. If None, then samples are equally weighted.

    Returns
    -------
    pipeline : _CalibratedClassifier instance
    """
    calibrators = []

    if method in ("isotonic", "sigmoid"):
        Y = label_binarize(y, classes=classes)
        label_encoder = LabelEncoder().fit(classes)
        pos_class_indices = label_encoder.transform(clf.classes_)
        for class_idx, this_pred in zip(pos_class_indices, predictions.T):
            if method == "isotonic":
                calibrator = IsotonicRegression(out_of_bounds="clip")
            else:  # "sigmoid"
                calibrator = _SigmoidCalibration()
            calibrator.fit(this_pred, Y[:, class_idx], sample_weight)
            calibrators.append(calibrator)
    elif method == "temperature":
        if classes.shape[0] == 2 and predictions.shape[-1] == 1:
            response_method_name = _check_response_method(
                clf,
                ["decision_function", "predict_proba"],
            ).__name__
            if response_method_name == "predict_proba":
                predictions = xp.concat([1 - predictions, predictions], axis=1)
        calibrator = _TemperatureScaling()
        calibrator.fit(predictions, y, sample_weight)
        calibrators.append(calibrator)

    pipeline = _CalibratedClassifier(clf, calibrators, method=method, classes=classes)
    return pipeline


class _CalibratedClassifier:
    """Pipeline-like chaining a fitted classifier and its fitted calibrators.

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier.

    calibrators : list of fitted estimator instances
        List of fitted calibrators (either 'IsotonicRegression' or
        '_SigmoidCalibration'). The number of calibrators equals the number of
        classes. However, if there are 2 classes, the list contains only one
        fitted calibrator.

    classes : array-like of shape (n_classes,)
        All the prediction classes.

    method : {'sigmoid', 'isotonic'}, default='sigmoid'
        The method to use for calibration. Can be 'sigmoid' which
        corresponds to Platt's method or 'isotonic' which is a
        non-parametric approach based on isotonic regression.
    """

    def __init__(self, estimator, calibrators, *, classes, method="sigmoid"):
        self.estimator = estimator
        self.calibrators = calibrators
        self.classes = classes
        self.method = method

    def predict_proba(self, X):
        """Calculate calibrated probabilities.

        Calculates classification calibrated probabilities
        for each class, in a one-vs-all manner, for `X`.

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_features)
            The sample data.

        Returns
        -------
        proba : array, shape (n_samples, n_classes)
            The predicted probabilities. Can be exact zeros.
        """
        predictions, _ = _get_response_values(
            self.estimator,
            X,
            response_method=["decision_function", "predict_proba"],
        )
        if predictions.ndim == 1:
            # Reshape binary output from `(n_samples,)` to `(n_samples, 1)`
            predictions = predictions.reshape(-1, 1)

        n_classes = self.classes.shape[0]

        proba = np.zeros((_num_samples(X), n_classes))

        if self.method in ("sigmoid", "isotonic"):
            label_encoder = LabelEncoder().fit(self.classes)
            pos_class_indices = label_encoder.transform(self.estimator.classes_)
            for class_idx, this_pred, calibrator in zip(
                pos_class_indices, predictions.T, self.calibrators
            ):
                if n_classes == 2:
                    # When binary, `predictions` consists only of predictions for
                    # clf.classes_[1] but `pos_class_indices` = 0
                    class_idx += 1
                proba[:, class_idx] = calibrator.predict(this_pred)
            # Normalize the probabilities
            if n_classes == 2:
                proba[:, 0] = 1.0 - proba[:, 1]
            else:
                denominator = np.sum(proba, axis=1)[:, np.newaxis]
                # In the edge case where for each class calibrator returns a zero
                # probability for a given sample, use the uniform distribution
                # instead.
                uniform_proba = np.full_like(proba, 1 / n_classes)
                proba = np.divide(
                    proba, denominator, out=uniform_proba, where=denominator != 0
                )
        elif self.method == "temperature":
            xp, _ = get_namespace(predictions)
            if n_classes == 2 and predictions.shape[-1] == 1:
                response_method_name = _check_response_method(
                    self.estimator,
                    ["decision_function", "predict_proba"],
                ).__name__
                if response_method_name == "predict_proba":
                    predictions = xp.concat([1 - predictions, predictions], axis=1)
            proba = self.calibrators[0].predict(predictions)

        # Deal with cases where the predicted probability minimally exceeds 1.0
        proba[(1.0 < proba) & (proba <= 1.0 + 1e-5)] = 1.0

        return proba


# The max_abs_prediction_threshold was approximated using
# logit(np.finfo(np.float64).eps) which is about -36
def _sigmoid_calibration(
    predictions, y, sample_weight=None, max_abs_prediction_threshold=30
):
    """Probability Calibration with sigmoid method (Platt 2000)

    Parameters
    ----------
    predictions : ndarray of shape (n_samples,)
        The decision function or predict proba for the samples.

    y : ndarray of shape (n_samples,)
        The targets.

    sample_weight : array-like of shape (n_samples,), default=None
        Sample weights. If None, then samples are equally weighted.

    Returns
    -------
    a : float
        The slope.

    b : float
        The intercept.

    References
    ----------
    Platt, "Probabilistic Outputs for Support Vector Machines"
    """
    predictions = column_or_1d(predictions)
    y = column_or_1d(y)

    F = predictions  # F follows Platt's notations

    scale_constant = 1.0
    max_prediction = np.max(np.abs(F))

    # If the predictions have large values we scale them in order to bring
    # them within a suitable range. This has no effect on the final
    # (prediction) result because linear models like Logisitic Regression
    # without a penalty are invariant to multiplying the features by a
    # constant.
    if max_prediction >= max_abs_prediction_threshold:
        scale_constant = max_prediction
        # We rescale the features in a copy: inplace rescaling could confuse
        # the caller and make the code harder to reason about.
        F = F / scale_constant

    # Bayesian priors (see Platt end of section 2.2):
    # It corresponds to the number of samples, taking into account the
    # `sample_weight`.
    mask_negative_samples = y <= 0
    if sample_weight is not None:
        prior0 = (sample_weight[mask_negative_samples]).sum()
        prior1 = (sample_weight[~mask_negative_samples]).sum()
    else:
        prior0 = float(np.sum(mask_negative_samples))
        prior1 = y.shape[0] - prior0
    T = np.zeros_like(y, dtype=predictions.dtype)
    T[y > 0] = (prior1 + 1.0) / (prior1 + 2.0)
    T[y <= 0] = 1.0 / (prior0 + 2.0)

    bin_loss = HalfBinomialLoss()

    def loss_grad(AB):
        # .astype below is needed to ensure y_true and raw_prediction have the
        # same dtype. With result = np.float64(0) * np.array([1, 2], dtype=np.float32)
        # - in Numpy 2, result.dtype is float64
        # - in Numpy<2, result.dtype is float32
        raw_prediction = -(AB[0] * F + AB[1]).astype(dtype=predictions.dtype)
        l, g = bin_loss.loss_gradient(
            y_true=T,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
        )
        loss = l.sum()
        # TODO: Remove casting to np.float64 when minimum supported SciPy is 1.11.2
        # With SciPy >= 1.11.2, the LBFGS implementation will cast to float64
        # https://github.com/scipy/scipy/pull/18825.
        # Here we cast to float64 to support SciPy < 1.11.2
        grad = np.asarray([-g @ F, -g.sum()], dtype=np.float64)
        return loss, grad

    AB0 = np.array([0.0, log((prior0 + 1.0) / (prior1 + 1.0))])

    opt_result = minimize(
        loss_grad,
        AB0,
        method="L-BFGS-B",
        jac=True,
        options={
            "gtol": 1e-6,
            "ftol": 64 * np.finfo(float).eps,
        },
    )
    AB_ = opt_result.x

    # The tuned multiplicative parameter is converted back to the original
    # input feature scale. The offset parameter does not need rescaling since
    # we did not rescale the outcome variable.
    return AB_[0] / scale_constant, AB_[1]


def _convert_to_logits(decision_values, eps=1e-12, xp=None):
    """Convert decision_function values to 2D and predict_proba values to logits.

    This function ensures that the output of `decision_function` is
    converted into a (n_samples, n_classes) array. For binary classification,
    each row contains logits for the negative and positive classes as (-x, x).

    If `predict_proba` is provided instead, it is converted into
    log-probabilities using `numpy.log`.

    Parameters
    ----------
    decision_values : array-like of shape (n_samples,) or (n_samples, 1) \
        or (n_samples, n_classes).

        The decision function values or probability estimates.
        - If shape is (n_samples,), converts to (n_samples, 2) with (-x, x).
        - If shape is (n_samples, 1), converts to (n_samples, 2) with (-x, x).
        - If shape is (n_samples, n_classes), returns unchanged.
        - For probability estimates, returns `numpy.log(decision_values + eps)`.

    eps : float
        Small positive value added to avoid log(0).

    Returns
    -------
    logits : ndarray of shape (n_samples, n_classes)
    """
    xp, _, device_ = get_namespace_and_device(decision_values, xp=xp)
    decision_values = check_array(
        decision_values, dtype=[xp.float64, xp.float32], ensure_2d=False
    )
    if (decision_values.ndim == 2) and (decision_values.shape[1] > 1):
        # Check if it is the output of predict_proba
        entries_zero_to_one = xp.all((decision_values >= 0) & (decision_values <= 1))
        # TODO: simplify once upstream issue is addressed
        # https://github.com/data-apis/array-api-extra/issues/478
        row_sums_to_one = xp.all(
            xpx.isclose(
                xp.sum(decision_values, axis=1),
                xp.asarray(1.0, device=device_, dtype=decision_values.dtype),
            )
        )

        if entries_zero_to_one and row_sums_to_one:
            logits = xp.log(decision_values + eps)
        else:
            logits = decision_values

    elif (decision_values.ndim == 2) and (decision_values.shape[1] == 1):
        logits = xp.concat([-decision_values, decision_values], axis=1)

    elif decision_values.ndim == 1:
        decision_values = xp.reshape(decision_values, (-1, 1))
        logits = xp.concat([-decision_values, decision_values], axis=1)

    return logits


class _SigmoidCalibration(RegressorMixin, BaseEstimator):
    """Sigmoid regression model.

    Attributes
    ----------
    a_ : float
        The slope.

    b_ : float
        The intercept.
    """

    def fit(self, X, y, sample_weight=None):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : array-like of shape (n_samples,)
            Training data.

        y : array-like of shape (n_samples,)
            Training target.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        X = column_or_1d(X)
        y = column_or_1d(y)
        X, y = indexable(X, y)

        self.a_, self.b_ = _sigmoid_calibration(X, y, sample_weight)
        return self

    def predict(self, T):
        """Predict new data by linear interpolation.

        Parameters
        ----------
        T : array-like of shape (n_samples,)
            Data to predict from.

        Returns
        -------
        T_ : ndarray of shape (n_samples,)
            The predicted data.
        """
        T = column_or_1d(T)
        return expit(-(self.a_ * T + self.b_))


class _TemperatureScaling(RegressorMixin, BaseEstimator):
    """Temperature scaling model.

    Attributes
    ----------
    beta_ : float
        The optimized inverse temperature.
    """

    def fit(self, X, y, sample_weight=None):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : ndarray of shape (n_samples,) or (n_samples, n_classes)
            Training data.

            This should be the output of `decision_function` or `predict_proba`.
            If the input appears to be probabilities (i.e., values between 0 and 1
            that sum to 1 across classes), it will be converted to logits using
            `np.log(p + eps)`.

            Binary decision function outputs (1D) will be converted to two-class
            logits of the form (-x, x). For shapes of the form (n_samples, 1), the
            same process applies.

        y : array-like of shape (n_samples,)
            Training target.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        xp, _, xp_device = get_namespace_and_device(X, y)
        X, y = indexable(X, y)
        check_consistent_length(X, y)
        logits = _convert_to_logits(X)  # guarantees xp.float64 or xp.float32

        dtype_ = logits.dtype
        labels = column_or_1d(y, dtype=dtype_)

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, labels, dtype=dtype_)

        if _is_numpy_namespace(xp):
            multinomial_loss = HalfMultinomialLoss(n_classes=logits.shape[1])
        else:
            multinomial_loss = partial(_half_multinomial_loss, xp=xp)

        def log_loss(log_beta=0.0):
            """Compute the log loss as a parameter of the inverse temperature
            (beta).

            Parameters
            ----------
            log_beta : float
                The current logarithm of the inverse temperature value during
                optimisation.

            Returns
            -------
            negative_log_likelihood_loss : float
                The negative log likelihood loss.

            """
            # TODO: numpy 2.0
            # Ensure raw_prediction has the same dtype as labels using .astype().
            # Without this, dtype promotion rules differ across NumPy versions:
            #
            #   beta = np.float64(0)
            #   logits = np.array([1, 2], dtype=np.float32)
            #
            #   result = beta * logits
            #   - NumPy < 2: result.dtype is float32
            #   - NumPy 2+:  result.dtype is float64
            #
            #  This can cause dtype mismatch errors downstream (e.g., buffer dtype).
            log_beta = xp.asarray(log_beta, dtype=dtype_, device=xp_device)
            raw_prediction = xp.exp(log_beta) * logits
            return multinomial_loss(labels, raw_prediction, sample_weight)

        xatol = 64 * xp.finfo(dtype_).eps
        log_beta_minimizer = minimize_scalar(
            log_loss,
            bounds=(-10.0, 10.0),
            options={
                "xatol": xatol,
            },
        )

        if not log_beta_minimizer.success:  # pragma: no cover
            raise RuntimeError(
                "Temperature scaling fails to optimize during calibration. "
                "Reason from `scipy.optimize.minimize_scalar`: "
                f"{log_beta_minimizer.message}"
            )

        self.beta_ = xp.exp(
            xp.asarray(log_beta_minimizer.x, dtype=dtype_, device=xp_device)
        )

        return self

    def predict(self, X):
        """Predict new data by linear interpolation.

        Parameters
        ----------
        X : ndarray of shape (n_samples,) or (n_samples, n_classes)
            Data to predict from.

            This should be the output of `decision_function` or `predict_proba`.
            If the input appears to be probabilities (i.e., values between 0 and 1
            that sum to 1 across classes), it will be converted to logits using
            `np.log(p + eps)`.

            Binary decision function outputs (1D) will be converted to two-class
            logits of the form (-x, x). For shapes of the form (n_samples, 1), the
            same process applies.

        Returns
        -------
        X_ : ndarray of shape (n_samples, n_classes)
             The predicted data.
        """
        logits = _convert_to_logits(X)
        return softmax(self.beta_ * logits)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.one_d_array = True
        tags.input_tags.two_d_array = False
        return tags


@validate_params(
    {
        "y_true": ["array-like"],
        "y_prob": ["array-like"],
        "pos_label": [Real, str, "boolean", None],
        "n_bins": [Interval(Integral, 1, None, closed="left")],
        "strategy": [StrOptions({"uniform", "quantile"})],
    },
    prefer_skip_nested_validation=True,
)
def calibration_curve(
    y_true,
    y_prob,
    *,
    pos_label=None,
    n_bins=5,
    strategy="uniform",
):
    """Compute true and predicted probabilities for a calibration curve.

    The method assumes the inputs come from a binary classifier, and
    discretize the [0, 1] interval into bins.

    Calibration curves may also be referred to as reliability diagrams.

    Read more in the :ref:`User Guide <calibration>`.

    Parameters
    ----------
    y_true : array-like of shape (n_samples,)
        True targets.

    y_prob : array-like of shape (n_samples,)
        Probabilities of the positive class.

    pos_label : int, float, bool or str, default=None
        The label of the positive class.

        .. versionadded:: 1.1

    n_bins : int, default=5
        Number of bins to discretize the [0, 1] interval. A bigger number
        requires more data. Bins with no samples (i.e. without
        corresponding values in `y_prob`) will not be returned, thus the
        returned arrays may have less than `n_bins` values.

    strategy : {'uniform', 'quantile'}, default='uniform'
        Strategy used to define the widths of the bins.

        uniform
            The bins have identical widths.
        quantile
            The bins have the same number of samples and depend on `y_prob`.

    Returns
    -------
    prob_true : ndarray of shape (n_bins,) or smaller
        The proportion of samples whose class is the positive class, in each
        bin (fraction of positives).

    prob_pred : ndarray of shape (n_bins,) or smaller
        The mean predicted probability in each bin.

    See Also
    --------
    CalibrationDisplay.from_predictions : Plot calibration curve using true
        and predicted labels.
    CalibrationDisplay.from_estimator : Plot calibration curve using an
        estimator and data.

    References
    ----------
    Alexandru Niculescu-Mizil and Rich Caruana (2005) Predicting Good
    Probabilities With Supervised Learning, in Proceedings of the 22nd
    International Conference on Machine Learning (ICML).
    See section 4 (Qualitative Analysis of Predictions).

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.calibration import calibration_curve
    >>> y_true = np.array([0, 0, 0, 0, 1, 1, 1, 1, 1])
    >>> y_pred = np.array([0.1, 0.2, 0.3, 0.4, 0.65, 0.7, 0.8, 0.9,  1.])
    >>> prob_true, prob_pred = calibration_curve(y_true, y_pred, n_bins=3)
    >>> prob_true
    array([0. , 0.5, 1. ])
    >>> prob_pred
    array([0.2  , 0.525, 0.85 ])
    """
    y_true = column_or_1d(y_true)
    y_prob = column_or_1d(y_prob)
    check_consistent_length(y_true, y_prob)
    pos_label = _check_pos_label_consistency(pos_label, y_true)

    if y_prob.min() < 0 or y_prob.max() > 1:
        raise ValueError("y_prob has values outside [0, 1].")

    labels = np.unique(y_true)
    if len(labels) > 2:
        raise ValueError(
            f"Only binary classification is supported. Provided labels {labels}."
        )
    y_true = y_true == pos_label

    if strategy == "quantile":  # Determine bin edges by distribution of data
        quantiles = np.linspace(0, 1, n_bins + 1)
        bins = np.percentile(y_prob, quantiles * 100)
    elif strategy == "uniform":
        bins = np.linspace(0.0, 1.0, n_bins + 1)
    else:
        raise ValueError(
            "Invalid entry to 'strategy' input. Strategy "
            "must be either 'quantile' or 'uniform'."
        )

    binids = np.searchsorted(bins[1:-1], y_prob)

    bin_sums = np.bincount(binids, weights=y_prob, minlength=len(bins))
    bin_true = np.bincount(binids, weights=y_true, minlength=len(bins))
    bin_total = np.bincount(binids, minlength=len(bins))

    nonzero = bin_total != 0
    prob_true = bin_true[nonzero] / bin_total[nonzero]
    prob_pred = bin_sums[nonzero] / bin_total[nonzero]

    return prob_true, prob_pred


class CalibrationDisplay(_BinaryClassifierCurveDisplayMixin):
    """Calibration curve (also known as reliability diagram) visualization.

    It is recommended to use
    :func:`~sklearn.calibration.CalibrationDisplay.from_estimator` or
    :func:`~sklearn.calibration.CalibrationDisplay.from_predictions`
    to create a `CalibrationDisplay`. All parameters are stored as attributes.

    Read more about calibration in the :ref:`User Guide <calibration>` and
    more about the scikit-learn visualization API in :ref:`visualizations`.

    For an example on how to use the visualization, see
    :ref:`sphx_glr_auto_examples_calibration_plot_calibration_curve.py`.

    .. versionadded:: 1.0

    Parameters
    ----------
    prob_true : ndarray of shape (n_bins,)
        The proportion of samples whose class is the positive class (fraction
        of positives), in each bin.

    prob_pred : ndarray of shape (n_bins,)
        The mean predicted probability in each bin.

    y_prob : ndarray of shape (n_samples,)
        Probability estimates for the positive class, for each sample.

    estimator_name : str, default=None
        Name of estimator. If None, the estimator name is not shown.

    pos_label : int, float, bool or str, default=None
        The positive class when calibration curve computed.
        If not `None`, this value is displayed in the x- and y-axes labels.

        .. versionadded:: 1.1

    Attributes
    ----------
    line_ : matplotlib Artist
        Calibration curve.

    ax_ : matplotlib Axes
        Axes with calibration curve.

    figure_ : matplotlib Figure
        Figure containing the curve.

    See Also
    --------
    calibration_curve : Compute true and predicted probabilities for a
        calibration curve.
    CalibrationDisplay.from_predictions : Plot calibration curve using true
        and predicted labels.
    CalibrationDisplay.from_estimator : Plot calibration curve using an
        estimator and data.

    Examples
    --------
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.calibration import calibration_curve, CalibrationDisplay
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, random_state=0)
    >>> clf = LogisticRegression(random_state=0)
    >>> clf.fit(X_train, y_train)
    LogisticRegression(random_state=0)
    >>> y_prob = clf.predict_proba(X_test)[:, 1]
    >>> prob_true, prob_pred = calibration_curve(y_test, y_prob, n_bins=10)
    >>> disp = CalibrationDisplay(prob_true, prob_pred, y_prob)
    >>> disp.plot()
    <...>
    """

    def __init__(
        self, prob_true, prob_pred, y_prob, *, estimator_name=None, pos_label=None
    ):
        self.prob_true = prob_true
        self.prob_pred = prob_pred
        self.y_prob = y_prob
        self.estimator_name = estimator_name
        self.pos_label = pos_label

    def plot(self, *, ax=None, name=None, ref_line=True, **kwargs):
        """Plot visualization.

        Extra keyword arguments will be passed to
        :func:`matplotlib.pyplot.plot`.

        Parameters
        ----------
        ax : Matplotlib Axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        name : str, default=None
            Name for labeling curve. If `None`, use `estimator_name` if
            not `None`, otherwise no labeling is shown.

        ref_line : bool, default=True
            If `True`, plots a reference line representing a perfectly
            calibrated classifier.

        **kwargs : dict
            Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        display : :class:`~sklearn.calibration.CalibrationDisplay`
            Object that stores computed values.
        """
        self.ax_, self.figure_, name = self._validate_plot_params(ax=ax, name=name)

        info_pos_label = (
            f"(Positive class: {self.pos_label})" if self.pos_label is not None else ""
        )

        default_line_kwargs = {"marker": "s", "linestyle": "-"}
        if name is not None:
            default_line_kwargs["label"] = name
        line_kwargs = _validate_style_kwargs(default_line_kwargs, kwargs)

        ref_line_label = "Perfectly calibrated"
        existing_ref_line = ref_line_label in self.ax_.get_legend_handles_labels()[1]
        if ref_line and not existing_ref_line:
            self.ax_.plot([0, 1], [0, 1], "k:", label=ref_line_label)
        self.line_ = self.ax_.plot(self.prob_pred, self.prob_true, **line_kwargs)[0]

        # We always have to show the legend for at least the reference line
        self.ax_.legend(loc="lower right")

        xlabel = f"Mean predicted probability {info_pos_label}"
        ylabel = f"Fraction of positives {info_pos_label}"
        self.ax_.set(xlabel=xlabel, ylabel=ylabel)

        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        y,
        *,
        n_bins=5,
        strategy="uniform",
        pos_label=None,
        name=None,
        ax=None,
        ref_line=True,
        **kwargs,
    ):
        """Plot calibration curve using a binary classifier and data.

        A calibration curve, also known as a reliability diagram, uses inputs
        from a binary classifier and plots the average predicted probability
        for each bin against the fraction of positive classes, on the
        y-axis.

        Extra keyword arguments will be passed to
        :func:`matplotlib.pyplot.plot`.

        Read more about calibration in the :ref:`User Guide <calibration>` and
        more about the scikit-learn visualization API in :ref:`visualizations`.

        .. versionadded:: 1.0

        Parameters
        ----------
        estimator : estimator instance
            Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
            in which the last estimator is a classifier. The classifier must
            have a :term:`predict_proba` method.

        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input values.

        y : array-like of shape (n_samples,)
            Binary target values.

        n_bins : int, default=5
            Number of bins to discretize the [0, 1] interval into when
            calculating the calibration curve. A bigger number requires more
            data.

        strategy : {'uniform', 'quantile'}, default='uniform'
            Strategy used to define the widths of the bins.

            - `'uniform'`: The bins have identical widths.
            - `'quantile'`: The bins have the same number of samples and depend
              on predicted probabilities.

        pos_label : int, float, bool or str, default=None
            The positive class when computing the calibration curve.
            By default, `estimators.classes_[1]` is considered as the
            positive class.

            .. versionadded:: 1.1

        name : str, default=None
            Name for labeling curve. If `None`, the name of the estimator is
            used.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        ref_line : bool, default=True
            If `True`, plots a reference line representing a perfectly
            calibrated classifier.

        **kwargs : dict
            Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        display : :class:`~sklearn.calibration.CalibrationDisplay`.
            Object that stores computed values.

        See Also
        --------
        CalibrationDisplay.from_predictions : Plot calibration curve using true
            and predicted labels.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.linear_model import LogisticRegression
        >>> from sklearn.calibration import CalibrationDisplay
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = LogisticRegression(random_state=0)
        >>> clf.fit(X_train, y_train)
        LogisticRegression(random_state=0)
        >>> disp = CalibrationDisplay.from_estimator(clf, X_test, y_test)
        >>> plt.show()
        """
        y_prob, pos_label, name = cls._validate_and_get_response_values(
            estimator,
            X,
            y,
            response_method="predict_proba",
            pos_label=pos_label,
            name=name,
        )

        return cls.from_predictions(
            y,
            y_prob,
            n_bins=n_bins,
            strategy=strategy,
            pos_label=pos_label,
            name=name,
            ref_line=ref_line,
            ax=ax,
            **kwargs,
        )

    @classmethod
    def from_predictions(
        cls,
        y_true,
        y_prob,
        *,
        n_bins=5,
        strategy="uniform",
        pos_label=None,
        name=None,
        ax=None,
        ref_line=True,
        **kwargs,
    ):
        """Plot calibration curve using true labels and predicted probabilities.

        Calibration curve, also known as reliability diagram, uses inputs
        from a binary classifier and plots the average predicted probability
        for each bin against the fraction of positive classes, on the
        y-axis.

        Extra keyword arguments will be passed to
        :func:`matplotlib.pyplot.plot`.

        Read more about calibration in the :ref:`User Guide <calibration>` and
        more about the scikit-learn visualization API in :ref:`visualizations`.

        .. versionadded:: 1.0

        Parameters
        ----------
        y_true : array-like of shape (n_samples,)
            True labels.

        y_prob : array-like of shape (n_samples,)
            The predicted probabilities of the positive class.

        n_bins : int, default=5
            Number of bins to discretize the [0, 1] interval into when
            calculating the calibration curve. A bigger number requires more
            data.

        strategy : {'uniform', 'quantile'}, default='uniform'
            Strategy used to define the widths of the bins.

            - `'uniform'`: The bins have identical widths.
            - `'quantile'`: The bins have the same number of samples and depend
              on predicted probabilities.

        pos_label : int, float, bool or str, default=None
            The positive class when computing the calibration curve.
            When `pos_label=None`, if `y_true` is in {-1, 1} or {0, 1},
            `pos_label` is set to 1, otherwise an error will be raised.

            .. versionadded:: 1.1

        name : str, default=None
            Name for labeling curve.

        ax : matplotlib axes, default=None
            Axes object to plot on. If `None`, a new figure and axes is
            created.

        ref_line : bool, default=True
            If `True`, plots a reference line representing a perfectly
            calibrated classifier.

        **kwargs : dict
            Keyword arguments to be passed to :func:`matplotlib.pyplot.plot`.

        Returns
        -------
        display : :class:`~sklearn.calibration.CalibrationDisplay`.
            Object that stores computed values.

        See Also
        --------
        CalibrationDisplay.from_estimator : Plot calibration curve using an
            estimator and data.

        Examples
        --------
        >>> import matplotlib.pyplot as plt
        >>> from sklearn.datasets import make_classification
        >>> from sklearn.model_selection import train_test_split
        >>> from sklearn.linear_model import LogisticRegression
        >>> from sklearn.calibration import CalibrationDisplay
        >>> X, y = make_classification(random_state=0)
        >>> X_train, X_test, y_train, y_test = train_test_split(
        ...     X, y, random_state=0)
        >>> clf = LogisticRegression(random_state=0)
        >>> clf.fit(X_train, y_train)
        LogisticRegression(random_state=0)
        >>> y_prob = clf.predict_proba(X_test)[:, 1]
        >>> disp = CalibrationDisplay.from_predictions(y_test, y_prob)
        >>> plt.show()
        """
        pos_label_validated, name = cls._validate_from_predictions_params(
            y_true, y_prob, sample_weight=None, pos_label=pos_label, name=name
        )

        prob_true, prob_pred = calibration_curve(
            y_true, y_prob, n_bins=n_bins, strategy=strategy, pos_label=pos_label
        )

        disp = cls(
            prob_true=prob_true,
            prob_pred=prob_pred,
            y_prob=y_prob,
            estimator_name=name,
            pos_label=pos_label_validated,
        )
        return disp.plot(ax=ax, ref_line=ref_line, **kwargs)
