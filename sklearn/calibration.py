"""Calibration of predicted probabilities."""

# Author: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#         Balazs Kegl <balazs.kegl@gmail.com>
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         Mathieu Blondel <mathieu@mblondel.org>
#
# License: BSD 3 clause

import warnings
from inspect import signature
from contextlib import suppress
from functools import partial

from math import log
import numpy as np
from joblib import Parallel

from scipy.special import expit
from scipy.special import xlogy
from scipy.optimize import fmin_bfgs

from .base import (BaseEstimator, ClassifierMixin, RegressorMixin, clone,
                   MetaEstimatorMixin)
from .preprocessing import label_binarize, LabelEncoder
from .utils import (
    check_array,
    column_or_1d,
    deprecated,
    indexable,
)
from .utils.multiclass import check_classification_targets
from .utils.fixes import delayed
from .utils.validation import check_is_fitted, check_consistent_length
from .utils.validation import _check_sample_weight
from .pipeline import Pipeline
from .isotonic import IsotonicRegression
from .svm import LinearSVC
from .model_selection import check_cv, cross_val_predict
from .utils.validation import _deprecate_positional_args


class CalibratedClassifierCV(ClassifierMixin,
                             MetaEstimatorMixin,
                             BaseEstimator):
    """Probability calibration with isotonic regression or logistic regression.

    This class uses cross-validation to both estimate the parameters of a
    classifier and subsequently calibrate a classifier. With default
    `ensemble=True`, for each cv split it
    fits a copy of the base estimator to the training subset, and calibrates it
    using the testing subset. For prediction, predicted probabilities are
    averaged across these individual calibrated classifiers. When
    `ensemble=False`, cross-validation is used to obtain unbiased predictions,
    via :func:`~sklearn.model_selection.cross_val_predict`, which are then
    used for calibration. For prediction, the base estimator, trained using all
    the data, is used. This is the method implemented when `probabilities=True`
    for :mod:`sklearn.svm` estimators.

    Already fitted classifiers can be calibrated via the parameter
    `cv="prefit"`. In this case, no cross-validation is used and all provided
    data is used for calibration. The user has to take care manually that data
    for model fitting and calibration are disjoint.

    The calibration is based on the :term:`decision_function` method of the
    `base_estimator` if it exists, else on :term:`predict_proba`.

    Read more in the :ref:`User Guide <calibration>`.

    Parameters
    ----------
    base_estimator : estimator instance, default=None
        The classifier whose output need to be calibrated to provide more
        accurate `predict_proba` outputs. The default classifier is
        a :class:`~sklearn.svm.LinearSVC`.

    method : {'sigmoid', 'isotonic'}, default='sigmoid'
        The method to use for calibration. Can be 'sigmoid' which
        corresponds to Platt's method (i.e. a logistic regression model) or
        'isotonic' which is a non-parametric approach. It is not advised to
        use isotonic calibration with too few calibration samples
        ``(<<1000)`` since it tends to overfit.

    cv : int, cross-validation generator, iterable or "prefit", \
            default=None
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

        If "prefit" is passed, it is assumed that `base_estimator` has been
        fitted already and all data is used for calibration.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    n_jobs : int, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors.

        Base estimator clones are fitted in parallel across cross-validation
        iterations. Therefore parallelism happens only when `cv != "prefit"`.

        See :term:`Glossary <n_jobs>` for more details.

        .. versionadded:: 0.24

    ensemble : bool, default=True
        Determines how the calibrator is fitted when `cv` is not `'prefit'`.
        Ignored if `cv='prefit'`.

        If `True`, the `base_estimator` is fitted using training data and
        calibrated using testing data, for each `cv` fold. The final estimator
        is an ensemble of `n_cv` fitted classifer and calibrator pairs, where
        `n_cv` is the number of cross-validation folds. The output is the
        average predicted probabilities of all pairs.

        If `False`, `cv` is used to compute unbiased predictions, via
        :func:`~sklearn.model_selection.cross_val_predict`, which are then
        used for calibration. At prediction time, the classifier used is the
        `base_estimator` trained on all the data.
        Note that this method is also internally implemented  in
        :mod:`sklearn.svm` estimators with the `probabilities=True` parameter.

        .. versionadded:: 0.24

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes,)
        The class labels.

    calibrated_classifiers_ : list (len() equal to cv or 1 if `cv="prefit"` \
            or `ensemble=False`)
        The list of classifier and calibrator pairs.

        - When `cv="prefit"`, the fitted `base_estimator` and fitted
          calibrator.
        - When `cv` is not "prefit" and `ensemble=True`, `n_cv` fitted
          `base_estimator` and calibrator pairs. `n_cv` is the number of
          cross-validation folds.
        - When `cv` is not "prefit" and `ensemble=False`, the `base_estimator`,
          fitted on all the data, and fitted calibrator.

        .. versionchanged:: 0.24
            Single calibrated classifier case when `ensemble=False`.

    Examples
    --------
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.naive_bayes import GaussianNB
    >>> from sklearn.calibration import CalibratedClassifierCV
    >>> X, y = make_classification(n_samples=100, n_features=2,
    ...                            n_redundant=0, random_state=42)
    >>> base_clf = GaussianNB()
    >>> calibrated_clf = CalibratedClassifierCV(base_estimator=base_clf, cv=3)
    >>> calibrated_clf.fit(X, y)
    CalibratedClassifierCV(base_estimator=GaussianNB(), cv=3)
    >>> len(calibrated_clf.calibrated_classifiers_)
    3
    >>> calibrated_clf.predict_proba(X)[:5, :]
    array([[0.110..., 0.889...],
           [0.072..., 0.927...],
           [0.928..., 0.071...],
           [0.928..., 0.071...],
           [0.071..., 0.928...]])

    >>> from sklearn.model_selection import train_test_split
    >>> X, y = make_classification(n_samples=100, n_features=2,
    ...                            n_redundant=0, random_state=42)
    >>> X_train, X_calib, y_train, y_calib = train_test_split(
    ...        X, y, random_state=42
    ... )
    >>> base_clf = GaussianNB()
    >>> base_clf.fit(X_train, y_train)
    GaussianNB()
    >>> calibrated_clf = CalibratedClassifierCV(
    ...     base_estimator=base_clf,
    ...     cv="prefit"
    ... )
    >>> calibrated_clf.fit(X_calib, y_calib)
    CalibratedClassifierCV(base_estimator=GaussianNB(), cv='prefit')
    >>> len(calibrated_clf.calibrated_classifiers_)
    1
    >>> calibrated_clf.predict_proba([[-0.5, 0.5]])
    array([[0.936..., 0.063...]])

    References
    ----------
    .. [1] Obtaining calibrated probability estimates from decision trees
           and naive Bayesian classifiers, B. Zadrozny & C. Elkan, ICML 2001

    .. [2] Transforming Classifier Scores into Accurate Multiclass
           Probability Estimates, B. Zadrozny & C. Elkan, (KDD 2002)

    .. [3] Probabilistic Outputs for Support Vector Machines and Comparisons to
           Regularized Likelihood Methods, J. Platt, (1999)

    .. [4] Predicting Good Probabilities with Supervised Learning,
           A. Niculescu-Mizil & R. Caruana, ICML 2005
    """
    @_deprecate_positional_args
    def __init__(self, base_estimator=None, *, method='sigmoid',
                 cv=None, n_jobs=None, ensemble=True):
        self.base_estimator = base_estimator
        self.method = method
        self.cv = cv
        self.n_jobs = n_jobs
        self.ensemble = ensemble

    def fit(self, X, y, sample_weight=None):
        """Fit the calibrated model.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        check_classification_targets(y)
        X, y = indexable(X, y)

        if self.base_estimator is None:
            # we want all classifiers that don't expose a random_state
            # to be deterministic (and we don't want to expose this one).
            base_estimator = LinearSVC(random_state=0)
        else:
            base_estimator = self.base_estimator

        self.calibrated_classifiers_ = []
        if self.cv == "prefit":
            # `classes_` and `n_features_in_` should be consistent with that
            # of base_estimator
            if isinstance(self.base_estimator, Pipeline):
                check_is_fitted(self.base_estimator[-1])
            else:
                check_is_fitted(self.base_estimator)
            with suppress(AttributeError):
                self.n_features_in_ = base_estimator.n_features_in_
            self.classes_ = self.base_estimator.classes_

            pred_method = _get_prediction_method(base_estimator)
            n_classes = len(self.classes_)
            predictions = _compute_predictions(pred_method, X, n_classes)

            calibrated_classifier = _fit_calibrator(
                base_estimator, predictions, y, self.classes_, self.method,
                sample_weight
            )
            self.calibrated_classifiers_.append(calibrated_classifier)
        else:
            X, y = self._validate_data(
                X, y, accept_sparse=['csc', 'csr', 'coo'],
                force_all_finite=False, allow_nd=True
            )
            # Set `classes_` using all `y`
            label_encoder_ = LabelEncoder().fit(y)
            self.classes_ = label_encoder_.classes_
            n_classes = len(self.classes_)

            # sample_weight checks
            fit_parameters = signature(base_estimator.fit).parameters
            supports_sw = "sample_weight" in fit_parameters
            if sample_weight is not None:
                sample_weight = _check_sample_weight(sample_weight, X)
                if not supports_sw:
                    estimator_name = type(base_estimator).__name__
                    warnings.warn(f"Since {estimator_name} does not support "
                                  "sample_weights, sample weights will only be"
                                  " used for the calibration itself.")

            # Check that each cross-validation fold can have at least one
            # example per class
            if isinstance(self.cv, int):
                n_folds = self.cv
            elif hasattr(self.cv, "n_splits"):
                n_folds = self.cv.n_splits
            else:
                n_folds = None
            if n_folds and np.any([np.sum(y == class_) < n_folds
                                   for class_ in self.classes_]):
                raise ValueError(f"Requesting {n_folds}-fold "
                                 "cross-validation but provided less than "
                                 f"{n_folds} examples for at least one class.")
            cv = check_cv(self.cv, y, classifier=True)

            if self.ensemble:
                parallel = Parallel(n_jobs=self.n_jobs)

                self.calibrated_classifiers_ = parallel(
                    delayed(_fit_classifier_calibrator_pair)(
                        clone(base_estimator), X, y, train=train, test=test,
                        method=self.method, classes=self.classes_,
                        supports_sw=supports_sw, sample_weight=sample_weight)
                    for train, test in cv.split(X, y)
                )
            else:
                this_estimator = clone(base_estimator)
                method_name = _get_prediction_method(this_estimator).__name__
                pred_method = partial(
                    cross_val_predict, estimator=this_estimator, X=X, y=y,
                    cv=cv, method=method_name, n_jobs=self.n_jobs
                )
                predictions = _compute_predictions(pred_method, X, n_classes)

                if sample_weight is not None and supports_sw:
                    this_estimator.fit(X, y, sample_weight)
                else:
                    this_estimator.fit(X, y)
                calibrated_classifier = _fit_calibrator(
                    this_estimator, predictions, y, self.classes_, self.method,
                    sample_weight
                )
                self.calibrated_classifiers_.append(calibrated_classifier)

        return self

    def predict_proba(self, X):
        """Calibrated probabilities of classification.

        This function returns calibrated probabilities of classification
        according to each class on an array of test vectors X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The samples.

        Returns
        -------
        C : ndarray of shape (n_samples, n_classes)
            The predicted probas.
        """
        check_is_fitted(self)
        X = check_array(X, accept_sparse=['csc', 'csr', 'coo'],
                        force_all_finite=False)
        # Compute the arithmetic mean of the predictions of the calibrated
        # classifiers
        mean_proba = np.zeros((X.shape[0], len(self.classes_)))
        for calibrated_classifier in self.calibrated_classifiers_:
            proba = calibrated_classifier.predict_proba(X)
            mean_proba += proba

        mean_proba /= len(self.calibrated_classifiers_)

        return mean_proba

    def predict(self, X):
        """Predict the target of new samples. The predicted class is the
        class that has the highest probability, and can thus be different
        from the prediction of the uncalibrated classifier.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The samples.

        Returns
        -------
        C : ndarray of shape (n_samples,)
            The predicted class.
        """
        check_is_fitted(self)
        return self.classes_[np.argmax(self.predict_proba(X), axis=1)]

    def _more_tags(self):
        return {
            '_xfail_checks': {
                'check_sample_weights_invariance':
                'zero sample_weight is not equivalent to removing samples',
            }
        }


def _fit_classifier_calibrator_pair(estimator, X, y, train, test, supports_sw,
                                    method, classes, sample_weight=None):
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

    train : ndarray, shape (n_train_indicies,)
        Indices of the training subset.

    test : ndarray, shape (n_test_indicies,)
        Indices of the testing subset.

    supports_sw : bool
        Whether or not the `estimator` supports sample weights.

    method : {'sigmoid', 'isotonic'}
        Method to use for calibration.

    classes : ndarray, shape (n_classes,)
        The target classes.

    sample_weight : array-like, default=None
        Sample weights for `X`.

    Returns
    -------
    calibrated_classifier : _CalibratedClassifier instance
    """
    if sample_weight is not None and supports_sw:
        estimator.fit(X[train], y[train],
                      sample_weight=sample_weight[train])
    else:
        estimator.fit(X[train], y[train])

    n_classes = len(classes)
    pred_method = _get_prediction_method(estimator)
    predictions = _compute_predictions(pred_method, X[test], n_classes)

    sw = None if sample_weight is None else sample_weight[test]
    calibrated_classifier = _fit_calibrator(
        estimator, predictions, y[test], classes, method, sample_weight=sw
    )
    return calibrated_classifier


def _get_prediction_method(clf):
    """Return prediction method.

    `decision_function` method of `clf` returned, if it
    exists, otherwise `predict_proba` method returned.

    Parameters
    ----------
    clf : Estimator instance
        Fitted classifier to obtain the prediction method from.

    Returns
    -------
    prediction_method : callable
        The prediction method.
    """
    if hasattr(clf, 'decision_function'):
        method = getattr(clf, 'decision_function')
    elif hasattr(clf, 'predict_proba'):
        method = getattr(clf, 'predict_proba')
    else:
        raise RuntimeError("'base_estimator' has no 'decision_function' or "
                           "'predict_proba' method.")
    return method


def _compute_predictions(pred_method, X, n_classes):
    """Return predictions for `X` and reshape binary outputs to shape
    (n_samples, 1).

    Parameters
    ----------
    pred_method : callable
        Prediction method.

    X : array-like or None
        Data used to obtain predictions.

    n_classes : int
        Number of classes present.

    Returns
    -------
    predictions : array-like, shape (X.shape[0], len(clf.classes_))
        The predictions. Note if there are 2 classes, array is of shape
        (X.shape[0], 1).
    """
    predictions = pred_method(X=X)
    if hasattr(pred_method, '__name__'):
        method_name = pred_method.__name__
    else:
        method_name = signature(pred_method).parameters['method'].default

    if method_name == 'decision_function':
        if predictions.ndim == 1:
            predictions = predictions[:, np.newaxis]
    elif method_name == 'predict_proba':
        if n_classes == 2:
            predictions = predictions[:, 1:]
    else:  # pragma: no cover
        # this branch should be unreachable.
        raise ValueError(f"Invalid prediction method: {method_name}")
    return predictions


def _fit_calibrator(clf, predictions, y, classes, method, sample_weight=None):
    """Fit calibrator(s) and return a `_CalibratedClassifier`
    instance.

    `n_classes` (i.e. `len(clf.classes_)`) calibrators are fitted.
    However, if `n_classes` equals 2, one calibrator is fitted.

    Parameters
    ----------
    clf : estimator instance
        Fitted classifier.

    predictions : array-like, shape (n_samples, n_classes) or (n_samples, 1) \
                    when binary.
        Raw predictions returned by the un-calibrated base classifier.

    y : array-like, shape (n_samples,)
        The targets.

    classes : ndarray, shape (n_classes,)
        All the prediction classes.

    method : {'sigmoid', 'isotonic'}
        The method to use for calibration.

    sample_weight : ndarray, shape (n_samples,), default=None
        Sample weights. If None, then samples are equally weighted.

    Returns
    -------
    pipeline : _CalibratedClassifier instance
    """
    Y = label_binarize(y, classes=classes)
    label_encoder = LabelEncoder().fit(classes)
    pos_class_indices = label_encoder.transform(clf.classes_)
    calibrators = []
    for class_idx, this_pred in zip(pos_class_indices, predictions.T):
        if method == 'isotonic':
            calibrator = IsotonicRegression(out_of_bounds='clip')
        elif method == 'sigmoid':
            calibrator = _SigmoidCalibration()
        else:
            raise ValueError("'method' should be one of: 'sigmoid' or "
                             f"'isotonic'. Got {method}.")
        calibrator.fit(this_pred, Y[:, class_idx], sample_weight)
        calibrators.append(calibrator)

    pipeline = _CalibratedClassifier(
        clf, calibrators, method=method, classes=classes
    )
    return pipeline


class _CalibratedClassifier:
    """Pipeline-like chaining a fitted classifier and its fitted calibrators.

    Parameters
    ----------
    base_estimator : estimator instance
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

    Attributes
    ----------
    calibrators_ : list of fitted estimator instances
        Same as `calibrators`. Exposed for backward-compatibility. Use
        `calibrators` instead.

        .. deprecated:: 0.24
           `calibrators_` is deprecated from 0.24 and will be removed in
           1.1 (renaming of 0.26). Use `calibrators` instead.
    """
    def __init__(self, base_estimator, calibrators, *, classes,
                 method='sigmoid'):
        self.base_estimator = base_estimator
        self.calibrators = calibrators
        self.classes = classes
        self.method = method

    # TODO: Remove in 1.1
    # mypy error: Decorated property not supported
    @deprecated(  # type: ignore
        "calibrators_ is deprecated in 0.24 and will be removed in 1.1"
        "(renaming of 0.26). Use calibrators instead."
    )
    @property
    def calibrators_(self):
        return self.calibrators

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
        n_classes = len(self.classes)
        pred_method = _get_prediction_method(self.base_estimator)
        predictions = _compute_predictions(pred_method, X, n_classes)

        label_encoder = LabelEncoder().fit(self.classes)
        pos_class_indices = label_encoder.transform(
            self.base_estimator.classes_
        )

        proba = np.zeros((X.shape[0], n_classes))
        for class_idx, this_pred, calibrator in \
                zip(pos_class_indices, predictions.T, self.calibrators):
            if n_classes == 2:
                # When binary, `predictions` consists only of predictions for
                # clf.classes_[1] but `pos_class_indices` = 0
                class_idx += 1
            proba[:, class_idx] = calibrator.predict(this_pred)

        # Normalize the probabilities
        if n_classes == 2:
            proba[:, 0] = 1. - proba[:, 1]
        else:
            denominator = np.sum(proba, axis=1)[:, np.newaxis]
            # In the edge case where for each class calibrator returns a null
            # probability for a given sample, use the uniform distribution
            # instead.
            uniform_proba = np.full_like(proba, 1 / n_classes)
            proba = np.divide(proba, denominator, out=uniform_proba,
                              where=denominator != 0)

        # Deal with cases where the predicted probability minimally exceeds 1.0
        proba[(1.0 < proba) & (proba <= 1.0 + 1e-5)] = 1.0

        return proba


def _sigmoid_calibration(predictions, y, sample_weight=None):
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

    # Bayesian priors (see Platt end of section 2.2)
    prior0 = float(np.sum(y <= 0))
    prior1 = y.shape[0] - prior0
    T = np.zeros(y.shape)
    T[y > 0] = (prior1 + 1.) / (prior1 + 2.)
    T[y <= 0] = 1. / (prior0 + 2.)
    T1 = 1. - T

    def objective(AB):
        # From Platt (beginning of Section 2.2)
        P = expit(-(AB[0] * F + AB[1]))
        loss = -(xlogy(T, P) + xlogy(T1, 1. - P))
        if sample_weight is not None:
            return (sample_weight * loss).sum()
        else:
            return loss.sum()

    def grad(AB):
        # gradient of the objective function
        P = expit(-(AB[0] * F + AB[1]))
        TEP_minus_T1P = T - P
        if sample_weight is not None:
            TEP_minus_T1P *= sample_weight
        dA = np.dot(TEP_minus_T1P, F)
        dB = np.sum(TEP_minus_T1P)
        return np.array([dA, dB])

    AB0 = np.array([0., log((prior0 + 1.) / (prior1 + 1.))])
    AB_ = fmin_bfgs(objective, AB0, fprime=grad, disp=False)
    return AB_[0], AB_[1]


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


@_deprecate_positional_args
def calibration_curve(y_true, y_prob, *, normalize=False, n_bins=5,
                      strategy='uniform'):
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

    normalize : bool, default=False
        Whether y_prob needs to be normalized into the [0, 1] interval, i.e.
        is not a proper probability. If True, the smallest value in y_prob
        is linearly mapped onto 0 and the largest one onto 1.

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

    if normalize:  # Normalize predicted values into interval [0, 1]
        y_prob = (y_prob - y_prob.min()) / (y_prob.max() - y_prob.min())
    elif y_prob.min() < 0 or y_prob.max() > 1:
        raise ValueError("y_prob has values outside [0, 1] and normalize is "
                         "set to False.")

    labels = np.unique(y_true)
    if len(labels) > 2:
        raise ValueError("Only binary classification is supported. "
                         "Provided labels %s." % labels)
    y_true = label_binarize(y_true, classes=labels)[:, 0]

    if strategy == 'quantile':  # Determine bin edges by distribution of data
        quantiles = np.linspace(0, 1, n_bins + 1)
        bins = np.percentile(y_prob, quantiles * 100)
        bins[-1] = bins[-1] + 1e-8
    elif strategy == 'uniform':
        bins = np.linspace(0., 1. + 1e-8, n_bins + 1)
    else:
        raise ValueError("Invalid entry to 'strategy' input. Strategy "
                         "must be either 'quantile' or 'uniform'.")

    binids = np.digitize(y_prob, bins) - 1

    bin_sums = np.bincount(binids, weights=y_prob, minlength=len(bins))
    bin_true = np.bincount(binids, weights=y_true, minlength=len(bins))
    bin_total = np.bincount(binids, minlength=len(bins))

    nonzero = bin_total != 0
    prob_true = bin_true[nonzero] / bin_total[nonzero]
    prob_pred = bin_sums[nonzero] / bin_total[nonzero]

    return prob_true, prob_pred
