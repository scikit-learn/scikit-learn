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

from math import log
import numpy as np

from scipy.special import expit
from scipy.special import xlogy
from scipy.optimize import fmin_bfgs

from .base import (BaseEstimator, ClassifierMixin, RegressorMixin, clone,
                   MetaEstimatorMixin)
from .preprocessing import label_binarize, LabelEncoder
from .utils import check_array, indexable, column_or_1d
from .utils.multiclass import check_classification_targets
from .utils.validation import check_is_fitted, check_consistent_length
from .utils.validation import _check_sample_weight
from .pipeline import Pipeline
from .isotonic import IsotonicRegression
from .svm import LinearSVC
from .model_selection import check_cv, cross_val_predict
from .utils.validation import _deprecate_positional_args


class CalibratedClassifierCV(BaseEstimator, ClassifierMixin,
                             MetaEstimatorMixin):
    """Probability calibration with isotonic regression or logistic regression.

    This class uses cross-validation to both estimate the parameters of a
    classifier and subsequently calibrate a classifier. For each cv split it
    fits a copy of the base estimator to the training folds, and calibrates it
    using the testing fold. For prediction, predicted probabilities are
    averaged across these individual calibrated classifiers.

    Already fitted classifiers can be calibrated via the parameter cv="prefit".
    In this case, no cross-validation is used and all provided data is used
    for calibration. The user has to take care manually that data for model
    fitting and calibration are disjoint.

    The calibration is based on the :term:`decision_function` method of the
    `base_estimator` if it exists, else on :term:`predict_proba`.

    Read more in the :ref:`User Guide <calibration>`.

    Parameters
    ----------
    base_estimator : instance BaseEstimator
        The classifier whose output need to be calibrated to provide more
        accurate `predict_proba` outputs.

    method : {'sigmoid', 'isotonic'}, default='sigmoid'
        The method to use for calibration. Can be 'sigmoid' which
        corresponds to Platt's method (i.e. a logistic regression model) or
        'isotonic' which is a non-parametric approach. It is not advised to
        use isotonic calibration with too few calibration samples
        ``(<<1000)`` since it tends to overfit.

    cv : integer, cross-validation generator, iterable or "prefit", \
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

    ensemble : bool, default=True
        Determines how the calibrator is fit, if `cv` is not `'prefit'`.
        Ignored if `cv='prefit'`.

        If `True`, the `base_estimator` is fit and calibrated on each
        `cv` fold. The final estimator is an ensemble that outputs the
        average predicted probabilities of all fitted classifier and calibrator
        pairs.

        If `False`, `cv` is used to compute unbiased predictions, which
        are concatenated and used to train the calibrator (sigmoid or isotonic
        model). The `base_estimator` trained on all the data is used at
        prediction time.
        Note this method is implemented when `probabilities=True` for
        :mod:`sklearn.svm` estimators.

        .. versionadded:: 0.24

    Attributes
    ----------
    classes_ : array, shape (n_classes)
        The class labels.

    calibrated_classifiers_ : list (len() equal to cv or 1 if cv == "prefit")
        The list of calibrated classifiers, one for each cross-validation
        split, which has been fitted on training folds and
        calibrated on the testing fold.

    n_features_in_ : int
        The number of features in `X`. If `cv='prefit'`, number of features
        in the data used to fit `base_estimator`.

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
    def __init__(self, base_estimator=None, *, method='sigmoid', cv=None,
                 ensemble=True):
        self.base_estimator = base_estimator
        self.method = method
        self.cv = cv
        self.ensemble = ensemble

    def fit(self, X, y, sample_weight=None):
        """Fit the calibrated model

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        y : array-like, shape (n_samples,)
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
            label_encoder_ = LabelEncoder().fit(self.classes_)

            preds = _get_predictions(
                base_estimator, X, label_encoder_
            )

            calibrated_classifier = _fit_calibrator(
                base_estimator, preds, y, label_encoder_, self.method,
                sample_weight
            )
            self.calibrated_classifiers_.append(calibrated_classifier)
        else:
            X, y = self._validate_data(
                X, y, accept_sparse=['csc', 'csr', 'coo'],
                force_all_finite=False, allow_nd=True
            )
            # Set attributes using all `y`
            label_encoder_ = LabelEncoder().fit(y)
            self.classes_ = label_encoder_.classes_

            fit_parameters = signature(base_estimator.fit).parameters
            base_estimator_supports_sw = "sample_weight" in fit_parameters

            if sample_weight is not None:
                sample_weight = _check_sample_weight(sample_weight, X)

                if not base_estimator_supports_sw:
                    estimator_name = type(base_estimator).__name__
                    warnings.warn("Since %s does not support sample_weights, "
                                  "sample weights will only be used for the "
                                  "calibration itself." % estimator_name)
            if self.ensemble:
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
                                     f"{n_folds} examples for at least one "
                                     "class.")
                cv = check_cv(self.cv, y, classifier=True)

                for train, test in cv.split(X, y):
                    this_estimator = clone(base_estimator)

                    if (sample_weight is not None
                            and base_estimator_supports_sw):
                        this_estimator.fit(X[train], y[train],
                                           sample_weight=sample_weight[train])
                    else:
                        this_estimator.fit(X[train], y[train])

                    preds = _get_predictions(
                        this_estimator, X[test], label_encoder_
                    )

                    sw = None if sample_weight is None else sample_weight[test]
                    calibrated_classifier = _fit_calibrator(
                        this_estimator, preds, y[test], label_encoder_,
                        self.method, sample_weight=sw
                    )
                    self.calibrated_classifiers_.append(calibrated_classifier)
            else:
                pred_method = get_prediction_method(
                    base_estimator, return_string=True
                )
                preds = cross_val_predict(base_estimator, X, y, cv=cv,
                                          method=pred_method)
                preds = _reshape_preds(
                    preds, pred_method, len(label_encoder_.classes_)
                )

                this_estimator = clone(base_estimator)
                if sample_weight is not None and base_estimator_supports_sw:
                    this_estimator.fit(X, y, sample_weight)
                else:
                    this_estimator.fit(X, y)
                calibrated_classifier = _fit_calibrator(
                    this_estimator, preds, y, label_encoder_, self.method,
                    sample_weight
                )
                self.calibrated_classifiers_.append(calibrated_classifier)
        return self

    def predict_proba(self, X):
        """Calibrated probabilities of classification

        This function returns calibrated probabilities of classification
        according to each class on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The samples.

        Returns
        -------
        C : array, shape (n_samples, n_classes)
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
        X : array-like, shape (n_samples, n_features)
            The samples.

        Returns
        -------
        C : array, shape (n_samples,)
            The predicted class.
        """
        check_is_fitted(self)
        return self.classes_[np.argmax(self.predict_proba(X), axis=1)]

    def _more_tags(self):
        return {
            '_xfail_checks': {
                'check_sample_weights_invariance(kind=zeros)':
                'zero sample_weight is not equivalent to removing samples',
            }
        }


def get_prediction_method(clf_fitted, return_string=False):
    """Return prediction method or their corresponding name as string.

    `decision_function` method of `clf_fitted` returned, if it
    exists, otherwise `predict_proba` method returned.

    Parameters
    ----------
    clf_fitted : Estimator instance
        Classifier to obtain the prediction method from.

    return_string : bool, default=False
        Whether to return the method name as string instead of the prediction
        method.

    Returns
    -------
    prediction_method : callable or str
        If `return_string=True`, name of the prediction method as string.
        If `return_string=False`, the prediction method.
    """
    if hasattr(clf_fitted, "decision_function"):
        method_name = "decision_function"
        method = getattr(clf_fitted, "decision_function")
    elif hasattr(clf_fitted, "predict_proba"):
        method_name = "predict_proba"
        method = getattr(clf_fitted, "predict_proba")
    else:
        raise RuntimeError("'base_estimator' has no 'decision_function' or "
                           "'predict_proba' method.")
    if return_string:
        return method_name
    else:
        return method


def _reshape_preds(preds, method, n_classes):
    """Reshape predictions when classification binary.

    Parameters
    ----------
    preds : array-like
        Predictions.

    method : {'decision_function', 'predict_proba'}
        Method used to obtain the predictions.

    n_classes : int
        Number of classes.

    Returns
    -------
    preds : array-like, shape (n_samples, 1)
        Reshaped predictions
    """
    if method == 'decision_function':
        if preds.ndim == 1:
            preds = preds[:, np.newaxis]
    elif method == 'predict_proba':
        if n_classes == 2:
            preds = preds[:, 1:]
    else:
        raise RuntimeError("'method' needs to be one of 'decision_function' "
                           "or 'predict_proba'.")
    return preds


def _get_predictions(clf_fitted, X, label_encoder_):
    """Returns predictions for `X` and the index of classes present.

    Parameters
    ----------
    clf_fitted : Estimator instance
        Fitted classifier.

    X : array-like
        Data used to obtain predictions.

    label_encoder_ : LabelEncoder instance
        LabelEncoder instance fitted on all the targets.

    Returns
    -------
    preds : array-like, shape (X.shape[0], len(clf_fitted.classes_))
        The predictions. Note if there are 2 classes, array is of shape
        (X.shape[0], 1).
    """
    pred_method = get_prediction_method(clf_fitted)
    preds = pred_method(X)
    n_classes = len(clf_fitted.classes_)
    preds = _reshape_preds(preds, pred_method.__name__, n_classes)

    return preds


def _fit_calibrator(clf_fitted, preds, y, label_encoder_, method,
                    sample_weight=None):
    """Fit calibrator(s) and return a `_CalibratedClassiferPipeline`
    instance.

    `n_classes` (i.e. `len(clf_fitted.classes_)`) calibrators are fitted.
    However, if `n_classes` equals 2, one calibrator is fit.

    Parameters
    ----------
    clf_fitted : Estimator instance
        Fitted classifier.

    preds :  array-like, shape (n_samples, n_classes)
        Predictions for calibrating the predictions.
        If binary, shape (n_samples, 1).

    y : ndarray, shape (n_samples,)
        The targets.

    label_encoder_ : LabelEncoder instance
        LabelEncoder instance fitted on all the targets.

    method : {'sigmoid', 'isotonic'}
        The method to use for calibration.

    sample_weight : ndarray, shape (n_samples,), default=None
        Sample weights. If None, then samples are equally weighted.

    Returns
    -------
    pipeline : _CalibratedClassiferPipeline instance
    """
    Y = label_binarize(y, classes=label_encoder_.classes_)
    pos_class_indices = label_encoder_.transform(clf_fitted.classes_)
    calibrated_classifiers = []
    for class_idx, this_pred in zip(pos_class_indices, preds.T):
        if method == 'isotonic':
            calibrator = IsotonicRegression(out_of_bounds='clip')
        elif method == 'sigmoid':
            calibrator = _SigmoidCalibration()
        else:
            raise ValueError("'method' should be one of: 'sigmoid' or "
                             f"'isotonic'. Got {method}.")
        calibrator.fit(this_pred, Y[:, class_idx], sample_weight)
        calibrated_classifiers.append(calibrator)

    pipeline = _CalibratedClassiferPipeline(
        clf_fitted, calibrated_classifiers, label_encoder_
    )
    return pipeline


class _CalibratedClassiferPipeline:
    """Pipeline-like chaining a fitted classifier and its fitted calibrators.

    Parameters
    ----------
    clf_fitted : Estimator instance
        Fitted classifier.

    calibrators_fitted : List of fitted estimator instances
        List of fitted calibrators (either 'IsotonicRegression' or
        '_SigmoidCalibration'). The number of calibrators equals the number of
        classes. However, if there are 2 classes, the list contains only one
        fitted calibrator.
    """
    def __init__(self, clf_fitted, calibrators_fitted, label_encoder_):
        self.clf_fitted = clf_fitted
        self.calibrators_fitted = calibrators_fitted
        self.label_encoder_ = label_encoder_

    def predict_proba(self, X):
        """Calculate calibrated probabilities.

        Calculates classification calibrated probabilities
        for each class, in a one-vs-all manner, for `X`.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The sample data.

        Returns
        -------
        proba : array, shape (n_samples, n_classes)
            The predicted probabilities. Can be exact zeros.
        """
        preds = _get_predictions(self.clf_fitted, X, self.label_encoder_)
        pos_class_indices = self.label_encoder_.transform(
            self.clf_fitted.classes_
        )

        n_classes = len(self.label_encoder_.classes_)
        proba = np.zeros((X.shape[0], n_classes))
        for class_idx, this_pred, calibrator in \
                zip(pos_class_indices, preds.T, self.calibrators_fitted):
            if n_classes == 2:
                # When binary, `preds` consists only of predictions for
                # clf_fitted.classes_[1] but `pos_class_indices` = 0
                class_idx += 1
            proba[:, class_idx] = calibrator.predict(this_pred)

        # Normalize the probabilities
        if n_classes == 2:
            proba[:, 0] = 1. - proba[:, 1]
        else:
            proba /= np.sum(proba, axis=1)[:, np.newaxis]

        # XXX : for some reason all probas can be 0
        proba[np.isnan(proba)] = 1. / n_classes

        # Deal with cases where the predicted probability minimally exceeds 1.0
        proba[(1.0 < proba) & (proba <= 1.0 + 1e-5)] = 1.0

        return proba


def _sigmoid_calibration(pred, y, sample_weight=None):
    """Probability Calibration with sigmoid method (Platt 2000)

    Parameters
    ----------
    pred : ndarray, shape (n_samples,)
        The decision function or predict proba for the samples.

    y : ndarray, shape (n_samples,)
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
    pred = column_or_1d(pred)
    y = column_or_1d(y)

    F = pred  # F follows Platt's notations

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
        X : array-like, shape (n_samples,)
            Training data.

        y : array-like, shape (n_samples,)
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
        T : array-like, shape (n_samples,)
            Data to predict from.

        Returns
        -------
        T_ : array, shape (n_samples,)
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
