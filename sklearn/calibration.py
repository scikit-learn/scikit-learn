"""Calibration of predicted probabilities."""

# Author: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#         Balazs Kegl <balazs.kegl@gmail.com>
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         Mathieu Blondel <mathieu@mblondel.org>
#
# License: BSD 3 clause

from copy import deepcopy
from inspect import signature
import numbers
import warnings

from math import log
import numpy as np

from scipy.special import expit
from scipy.special import xlogy
from scipy.optimize import fmin_bfgs
from .preprocessing import LabelEncoder

from .base import BaseEstimator
from .base import ClassifierMixin
from .base import MetaEstimatorMixin
from .base import RegressorMixin
from .base import clone
from .exceptions import NotFittedError
from .isotonic import IsotonicRegression
from .metrics import precision_recall_curve
from .metrics import roc_curve
from .model_selection import check_cv
from .preprocessing import label_binarize
from .preprocessing import LabelBinarizer
from .svm import LinearSVC
from .utils import check_X_y
from .utils import check_array
from .utils import column_or_1d
from .utils import indexable
from .utils import _safe_indexing
from .utils.multiclass import check_classification_targets
from .utils.multiclass import type_of_target
from .utils.validation import check_is_fitted
from .utils.validation import check_consistent_length
from .utils.validation import _check_sample_weight
from .utils.validation import _deprecate_positional_args


class CalibratedClassifierCV(BaseEstimator, ClassifierMixin,
                             MetaEstimatorMixin):
    """Probability calibration with isotonic regression or logistic regression.

    The calibration is based on the :term:`decision_function` method of the
    `base_estimator` if it exists, else on :term:`predict_proba`.

    Read more in the :ref:`User Guide <calibration>`.

    Parameters
    ----------
    base_estimator : instance BaseEstimator
        The classifier whose output need to be calibrated to provide more
        accurate `predict_proba` outputs.

    method : 'sigmoid' or 'isotonic'
        The method to use for calibration. Can be 'sigmoid' which
        corresponds to Platt's method (i.e. a logistic regression model) or
        'isotonic' which is a non-parametric approach. It is not advised to
        use isotonic calibration with too few calibration samples
        ``(<<1000)`` since it tends to overfit.

    cv : integer, cross-validation generator, iterable or "prefit", optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`sklearn.model_selection.StratifiedKFold` is used. If ``y`` is
        neither binary nor multiclass, :class:`sklearn.model_selection.KFold`
        is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        If "prefit" is passed, it is assumed that `base_estimator` has been
        fitted already and all data is used for calibration.

        .. versionchanged:: 0.22
            ``cv`` default value if None changed from 3-fold to 5-fold.

    Attributes
    ----------
    classes_ : array, shape (n_classes)
        The class labels.

    calibrated_classifiers_ : list (len() equal to cv or 1 if cv == "prefit")
        The list of calibrated classifiers, one for each cross-validation fold,
        which has been fitted on all but the validation fold and calibrated
        on the validation fold.

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
    def __init__(self, base_estimator=None, *, method='sigmoid', cv=None):
        self.base_estimator = base_estimator
        self.method = method
        self.cv = cv

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
        X, y = check_X_y(X, y, accept_sparse=['csc', 'csr', 'coo'],
                         force_all_finite=False, allow_nd=True)
        X, y = indexable(X, y)
        le = LabelBinarizer().fit(y)
        self.classes_ = le.classes_

        # Check that each cross-validation fold can have at least one
        # example per class
        n_folds = self.cv if isinstance(self.cv, int) \
            else self.cv.n_folds if hasattr(self.cv, "n_folds") else None
        if n_folds and \
                np.any([np.sum(y == class_) < n_folds for class_ in
                        self.classes_]):
            raise ValueError("Requesting %d-fold cross-validation but provided"
                             " less than %d examples for at least one class."
                             % (n_folds, n_folds))

        self.calibrated_classifiers_ = []
        if self.base_estimator is None:
            # we want all classifiers that don't expose a random_state
            # to be deterministic (and we don't want to expose this one).
            base_estimator = LinearSVC(random_state=0)
        else:
            base_estimator = self.base_estimator

        if self.cv == "prefit":
            calibrated_classifier = _CalibratedClassifier(
                base_estimator, method=self.method)
            calibrated_classifier.fit(X, y, sample_weight)
            self.calibrated_classifiers_.append(calibrated_classifier)
        else:
            cv = check_cv(self.cv, y, classifier=True)
            fit_parameters = signature(base_estimator.fit).parameters
            base_estimator_supports_sw = "sample_weight" in fit_parameters

            if sample_weight is not None:
                sample_weight = _check_sample_weight(sample_weight, X)

                if not base_estimator_supports_sw:
                    estimator_name = type(base_estimator).__name__
                    warnings.warn("Since %s does not support sample_weights, "
                                  "sample weights will only be used for the "
                                  "calibration itself." % estimator_name)

            for train, test in cv.split(X, y):
                this_estimator = clone(base_estimator)

                if sample_weight is not None and base_estimator_supports_sw:
                    this_estimator.fit(X[train], y[train],
                                       sample_weight=sample_weight[train])
                else:
                    this_estimator.fit(X[train], y[train])

                calibrated_classifier = _CalibratedClassifier(
                    this_estimator, method=self.method, classes=self.classes_)
                sw = None if sample_weight is None else sample_weight[test]
                calibrated_classifier.fit(X[test], y[test], sample_weight=sw)
                self.calibrated_classifiers_.append(calibrated_classifier)

        return self

    def predict_proba(self, X):
        """Posterior probabilities of classification

        This function returns posterior probabilities of classification
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


class _CalibratedClassifier:
    """Probability calibration with isotonic regression or sigmoid.

    It assumes that base_estimator has already been fit, and trains the
    calibration on the input set of the fit function. Note that this class
    should not be used as an estimator directly. Use CalibratedClassifierCV
    with cv="prefit" instead.

    Parameters
    ----------
    base_estimator : instance BaseEstimator
        The classifier whose output decision function needs to be calibrated
        to offer more accurate predict_proba outputs. No default value since
        it has to be an already fitted estimator.

    method : 'sigmoid' | 'isotonic'
        The method to use for calibration. Can be 'sigmoid' which
        corresponds to Platt's method or 'isotonic' which is a
        non-parametric approach based on isotonic regression.

    classes : array-like, shape (n_classes,), optional
            Contains unique classes used to fit the base estimator.
            if None, then classes is extracted from the given target values
            in fit().

    See also
    --------
    CalibratedClassifierCV

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
    def __init__(self, base_estimator, *, method='sigmoid', classes=None):
        self.base_estimator = base_estimator
        self.method = method
        self.classes = classes

    def _preproc(self, X):
        n_classes = len(self.classes_)
        if hasattr(self.base_estimator, "decision_function"):
            df = self.base_estimator.decision_function(X)
            if df.ndim == 1:
                df = df[:, np.newaxis]
        elif hasattr(self.base_estimator, "predict_proba"):
            df = self.base_estimator.predict_proba(X)
            if n_classes == 2:
                df = df[:, 1:]
        else:
            raise RuntimeError('classifier has no decision_function or '
                               'predict_proba method.')

        idx_pos_class = self.label_encoder_.\
            transform(self.base_estimator.classes_)

        return df, idx_pos_class

    def fit(self, X, y, sample_weight=None):
        """Calibrate the fitted model

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

        self.label_encoder_ = LabelEncoder()
        if self.classes is None:
            self.label_encoder_.fit(y)
        else:
            self.label_encoder_.fit(self.classes)

        self.classes_ = self.label_encoder_.classes_
        Y = label_binarize(y, self.classes_)

        df, idx_pos_class = self._preproc(X)
        self.calibrators_ = []

        for k, this_df in zip(idx_pos_class, df.T):
            if self.method == 'isotonic':
                calibrator = IsotonicRegression(out_of_bounds='clip')
            elif self.method == 'sigmoid':
                calibrator = _SigmoidCalibration()
            else:
                raise ValueError('method should be "sigmoid" or '
                                 '"isotonic". Got %s.' % self.method)
            calibrator.fit(this_df, Y[:, k], sample_weight)
            self.calibrators_.append(calibrator)

        return self

    def predict_proba(self, X):
        """Posterior probabilities of classification

        This function returns posterior probabilities of classification
        according to each class on an array of test vectors X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The samples.

        Returns
        -------
        C : array, shape (n_samples, n_classes)
            The predicted probas. Can be exact zeros.
        """
        n_classes = len(self.classes_)
        proba = np.zeros((X.shape[0], n_classes))

        df, idx_pos_class = self._preproc(X)

        for k, this_df, calibrator in \
                zip(idx_pos_class, df.T, self.calibrators_):
            if n_classes == 2:
                k += 1
            proba[:, k] = calibrator.predict(this_df)

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


def _sigmoid_calibration(df, y, sample_weight=None):
    """Probability Calibration with sigmoid method (Platt 2000)

    Parameters
    ----------
    df : ndarray, shape (n_samples,)
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
    df = column_or_1d(df)
    y = column_or_1d(y)

    F = df  # F follows Platt's notations

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


def calibration_curve(y_true, y_prob, normalize=False, n_bins=5,
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
    y_true = label_binarize(y_true, labels)[:, 0]

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


class CutoffClassifier(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
    """Decision threshold calibration for binary classification.

    Estimator that calibrates the decision threshold (cutoff point) that is
    used for prediction. The methods for picking cutoff points make use of
    traditional binary classification evaluation statistics such as the true
    positive and true negative rates and F-scores.

    If cv="prefit" the base estimator is assumed to be fitted and all data will
    be used for the selection of the cutoff point. Otherwise the decision
    threshold is calculated as the average of the thresholds resulting from the
    cross-validation loop.

    Parameters
    ----------
    base_estimator : estimator instance
        The binary classifier whose decision threshold will be adapted
        according to the acquired cutoff point. The estimator must implement
        `decision_function` or `predict_proba` function.

    strategy : {"roc", "f_beta", "max_tpr", "max_tnr", "constant"}, \
            default="roc" The strategy to use for choosing the cutoff point.

        - "roc" selects the point on the ROC curve that is closest to the ideal
          corner (0, 1).

        - "f_beta" selects a decision threshold that maximizes the `f_beta`
          score.

        - "max_tpr" selects the point that yields the highest true positive
          rate (TPR) with true negative rate (TNR) at least equal to the value
          of the parameter threshold.

        - "max_tnr" selects the point that yields the highest true negative
          rate (TNR) with true positive rate (TPR) at least equal to the value
          of the parameter threshold.

        - "constant" will use the threshold specified by the parameter
          `decision_threshold`.

    method : {"auto", "decision_function", "predict_proba"}, default="auto" The
        method to be used to get the predictions. If `"auto"` (default), the
        base estimator will try to invoke `decision_function` or
        `predict_proba`, in that order.

    beta : float in [0, 1], optional (default=None) beta value to be used in
        case strategy == 'f_beta'

    threshold : float in [0, 1] or None, (default=None) In case strategy is
        'max_tpr' or 'max_tnr' this parameter must be set to specify the
        threshold for the true negative rate or true positive rate respectively
        that needs to be achieved

    decision_threshold : float, default=0.5
        When `strategy="constant"`, decision threshold used as cutoff point.

    pos_label : object, optional (default=1) Object representing the positive
        label

    cv : int, cross-validation generator, iterable or 'prefit', optional
        (default=3). Determines the cross-validation splitting strategy. If
        cv='prefit' the base estimator is assumed to be fitted and all data
        will be used for the calibration of the probability threshold

    Attributes
    ----------
    decision_threshold_ : float Decision threshold for the positive class.
        Determines the output of predict

    std_ : float Standard deviation of the obtained decision thresholds for
        when the provided base estimator is not pre-trained and the
        decision_threshold_ is computed as the mean of the decision threshold
        of each cross-validation iteration. If the base estimator is
        pre-trained then std_ = None

    classes_ : array, shape (n_classes) The class labels.

    References
    ----------
    .. [1] Receiver-operating characteristic (ROC) plots: a fundamental
           evaluation tool in clinical medicine, MH Zweig, G Campbell -
           Clinical chemistry, 1993

    """
    def __init__(self, base_estimator, strategy="roc", method="auto",
                 beta=None, threshold=None, decision_threshold=0.5,
                 pos_label=1, cv=3):
        self.base_estimator = base_estimator
        self.strategy = strategy
        self.method = method
        self.beta = beta
        self.threshold = threshold
        self.decision_threshold = decision_threshold
        self.pos_label = pos_label
        self.cv = cv

    def _validate_parameters(self):
        """Validate the input parameters."""
        supported_methods = ("decision_function", "predict_proba")
        if self.method == "auto":
            has_methods = [
                hasattr(self.base_estimator, m) for m in supported_methods
            ]
            if not any(has_methods):
                raise TypeError(
                    f"'base_estimator' must implement one of the "
                    f"{', '.join(supported_methods)} methods."
                )
            self._method = next(
                (m for m, i in zip(supported_methods, has_methods) if i), None
            )
        else:
            if self.method not in supported_methods:
                raise ValueError(
                    f"'method' should be one of {', '.join(supported_methods)}"
                    f". Got {self.method} instead."
                )
            elif not hasattr(self.base_estimator, self.method):
                raise TypeError(
                    f"'base_estimator' does not implement {self.method}."
                )
            self._method = self.method

        strategies = ("roc", "f_beta", "max_tpr", "max_tnr", "constant")
        if self.strategy not in strategies:
            raise ValueError(
                f"'strategy' must be of {', '.join(strategies)}. "
                f"Got {self.strategy} instead."
            )
        elif self.strategy in ("max_tpr", "max_tnr"):
            if not isinstance(self.threshold, numbers.Real):
                raise TypeError(
                    "When strategy is max_tpr or max_tnr, threshold should be "
                    f"a real in [0, 1]. Got {type(self.threshold)} instead."
                )
            elif not (0 < self.threshold < 1):
                raise ValueError(
                    f"threshold should be in the range [0, 1]. "
                    f"Got {self.threshold} instead."
                )
        elif self.strategy == "f_beta":
            if not isinstance(self.beta, numbers.Real):
                raise TypeError(
                    "When strategy is f_beta, beta should be a real. "
                    f"Got {type(self.beta)} instead."
                )
        elif self.strategy == "constant":
            if (self.method == "predict_proba" and
                    not (0 < self.decision_threshold < 1)):
                raise ValueError(
                    f"decision_threshold should be in the range [0, 1] when "
                    f"using 'predict_proba'. Got {self.decision_threshold} "
                    "instead."
                )

    def _validate_data(self, X, y):
        X = check_array(
            X, accept_sparse=['csc', 'csr'], force_all_finite=False,
            allow_nd=True,
        )
        y = check_array(y, ensure_2d=False, dtype=None)
        check_classification_targets(y)
        y_type = type_of_target(y)
        if y_type != 'binary':
            raise ValueError(f'Expected target of binary type. Got {y_type}.')
        return X, y

    def fit(self, X, y):
        """Fit model

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data

        y : array-like, shape (n_samples,)
            Target values. There must be two 2 distinct values

        Returns
        -------
        self : object
            Instance of self
        """
        self._validate_parameters()
        X, y = self._validate_data(X, y)

        self._label_encoder = LabelEncoder().fit(y)
        self.classes_ = self._label_encoder.classes_

        try:
            check_is_fitted(self.base_estimator)
            self._base_estimator = deepcopy(self.base_estimator)
        except NotFittedError:
            self._base_estimator = clone(self.base_estimator).fit(X, y)

        if self.strategy == "constant":
            self.decision_threshold_ = self.decision_threshold
        elif self.cv == 'prefit':
            self.decision_threshold_ = _find_optimal_decision_threshold(
                self._base_estimator, X, y, self.strategy, self._method,
                self.beta, self.threshold, self.pos_label
            )
            self.std_ = None
        else:
            cv = check_cv(self.cv, y, classifier=True)
            decision_thresholds = []

            for train, test in cv.split(X, y):
                estimator = clone(self._base_estimator).fit(
                    _safe_indexing(X, train), _safe_indexing(y, train)
                )
                decision_thresholds.append(
                    _find_optimal_decision_threshold(
                        estimator,
                        _safe_indexing(X, test), _safe_indexing(y, test),
                        self.strategy, self._method, self.beta, self.threshold,
                        self.pos_label
                    )
                )
            self.decision_threshold_ = np.mean(decision_thresholds)
            self.std_ = np.std(decision_thresholds)
        return self

    def predict(self, X):
        """Predict using the calibrated decision threshold

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The samples

        Returns
        -------
        C : array, shape (n_samples,)
            The predicted class
        """
        check_is_fitted(self)

        y_score = _get_binary_score(
            self._base_estimator, X, self._method, self.pos_label
        )
        return self._label_encoder.inverse_transform(
            (y_score > self.decision_threshold_).astype(int)
        )

    def _more_tags(self):
        return {"binary_only": True}


def _find_optimal_decision_threshold(estimator, X, y, strategy, method, beta,
                                     threshold, pos_label):
    y_score = _get_binary_score(
        estimator, X, method=method, pos_label=pos_label
    )
    if strategy == 'f_beta':
        precision, recall, thresholds = precision_recall_curve(
            y, y_score, pos_label=pos_label
        )
        f_beta = ((1 + beta ** 2) * (precision * recall) /
                  (beta ** 2 * precision + recall))
        return thresholds[np.argmax(f_beta)]

    fpr, tpr, thresholds = roc_curve(y, y_score, pos_label=pos_label)

    if strategy == 'roc':
        # we find the threshold of the point (fpr, tpr) with the smallest
        # euclidean distance from the "ideal" corner (0, 1)
        return thresholds[np.argmin(fpr ** 2 + (tpr - 1) ** 2)]
    elif strategy == 'max_tpr':
        indices = np.where(1 - fpr >= threshold)[0]
        max_tpr_index = np.argmax(tpr[indices])
        return thresholds[indices[max_tpr_index]]
    indices = np.where(tpr >= threshold)[0]
    max_tnr_index = np.argmax(1 - fpr[indices])
    return thresholds[indices[max_tnr_index]]


def _get_binary_score(estimator, X, method, pos_label):
    """Binary classification score for the positive label (0 or 1)

    Returns the score that a binary classifier outputs for the positive label
    acquired either from decision_function or predict_proba

    Parameters
    ----------
    estimator : estimator object
        Fitted estimator to get prediction from.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        The data matrix.

    pos_label : int or str
        The positive label.

    method : str or None, optional (default=None)
        The method to be used for acquiring the score. Can either be
        "decision_function" or "predict_proba" or None. If None then
        decision_function will be used first and if not available
        predict_proba

    Returns
    -------
    y_score : array-like, shape (n_samples,)
        The return value of the provided classifier's decision_function or
        predict_proba depending on the method used.
    """
    # FIXME: what if estimator was fitted on encoded label??
    y_score = getattr(estimator, method)(X)
    if y_score.ndim == 2:
        # probabilities
        y_score = y_score[:, pos_label]
    elif pos_label == estimator.classes_[0]:
        y_score = -y_score
    return y_score
