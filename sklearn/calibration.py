"""Calibration of predicted probabilities."""

# Author: Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#         Balazs Kegl <balazs.kegl@gmail.com>
#         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#         Mathieu Blondel <mathieu@mblondel.org>
#
# License: BSD 3 clause

import warnings
from inspect import signature

from math import log
import numpy as np

from scipy.special import expit
from scipy.special import xlogy
from scipy.optimize import fmin_bfgs
from .preprocessing import LabelEncoder

from .base import (BaseEstimator, ClassifierMixin, RegressorMixin, clone,
                   MetaEstimatorMixin)
from .preprocessing import label_binarize, LabelBinarizer
from .utils import check_array, indexable, column_or_1d
from .utils.validation import check_is_fitted, check_consistent_length
from .utils.validation import _check_sample_weight
from .isotonic import IsotonicRegression
from .svm import LinearSVC
from .model_selection import check_cv
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
        X, y = self._validate_data(X, y, accept_sparse=['csc', 'csr', 'coo'],
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
        Y = label_binarize(y, classes=self.classes_)

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
