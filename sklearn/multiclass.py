"""
Multiclass and multilabel classification strategies
===================================================

This module implements multiclass learning algorithms:
    - one-vs-the-rest / one-vs-all
    - one-vs-one
    - error correcting output codes

The estimators provided in this module are meta-estimators: they require a base
estimator to be provided in their constructor. For example, it is possible to
use these estimators to turn a binary classifier or a regressor into a
multiclass classifier. It is also possible to use these estimators with
multiclass estimators in the hope that their accuracy or runtime performance
improves.

All classifiers in scikit-learn implement multiclass classification; you
only need to use this module if you want to experiment with custom multiclass
strategies.

The one-vs-the-rest meta-classifier also implements a `predict_proba` method,
so long as such a method is implemented by the base classifier. This method
returns probabilities of class membership in both the single label and
multilabel case.  Note that in the multilabel case, probabilities are the
marginal probability that a given sample falls in the given class. As such, in
the multilabel case the sum of these probabilities over all possible labels
for a given sample *will not* sum to unity, as they do in the single label
case.
"""

# Author: Mathieu Blondel <mathieu@mblondel.org>
#
# License: BSD 3 clause

import numpy as np
import warnings

from .base import BaseEstimator, ClassifierMixin, clone, is_classifier
from .base import MetaEstimatorMixin
from .preprocessing import LabelBinarizer
from .metrics.pairwise import euclidean_distances
from .utils import check_random_state
from .externals.joblib import Parallel
from .externals.joblib import delayed


def _fit_binary(estimator, X, y, classes=None):
    """Fit a single binary estimator."""
    unique_y = np.unique(y)
    if len(unique_y) == 1:
        if classes is not None:
            if y[0] == -1:
                c = 0
            else:
                c = y[0]
            warnings.warn("Label %s is present in all training examples." %
                          str(classes[c]))
        estimator = _ConstantPredictor().fit(X, unique_y)
    else:
        estimator = clone(estimator)
        estimator.fit(X, y)
    return estimator


def _predict_binary(estimator, X):
    """Make predictions using a single binary estimator."""
    try:
        score = np.ravel(estimator.decision_function(X))
    except (AttributeError, NotImplementedError):
        # probabilities of the positive class
        score = estimator.predict_proba(X)[:, 1]
    return score


def _check_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""
    if (not hasattr(estimator, "decision_function") and
            not hasattr(estimator, "predict_proba")):
        raise ValueError("The base estimator should implement "
                         "decision_function or predict_proba!")


def fit_ovr(estimator, X, y, n_jobs=1):
    """Fit a one-vs-the-rest strategy."""
    _check_estimator(estimator)

    lb = LabelBinarizer()
    Y = lb.fit_transform(y)

    estimators = Parallel(n_jobs=n_jobs)(
        delayed(_fit_binary)(estimator, X, Y[:, i], classes=["not %s" % i, i])
        for i in range(Y.shape[1]))
    return estimators, lb


def predict_ovr(estimators, label_binarizer, X):
    """Make predictions using the one-vs-the-rest strategy."""
    Y = np.array([_predict_binary(e, X) for e in estimators])
    e = estimators[0]
    thresh = 0 if hasattr(e, "decision_function") and is_classifier(e) else .5
    return label_binarizer.inverse_transform(Y.T, threshold=thresh)


def predict_proba_ovr(estimators, X, is_multilabel):
    """Estimate probabilities using the one-vs-the-rest strategy.

    If multilabel is true, returned matrix will not sum to one.  Estimators
    must have a predict_proba method."""

    # Y[i,j] gives the probability that sample i has the label j.
    # In the multi-label case, these are not disjoint.
    Y = np.array([est.predict_proba(X)[:, 1] for est in estimators]).T

    if not is_multilabel:
        # Then, probabilities should be normalized to 1.
        Y /= np.sum(Y, axis=1)[:, np.newaxis]
    return Y


class _ConstantPredictor(BaseEstimator):
    def fit(self, X, y):
        self.y_ = y
        return self

    def predict(self, X):
        return np.repeat(self.y_, X.shape[0])

    def decision_function(self, X):
        return np.repeat(self.y_, X.shape[0])


class OneVsRestClassifier(BaseEstimator, ClassifierMixin, MetaEstimatorMixin):
    """One-vs-the-rest (OvR) multiclass/multilabel strategy

    Also known as one-vs-all, this strategy consists in fitting one classifier
    per class. For each classifier, the class is fitted against all the other
    classes. In addition to its computational efficiency (only `n_classes`
    classifiers are needed), one advantage of this approach is its
    interpretability. Since each class is represented by one and one classifier
    only, it is possible to gain knowledge about the class by inspecting its
    corresponding classifier. This is the most commonly used strategy for
    multiclass classification and is a fair default choice.

    This strategy can also be used for multilabel learning, where a classifier
    is used to predict multiple labels for instance, by fitting on a sequence
    of sequences of labels (e.g., a list of tuples) rather than a single
    target vector. For multilabel learning, the number of classes must be at
    least three, since otherwise OvR reduces to binary classification.

    In the multilabel learning literature, OvR is also known as the binary
    relevance method.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit` and one of `decision_function`
        or `predict_proba`.

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. If -1 all CPUs are used.
        If 1 is given, no parallel computing code is used at all, which is
        useful for debugging. For n_jobs below -1, (n_cpus + 1 + n_jobs) are
        used. Thus for n_jobs = -2, all CPUs but one are used.

    Attributes
    ----------
    `estimators_` : list of `n_classes` estimators
        Estimators used for predictions.

    `classes_` : array, shape = [`n_classes`]
        Class labels.
    `label_binarizer_` : LabelBinarizer object
        Object used to transform multiclass labels to binary labels and
        vice-versa.
    `multilabel_` : boolean
        Whether a OneVsRestClassifier is a multilabel classifier.
    """

    def __init__(self, estimator, n_jobs=1):
        self.estimator = estimator
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data.

        y : array-like, shape = [n_samples]
         or sequence of sequences, len = n_samples
            Multi-class targets. A sequence of sequences turns on multilabel
            classification.

        Returns
        -------
        self
        """
        self.estimators_, self.label_binarizer_ = fit_ovr(self.estimator, X, y,
                                                          n_jobs=self.n_jobs)
        return self

    def _check_is_fitted(self):
        if not hasattr(self, "estimators_"):
            raise ValueError("The object hasn't been fitted yet!")

    def predict(self, X):
        """Predict multi-class targets using underlying estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : array-like, shape = [n_samples]
            Predicted multi-class targets.
        """
        self._check_is_fitted()

        return predict_ovr(self.estimators_, self.label_binarizer_, X)

    def predict_proba(self, X):
        """Probability estimates.

        The returned estimates for all classes are ordered by label of classes.

        Note that in the multilabel case, each sample can have any number of
        labels. This returns the marginal probability that the given sample has
        the label in question. For example, it is entirely consistent that two
        labels both have a 90% probability of applying to a given sample.

        In the single label multiclass case, the rows of the returned matrix
        sum to 1.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in `self.classes_`.
        """
        return predict_proba_ovr(self.estimators_, X,
                                 is_multilabel=self.multilabel_)

    def decision_function(self, X):
        """Returns the distance of each sample from the decision boundary for
        each class. This can only be used with estimators which implement the
        decision_function method.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = [n_samples, n_classes]
        """
        if not hasattr(self.estimators_[0], "decision_function"):
            raise AttributeError(
                "Base estimator doesn't have a decision_function attribute.")
        return np.array([est.decision_function(X).ravel()
                         for est in self.estimators_]).T

    @property
    def multilabel_(self):
        """Whether this is a multilabel classifier"""
        return self.label_binarizer_.multilabel_

    def score(self, X, y):
        if self.multilabel_:
            raise NotImplementedError(
                "score is not supported for multilabel classifiers")
        else:
            return super(OneVsRestClassifier, self).score(X, y)

    @property
    def classes_(self):
        return self.label_binarizer_.classes_

    @property
    def coef_(self):
        self._check_is_fitted()
        if not hasattr(self.estimators_[0], "coef_"):
            raise AttributeError(
                "Base estimator doesn't have a coef_ attribute.")
        return np.array([e.coef_.ravel() for e in self.estimators_])

    @property
    def intercept_(self):
        self._check_is_fitted()
        if not hasattr(self.estimators_[0], "intercept_"):
            raise AttributeError(
                "Base estimator doesn't have an intercept_ attribute.")
        return np.array([e.intercept_.ravel() for e in self.estimators_])


def _fit_ovo_binary(estimator, X, y, i, j):
    """Fit a single binary estimator (one-vs-one)."""
    cond = np.logical_or(y == i, y == j)
    y = y[cond]
    y[y == i] = 0
    y[y == j] = 1
    ind = np.arange(X.shape[0])
    return _fit_binary(estimator, X[ind[cond]], y, classes=[i, j])


def fit_ovo(estimator, X, y, n_jobs=1):
    """Fit a one-vs-one strategy."""
    classes = np.unique(y)
    n_classes = classes.shape[0]
    estimators = Parallel(n_jobs=n_jobs)(
        delayed(_fit_ovo_binary)(
            estimator, X, y, classes[i], classes[j])
        for i in range(n_classes) for j in range(i + 1, n_classes))

    return estimators, classes


def predict_ovo(estimators, classes, X):
    """Make predictions using the one-vs-one strategy."""
    n_samples = X.shape[0]
    n_classes = classes.shape[0]
    votes = np.zeros((n_samples, n_classes))
    scores = np.zeros((n_samples, n_classes))

    k = 0
    for i in range(n_classes):
        for j in range(i + 1, n_classes):
            pred = estimators[k].predict(X)
            score = _predict_binary(estimators[k], X)
            scores[:, i] += score
            scores[:, j] -= score
            votes[pred == 0, i] += 1
            votes[pred == 1, j] += 1
            k += 1

    # find all places with maximum votes per sample
    maxima = votes == np.max(votes, axis=1)[:, np.newaxis]

    # if there are ties, use scores to break them
    if np.any(maxima.sum(axis=1) > 1):
        scores[~maxima] = -np.inf
        prediction = scores.argmax(axis=1)
    else:
        prediction = votes.argmax(axis=1)
    return classes[prediction]


class OneVsOneClassifier(BaseEstimator, ClassifierMixin, MetaEstimatorMixin):
    """One-vs-one multiclass strategy

    This strategy consists in fitting one classifier per class pair.
    At prediction time, the class which received the most votes is selected.
    Since it requires to fit `n_classes * (n_classes - 1) / 2` classifiers,
    this method is usually slower than one-vs-the-rest, due to its
    O(n_classes^2) complexity. However, this method may be advantageous for
    algorithms such as kernel algorithms which don't scale well with
    `n_samples`. This is because each individual learning problem only involves
    a small subset of the data whereas, with one-vs-the-rest, the complete
    dataset is used `n_classes` times.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit` and `predict`.

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. If -1 all CPUs are used.
        If 1 is given, no parallel computing code is used at all, which is
        useful for debugging. For n_jobs below -1, (n_cpus + 1 + n_jobs) are
        used. Thus for n_jobs = -2, all CPUs but one are used.

    Attributes
    ----------
    `estimators_` : list of `n_classes * (n_classes - 1) / 2` estimators
        Estimators used for predictions.

    `classes_` : numpy array of shape [n_classes]
        Array containing labels.
    """

    def __init__(self, estimator, n_jobs=1):
        self.estimator = estimator
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data.

        y : numpy array of shape [n_samples]
            Multi-class targets.

        Returns
        -------
        self
        """
        self.estimators_, self.classes_ = fit_ovo(self.estimator, X, y,
                                                  self.n_jobs)
        return self

    def predict(self, X):
        """Predict multi-class targets using underlying estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : numpy array of shape [n_samples]
            Predicted multi-class targets.
        """
        if not hasattr(self, "estimators_"):
            raise ValueError("The object hasn't been fitted yet!")

        return predict_ovo(self.estimators_, self.classes_, X)


def fit_ecoc(estimator, X, y, code_size=1.5, random_state=None, n_jobs=1):
    """
    Fit an error-correcting output-code strategy.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit` and one of `decision_function`
        or `predict_proba`.

    code_size : float, optional
        Percentage of the number of classes to be used to create the code book.

    random_state : numpy.RandomState, optional
        The generator used to initialize the codebook. Defaults to
        numpy.random.


    Returns
    --------
    estimators : list of `int(n_classes * code_size)` estimators
        Estimators used for predictions.

    classes : numpy array of shape [n_classes]
        Array containing labels.

    `code_book_`: numpy array of shape [n_classes, code_size]
        Binary array containing the code of each class.
    """
    _check_estimator(estimator)
    random_state = check_random_state(random_state)

    classes = np.unique(y)
    n_classes = classes.shape[0]
    code_size = int(n_classes * code_size)

    # FIXME: there are more elaborate methods than generating the codebook
    # randomly.
    code_book = random_state.random_sample((n_classes, code_size))
    code_book[code_book > 0.5] = 1

    if hasattr(estimator, "decision_function"):
        code_book[code_book != 1] = -1
    else:
        code_book[code_book != 1] = 0

    cls_idx = dict((c, i) for i, c in enumerate(classes))

    Y = np.array([code_book[cls_idx[y[i]]] for i in range(X.shape[0])],
                 dtype=np.int)

    estimators = Parallel(n_jobs=n_jobs)(
        delayed(_fit_binary)(estimator, X, Y[:, i])
        for i in range(Y.shape[1]))

    return estimators, classes, code_book


def predict_ecoc(estimators, classes, code_book, X):
    """Make predictions using the error-correcting output-code strategy."""
    Y = np.array([_predict_binary(e, X) for e in estimators]).T
    pred = euclidean_distances(Y, code_book).argmin(axis=1)
    return classes[pred]


class OutputCodeClassifier(BaseEstimator, ClassifierMixin, MetaEstimatorMixin):
    """(Error-Correcting) Output-Code multiclass strategy

    Output-code based strategies consist in representing each class with a
    binary code (an array of 0s and 1s). At fitting time, one binary
    classifier per bit in the code book is fitted.  At prediction time, the
    classifiers are used to project new points in the class space and the class
    closest to the points is chosen. The main advantage of these strategies is
    that the number of classifiers used can be controlled by the user, either
    for compressing the model (0 < code_size < 1) or for making the model more
    robust to errors (code_size > 1). See the documentation for more details.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing `fit` and one of `decision_function`
        or `predict_proba`.

    code_size : float
        Percentage of the number of classes to be used to create the code book.
        A number between 0 and 1 will require fewer classifiers than
        one-vs-the-rest. A number greater than 1 will require more classifiers
        than one-vs-the-rest.

    random_state : numpy.RandomState, optional
        The generator used to initialize the codebook. Defaults to
        numpy.random.

    n_jobs : int, optional, default: 1
        The number of jobs to use for the computation. If -1 all CPUs are used.
        If 1 is given, no parallel computing code is used at all, which is
        useful for debugging. For n_jobs below -1, (n_cpus + 1 + n_jobs) are
        used. Thus for n_jobs = -2, all CPUs but one are used.

    Attributes
    ----------
    `estimators_` : list of `int(n_classes * code_size)` estimators
        Estimators used for predictions.

    `classes_` : numpy array of shape [n_classes]
        Array containing labels.

    `code_book_` : numpy array of shape [n_classes, code_size]
        Binary array containing the code of each class.

    References
    ----------

    .. [1] "Solving multiclass learning problems via error-correcting output
       codes",
       Dietterich T., Bakiri G.,
       Journal of Artificial Intelligence Research 2,
       1995.

    .. [2] "The error coding method and PICTs",
       James G., Hastie T.,
       Journal of Computational and Graphical statistics 7,
       1998.

    .. [3] "The Elements of Statistical Learning",
       Hastie T., Tibshirani R., Friedman J., page 606 (second-edition)
       2008.
    """

    def __init__(self, estimator, code_size=1.5, random_state=None, n_jobs=1):
        if (code_size <= 0):
            raise ValueError("code_size should be greater than 0!")

        self.estimator = estimator
        self.code_size = code_size
        self.random_state = random_state
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data.

        y : numpy array of shape [n_samples]
            Multi-class targets.

        Returns
        -------
        self
        """
        self.estimators_, self.classes_, self.code_book_ = \
            fit_ecoc(self.estimator, X, y, self.code_size, self.random_state,
                     self.n_jobs)
        return self

    def predict(self, X):
        """Predict multi-class targets using underlying estimators.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : numpy array of shape [n_samples]
            Predicted multi-class targets.
        """
        if not hasattr(self, "estimators_"):
            raise ValueError("The object hasn't been fitted yet!")

        return predict_ecoc(self.estimators_, self.classes_,
                            self.code_book_, X)
