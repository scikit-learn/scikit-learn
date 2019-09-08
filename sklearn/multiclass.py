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
# Author: Hamzeh Alsalhi <93hamsal@gmail.com>
#
# License: BSD 3 clause

import array
import numpy as np
import warnings
import scipy.sparse as sp
import itertools

from .base import BaseEstimator, ClassifierMixin, clone, is_classifier
from .base import MultiOutputMixin
from .base import MetaEstimatorMixin, is_regressor
from .preprocessing import LabelBinarizer
from .metrics.pairwise import euclidean_distances
from .utils import check_random_state
from .utils.validation import _num_samples
from .utils.validation import check_is_fitted
from .utils.validation import check_X_y, check_array
from .utils.multiclass import (_check_partial_fit_first_call,
                               check_classification_targets,
                               _ovr_decision_function)
from .utils.metaestimators import _safe_split, if_delegate_has_method

from joblib import Parallel, delayed

__all__ = [
    "OneVsRestClassifier",
    "OneVsOneClassifier",
    "OutputCodeClassifier",
]


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


def _partial_fit_binary(estimator, X, y):
    """Partially fit a single binary estimator."""
    estimator.partial_fit(X, y, np.array((0, 1)))
    return estimator


def _predict_binary(estimator, X):
    """Make predictions using a single binary estimator."""
    if is_regressor(estimator):
        return estimator.predict(X)
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


class _ConstantPredictor(BaseEstimator):

    def fit(self, X, y):
        self.y_ = y
        return self

    def predict(self, X):
        check_is_fitted(self)

        return np.repeat(self.y_, X.shape[0])

    def decision_function(self, X):
        check_is_fitted(self)

        return np.repeat(self.y_, X.shape[0])

    def predict_proba(self, X):
        check_is_fitted(self)

        return np.repeat([np.hstack([1 - self.y_, self.y_])],
                         X.shape[0], axis=0)


class OneVsRestClassifier(MultiOutputMixin, ClassifierMixin,
                          MetaEstimatorMixin, BaseEstimator):
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
    is used to predict multiple labels for instance, by fitting on a 2-d matrix
    in which cell [i, j] is 1 if sample i has label j and 0 otherwise.

    In the multilabel learning literature, OvR is also known as the binary
    relevance method.

    Read more in the :ref:`User Guide <ovr_classification>`.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing :term:`fit` and one of
        :term:`decision_function` or :term:`predict_proba`.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    estimators_ : list of `n_classes` estimators
        Estimators used for predictions.

    classes_ : array, shape = [`n_classes`]
        Class labels.
    label_binarizer_ : LabelBinarizer object
        Object used to transform multiclass labels to binary labels and
        vice-versa.
    multilabel_ : boolean
        Whether a OneVsRestClassifier is a multilabel classifier.
    """
    def __init__(self, estimator, n_jobs=None):
        self.estimator = estimator
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying estimators.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        y : (sparse) array-like, shape = [n_samples, ], [n_samples, n_classes]
            Multi-class targets. An indicator matrix turns on multilabel
            classification.

        Returns
        -------
        self
        """
        # A sparse LabelBinarizer, with sparse_output=True, has been shown to
        # outperform or match a dense label binarizer in all cases and has also
        # resulted in less or equal memory consumption in the fit_ovr function
        # overall.
        self.label_binarizer_ = LabelBinarizer(sparse_output=True)
        Y = self.label_binarizer_.fit_transform(y)
        Y = Y.tocsc()
        self.classes_ = self.label_binarizer_.classes_
        columns = (col.toarray().ravel() for col in Y.T)
        # In cases where individual estimators are very fast to train setting
        # n_jobs > 1 in can results in slower performance due to the overhead
        # of spawning threads.  See joblib issue #112.
        self.estimators_ = Parallel(n_jobs=self.n_jobs)(delayed(_fit_binary)(
            self.estimator, X, column, classes=[
                "not %s" % self.label_binarizer_.classes_[i],
                self.label_binarizer_.classes_[i]])
            for i, column in enumerate(columns))

        return self

    @if_delegate_has_method('estimator')
    def partial_fit(self, X, y, classes=None):
        """Partially fit underlying estimators

        Should be used when memory is inefficient to train all data.
        Chunks of data can be passed in several iteration.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        y : (sparse) array-like, shape = [n_samples, ], [n_samples, n_classes]
            Multi-class targets. An indicator matrix turns on multilabel
            classification.

        classes : array, shape (n_classes, )
            Classes across all calls to partial_fit.
            Can be obtained via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is only required in the first call of partial_fit
            and can be omitted in the subsequent calls.

        Returns
        -------
        self
        """
        if _check_partial_fit_first_call(self, classes):
            if not hasattr(self.estimator, "partial_fit"):
                raise ValueError(("Base estimator {0}, doesn't have "
                                 "partial_fit method").format(self.estimator))
            self.estimators_ = [clone(self.estimator) for _ in range
                                (self.n_classes_)]

            # A sparse LabelBinarizer, with sparse_output=True, has been
            # shown to outperform or match a dense label binarizer in all
            # cases and has also resulted in less or equal memory consumption
            # in the fit_ovr function overall.
            self.label_binarizer_ = LabelBinarizer(sparse_output=True)
            self.label_binarizer_.fit(self.classes_)

        if len(np.setdiff1d(y, self.classes_)):
            raise ValueError(("Mini-batch contains {0} while classes " +
                             "must be subset of {1}").format(np.unique(y),
                                                             self.classes_))

        Y = self.label_binarizer_.transform(y)
        Y = Y.tocsc()
        columns = (col.toarray().ravel() for col in Y.T)

        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_partial_fit_binary)(estimator, X, column)
            for estimator, column in zip(self.estimators_, columns))

        return self

    def predict(self, X):
        """Predict multi-class targets using underlying estimators.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : (sparse) array-like, shape = [n_samples, ], [n_samples, n_classes].
            Predicted multi-class targets.
        """
        check_is_fitted(self)

        n_samples = _num_samples(X)
        if self.label_binarizer_.y_type_ == "multiclass":
            maxima = np.empty(n_samples, dtype=float)
            maxima.fill(-np.inf)
            argmaxima = np.zeros(n_samples, dtype=int)
            for i, e in enumerate(self.estimators_):
                pred = _predict_binary(e, X)
                np.maximum(maxima, pred, out=maxima)
                argmaxima[maxima == pred] = i
            return self.classes_[argmaxima]
        else:
            if (hasattr(self.estimators_[0], "decision_function") and
                    is_classifier(self.estimators_[0])):
                thresh = 0
            else:
                thresh = .5
            indices = array.array('i')
            indptr = array.array('i', [0])
            for e in self.estimators_:
                indices.extend(np.where(_predict_binary(e, X) > thresh)[0])
                indptr.append(len(indices))
            data = np.ones(len(indices), dtype=int)
            indicator = sp.csc_matrix((data, indices, indptr),
                                      shape=(n_samples, len(self.estimators_)))
            return self.label_binarizer_.inverse_transform(indicator)

    @if_delegate_has_method(['_first_estimator', 'estimator'])
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
        T : (sparse) array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in the model,
            where classes are ordered as they are in `self.classes_`.
        """
        check_is_fitted(self)
        # Y[i, j] gives the probability that sample i has the label j.
        # In the multi-label case, these are not disjoint.
        Y = np.array([e.predict_proba(X)[:, 1] for e in self.estimators_]).T

        if len(self.estimators_) == 1:
            # Only one estimator, but we still want to return probabilities
            # for two classes.
            Y = np.concatenate(((1 - Y), Y), axis=1)

        if not self.multilabel_:
            # Then, probabilities should be normalized to 1.
            Y /= np.sum(Y, axis=1)[:, np.newaxis]
        return Y

    @if_delegate_has_method(['_first_estimator', 'estimator'])
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
        check_is_fitted(self)
        if len(self.estimators_) == 1:
            return self.estimators_[0].decision_function(X)
        return np.array([est.decision_function(X).ravel()
                         for est in self.estimators_]).T

    @property
    def multilabel_(self):
        """Whether this is a multilabel classifier"""
        return self.label_binarizer_.y_type_.startswith('multilabel')

    @property
    def n_classes_(self):
        return len(self.classes_)

    @property
    def coef_(self):
        check_is_fitted(self)
        if not hasattr(self.estimators_[0], "coef_"):
            raise AttributeError(
                "Base estimator doesn't have a coef_ attribute.")
        coefs = [e.coef_ for e in self.estimators_]
        if sp.issparse(coefs[0]):
            return sp.vstack(coefs)
        return np.vstack(coefs)

    @property
    def intercept_(self):
        check_is_fitted(self)
        if not hasattr(self.estimators_[0], "intercept_"):
            raise AttributeError(
                "Base estimator doesn't have an intercept_ attribute.")
        return np.array([e.intercept_.ravel() for e in self.estimators_])

    @property
    def _pairwise(self):
        """Indicate if wrapped estimator is using a precomputed Gram matrix"""
        return getattr(self.estimator, "_pairwise", False)

    @property
    def _first_estimator(self):
        return self.estimators_[0]


def _fit_ovo_binary(estimator, X, y, i, j):
    """Fit a single binary estimator (one-vs-one)."""
    cond = np.logical_or(y == i, y == j)
    y = y[cond]
    y_binary = np.empty(y.shape, np.int)
    y_binary[y == i] = 0
    y_binary[y == j] = 1
    indcond = np.arange(X.shape[0])[cond]
    return _fit_binary(estimator,
                       _safe_split(estimator, X, None, indices=indcond)[0],
                       y_binary, classes=[i, j]), indcond


def _partial_fit_ovo_binary(estimator, X, y, i, j):
    """Partially fit a single binary estimator(one-vs-one)."""

    cond = np.logical_or(y == i, y == j)
    y = y[cond]
    if len(y) != 0:
        y_binary = np.zeros_like(y)
        y_binary[y == j] = 1
        return _partial_fit_binary(estimator, X[cond], y_binary)
    return estimator


class OneVsOneClassifier(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
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

    Read more in the :ref:`User Guide <ovo_classification>`.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing :term:`fit` and one of
        :term:`decision_function` or :term:`predict_proba`.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    estimators_ : list of ``n_classes * (n_classes - 1) / 2`` estimators
        Estimators used for predictions.

    classes_ : numpy array of shape [n_classes]
        Array containing labels.

    n_classes_ : int
        Number of classes

    pairwise_indices_ : list, length = ``len(estimators_)``, or ``None``
        Indices of samples used when training the estimators.
        ``None`` when ``estimator`` does not have ``_pairwise`` attribute.
    """

    def __init__(self, estimator, n_jobs=None):
        self.estimator = estimator
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying estimators.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        y : array-like, shape = [n_samples]
            Multi-class targets.

        Returns
        -------
        self
        """
        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc'])
        check_classification_targets(y)

        self.classes_ = np.unique(y)
        if len(self.classes_) == 1:
            raise ValueError("OneVsOneClassifier can not be fit when only one"
                             " class is present.")
        n_classes = self.classes_.shape[0]
        estimators_indices = list(zip(*(Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_ovo_binary)
            (self.estimator, X, y, self.classes_[i], self.classes_[j])
            for i in range(n_classes) for j in range(i + 1, n_classes)))))

        self.estimators_ = estimators_indices[0]
        self.pairwise_indices_ = (
            estimators_indices[1] if self._pairwise else None)

        return self

    @if_delegate_has_method(delegate='estimator')
    def partial_fit(self, X, y, classes=None):
        """Partially fit underlying estimators

        Should be used when memory is inefficient to train all data. Chunks
        of data can be passed in several iteration, where the first call
        should have an array of all target variables.


        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        y : array-like, shape = [n_samples]
            Multi-class targets.

        classes : array, shape (n_classes, )
            Classes across all calls to partial_fit.
            Can be obtained via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is only required in the first call of partial_fit
            and can be omitted in the subsequent calls.

        Returns
        -------
        self
        """
        if _check_partial_fit_first_call(self, classes):
            self.estimators_ = [clone(self.estimator) for _ in
                                range(self.n_classes_ *
                                      (self.n_classes_ - 1) // 2)]

        if len(np.setdiff1d(y, self.classes_)):
            raise ValueError("Mini-batch contains {0} while it "
                             "must be subset of {1}".format(np.unique(y),
                                                            self.classes_))

        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc'])
        check_classification_targets(y)
        combinations = itertools.combinations(range(self.n_classes_), 2)
        self.estimators_ = Parallel(
            n_jobs=self.n_jobs)(
                delayed(_partial_fit_ovo_binary)(
                    estimator, X, y, self.classes_[i], self.classes_[j])
                for estimator, (i, j) in zip(self.estimators_,
                                             (combinations)))

        self.pairwise_indices_ = None

        return self

    def predict(self, X):
        """Estimate the best class label for each sample in X.

        This is implemented as ``argmax(decision_function(X), axis=1)`` which
        will return the label of the class with most votes by estimators
        predicting the outcome of a decision for each possible class pair.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : numpy array of shape [n_samples]
            Predicted multi-class targets.
        """
        Y = self.decision_function(X)
        if self.n_classes_ == 2:
            return self.classes_[(Y > 0).astype(np.int)]
        return self.classes_[Y.argmax(axis=1)]

    def decision_function(self, X):
        """Decision function for the OneVsOneClassifier.

        The decision values for the samples are computed by adding the
        normalized sum of pair-wise classification confidence levels to the
        votes in order to disambiguate between the decision values when the
        votes for all the classes are equal leading to a tie.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        Y : array-like, shape = [n_samples, n_classes]
        """
        check_is_fitted(self)

        indices = self.pairwise_indices_
        if indices is None:
            Xs = [X] * len(self.estimators_)
        else:
            Xs = [X[:, idx] for idx in indices]

        predictions = np.vstack([est.predict(Xi)
                                 for est, Xi in zip(self.estimators_, Xs)]).T
        confidences = np.vstack([_predict_binary(est, Xi)
                                 for est, Xi in zip(self.estimators_, Xs)]).T
        Y = _ovr_decision_function(predictions,
                                   confidences, len(self.classes_))
        if self.n_classes_ == 2:
            return Y[:, 1]
        return Y

    @property
    def n_classes_(self):
        return len(self.classes_)

    @property
    def _pairwise(self):
        """Indicate if wrapped estimator is using a precomputed Gram matrix"""
        return getattr(self.estimator, "_pairwise", False)


class OutputCodeClassifier(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
    """(Error-Correcting) Output-Code multiclass strategy

    Output-code based strategies consist in representing each class with a
    binary code (an array of 0s and 1s). At fitting time, one binary
    classifier per bit in the code book is fitted.  At prediction time, the
    classifiers are used to project new points in the class space and the class
    closest to the points is chosen. The main advantage of these strategies is
    that the number of classifiers used can be controlled by the user, either
    for compressing the model (0 < code_size < 1) or for making the model more
    robust to errors (code_size > 1). See the documentation for more details.

    Read more in the :ref:`User Guide <ecoc>`.

    Parameters
    ----------
    estimator : estimator object
        An estimator object implementing :term:`fit` and one of
        :term:`decision_function` or :term:`predict_proba`.

    code_size : float
        Percentage of the number of classes to be used to create the code book.
        A number between 0 and 1 will require fewer classifiers than
        one-vs-the-rest. A number greater than 1 will require more classifiers
        than one-vs-the-rest.

    random_state : int, RandomState instance or None, optional, default: None
        The generator used to initialize the codebook.  If int, random_state is
        the seed used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    estimators_ : list of `int(n_classes * code_size)` estimators
        Estimators used for predictions.

    classes_ : numpy array of shape [n_classes]
        Array containing labels.

    code_book_ : numpy array of shape [n_classes, code_size]
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

    def __init__(self, estimator, code_size=1.5, random_state=None,
                 n_jobs=None):
        self.estimator = estimator
        self.code_size = code_size
        self.random_state = random_state
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying estimators.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        y : numpy array of shape [n_samples]
            Multi-class targets.

        Returns
        -------
        self
        """
        X, y = check_X_y(X, y)
        if self.code_size <= 0:
            raise ValueError("code_size should be greater than 0, got {0}"
                             "".format(self.code_size))

        _check_estimator(self.estimator)
        random_state = check_random_state(self.random_state)
        check_classification_targets(y)

        self.classes_ = np.unique(y)
        n_classes = self.classes_.shape[0]
        code_size_ = int(n_classes * self.code_size)

        # FIXME: there are more elaborate methods than generating the codebook
        # randomly.
        self.code_book_ = random_state.random_sample((n_classes, code_size_))
        self.code_book_[self.code_book_ > 0.5] = 1

        if hasattr(self.estimator, "decision_function"):
            self.code_book_[self.code_book_ != 1] = -1
        else:
            self.code_book_[self.code_book_ != 1] = 0

        classes_index = {c: i for i, c in enumerate(self.classes_)}

        Y = np.array([self.code_book_[classes_index[y[i]]]
                      for i in range(X.shape[0])], dtype=np.int)

        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_binary)(self.estimator, X, Y[:, i])
            for i in range(Y.shape[1]))

        return self

    def predict(self, X):
        """Predict multi-class targets using underlying estimators.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : numpy array of shape [n_samples]
            Predicted multi-class targets.
        """
        check_is_fitted(self)
        X = check_array(X)
        Y = np.array([_predict_binary(e, X) for e in self.estimators_]).T
        pred = euclidean_distances(Y, self.code_book_).argmin(axis=1)
        return self.classes_[pred]
