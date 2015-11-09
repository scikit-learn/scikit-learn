"""
Multiclass and multilabel classification strategies
===================================================

This module implements multiclass learning algorithms:
    - one-vs-the-rest / one-vs-all
    - one-vs-one
    - error correcting output codes
    - rakel (multilabel)

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
# Author: Alain Pena <alain.pena.r@gmail.com>
#
# License: BSD 3 clause

import array
import numpy as np
import warnings
import scipy.sparse as sp
import copy as cp

from itertools import izip

from .base import BaseEstimator, ClassifierMixin, clone, is_classifier
from .base import MetaEstimatorMixin, is_regressor
from .preprocessing import LabelBinarizer
from .metrics.pairwise import euclidean_distances
from .utils import check_random_state, check_X_y
from .utils.validation import _num_samples
from .utils.validation import check_consistent_length
from .utils.validation import check_is_fitted
from .externals.joblib import Parallel
from .externals.joblib import delayed
from .utilities import get_random_numbers, cycle_permutations
from .basic_checks import check_is_in_range, check_is_binary_01
from .powerset import Powerset

__all__ = [
    "OneVsRestClassifier",
    "OneVsOneClassifier",
    "OutputCodeClassifier",
    "RakelClassifier",
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
        check_is_fitted(self, 'y_')

        return np.repeat(self.y_, X.shape[0])

    def decision_function(self, X):
        check_is_fitted(self, 'y_')

        return np.repeat(self.y_, X.shape[0])

    def predict_proba(self, X):
        check_is_fitted(self, 'y_')

        return np.repeat([np.hstack([1 - self.y_, self.y_])],
                         X.shape[0], axis=0)


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
    is used to predict multiple labels for instance, by fitting on a 2-d matrix
    in which cell [i, j] is 1 if sample i has label j and 0 otherwise.

    In the multilabel learning literature, OvR is also known as the binary
    relevance method.

    Read more in the :ref:`User Guide <ovr_classification>`.

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

    def __init__(self, estimator, n_jobs=1):
        self.estimator = estimator
        self.n_jobs = n_jobs

    def fit(self, X, y):
        """Fit underlying estimators.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        y : (sparse) array-like, shape = [n_samples] or [n_samples, n_classes]
            Multi-class targets. An indicator matrix turns on multilabel
            classification.

        Returns
        -------
        self
        """
        # A sparse LabelBinarizer, with sparse_output=True, has been shown to
        # outpreform or match a dense label binarizer in all cases and has also
        # resulted in less or equal memory consumption in the fit_ovr function
        # overall.
        self.label_binarizer_ = LabelBinarizer(sparse_output=True)
        Y = self.label_binarizer_.fit_transform(y)
        Y = Y.tocsc()
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

    def predict(self, X):
        """Predict multi-class targets using underlying estimators.

        Parameters
        ----------
        X : (sparse) array-like, shape = [n_samples, n_features]
            Data.

        Returns
        -------
        y : (sparse) array-like, shape = [n_samples] or [n_samples, n_classes].
            Predicted multi-class targets.
        """
        check_is_fitted(self, 'estimators_')
        if (hasattr(self.estimators_[0], "decision_function") and
                is_classifier(self.estimators_[0])):
            thresh = 0
        else:
            thresh = .5

        n_samples = _num_samples(X)
        if self.label_binarizer_.y_type_ == "multiclass":
            maxima = np.empty(n_samples, dtype=float)
            maxima.fill(-np.inf)
            argmaxima = np.zeros(n_samples, dtype=int)
            for i, e in enumerate(self.estimators_):
                pred = _predict_binary(e, X)
                np.maximum(maxima, pred, out=maxima)
                argmaxima[maxima == pred] = i
            return self.label_binarizer_.classes_[np.array(argmaxima.T)]
        else:
            indices = array.array('i')
            indptr = array.array('i', [0])
            for e in self.estimators_:
                indices.extend(np.where(_predict_binary(e, X) > thresh)[0])
                indptr.append(len(indices))
            data = np.ones(len(indices), dtype=int)
            indicator = sp.csc_matrix((data, indices, indptr),
                                      shape=(n_samples, len(self.estimators_)))
            return self.label_binarizer_.inverse_transform(indicator)

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
        check_is_fitted(self, 'estimators_')
        # Y[i,j] gives the probability that sample i has the label j.
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
        check_is_fitted(self, 'estimators_')
        if not hasattr(self.estimators_[0], "decision_function"):
            raise AttributeError(
                "Base estimator doesn't have a decision_function attribute.")
        return np.array([est.decision_function(X).ravel()
                         for est in self.estimators_]).T

    @property
    def multilabel_(self):
        """Whether this is a multilabel classifier"""
        return self.label_binarizer_.y_type_.startswith('multilabel')

    @property
    def classes_(self):
        return self.label_binarizer_.classes_

    @property
    def coef_(self):
        check_is_fitted(self, 'estimators_')
        if not hasattr(self.estimators_[0], "coef_"):
            raise AttributeError(
                "Base estimator doesn't have a coef_ attribute.")
        coefs = [e.coef_ for e in self.estimators_]
        if sp.issparse(coefs[0]):
            return sp.vstack(coefs)
        return np.vstack(coefs)

    @property
    def intercept_(self):
        check_is_fitted(self, 'estimators_')
        if not hasattr(self.estimators_[0], "intercept_"):
            raise AttributeError(
                "Base estimator doesn't have an intercept_ attribute.")
        return np.array([e.intercept_.ravel() for e in self.estimators_])


def _fit_ovo_binary(estimator, X, y, i, j):
    """Fit a single binary estimator (one-vs-one)."""
    cond = np.logical_or(y == i, y == j)
    y = y[cond]
    y_binary = np.empty(y.shape, np.int)
    y_binary[y == i] = 0
    y_binary[y == j] = 1
    ind = np.arange(X.shape[0])
    return _fit_binary(estimator, X[ind[cond]], y_binary, classes=[i, j])


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

    Read more in the :ref:`User Guide <ovo_classification>`.

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
    estimators_ : list of `n_classes * (n_classes - 1) / 2` estimators
        Estimators used for predictions.

    classes_ : numpy array of shape [n_classes]
        Array containing labels.
    """

    def __init__(self, estimator, n_jobs=1):
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
        y = np.asarray(y)
        check_consistent_length(X, y)

        self.classes_ = np.unique(y)
        n_classes = self.classes_.shape[0]
        self.estimators_ = Parallel(n_jobs=self.n_jobs)(
            delayed(_fit_ovo_binary)(
                self.estimator, X, y, self.classes_[i], self.classes_[j])
            for i in range(n_classes) for j in range(i + 1, n_classes))

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
        check_is_fitted(self, 'estimators_')

        predictions = np.vstack([est.predict(X) for est in self.estimators_]).T
        confidences = np.vstack([_predict_binary(est, X) for est in self.estimators_]).T
        return _ovr_decision_function(predictions, confidences,
                                      len(self.classes_))


def _ovr_decision_function(predictions, confidences, n_classes):
    """Compute a continuous, tie-breaking ovr decision function.

    It is important to include a continuous value, not only votes,
    to make computing AUC or calibration meaningful.

    Parameters
    ----------
    predictions : array-like, shape (n_samples, n_classifiers)
        Predicted classes for each binary classifier.

    confidences : array-like, shape (n_samples, n_classifiers)
        Decision functions or predicted probabilities for positive class
        for each binary classifier.

    n_classes : int
        Number of classes. n_classifiers must be
        ``n_classes * (n_classes - 1 ) / 2``
    """
    n_samples = predictions.shape[0]
    votes = np.zeros((n_samples, n_classes))
    sum_of_confidences = np.zeros((n_samples, n_classes))

    k = 0
    for i in range(n_classes):
        for j in range(i + 1, n_classes):
            sum_of_confidences[:, i] -= confidences[:, k]
            sum_of_confidences[:, j] += confidences[:, k]
            votes[predictions[:, k] == 0, i] += 1
            votes[predictions[:, k] == 1, j] += 1
            k += 1

    max_confidences = sum_of_confidences.max()
    min_confidences = sum_of_confidences.min()

    if max_confidences == min_confidences:
        return votes

    # Scale the sum_of_confidences to (-0.5, 0.5) and add it with votes.
    # The motivation is to use confidence levels as a way to break ties in
    # the votes without switching any decision made based on a difference
    # of 1 vote.
    eps = np.finfo(sum_of_confidences.dtype).eps
    max_abs_confidence = max(abs(max_confidences), abs(min_confidences))
    scale = (0.5 - eps) / max_abs_confidence
    return votes + sum_of_confidences * scale


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

    Read more in the :ref:`User Guide <ecoc>`.

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

    def __init__(self, estimator, code_size=1.5, random_state=None, n_jobs=1):
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
        if self.code_size <= 0:
            raise ValueError("code_size should be greater than 0, got {1}"
                             "".format(self.code_size))

        _check_estimator(self.estimator)
        random_state = check_random_state(self.random_state)

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

        classes_index = dict((c, i) for i, c in enumerate(self.classes_))

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
        check_is_fitted(self, 'estimators_')
        Y = np.array([_predict_binary(e, X) for e in self.estimators_]).T
        pred = euclidean_distances(Y, self.code_book_).argmin(axis=1)
        return self.classes_[pred]


def _binomialCoeff(n, k):
    # see http://rosettacode.org/wiki/Evaluate_binomial_coefficients
    result = 1
    if (n <= k):
        return result
    for i in range(1, k+1):
        result = result * (n-i+1) / i
    return result


def _cmp_fct(val1, val2):
    """
    Comparing function between 2 arrays of same size.
    """
    retv_ = 0
    for i1, i2 in izip(val1, val2):
        retv_ = cmp(i1, i2)
        if (retv_ != 0):
            return retv_
    return 0


def _valid_possibility(possibility, sorted_array):
    """
    Tests if the possibility is absent from the sorted array.
    """
    cmped = 0
    for elem in sorted_array:
        cmped = _cmp_fct(possibility, elem)
        if (cmped == 0):
            return False
        elif (cmped < 0):
            return True
    return True


def _make_possibility(poss, n_labels):
    """
    Converts an array of indexes into a "boolean" array.
    """
    possibility = np.zeros((n_labels,), dtype=np.int8)
    possibility[poss] = 1
    return possibility


def _get_possibility(start, n_max, sorted_array):
    """
    Returns a possible permutation from a start which does not appear
    in a sorted array or raises a `ValueError` if cannot be found.
    """
    for p in cycle_permutations(start, n_iter=n_max):
        if _valid_possibility(p, sorted_array):
            return p
    raise ValueError("Error: no possibility found.")


class RakelClassifier(BaseEstimator, ClassifierMixin, MetaEstimatorMixin):

    """
    Implements a rakel meta-classifier.

    Parameters
    ----------
    estimator : sklearn.base.BaseEstimator
        The underlying estimator.

    k : Number
        Size of the labelset to generate.
        Integer value means its absolute size,
        while float value in range [0, 1[ means the relative size.

    n_estimators : Number
        Number of labelsets to generate. It is the "m" parameter.

    threshold : float
        Threshold for the meta-estimator to consider as a true or false value.

    no_hole : bool
        Whether or not the union of the labelsets must totally cover the order.

    random_state : int or RandomState
        RNG seed.

    fill_method : float or string
        The way the classifier fills predictions if a label is not
        in any labelset.
            "most_frequent" for using the most frequent value for
            this label (using percentage)

            float in [0, 1] for a constant value

    powerset : str
        The method to decompress the predictions.

        `majority` for taking the most dominant class, then split it back
        into multilabel

        `probabilities` for decompressing each class following its estimated
        probabilities.

    verbose : int
        Verbosity level.

    n_jobs : int
        Number of jobs to run concurrently, currently not used.

    Reference
    ---------
    Tsoumakas, G., Vlahavas, I.:

        Random k-labelsets: An ensemble method for multilabel
        classification.

        In: Proceedings of the 18th European Conference
        on Machine Learning (ECML 2007),

        Warsaw, Poland (2007) 406-417
    """

    def __init__(self, estimator, k=1, n_estimators=1, threshold=0.5,
                 no_hole=True, random_state=None,
                 fill_method="most_frequent",
                 uniqueness=True,
                 powerset="majority",
                 verbose=0, n_jobs=1):
        self.estimator = estimator
        self.k = k
        self.n_estimators = n_estimators
        self.threshold = threshold
        self.no_hole = no_hole
        self.random_state = random_state
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.n_labels_ = 0
        self.fill_method = fill_method
        self.uniqueness = uniqueness
        self.powerset = powerset
        self.multi_label_ = True
        self.X_ = None
        self.y_ = None
        self.view_ = None
        self.labelsets_ = None
        self._ACCEPTABLE_COLLISION = 0.5

    def _get_labelsets(self, n_labelsets, n_labels, size_labelsets,
                       random_state=None, uniqueness=False):
        # n_max =                   factorial(n_labels)
        #    -----------------------------------------------------------------
        #    (factorial(size_labelsets) * factorial(n_labels - size_labelsets)
        if ((n_labels < size_labelsets) or
                (size_labelsets < 0) or
                (n_labels < 1)):
            return []
        if (uniqueness):
            n_max = _binomialCoeff(n_labels, size_labelsets)
            if (n_max < n_labelsets):
                raise ValueError("The combination of the arguments will "
                                 "force the generation of duplicates.")

        retv = [np.zeros((n_labels,), dtype=np.int8)
                for _ in xrange(n_labelsets)]
        if (size_labelsets == 0):
            return retv

        rnd = check_random_state(seed=random_state)

        r = [rnd.choice(n_labels, size=size_labelsets, replace=False)
             for _ in xrange(n_labelsets)]
        for idx, arr in izip(r, retv):
            arr[idx] = 1

        retv.sort(cmp=_cmp_fct)
        if not uniqueness:
            return retv

        possibility = None
        idx = 0
        collide_chance = n_max * self._ACCEPTABLE_COLLISION > n_labelsets
        while idx < (len(retv) - 1):
            # Iterates through the return value to be sure that each element
            # is unique.
            if (_cmp_fct(retv[idx], retv[idx+1]) == 0):
                if (collide_chance):
                    # If enough possibilities, let randomness decide
                    possibility = rnd.choice(n_labels, size=size_labelsets,
                                             replace=False)
                    possibility = _make_possibility(possibility, n_labels)
                    while not _valid_possibility(possibility, retv):
                        possibility = rnd.choice(n_labels, size=size_labelsets,
                                                 replace=False)
                        possibility = _make_possibility(possibility, n_labels)
                else:
                    # Good chances of collisions, so use deterministic way
                    possibility = _get_possibility(retv[idx+1], n_max, retv)
                retv[idx+1] = np.asarray(possibility, dtype=np.int8)
                retv.sort(cmp=_cmp_fct)
            else:
                idx = idx + 1
        return retv

    def _force_no_hole(self, labelsets):
        labs = np.array([sum(i) for i in izip(*labelsets)])
        if (np.sum(labs) < labs.size):
            raise ValueError("Impossible to get no hole.")

        indexes = np.count_nonzero(labs)
        if indexes == labs.size:
            return labelsets
        amin = np.argmin(labs)
        while not labs[amin]:
            # 1- select the most appearing element
            amax = np.argmax(labs)
            # 2- remove one and resort the array of labs
            labs[amax] = labs[amax] - 1
            labs[amin] = labs[amin] + 1
            idx_ = labs.size - 1
            while ((idx_ > 0) and (labs[idx_] < labs[idx_ - 1])):
                labs[idx_], labs[idx_ - 1] = labs[idx_ - 1], labs[idx_]
                idx_ = idx_ - 1
            # 3- find a labelset containing `most`th label and swap them
            idx_ = 0
            while not labelsets[idx_][amax]:
                idx_ = idx_ + 1
            labelsets[idx_][amax], labelsets[idx_][amin] = \
                labelsets[idx_][amin], labelsets[idx_][amax]
            amin = np.argmin(labs)
        return labelsets

    def _create_labelsets(self):
        safe_k = int(self.k)
        if check_is_in_range(self.k, lower_bound=0,
                             low_exclusive=True,
                             higher_bound=1,
                             high_exclusive=True,
                             launch_exception=False):
            safe_k = int(np.ceil(self.k * self.n_labels_))
        else:
            check_is_in_range(safe_k, lower_bound=1,
                              low_exclusive=False,
                              higher_bound=self.n_labels_,
                              high_exclusive=False,
                              launch_exception=True,
                              msg=("Invalid size of labelset ({0})."
                                   ).format(str(safe_k)),
                              exception_type=ValueError)

        def get_m():
            safe_m_ = int(self.n_estimators)
            if not self.uniqueness:
                return safe_m_
            return min(safe_m_, _binomialCoeff(self.n_labels_, safe_k))
        safe_m = get_m()

        check_is_in_range(safe_m, lower_bound=1,
                          low_exclusive=False,
                          launch_exception=True,
                          msg=("Invalid number of labelsets ({0})."
                               ).format(str(safe_m)),
                          exception_type=ValueError)
        if (self.no_hole and ((safe_k * safe_m) < self.n_labels_)):
            raise ValueError(("Impossible to get no hole for {0}"
                              " labels and k = {1}, m = {2}.").format(
                              str(self.n_labels_), str(safe_k), str(safe_m)))
        retv = self._get_labelsets(safe_m, self.n_labels_, safe_k,
                                   random_state=cp.copy(self.random_state),
                                   uniqueness=self.uniqueness)
        if self.no_hole:
            retv = self._force_no_hole(retv)
        retv.sort(cmp=_cmp_fct)
        return retv

    def fit(self, X, y, copy=False):
        """
        Fits the underlying estimator on the data.

        Parameters
        ----------
        X : array_like of shape [n_samples, n_features]
            The training input samples.

        y : array_like of shape [n_samples] or [n_samples, n_outputs]
            The target values.

        Returns
        -------
        RakelClassifier
            self.

        WARNING: IF COPY = FALSE:
        This estimator will use the X and y arrays for the
        prediction part, but those will not be copied. Thus, a
        modification in one of those will have as effect to fit and
        predict on the modified data.
        """
        self.X_, self.y_ = check_X_y(X, y, multi_output=True,
                                     accept_sparse=None, copy=copy,
                                     force_all_finite=True,
                                     allow_nd=False)

        if copy:
            self.y_ = cp.deepcopy(self.y_)
        try:
            n_labels = self.y_.shape[1]
        except IndexError:
            n_labels = 1
            self.view_ = np.reshape(self.y_, (len(self.y_), 1))
            self.multi_label_ = False
        else:
            self.view_ = self.y_
            self.multi_label_ = True

        self.n_labels_ = n_labels
        try:
            self.labelsets_ = self._create_labelsets()
        except ValueError:
            self.view_ = None
            self.n_labels_ = 0
            self.labelsets_ = None
            self.X_ = None
            self.y_ = None
            raise
        return self

    def __predict(self, X, predict_only=True):
        """
        Prediction helper.
        """
        if (self.labelsets_ is None):
            raise Exception(("Rakel not initialized. "
                             "Perform a fit first."))
        check_is_binary_01(
            self.view_, launch_exception=True,
            exception_type=ValueError)
        if self.fill_method == "most_frequent":
            def get_most_freqs():
                frequencies1 = np.array(
                    [float((self.view_[:, i] == 1).sum())
                        for i in xrange(self.n_labels_)])
                frequencies1 /= float(self.view_.shape[0])
                return frequencies1
            fill_vals = get_most_freqs()
        else:
            check_is_in_range(elem=self.fill_method,
                              lower_bound=0.0,
                              higher_bound=1.0, low_exclusive=False,
                              high_exclusive=False,
                              launch_exception=True,
                              msg="Error: invalid fill_method.",
                              exception_type=ValueError)
            fill_vals = [self.fill_method for _ in xrange(self.n_labels_)]

        if (self.verbose > 1):
            print("Number of estimators: {0}.".format(str(self.n_estimators)))
        votes = np.zeros((X.shape[0], self.n_labels_))
        if predict_only:
            avg = np.zeros((X.shape[0], self.n_labels_))
            retv = np.zeros((X.shape[0], self.n_labels_))
        else:
            retv = [np.zeros((X.shape[0], 2)) for _ in xrange(self.n_labels_)]

        rnds = get_random_numbers(shape=(self.n_estimators, ),
                                  random_state=cp.copy(self.random_state))
        for bin_labelset, rnd in izip(self.labelsets_, rnds):
            labelset = bin_labelset.nonzero()[0]
            if (self.verbose > 2):
                print("labelset: {0}".format(labelset))
            power = Powerset()
            y_train = power.compressed_powerize(self.view_[:, labelset])
            estim = clone(self.estimator)
            try:
                estim.set_params(random_state=rnd)
            except ValueError:
                pass
            estim.fit(self.X_, y_train)

            if predict_only:
                result = estim.predict(X)
                uncompressed = power.unpowerize(result)
                for sample, uncom in enumerate(uncompressed):
                    for label, unc in izip(labelset, uncom):
                        avg[sample, label] += unc
                        votes[sample, label] += 1
            else:
                r = estim.predict_proba(X)
                if self.powerset == "probabilities":
                    odd = power.probas_unpowerize(r)
                    for idx_label, label in enumerate(odd):
                        for i in xrange(label.shape[0]):
                            retv[labelset[idx_label]][i, :] += label[i, :]
                            votes[i, labelset[idx_label]] += 1
                else:
                    # odd = [[0 for _ in xrange(n_samples)]
                    #        for i in xrange(n_labels)]
                    odd = power.majority_unpowerize(r)
                    for idx_label, label in enumerate(odd):
                        retv[labelset[idx_label]][:, 1] += label
                        votes[:, labelset[idx_label]] += 1
                    for idx_label in xrange(self.n_labels_):
                        retv[idx_label][:, 0] = (votes[:, idx_label] -
                                                 retv[idx_label][:, 1])

            del estim

        if predict_only:
            for i in xrange(votes.shape[0]):
                for j in xrange(votes.shape[1]):
                    if votes[i, j] == 0:
                        if (self.verbose > 0):
                            print(("The label {0} of the sample {1} "
                                   "cannot be estimated. Set to fill "
                                   "value.").format(str(j), str(i)))
                        votes[i, j] = 1.0
                        avg[i, j] = fill_vals[j]
                    avg[i, j] /= votes[i, j]
                    if (avg[i, j] > self.threshold):
                        retv[i, j] = 1
        else:
            for i in xrange(votes.shape[0]):
                for j in xrange(votes.shape[1]):
                    if votes[i, j] == 0:
                        if (self.verbose > 0):
                            print(("The label {0} of the sample {1} "
                                   "cannot be estimated. Set to fill "
                                   "value.").format(str(j), str(i)))
                        votes[i, j] = 1.0
                        retv[j][i, 1] = fill_vals[j]
                    else:
                        retv[j][i, :] /= votes[i, j]
                    retv[j][i, 0] = 1.0 - retv[j][i, 1]

        if not self.multi_label_:
            if predict_only:
                return retv.ravel()
            else:
                return retv[0]

        return retv

    def predict(self, X):
        """
        Predict class or regression value for X using the underlying
            estimator. For a classification model, the predicted class for
        each sample in X is returned. For a regression model, the predicted
        value based on X is returned.

        Parameters
        ----------
        X : array_like of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted classes, or the predict values.
        """
        return self.__predict(X, predict_only=True)

    def predict_proba(self, X):
        """
        Predict class probabilities of the input samples X using
        the underlying estimator.

        Parameters
        ----------
        X : array_like of shape [n_samples, n_features]
            The input samples.

        Returns
        -------
        array of shape = [n_samples, n_classes], or a list of n_outputs
            such arrays if n_outputs > 1. The class probabilities of
            the input samples. The order of the classes corresponds to
            the one from the y used for fitting.
        """
        return self.__predict(X, predict_only=False)
