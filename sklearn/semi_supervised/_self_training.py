import numpy as np

from ..base import MetaEstimatorMixin, clone, BaseEstimator
from ..utils.validation import check_X_y, check_array, check_is_fitted
from ..utils.metaestimators import if_delegate_has_method
from ..utils import safe_mask
import warnings

__all__ = ["SelfTrainingClassifier"]

# Authors: Oliver Rausch   <rauscho@ethz.ch>
#          Patrice Becker  <beckerp@ethz.ch>
# License: BSD 3 clause


def _validate_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""
    if not hasattr(estimator, "predict_proba"):
        msg = "base_classifier ({}) should implement predict_proba!"
        raise ValueError(msg.format(type(estimator).__name__))


class SelfTrainingClassifier(MetaEstimatorMixin, BaseEstimator):
    """Self-training classifier

    This class allows a given supervised classifier to function as a
    semi-supervised classifier, allowing it to learn from unlabeled data. It
    does this by iteratively predicting pseudo-labels for the unlabeled data
    and adding them to the training set.

    The classifier will continue iterating until either max_iter is reached, or
    no pseudo-labels were added to the training set in the previous iteration.

    Read more in the :ref:`User Guide <self_training>`.

    Parameters
    ----------
    base_classifier : estimator object
        An estimator object implementing ``fit`` and ``predict_proba``.
        Invoking the ``fit`` method will fit a clone of the passed estimator,
        which will be stored in the ``base_classifier_`` attribute.

    selection_criterion : {'threshold', 'n_best'}, optional
    (default='threshold')
        The selection criterion used to select which labels to add to the
        training set. If 'threshold', pseudo-labels with prediction
        probabilities above `threshold` are added to the dataset. If 'n_best',
        the `n_best` pseudo-labels with highest prediction probabilities are
        added to the dataset.

    threshold : float, optional (default=0.75)
        The decision threshold for use with `selection_criterion`='threshold'.
        Should be in [0, 1).

    n_best : int, optional (default=10)
        The amount of samples to add in each iteration. Only used when
        `selection_criterion`='n_best'.

    max_iter : int or ``None``, optional (default=10)
        Maximum number of iterations allowed. Should be greater than or equal
        to 0. If it is ``None``, the classifier will continue to predict labels
        until no new pseudo-labels are added, or all unlabeled samples have
        been labeled.

    verbose: bool, (default=False)
        Enable verbose output.

    Attributes
    ----------
    base_classifier_ : estimator object
        The fitted estimator.

    transduction_ : array, shape=(n_samples,)
        The labels used for the final fit of the classifier, including
        pseudo-labels added during fit.

    labeled_iter_ : array, shape=(n_samples,)
        The iteration in which each sample was labeled. When a sample has
        iteration 0, the sample was already labeled in the original dataset.
        When a sample has iteration -1, the sample was not labeled in any
        iteration.

    n_iter_ : int
        The number of rounds of self-training, that is the number of times the
        base estimator is fitted on relabeled variants of the training set.

    termination_condition_ : {'max_iter', 'no_change', 'all_labeled'}
        The reason that fitting was stopped.

        - 'max_iter': `n_iter_` reached `max_iter`.
        - 'no_change': no new labels were predicted.
        - 'all_labeled': all unlabeled samples were labeled before `max_iter`
          was reached.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import datasets
    >>> from sklearn.semi_supervised import SelfTrainingClassifier
    >>> from sklearn.svm import SVC
    >>> rng = np.random.RandomState(42)
    >>> iris = datasets.load_iris()
    >>> random_unlabeled_points = rng.rand(iris.target.shape[0]) < 0.3
    >>> iris.target[random_unlabeled_points] = -1
    >>> svc = SVC(probability=True, gamma="auto")
    >>> self_training_model = SelfTrainingClassifier(svc)
    >>> self_training_model.fit(iris.data, iris.target)
    ... # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    SelfTrainingClassifier(...)

    References
    ----------
    David Yarowsky. 1995. Unsupervised word sense disambiguation rivaling
    supervised methods. In Proceedings of the 33rd annual meeting on
    Association for Computational Linguistics (ACL '95). Association for
    Computational Linguistics, Stroudsburg, PA, USA, 189-196. DOI:
    https://doi.org/10.3115/981658.981684
    """
    _estimator_type = "classifier"

    def __init__(self,
                 base_classifier,
                 threshold=0.75,
                 selection_criterion='threshold',
                 n_best=10,
                 max_iter=10,
                 verbose=False):
        self.base_classifier = base_classifier
        self.threshold = threshold
        self.selection_criterion = selection_criterion
        self.n_best = n_best
        self.max_iter = max_iter
        self.verbose = verbose

    def fit(self, X, y):
        """
        Fits this ``SelfTrainingClassifier`` to a dataset.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        y : array-like, shape = (n_samples,)
            array representing the labels. Unlabeled samples should have the
            label -1.

        Returns
        -------
        self : object
            returns an instance of self.
        """
        X, y = check_X_y(X, y)

        if self.base_classifier is None:
            raise ValueError("base_classifier cannot be None!")

        self.base_classifier_ = clone(self.base_classifier)

        if self.max_iter is not None and self.max_iter < 0:
            msg = "max_iter must be >= 0 or None, got {}".format(self.max_iter)
            raise ValueError(msg)

        if not 0 <= self.threshold < 1:
            msg = "threshold must be in [0,1), got {}".format(self.threshold)
            raise ValueError(msg)

        if self.selection_criterion not in ['threshold', 'n_best']:
            raise ValueError("selection_criterion must be either 'threshold' "
                             "or 'n_best', got "
                             "{}".format(self.selection_criterion))

        if y.dtype.kind in ['U', 'S']:
            raise ValueError("y has dtype string. If you wish to predict on "
                             "string targets, use dtype " "object, and use -1"
                             " as the label for unlabeled samples.")

        has_label = y != -1

        if np.all(has_label):
            warnings.warn("y contains no unlabeled samples", UserWarning)

        if self.selection_criterion == 'n_best' and (self.n_best > X.shape[0] -
                                                     np.sum(has_label)):
            warnings.warn("n_best is larger than the amount of unlabeled "
                          "samples. All unlabeled samples will be labeled in "
                          "the first iteration", UserWarning)

        self.transduction_ = np.copy(y)
        self.labeled_iter_ = np.full_like(y, -1)
        self.labeled_iter_[has_label] = 0

        self.n_iter_ = 0

        while not np.all(has_label) and (self.max_iter is None or
                                         self.n_iter_ < self.max_iter):
            self.n_iter_ += 1
            self.base_classifier_.fit(
                X[safe_mask(X, has_label)],
                self.transduction_[safe_mask(self.transduction_, has_label)])

            if self.n_iter_ == 1:
                _validate_estimator(self.base_classifier)

            # Predict on the unlabeled samples
            prob = self.base_classifier_.predict_proba(
                X[safe_mask(X, ~has_label)])
            pred = self.base_classifier_.classes_[np.argmax(prob, axis=1)]
            max_proba = np.max(prob, axis=1)

            # Select samples
            if self.selection_criterion == 'threshold':
                new_labels_unlabeled = max_proba > self.threshold
            else:
                n_to_select = min(self.n_best, max_proba.shape[0])
                if n_to_select == max_proba.shape[0]:
                    new_labels_unlabeled = np.ones_like(max_proba, dtype=bool)
                else:
                    # NB these are indicies, not a mask
                    new_labels_unlabeled = \
                        np.argpartition(max_proba, n_to_select)[:n_to_select]

            # new_labels_unlabeled indexes into only the unlabeled samples
            # new_labels_full indexes into the full X
            new_labels_full = np.nonzero(~has_label)[0][new_labels_unlabeled]

            # Add newly labeled confident predictions to the dataset
            self.transduction_[new_labels_full] = pred[new_labels_unlabeled]
            has_label[new_labels_full] = True
            self.labeled_iter_[new_labels_full] = self.n_iter_

            if new_labels_full.shape[0] == 0:
                # no changed labels
                self.termination_condition_ = "no_change"
                break

            if self.verbose:
                msg = "End of iteration {}, added {} new labels."
                print(msg.format(self.n_iter_, new_labels_full.shape[0]))

        if self.n_iter_ == self.max_iter:
            self.termination_condition_ = "max_iter"
        if np.all(has_label):
            self.termination_condition_ = "all_labeled"

        self.base_classifier_.fit(
            X[safe_mask(X, has_label)],
            self.transduction_[safe_mask(self.transduction_, has_label)])
        self.classes_ = self.base_classifier_.classes_
        return self

    @if_delegate_has_method(delegate='base_classifier')
    def predict(self, X):
        """Predict the classes of X.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples,)
            array with predicted labels
        """
        check_is_fitted(self, 'transduction_')
        X = check_array(X)
        return self.base_classifier_.predict(X)

    def predict_proba(self, X):
        """Predict probability for each possible outcome.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples, n_features)
            array with prediction probabilities
        """
        check_is_fitted(self, 'transduction_')
        return self.base_classifier_.predict_proba(X)

    @if_delegate_has_method(delegate='base_classifier')
    def decision_function(self, X):
        """Calls decision function of the base_classifier

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples, n_features)
            result of the decision function of the base_classifier
        """
        check_is_fitted(self, 'transduction_')
        return self.base_classifier_.decision_function(X)

    @if_delegate_has_method(delegate='base_classifier')
    def predict_log_proba(self, X):
        """Predict log probability for each possible outcome.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape=(n_samples, n_features)
            array with log prediction probabilities
        """
        check_is_fitted(self, 'transduction_')
        return self.base_classifier_.predict_log_proba(X)

    @if_delegate_has_method(delegate='base_classifier')
    def score(self, X, y):
        """Calls score on the base_classifier

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            array representing the data

        y : array-like, shape=(n_samples,)
            array representing the labels

        Returns
        -------
        score : float
            result of calling score on the base_classifier
        """
        check_is_fitted(self, 'transduction_')
        return self.base_classifier_.score(X, y)
