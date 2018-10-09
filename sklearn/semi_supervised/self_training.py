import numpy as np

from ..base import BaseEstimator, clone
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
    if estimator is None:
        raise ValueError("base_classifier cannot be None!")

    if not hasattr(estimator, "predict_proba"):
        name = type(estimator).__name__
        msg = "base_classifier ({}) should implement predict_proba!"
        msg = msg.format(name)
        raise ValueError(msg)


class SelfTrainingClassifier(BaseEstimator):
    """Self-training classifier

    Parameters
    ----------
    base_classifier : estimator object
        An estimator object implementing ``fit`` and ``predict_proba``.
        Invoking the ``fit`` method will fit a clone of the passed estimator,
        which will be stored in the ``self.base_classifier_`` attribute.

    threshold : float, optional (default=0.75)
        Threshold above which predictions are added to the labeled dataset.
        Should be in [0, 1).

    max_iter : integer, optional (default=20)
        Maximum number of iterations allowed. Should be greater than or equal
        to 0.

    Attributes
    ----------
    base_classifier_ : estimator object
        The fitted estimator.

    y_labeled_ : array, shape = (n_samples,)
        The labels used for the final fit of the classifier.

    y_labeled_iter_ : array, shape = (n_samples,)
        The iteration in which each sample was labeled. When a sample has
        iteration 0, the sample was already labeled in the given dataset. When
        a sample has iteration -1, the sample was not labeled in any iteration.

    n_iter_ : int
        The amount of iterations the classifier used during fitting.


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
                 max_iter=10):
        self.base_classifier = base_classifier
        self.threshold = threshold
        self.max_iter = max_iter

    def fit(self, X, y):
        """
        Fits this ``SelfTrainingClassifier`` to dataset, using the
        base_classifier passed. The fitted base_classifier is stored as
        ``self.base_classifier_``.

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
        _validate_estimator(self.base_classifier)
        self.base_classifier_ = clone(self.base_classifier)

        if self.max_iter < 0:
            msg = "max_iter must be >= 0, got {}".format(self.max_iter)
            raise ValueError(msg)

        if not 0 <= self.threshold < 1:
            msg = "threshold must be in [0,1), got {}".format(self.threshold)
            raise ValueError(msg)

        if y.dtype.kind == 'U' or y.dtype.kind == 'S':
            has_label = y != '-1'
        else:
            has_label = y != -1

        if np.all(has_label):
            warnings.warn(RuntimeWarning("y contains no unlabeled samples"))

        self.y_labels_ = np.copy(y)
        self.y_labeled_iter_ = np.full_like(y, -1)
        self.y_labeled_iter_[has_label] = 0

        self.n_iter_ = 0
        while not np.all(has_label) and self.n_iter_ < self.max_iter:
            self.n_iter_ += 1
            self.base_classifier_.fit(
                X[safe_mask(X, has_label)],
                self.y_labels_[safe_mask(self.y_labels_, has_label)])

            # Predict on the unlabeled samples
            prob = self.base_classifier_.predict_proba(
                X[safe_mask(X, ~has_label)])
            pred = self.base_classifier_.classes_[np.argmax(prob, axis=1)]
            max_proba = np.max(prob, axis=1)

            # Select samples where confidence is above the threshold
            confident_labels = max_proba > self.threshold

            new_labels_idx = np.flatnonzero(~has_label)[confident_labels]

            # Add newly labeled confident predictions to the dataset
            self.y_labels_[new_labels_idx] = pred[confident_labels]
            has_label[new_labels_idx] = True
            self.y_labeled_iter_[new_labels_idx] = self.n_iter_

        self.base_classifier_.fit(
            X[safe_mask(X, has_label)],
            self.y_labels_[safe_mask(self.y_labels_, has_label)])
        self.classes_ = self.base_classifier_.classes_
        return self

    @if_delegate_has_method(delegate='base_classifier')
    def predict(self, X):
        """Predict the classes of X.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples,)
            array with predicted labels
        """
        check_is_fitted(self, 'y_labeled_iter_')
        X = check_array(X)
        return self.base_classifier_.predict(X)

    def predict_proba(self, X):
        """Predict probability for each possible outcome.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples, n_features)
            array with prediction probabilities
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.predict_proba(X)

    @if_delegate_has_method(delegate='base_classifier')
    def decision_function(self, X):
        """Calls decision function of the base_classifier

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples, n_features)
            result of the decision function of the base_classifier
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.decision_function(X)

    @if_delegate_has_method(delegate='base_classifier')
    def predict_log_proba(self, X):
        """Predict log probability for each possible outcome.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples, n_features)
            array with log prediction probabilities
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.predict_log_proba(X)

    @if_delegate_has_method(delegate='base_classifier')
    def score(self, X, y):
        """Calls score on the base_classifier

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        y : array-like, shape = (n_samples, 1)
            array representing the labels

        Returns
        -------
        score : float
            result of calling score on the base_classifier
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_classifier_.score(X, y)
