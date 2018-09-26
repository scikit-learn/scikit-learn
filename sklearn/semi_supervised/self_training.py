import numpy as np

from ..base import BaseEstimator, clone
from ..utils.validation import check_X_y, check_array, check_is_fitted
from ..utils.metaestimators import if_delegate_has_method
from ..utils import safe_mask
from ..base import is_classifier

__all__ = ["SelfTrainingClassifier"]

# Authors: Oliver Rausch   <rauscho@ethz.ch>
#          Patrice Becker  <beckerp@ethz.ch>
# License: BSD 3 clause


def _check_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""
    if not is_classifier(estimator):
        raise ValueError("The base_estimator should be a classifier!")

    if not hasattr(estimator, "predict_proba"):
        raise ValueError("The base_estimator should implement predict_proba!")


class SelfTrainingClassifier(BaseEstimator):
    """Self-training classifier

    Parameters
    ----------
    base_estimator : estimator object
        An estimator object implementing ``fit`` and ``predict_proba``.
        Invoking the ``fit`` method will fit a clone of the passed estimator,
        which will be stored in the ``self.base_estimator_`` attribute.

    threshold : float
        Threshold above which predictions are added to the labeled dataset.
        Should be in [0, 1).

    max_iter : integer
        Maximum number of iterations allowed. Should be greater than or equal
        to 0.

    Attributes
    ----------
    base_estimator_ : estimator object
        The fitted estimator.

    y_labeled_ : array, shape = (n_samples,)
        The labels used for the final fit of the classifier.

    y_labeled_iter_ : array, shape = (n_samples,)
        The iteration in which each sample was labeled. When a sample has
        iteration 0, the sample was already labeled in the given dataset. When
        a sample has iteration -1, the sample was not labeled in any iteration.

    n_iter_ : int
        The amount of iterations the classifier during fitting.


    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import datasets
    >>> from sklearn.semi_supervised import SelfTrainingClassifier
    >>> from sklearn.svm import SVC
    >>> svc = SVC(probability=True, gamma="auto")
    >>> self_training_model = SelfTrainingClassifier(svc)
    >>> iris = datasets.load_iris()
    >>> rng = np.random.RandomState(42)
    >>> random_unlabeled_points = rng.rand(len(iris.target)) < 0.3
    >>> labels = np.copy(iris.target)
    >>> labels[random_unlabeled_points] = -1
    >>> self_training_model.fit(iris.data, labels)
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
                 base_estimator,
                 threshold=0.75,
                 max_iter=100):
        self.base_estimator = base_estimator
        self.threshold = threshold
        self.max_iter = max_iter

    def fit(self, X, y):
        """
        Fits this ``SelfTrainingClassifier`` to dataset, using the
        base_estimator passed. The fitted base_estimator is stored as
        ``self.base_estimator_``.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        y : array-like, shape = (n_samples,)
            array representing the labels. Unlabeled samples should have the
            label -1.

        Returns
        -------
        self: returns an instance of self.
        """
        X, y = check_X_y(X, y)
        _check_estimator(self.base_estimator)
        self.base_estimator_ = clone(self.base_estimator)

        if not 0 <= self.max_iter:
            raise ValueError("max_iter must be >= 0")

        if not 0 <= self.threshold < 1:
            raise ValueError("threshold must be in [0,1)")

        if y.dtype.kind == 'U' or y.dtype.kind == 'S':
            has_label = y != '-1'
        else:
            has_label = y != -1

        self.y_labels_ = np.copy(y)
        self.y_labeled_iter_ = np.full_like(y, -1)
        self.y_labeled_iter_[has_label] = 0

        self.n_iter_ = 0
        while np.sum(has_label) < len(X) and self.n_iter_ < self.max_iter:
            self.n_iter_ += 1
            self.base_estimator_.fit(X[safe_mask(X, has_label)],
                                     self.y_labels_[safe_mask(self.y_labels_,
                                                              has_label)])

            # Predict on the unlabeled samples
            prob = self.base_estimator_.predict_proba(X[safe_mask(X,
                                                                  ~has_label)])
            pred = self.base_estimator_.classes_[np.argmax(prob, axis=1)]
            max_proba = np.max(prob, axis=1)

            # Select samples where confidence is above the threshold
            confident_labels = max_proba > self.threshold

            new_labels_idx = np.flatnonzero(~has_label)[confident_labels]

            # Add newly labeled confident predictions to the dataset
            self.y_labels_[new_labels_idx] = pred[confident_labels]
            has_label[new_labels_idx] = True
            self.y_labeled_iter_[new_labels_idx] = self.n_iter_

        self.base_estimator_.fit(X[safe_mask(X, has_label)],
                                 self.y_labels_[safe_mask(self.y_labels_,
                                                          has_label)])
        return self

    @if_delegate_has_method(delegate='base_estimator')
    def predict(self, X):
        """Predict on a dataset.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples, 1)
            array with predicted labels
        """
        check_is_fitted(self, 'y_labeled_iter_')
        X = check_array(X)
        return self.base_estimator_.predict(X)

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
        return self.base_estimator_.predict_proba(X)

    @if_delegate_has_method(delegate='base_estimator')
    def decision_function(self, X):
        """Calls decision function of the base_estimator

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        Returns
        -------
        y : array-like, shape = (n_samples, n_features)
            result of the decision function of the base_estimator
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_estimator_.decision_function(X)

    @if_delegate_has_method(delegate='base_estimator')
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
        return self.base_estimator_.predict_log_proba(X)

    @if_delegate_has_method(delegate='base_estimator')
    def score(self, X, y):
        """Calls score on the base_estimator

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            array representing the data

        y : array-like, shape = (n_samples, 1)
            array representing the labels

        Returns
        -------
        score : float
            result of calling score on the base_estimator
        """
        check_is_fitted(self, 'y_labeled_iter_')
        return self.base_estimator_.score(X, y)
