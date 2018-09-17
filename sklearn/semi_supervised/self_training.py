import numpy as np

from ..base import BaseEstimator, clone
from ..utils.validation import check_X_y, check_array, check_is_fitted
from ..utils.metaestimators import if_delegate_has_method
from ..utils import safe_mask

__all__ = ["SelfTrainingClassifier"]

# Authors:  Oliver Rausch   <rauscho@ethz.ch>
#           Patrice Becker  <beckerp@ethz.ch>
# License: BSD 3 clause

def _check_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""

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
    base_estimator_: estimator object
        The fitted estimator.

    y_labeled_ : array, shape = (n_samples)
        The labels assigned to unlabeled datapoints during fitting.

    X_labeled_ : array, shape = (n_samples, n_features)
        The labeled samples used for the final fit.

    y_labeled_iter_ : array, shape = (n_samples)
        The iteration in which each sample was labeled. When a sample has
        iteration 0, the sample was labeled in the given dataset.

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

        y : array-like, shape = (n_samples, 1)
            array representing the labels

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

        # Data usable for supervised training
        self.X_labeled_ = X[safe_mask(X, np.where(y != -1))][0]
        self.y_labeled_ = y[safe_mask(y, np.where(y != -1))][0]
        self.y_labeled_iter_ = np.full_like(self.y_labeled_, 0)

        # Unlabeled data
        X_unlabeled = X[safe_mask(X, np.where(y == -1))][0]
        y_unlabeled = y[safe_mask(y, np.where(y == -1))][0]

        iter = 0
        while len(self.X_labeled_) < len(X) and iter < self.max_iter:
            iter += 1
            self.base_estimator_.fit(self.X_labeled_, self.y_labeled_)

            # Select predictions where confidence is above the threshold
            predict_proba = self.base_estimator_.predict_proba(X_unlabeled)

            pred = np.argmax(predict_proba, axis=1)
            max_proba = np.max(predict_proba, axis=1)

            confident = np.where(max_proba > self.threshold)[0]

            # Add newly labeled confident predictions to the dataset
            self.X_labeled_ = np.append(
                self.X_labeled_, X_unlabeled[confident], axis=0)
            self.y_labeled_ = np.append(
                self.y_labeled_, pred[confident], axis=0)
            self.y_labeled_iter_ = np.append(
                self.y_labeled_iter_, np.full_like(confident, iter))

            # Remove already labeled data from unlabeled dataset
            X_unlabeled = np.delete(X_unlabeled, confident, axis=0)
            y_unlabeled = np.delete(y_unlabeled, confident, axis=0)

        self.base_estimator_.fit(self.X_labeled_, self.y_labeled_)
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
