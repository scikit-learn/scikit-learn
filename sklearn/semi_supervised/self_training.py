import numpy as np

from ..base import BaseEstimator
from ..utils.validation import check_X_y, check_array, check_is_fitted
from ..utils import safe_mask


def _check_estimator(estimator):
    """Make sure that an estimator implements the necessary methods."""
    if not hasattr(estimator, "predict_proba"):
        raise ValueError("The base estimator should implement predict_proba!")


class SelfTraining(BaseEstimator):

    """Self-training classifier

    Parameters
    ----------
    base : estimator object
        An estimator object implementing `fit` and `predict_proba`.

    threshold : float
        Threshold above which predictions are added to the labeled dataset

    max_iter : integer
        Maximum number of iterations allowed

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn import datasets
    >>> from sklearn.semi_supervised import SelfTraining
    >>> from sklearn.svm import SVC
    >>> svc = SVC(probability=True)
    >>> self_training_model = SelfTraining(svc)
    >>> iris = datasets.load_iris()
    >>> rng = np.random.RandomState(42)
    >>> random_unlabeled_points = rng.rand(len(iris.target)) < 0.3
    >>> labels = np.copy(iris.target)
    >>> labels[random_unlabeled_points] = -1
    >>> self_training_model.fit(iris.data, labels)
    ... # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    SVC(...)

    References
    ----------
    David Yarowsky. 1995. Unsupervised word sense disambiguation rivaling
    supervised methods. In Proceedings of the 33rd annual meeting on
    Association for Computational Linguistics (ACL '95). Association for
    Computational Linguistics, Stroudsburg, PA, USA, 189-196. DOI:
    https://doi.org/10.3115/981658.981684
    """
    def __init__(self, base, threshold=0.75, max_iter=100):
        self.base = base
        self.threshold = threshold
        self.max_iter = max_iter

    def fit(self, X, y):
        """
        Fits SelfTraining Estimator to dataset

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
        _check_estimator(self.base)

        # Data usable for supervised training
        X_labeled = X[safe_mask(X, np.where(y != -1))][0]
        y_labeled = y[safe_mask(y, np.where(y != -1))][0]

        # Unlabeled data
        X_unlabeled = X[safe_mask(X, np.where(y == -1))][0]
        y_unlabeled = y[safe_mask(y, np.where(y == -1))][0]

        iter = 0
        while (len(X_labeled) < len(X) and iter < self.max_iter):
            iter += 1
            self.base.fit(X_labeled, y_labeled)

            # Select prediction where confidence is above the threshold
            pred = self.predict(X_unlabeled)
            max_proba = np.max(self.predict_proba(X_unlabeled), axis=1)
            confident = np.where(max_proba > self.threshold)[0]

            # Add newly labeled confident predictions to the dataset
            X_labeled = np.append(X_labeled, X_unlabeled[confident], axis=0)
            y_labeled = np.append(y_labeled, pred[confident], axis=0)

            # Remove already labeled data from unlabeled dataset
            X_unlabeled = np.delete(X_unlabeled, confident, axis=0)
            y_unlabeled = np.delete(y_unlabeled, confident, axis=0)

        self.base.fit(X_labeled, y_labeled)
        return self.base

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
        check_is_fitted(self, 'base')
        X = check_array(X)
        return self.base.predict(X)

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
        _check_estimator(self.base)
        check_is_fitted(self, 'base')
        return self.base.predict_proba(X)
