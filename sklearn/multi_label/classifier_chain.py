"""
Classifier Chain.

Use classifier chain as a meta-algorithm to
combine classifiers for each label.
"""

# Author: Chih-Wei Chang <jrweizhang@gmail.com>

import numpy as np

from ..base import BaseEstimator


class ClassifierChain(BaseEstimator):
    """Classifier Chain

    Classifier chain [1] fit classifiers for each label,
    and chains the classifiers to make prediction.

    Parameters
    ----------
    base_estimator : object
        The base estimator used for fitting each label.

    Attributes
    ----------
    n_labels : int
        How many labels are there in this model.

    classifiers_ : array
        List of classifiers, which will be used to chain prediction.

    References
    ----------
    .. [1] Jesse Read, Bernhard Pfahringer, Geoff Holmes, Eibe Frank,
           "Classifier Chains for Multi-label Classification", 2009.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.datasets import make_multilabel_classification
    >>> from sklearn.multi_label import ClassifierChain
    >>> from sklearn.svm import LinearSVC
    >>> X, Y = make_multilabel_classification(return_indicator=True, random_state=0)
    >>> cc = ClassifierChain(base_estimator=LinearSVC)
    >>> cc.fit(X, Y)
    ClassifierChain(base_estimator=<class 'sklearn.svm.classes.LinearSVC'>)
    """

    def __init__(self, base_estimator):
        self.base_estimator = base_estimator
        self.classifiers_ = []

    def fit(self, X, Y):
        """ Build a sequence of classifiers

        Parameters
        ----------
        X : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        Y : array_like, shape (n, n_labels)
            List of n_labels-dimensional data points.  Each row
            corresponds to a label indicators.
        """
        X = np.asarray(X, dtype=np.float)
        Y = np.asarray(Y, dtype=np.int)

        self.n_labels_ = Y.shape[1]
        self.classifiers_ = []

        for i in xrange(self.n_labels_):
            y = Y[:, i]

            clf = self.base_estimator()
            clf.fit(X, y)
            self.classifiers_.append(clf)

            X = self._predict_and_chain(clf, X)

        return self

    def predict(self, X):
        """Predict label for data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = (n_samples,)
        """

        for clf in self.classifiers_:
            X = self._predict_and_chain(clf, X)

        return X[:, -self.n_labels:]

    def _predict_and_chain(self, clf, X):
        return np.hstack((X, clf.predict(X).reshape(-1, 1)))
