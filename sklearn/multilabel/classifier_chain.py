import numpy as np
from ..base import BaseEstimator, clone
from ..utils import check_random_state
from scipy.sparse import issparse


class ClassifierChain(BaseEstimator):
    """
    A multi-label model that arranges binary classifiers into a chain. Each
    model makes a prediction in the order specified by the chain using all
    of the available features provided to the model plus the predictions of
    models that are earlier in the chain.

    Parameters
    ----------
    base_estimator : estimator
        The base estimator from which the classifier chain is built.

    random_state : int or RandomState, optional, default None
        State or seed for random number generator.

    order : list of integers, None, or 'random'
        The order of the chain can be explicitly set by providing a list of
        integers. For example, for a chain of length 5
            order = [1, 3, 2, 4, 0]
        means that the first model in the chain will make predictions for
        column 1 in the Y matrix, the second model will make predictions
        for column 3, etc.

        If order is None the order will be determined by the
        order of columns in the label matrix Y.
            order = [0, 1, 2, ..., Y.shape[1] - 1]

        If order is 'random' a random ordering will be used.

    Attributes
    ----------
    estimators_ : list
        A list of copies of base_estimator. Once fit the estimators in
        this list will be ordered as specified by the order attribute.
        The original order of labels (as specified by Y) is recovered in
        the predict method by indexing into the predictions with the
        list of indices in the order attribute.


    References
    ----------
    Jesse Read, Bernhard Pfahringer, Geoff Holmes, Eibe Frank, "Classifier
    Chains for Multi-label Classification", 2009.

    """

    def __init__(self, base_estimator, random_state=None, order=None):
        self.base_estimator = base_estimator
        self.random_state = random_state
        self.order = order

    def fit(self, X, Y):
        """Fit the model to data matrix X and targets Y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Y : {array-like, sparse matrix}, shape (n_samples, n_classes)

        Returns
        -------
        self : object
            Returns self.
        """

        if issparse(Y):
            Y = Y.toarray()

        if issparse(X):
            X = X.toarray()

        random_state = check_random_state(self.random_state)

        if self.order is None:
            self.order = list(range(Y.shape[1]))
        elif self.order == 'random':
            self.order = list(range(Y.shape[1]))
            random_state.shuffle(self.order)
        elif sorted(self.order) != list(range(Y.shape[1])):
                raise ValueError("invalid order")

        self.estimators_ = [clone(self.base_estimator)
                            for _ in range(Y.shape[1])]

        for chain_idx, estimator in enumerate(self.estimators_):
            previous_labels = Y[:, self.order[:chain_idx]]
            y = Y[:, self.order[chain_idx]]
            X_aug = np.hstack((X, previous_labels))
            estimator.fit(X_aug, y)

    def predict(self, X):
        """Predict on the data matrix X using the ClassifierChain model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)

        Returns
        -------
        Y_pred-like : array-like, shape (n_samples, n_classes)
        """

        if issparse(X):
            X = X.toarray()

        Y_pred_chain = np.zeros((X.shape[0], len(self.estimators_)))
        for chain_idx, estimator in enumerate(self.estimators_):
            previous_predictions = Y_pred_chain[:, :chain_idx]
            X_aug = np.hstack((X, previous_predictions))
            Y_pred_chain[:, chain_idx] = estimator.predict(X_aug)
        chain_key = [self.order.index(i) for i in range(len(self.order))]
        Y_pred = Y_pred_chain[:, chain_key]

        return Y_pred
