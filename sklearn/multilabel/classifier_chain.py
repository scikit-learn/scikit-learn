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

    By default the order of the chain is random although it can be specified
    at the time of fitting with the chain_order parameter. Since the optimal
    chain_order is not known a priori it is common to use the mean prediction
    of an ensemble randomly ordered classifier chains.

    References
    ----------
    Jesse Read, Bernhard Pfahringer, Geoff Holmes, Eibe Frank, "Classifier
    Chains for Multi-label Classification", 2009.

    """

    def __init__(self, base_estimator, random_state=None,
                 order=None, shuffle=True):
        """
        Parameters
        ----------
        base_estimator : estimator
            The base estimator from which the classifier chain is built.

        random_state : int or RandomState, optional, default None
            State or seed for random number generator.

        order : list
            A list of integers specifying the order of the classes in the
            chain (0 indexed). For example, for a chain of length 5
                order = [1, 3, 2, 4, 0]
            means that the first model in the chain will make predictions for
            column 1 in the Y matrix, the second model will make predictions
            for column 3, etc. If chain_order is not None it must have a length
            equal to the number of columns in Y.
            If chain_order is None an ordered list of integers will be used
                chain_order = [0, 1, 2, ..., Y.shape[1]]
            where Y is the label matrix passed to the fit method.

        shuffle : bool
            If true order is shuffled
        """
        self.base_estimator = base_estimator
        self.random_state = random_state
        if order is not None:
            if not isinstance(order, list):
                order = list(order)
            if not all([isinstance(i, int) for i in order]):
                raise ValueError("all elements of order parameter must be "
                                 "integers")

        self.order = order

        self.shuffle = shuffle

    def fit(self, X, Y):
        """Fit the ClassifierChain model according to the given training (X)
        and the true labels (Y) of classes that appear earlier in the chain as
        defined by the order parameter.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
        Y : array, shape (n_samples, n_classes)

        Attributes
        ----------
        classifiers_ : list
            A list of copies of base_estimator. Once fit the classifiers in
            this list will be ordered as specified by the order attribute.
            The original order of labels (as specified by Y) is recovered in
            the predict method by indexing into the predictions with the
            list of indices in the order attribute.

        Returns
        -------
        self : object
            Returns self.
        """

        if issparse(Y):
            Y = Y.toarray()

        random_state = check_random_state(self.random_state)
        self.classifiers_ = [clone(self.base_estimator)
                             for _ in range(Y.shape[1])]

        if self.order is not None:
            if not len(self.order) == Y.shape[1]:
                raise ValueError("chain_order length must equal n_labels")
        else:
            self.order = list(range(Y.shape[1]))

        if self.shuffle:
            random_state.shuffle(self.order)

        for chain_idx, classifier in enumerate(self.classifiers_):
            previous_labels = Y[:, self.order[:chain_idx]]
            previous_labels = previous_labels.toarray()
            y = Y[:, self.order[chain_idx]]
            y = y.toarray()[:, 0]

            X_aug = np.hstack((X, previous_labels))
            classifier.fit(X_aug, y)

    def predict(self, X):
        """For each row in X make a prediction for each of the models in
        the chain using the features available in X plus the predictions of
        previous models in the chain as defined by the order parameter.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)

        Returns
        -------
        Y_pred : array, shape (n_samples, n_classes)
        """

        Y_pred_chain = np.zeros((X.shape[0], len(self.classifiers_)))
        for chain_idx, classifier in enumerate(self.classifiers_):
            previous_predictions = Y_pred_chain[:, :chain_idx]
            X_aug = np.hstack((X, previous_predictions))
            Y_pred_chain[:, chain_idx] = classifier.predict(X_aug)
        chain_key = [self.order.index(i) for i in range(len(self.order))]
        Y_pred = Y_pred_chain[:, chain_key]

        return Y_pred
