import numpy as np
from ..base import BaseEstimator, clone
from sklearn.linear_model import LogisticRegression
import random


class ClassifierChain(BaseEstimator):
    """Classifier Chain
    A multi-label model that arranges classifiers into a chain. Predictions are made
    in the order specified by the chain. The predictions of earlier models are used
    as features by later models.

    By default the order of the chain is random although it can be specified at the time of fitting with the
    chain_order parameter. Since the optimal chain_order is not known a priori it is common to use the mean
    prediction of an ensemble randomly ordered classifier chains.

    Parameters
    ----------
    base_estimator : object
        The base estimator used for fitting the model for each label. Defaults


    Attributes
    ----------
    classifiers_ : array
        List of classifiers, which will be used to chain prediction.
    chain_order : list of ints
        A list of integers specifying the order of the classes in the chain. For example, for a chain of length 5
            chain_order = [1, 3, 2, 4, 0]
        means that the first model in the chain will make predictions for column 1 in the Y matrix, the second model
        will make predictions for column 3, etc.
        If chain_order is not None it must have a length equal to the number of columns in Y.
        If chain_order is None a random sequence of integers will be generated.
        To fit a ClassifierChain with the chain order corresponding to the order of classes in the Y matrix use
            chain_order = range(Y.shape[1])


    References
    ----------
    Jesse Read, Bernhard Pfahringer, Geoff Holmes, Eibe Frank,
        "Classifier Chains for Multi-label Classification", 2009.

    """

    def __init__(self, base_estimator=LogisticRegression(), random_state=None):
        self.base_estimator = base_estimator
        random.seed(random_state)
        self.classifiers = []
        self.chain_order = None

    def fit(self, X, Y, chain_order = None):
        self.classifiers = [clone(self.base_estimator) for _ in range(Y.shape[1])]

        if chain_order is not None:
            assert len(chain_order) == Y.shape[1], "chain_order length must equal n_labels"
            assert all([isinstance(i, int) for i in chain_order]), "chain_order must be a list of integers"
            self.chain_order = chain_order
        else:
            self.chain_order = range(Y.shape[1])
            random.shuffle(self.chain_order)

        for chain_idx, classifier in enumerate(self.classifiers):
            X_augmented = np.hstack((X, Y[:, self.chain_order[:chain_idx]]))
            y = np.squeeze(Y[:, self.chain_order[chain_idx]])
            classifier.fit(X_augmented, y)

    def predict(self, X):
        Y_pred_chain = np.zeros((X.shape[0], len(self.classifiers)))
        for chain_idx, classifier in enumerate(self.classifiers):
            X_augmented = np.hstack((X, Y_pred_chain[:, self.chain_order[:chain_idx]]))
            Y_pred_chain[:, chain_idx] = classifier.predict_proba(X_augmented)[:, 1]

        chain_key = [self.chain_order.index(i) for i in range(len(self.chain_order))]
        Y_pred = Y_pred_chain[:, chain_key]

        return Y_pred