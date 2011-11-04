# Authors: Gilles Louppe,
#          Brian Holt

import numpy as np

from ..base import clone
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor, \
                   ExtraTreeClassifier, ExtraTreeRegressor
from ..utils import check_random_state

__all__ = [
    "RandomForestClassifier",
    "RandomForestRegressor",
    "ExtraTreesClassifier",
    "ExtraTreesRegressor"
]


class Forest(BaseEstimator):
    def __init__(self, base_tree, n_trees=10, bootstrap=False, random_state=None):
        self.base_tree = base_tree
        self.n_trees = n_trees
        self.bootstrap = bootstrap
        self.random_state = check_random_state(random_state)
        self.forest = []

    def fit(self, X, y):
        X = np.atleast_2d(X)
        y = np.atleast_1d(y)

        if isinstance(self.base_tree, ClassifierMixin):
            self.classes = np.unique(y)
            self.n_classes = len(self.classes)
            y = np.searchsorted(self.classes, y)

        for i in xrange(self.n_trees):
            tree = clone(self.base_tree)

            if self.bootstrap:
                n_samples = X.shape[0]
                indices = np.random.randint(0, n_samples, n_samples)

                X = X[indices]
                y = y[indices]

            tree.fit(X, y)
            self.forest.append(tree)

class ForestClassifier(Forest, ClassifierMixin):
    def __init__(self, base_tree, n_trees=10, bootstrap=False, random_state=None):
        Forest.__init__(self, base_tree,
                              n_trees,
                              bootstrap=bootstrap,
                              random_state=random_state)

    def predict(self, X):
        return self.classes[np.argmax(self.predict_proba(X), axis=1)]

    def predict_proba(self, X):
        X = np.atleast_2d(X)
        p = np.zeros((X.shape[0], self.n_classes))

        for tree in self.forest:
            p += tree.predict_proba(X)

        p /= self.n_trees

        return p

class ForestRegressor(Forest, RegressorMixin):
    def __init__(self, base_tree, n_trees=10, bootstrap=False, random_state=None):
        Forest.__init__(self, base_tree,
                              n_trees,
                              bootstrap=bootstrap,
                              random_state=random_state)

    def predict(self, X):
        X = np.atleast_2d(X)
        y_hat = np.zeros(X.shape[0])

        for tree in self.forest:
            y_hat += tree.predict(X)

        y_hat /= self.n_trees

        return y_hat

class RandomForestClassifier(ForestClassifier):
    def __init__(self, n_trees=10, bootstrap=True, random_state=None, **tree_args):
        ForestClassifier.__init__(self, DecisionTreeClassifier(**tree_args),
                                        n_trees,
                                        bootstrap=bootstrap,
                                        random_state=random_state)

class RandomForestRegressor(ForestRegressor):
    def __init__(self, n_trees=10, bootstrap=True, random_state=None, **tree_args):
        ForestRegressor.__init__(self, DecisionTreeRegressor(**tree_args),
                                       n_trees,
                                       bootstrap=bootstrap,
                                       random_state=random_state)

class ExtraTreesClassifier(ForestClassifier):
    def __init__(self, n_trees=10, bootstrap=False, random_state=None, **tree_args):
        ForestClassifier.__init__(self, ExtraTreeClassifier(**tree_args),
                                        n_trees,
                                        bootstrap=bootstrap,
                                        random_state=random_state)

class ExtraTreesRegressor(ForestRegressor):
    def __init__(self, n_trees=10, bootstrap=False, random_state=None, **tree_args):
        ForestRegressor.__init__(self, ExtraTreeRegressor(**tree_args),
                                       n_trees,
                                       bootstrap=bootstrap,
                                       random_state=random_state)
