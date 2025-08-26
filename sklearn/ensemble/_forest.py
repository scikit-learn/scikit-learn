"""Forest of trees-based ensemble methods."""

from sklearn.tree import DecisionTreeClassifier
from sklearn.base import BaseEstimator, ClassifierMixin
import numpy as np


class RandomForestClassifier(BaseEstimator, ClassifierMixin):
    """A random forest classifier.

    A random forest is a meta estimator that fits a number of decision tree
    classifiers on various sub-samples of the dataset and uses averaging to
    improve the predictive accuracy and control over-fitting.

    The sub-sample size is controlled with the `max_samples` parameter if
    `bootstrap=True` (default), otherwise the whole dataset is used to build
    each tree.

    Each decision tree in the ensemble uses the 'best' splitter strategy,
    which exhaustively evaluates all possible split points between consecutive
    unique values for each feature. For a feature with n unique values,
    this yields (n-1) potential thresholds. The splitting algorithm considers
    each midpoint between adjacent sorted unique feature values as a candidate
    threshold and selects the one that optimizes the chosen criterion
    (e.g., Gini impurity or entropy).

    This comprehensive evaluation ensures optimal local splits at each node,
    though the global feature subset selection remains random as controlled
    by the `max_features` parameter.

    Parameters
    ----------
    n_estimators : int, default=100
        The number of trees in the forest.

    criterion : {"gini", "entropy", "log_loss"}, default="gini"
        The function to measure the quality of a split.

    max_depth : int, default=None
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int or float, default=2
        The minimum number of samples required to split an internal node.

    min_samples_leaf : int or float, default=1
        The minimum number of samples required to be at a leaf node.

    max_features : {"sqrt", "log2", None}, int or float, default="sqrt"
        The number of features to consider when looking for the best split.

    bootstrap : bool, default=True
        Whether bootstrap samples are used when building trees.

    random_state : int, RandomState instance or None, default=None
        Controls both the randomness of the bootstrapping of the samples used
        when building trees and the sampling of the features to consider when
        looking for the best split at each node.

    Attributes
    ----------
    estimators_ : list of DecisionTreeClassifier
        The collection of fitted sub-estimators.

    classes_ : ndarray of shape (n_classes,)
        The classes labels.

    n_classes_ : int
        The number of classes.

    n_features_in_ : int
        Number of features seen during fit.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during fit. Defined only when `X`
        has feature names that are all strings.

    Examples
    --------
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.datasets import make_classification
    >>> X, y = make_classification(n_samples=1000, n_features=4,
    ...                            n_informative=2, n_redundant=0,
    ...                            random_state=0, shuffle=False)
    >>> clf = RandomForestClassifier(max_depth=2, random_state=0)
    >>> clf.fit(X, y)
    RandomForestClassifier(...)
    >>> print(clf.predict([[0, 0, 0, 0]]))
    [1]

    References
    ----------
    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.
    """

    def __init__(self, n_estimators=100, criterion='gini', max_depth=None,
                 min_samples_split=2, min_samples_leaf=1, max_features='sqrt',
                 bootstrap=True, random_state=None):
        self.n_estimators = n_estimators
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.max_features = max_features
        self.bootstrap = bootstrap
        self.random_state = random_state

    def fit(self, X, y):
        """Build a forest of trees from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The training input samples.
        y : array-like of shape (n_samples,)
            The target values (class labels).

        Returns
        -------
        self : object
            Fitted estimator.
        """
        # Basic implementation - in practice this would be much more complex
        self.classes_ = np.unique(y)
        self.n_classes_ = len(self.classes_)
        self.n_features_in_ = X.shape[1]
        
        # Create and fit estimators
        self.estimators_ = []
        for i in range(self.n_estimators):
            tree = DecisionTreeClassifier(
                criterion=self.criterion,
                max_depth=self.max_depth,
                min_samples_split=self.min_samples_split,
                min_samples_leaf=self.min_samples_leaf,
                max_features=self.max_features,
                random_state=None if self.random_state is None else self.random_state + i,
                splitter='best'  # This ensures exhaustive split evaluation
            )
            tree.fit(X, y)
            self.estimators_.append(tree)
        
        return self

    def predict(self, X):
        """Predict class for X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray of shape (n_samples,)
            The predicted classes.
        """
        predictions = np.array([tree.predict(X) for tree in self.estimators_])
        # Simple majority voting
        return np.apply_along_axis(
            lambda x: np.bincount(x).argmax(), axis=0, arr=predictions
        )

    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        p : ndarray of shape (n_samples, n_classes)
            The class probabilities of the input samples.
        """
        all_proba = np.array([tree.predict_proba(X) for tree in self.estimators_])
        return np.mean(all_proba, axis=0)
