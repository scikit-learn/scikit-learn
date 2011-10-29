# Authors: Peter Prettenhofer
#
# License: BSD Style.
#
# TODO:
#      * Multi-class classification
#      * Huber loss for regression
#

from __future__ import division
import numpy as np

from ..base import BaseEstimator
from ..base import ClassifierMixin
from ..base import RegressorMixin
from ..utils import check_random_state

from ..tree.tree import _build_tree
from ..tree._tree import apply_tree
from ..tree._tree import MSE
from ..tree._tree import DTYPE


class VariableImportance(object):

    def __init__(self, n_features):
        self.variable_importance = np.zeros((n_features), dtype=np.float64)

    def visit_nonterminal_region(self, node):
        if node.is_leaf:
            return
        else:
            feature = node.feature
            error_improvement = (node.initial_error - node.best_error) \
                                / node.initial_error
            self.variable_importance[feature] += error_improvement ** 2.0

            # tail recursion
            self.visit_nonterminal_region(node.left)
            self.visit_nonterminal_region(node.right)


class MedianPredictor(object):
    """A simple initial estimator that predicts the median
    of the training targets.
    """

    median = None

    def fit(self, X, y):
        self.median = np.median(y)

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.median)
        return y


class MeanPredictor(object):
    """A simple initial estimator that predicts the mean
    of the training targets.
    """

    mean = None

    def fit(self, X, y):
        self.mean = np.mean(y)

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.mean)
        return y


class ClassPriorPredictor(object):
    """A simple initial estimator that predicts the mean
    of the training targets.
    """

    prior = None

    def fit(self, X, y):
        self.prior = np.log(np.sum(y) / np.sum(1.0 - y))

    def predict(self, X):
        y = np.empty((X.shape[0],), dtype=np.float64)
        y.fill(self.prior)
        return y


class LossFunction(object):
    """Abstract base class for various loss functions."""

    def init_estimator(self, X, y):
        pass

    def __call__(self, y, pred):
        pass

    def negative_gradient(self, y, pred):
        """Compute the negative gradient."""
        pass

    def update_terminal_regions(self, tree, X, y, residual, y_pred,
                                learn_rate=1.0):
        """Update the terminal regions (=leafs) of the given
        tree. Traverses tree and invokes template method
        `_update_terminal_region`. """
        if tree.is_leaf:
            self._update_terminal_region(tree, X, y, residual, y_pred)

            # update predictions
            y_pred[tree.terminal_region] += learn_rate * tree.value

            # Delete terminal_region to save memory
            del tree.terminal_region
            tree.terminal_region = None
        else:
            #print "%d: fx:%d, thres:%.8f" % (tree.id, tree.feature,
            #                                 tree.threshold)
            self.update_terminal_regions(tree.left, X, y, residual, y_pred,
                                         learn_rate)
            self.update_terminal_regions(tree.right, X, y, residual, y_pred,
                                         learn_rate)

    def _update_terminal_region(self, node, X, y, residual, pred):
        """Template method for updating terminal regions (=leafs). """
        pass


class LeastSquaresError(LossFunction):
    """Loss function for least squares (LS) estimation.
    Terminal regions need not to be updated for least squares. """

    def init_estimator(self):
        return MeanPredictor()

    def __call__(self, y, pred):
        return 0.5 * np.sum((y - pred) ** 2.0)

    def negative_gradient(self, y, pred):
        return y - pred


class LeastAbsoluteError(LossFunction):
    """Loss function for least absolute deviation (LAD) regression. """

    def init_estimator(self):
        return MedianPredictor()

    def __call__(self, y, pred):
        return np.abs(y - pred)

    def negative_gradient(self, y, pred):
        return np.sign(y - pred)

    def _update_terminal_region(self, node, X, y, residual, pred):
        """LAD updates terminal regions to median estimates. """
        node.value = np.asanyarray(
            np.median(y.take(node.terminal_region, axis=0) - \
                      pred.take(node.terminal_region, axis=0)))


## class HuberError(LossFunction):
##     """Loss function for least absolute deviation (LAD) regression."""

##     def init_estimator(self):
##         return MedianPredictor()

##     def __call__(self, y, pred):
##         return np.sum(np.abs(y - pred))

##     def negative_gradient(self, y, pred):
##         return np.sign(y - pred)

##     def _update_terminal_region(self, node, X, y, residual, pred):
##         """LAD updates terminal regions to median estimates. """
##         ## FIXME copied from LAD, still TODO
##         node.value = np.asanyarray(np.median(y.take(node.terminal_region, axis=0) - \
##                                              pred.take(node.terminal_region, axis=0)))


class BinomialDeviance(LossFunction):

    def init_estimator(self):
        return ClassPriorPredictor()

    def __call__(self, y, pred):
        """Compute the deviance (= negative log-likelihood). """
        return -2.0 * np.sum((y * pred) -
                             np.log(1.0 + np.exp(pred))) / y.shape[0]

    def negative_gradient(self, y, pred):
        return y - (1.0 / (1.0 + np.exp(-pred)))

    def _update_terminal_region(self, node, X, y, residual, pred):
        """Make a single Newton-Raphson step. """
        residual = residual.take(node.terminal_region, axis=0)
        y = y.take(node.terminal_region, axis=0)

        numerator = residual.sum()
        denominator = np.sum((y - residual) * (1.0 - y + residual))

        if denominator == 0.0:
            node.value = np.array(0.0, dtype=np.float64)
        else:
            node.value = np.asanyarray(numerator / denominator,
                                       dtype=np.float64)


LOSS_FUNCTIONS = {'ls': LeastSquaresError,
                  'lad': LeastAbsoluteError,
                  'deviance': BinomialDeviance}


class BaseGradientBoosting(BaseEstimator):

    trees = []

    def __init__(self, loss, learn_rate, n_iter, min_split, max_depth, init,
                 subsample, random_state):
        if n_iter <= 0:
            raise ValueError("n_iter must be greater than 0")
        self.n_iter = n_iter

        if learn_rate <= 0.0:
            raise ValueError("learn_rate must be greater than 0")
        self.learn_rate = learn_rate

        if loss not in LOSS_FUNCTIONS:
            raise ValueError("loss not supported")
        self.loss = LOSS_FUNCTIONS[loss]()

        if min_split <= 0:
            raise ValueError("min_split must be larger than 0.")
        self.min_split = min_split

        if subsample <= 0.0 or subsample > 1:
            raise ValueError("subsample must be in (0,1]")
        self.subsample = subsample

        if max_depth <= 0:
            raise ValueError("max_depth must be larger than 0.")
        self.max_depth = max_depth

        if init is None:
            self.init = self.loss.init_estimator()
        else:
            if not hasattr(init, 'fit') or not hasattr(init, 'predict'):
                raise ValueError("init must be valid estimator")
            self.init = init

        self.random_state = check_random_state(random_state)

    def fit(self, X, y, monitor=None):
        """Fit the gradient boosting model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes
            0, 1, ..., n_classes-1

        Returns
        -------
        self : object
            Returns self.
        """
        X = np.asfortranarray(X, dtype=DTYPE)
        y = np.ascontiguousarray(y)

        n_samples, n_features = X.shape
        if y.shape[0] != n_samples:
            raise ValueError("Number of labels does not match " \
                             "number of samples.")
        self.n_features = n_features

        # create argsorted X for fast tree induction
        X_argsorted = np.asfortranarray(
            np.argsort(X.T, axis=1).astype(np.int32).T)

        # fit initial model
        self.init.fit(X, y)

        # init predictions
        y_pred = self.init.predict(X)

        self.trees = []
        loss = self.loss
        #oob = np.zeros((self.n_iter), dtype=np.float64)
        #obj_func = np.zeros((self.n_iter), dtype=np.float64)

        # perform boosting iterations
        for i in xrange(self.n_iter):

            # subsampling
            sample_mask = np.random.rand(n_samples) > (1.0 - self.subsample)

            residual = loss.negative_gradient(y, y_pred)

            # induce regression tree on residuals
            tree = _build_tree(False, X, residual, MSE(), self.max_depth,
                               self.min_split, None, 1, self.random_state,
                               0.0, sample_mask, X_argsorted, True)

            assert tree.is_leaf != True

            # update tree leafs
            loss.update_terminal_regions(tree, X, y, residual, y_pred,
                                         self.learn_rate)

            # add tree to ensemble
            self.trees.append(tree)

            # update OOB predictions
            if self.subsample < 1.0:
                y_pred[~sample_mask] = self._predict(
                    X[~sample_mask], old_pred=y_pred[~sample_mask])

            if monitor:
                monitor(self, i)

            #if self.subsample < 1.0:
            #    oob[i] = loss.loss(y[~sample_mask], y_pred[~sample_mask])
            #print "Iteration %d - in %fs" % (i, time() - t0)
            #print "Obj_func: %.4f     OOB: %.4f" % (obj_func[i], oob[i])

        #self.obj_func = obj_func
        #self.oob = oob

        return self

    def _predict(self, X, old_pred=None):
        """Predict targets with current model. Re-uses predictions
        from previous iteration if available.
        """
        if old_pred is not None:
            return old_pred + self.learn_rate * apply_tree(self.trees[-1],
                                                           X, 1).ravel()
        else:
            y = self.init.predict(X)
            for tree in self.trees:
                y += self.learn_rate * apply_tree(tree, X, 1).ravel()
            return y

    @property
    def variable_importance(self):
        if not self.trees or len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `variable_importance`.")
        variable_importance = np.zeros((self.n_features,), dtype=np.float64)
        for tree in self.trees:
            vi = VariableImportance(self.n_features)
            vi.visit_nonterminal_region(tree)
            variable_importance += vi.variable_importance
        variable_importance /= len(self.trees)
        variable_importance = 100.0 * (variable_importance /
                                       variable_importance.max())
        return variable_importance


class GradientBoostingClassifier(BaseGradientBoosting, ClassifierMixin):
    """Gradient Boosting for classification. GB builds an additive model in a
    forward stage-wise fashion; it allows for the optimization of
    arbitrary differentiable loss functions. In each stage a regression
    tree is fit on the negative gradient of binomial or multinomial
    deviance.

    Parameters
    ----------
    n_iter : int
        The number of boosting stages to perform. Gradient boosting
        is fairly robust to over-fitting so a large number usually
        results in better performance.

    learn_rate : float, optional (default=0.1)
        learning rate shrinks the contribution of each tree by `learn_rate`.
        There is a trade-off between learn_rate and n_iter (see Discussion).

    max_depth : integer, optional (default=3)
        maximum depth of the individual regression trees. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.

    min_split : integer, optional (default=1)
        minimum number of samples required at any leaf node. Use larger
        value if number of samples is large.

    subsample : float, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. `subsample` interacts with the parameter `n_iter`
        (see Discussion).

    Examples
    --------
    >>> samples = [[0, 0, 2], [1, 0, 0]]
    >>> labels = [0, 1]
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> gb = GradientBoostingRegressor(10)
    >>> gb.fit(samples, labels)
    GradientBoostingRegressor()
    >>> print gb.predict([[0, 0, 0]])
    [1]

    See also
    --------
    DecisionTreeRegressor, RandomForestRegressor

    Discussion
    ----------
    The optimal algorithm for a given dataset is a complicated choice, and
    depends on a number of factors:
    * n_iter vs. learn_rate
        TODO
    * n_iter vs. subsample
        TODO

    References
    ----------
    J. Friedman, Greedy Function Approximation: A Gradient Boosting
    Machine, The Annals of Statistics, Vol. 29, No. 5, 2001.

    J. Friedman, Stochastic Gradient Boosting, 1999

    T. Hastie, R. Tibshirani and J. Friedman.
    Elements of Statistical Learning Ed. 2, Springer, 2009.
    """

    def __init__(self, loss='deviance', learn_rate=0.1, n_iter=100,
                 subsample=1.0, min_split=1, max_depth=3,
                 init=None, random_state=None):

        super(GradientBoostingClassifier, self).__init__(
            loss, learn_rate, n_iter, min_split, max_depth, init, subsample,
            random_state)

    def fit(self, X, y, monitor=None):
        self.classes = np.unique(y)
        if self.classes.shape[0] != 2:
            raise ValueError("only binary classification supported")
        y = np.searchsorted(self.classes, y)
        super(GradientBoostingClassifier, self).fit(X, y, monitor=monitor)

    def predict(self, X):
        X = np.atleast_2d(X)
        X = X.astype(DTYPE)
        if len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `predict`.")
        f = self._predict(X)
        return self.classes.take(f >= 0.0, axis=0)

    def predict_proba(self, X):
        X = np.atleast_2d(X)
        X = X.astype(DTYPE)
        if len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `predict_proba`.")
        f = self._predict(X)
        P = np.ones((X.shape[0], 2), dtype=np.float64)
        P[:, 1] = 1.0 / (1.0 + np.exp(-f))
        P[:, 0] -= P[:, 1]
        return P


class GradientBoostingRegressor(BaseGradientBoosting, RegressorMixin):
    """Gradient Boosting for regression. GB builds an additive model in a
    forward stage-wise fashion; it allows for the optimization of
    arbitrary differentiable loss functions. In each stage a regression
    tree is fit on the negative gradient of the given loss function.

    Parameters
    ----------
    n_iter : int
        The number of boosting stages to perform. Gradient boosting
        is fairly robust to over-fitting so a large number usually
        results in better performance.

    learn_rate : float, optional (default=0.1)
        learning rate shrinks the contribution of each tree by `learn_rate`.
        There is a trade-off between learn_rate and n_iter (see Discussion).

    loss : {'ls', 'lad'}, optional
        loss function to be optimized. 'ls' refers to least squares
        regression. 'lad' (least absolute deviation) is a highly robust
        loss function soley based on order information of the input
        variables.

    max_depth : integer, optional (default=3)
        maximum depth of the individual regression trees. The maximum
        depth limits the number of nodes in the tree. Tune this parameter
        for best performance; the best value depends on the interaction
        of the input variables.

    min_split : integer, optional (default=1)
        minimum number of samples required at any leaf node. Use larger
        value if number of samples is large.

    subsample : float, optional (default=1.0)
        The fraction of samples to be used for fitting the individual base
        learners. If smaller than 1.0 this results in Stochastic Gradient
        Boosting. `subsample` interacts with the parameter `n_iter`
        (see Discussion).

    Examples
    --------
    >>> samples = [[0, 0, 2], [1, 0, 0]]
    >>> labels = [0, 1]
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> gb = GradientBoostingRegressor(10)
    >>> gb.fit(samples, labels)
    GradientBoostingRegressor()
    >>> print gb.predict([[0, 0, 0]])
    [1]

    See also
    --------
    DecisionTreeRegressor, RandomForestRegressor

    Discussion
    ----------
    The optimal algorithm for a given dataset is a complicated choice, and
    depends on a number of factors:
    * n_iter vs. learn_rate
        TODO
    * n_iter vs. subsample
        TODO

    References
    ----------
    J. Friedman, Greedy Function Approximation: A Gradient Boosting
    Machine, The Annals of Statistics, Vol. 29, No. 5, 2001.

    J. Friedman, Stochastic Gradient Boosting, 1999

    T. Hastie, R. Tibshirani and J. Friedman.
    Elements of Statistical Learning Ed. 2, Springer, 2009.
    """

    def __init__(self, loss='ls', learn_rate=0.1, n_iter=100, subsample=1.0,
                 min_split=1, max_depth=3, init=None, random_state=None):

        super(GradientBoostingRegressor, self).__init__(
            loss, learn_rate, n_iter, min_split, max_depth, init, subsample,
            random_state)

    def predict(self, X):
        X = np.atleast_2d(X)
        X = X.astype(DTYPE)
        if len(self.trees) == 0:
            raise ValueError("Estimator not fitted, " \
                             "call `fit` before `predict`.")
        y = self._predict(X)
        return y
