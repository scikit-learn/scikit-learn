"""
Tune the parameters of an estimator by cross-validation.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux    <gael.varoquaux@normalesup.org>
# License: BSD Style.

import copy

from .externals.joblib import Parallel, delayed
from .cross_val import KFold, StratifiedKFold
from .base import BaseEstimator, is_classifier, clone

try:
    from itertools import product
except:
    def product(*args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)


def iter_grid(param_grid):
    """ Generators on the combination of the various parameter lists given.

        Parameters
        -----------
        kwargs: keyword arguments, lists
            Each keyword argument must be a list of values that should
            be explored.

        Returns
        --------
        params: dictionary
            Dictionnary with the input parameters taking the various
            values succesively.

        Examples
        ---------
        >>> from scikits.learn.grid_search import iter_grid
        >>> param_grid = {'a':[1, 2], 'b':[True, False]}
        >>> list(iter_grid(param_grid))
        [{'a': 1, 'b': True}, {'a': 1, 'b': False}, {'a': 2, 'b': True}, {'a': 2, 'b': False}]

    """
    if hasattr(param_grid, 'has_key'):
        param_grid = [param_grid]
    for p in param_grid:
        # Always sort the keys of a dictionary, for reproducibility
        items = sorted(p.items())
        keys, values = zip(*items)
        for v in product(*values):
            params = dict(zip(keys, v))
            yield params


def fit_grid_point(X, y, base_clf, clf_params, cv, loss_func, iid,
                   **fit_params):
    """Run fit on one set of parameters
    Returns the score and the instance of the classifier
    """
    # update parameters of the classifier after a copy of its base structure
    clf = copy.deepcopy(base_clf)
    clf._set_params(**clf_params)

    score = 0.
    n_test_samples = 0.
    for train, test in cv:
        clf.fit(X[train], y[train], **fit_params)
        y_test = y[test]
        if loss_func is not None:
            y_pred = clf.predict(X[test])
            this_score = -loss_func(y_test, y_pred)
        else:
            this_score = clf.score(X[test], y_test)
        if iid:
            this_score *= len(y_test)
            n_test_samples += len(y_test)
        score += this_score
    if iid:
        score /= n_test_samples

    return score, clf


################################################################################
class GridSearchCV(BaseEstimator):
    """
    Grid search on the parameters of a classifier.

    Important members are fit, predict.

    GridSearchCV implements a "fit" method and a "predict" method like
    any classifier except that the parameters of the classifier
    used to predict is optimized by cross-validation

    Parameters
    ----------
    estimator: object type that implements the "fit" and "predict" methods
        A object of that type is instanciated for each grid point

    param_grid: dict
        a dictionary of parameters that are used the generate the grid

    loss_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediciton (small is good)
        if None is passed, the score of the estimator is maximized

    fit_params : dict, optional
        parameters to pass to the fit method

    n_jobs: int, optional
        number of jobs to run in parallel (default 1)

    iid: boolean, optional
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    Methods
    -------
    fit(X, Y) : self
        Fit the model

    predict(X) : array
        Predict using the model.

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn.cross_val import LeaveOneOut
    >>> from scikits.learn.svm import SVR
    >>> from scikits.learn.grid_search import GridSearchCV
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> parameters = {'kernel':('linear', 'rbf'), 'C':[1, 10]}
    >>> svr = SVR()
    >>> clf = GridSearchCV(svr, parameters, n_jobs=1)
    >>> clf.fit(X, y).predict([[-0.8, -1]])
    array([ 1.])
    """

    def __init__(self, estimator, param_grid, loss_func=None,
                        fit_params={}, n_jobs=1, iid=True):
        assert hasattr(estimator, 'fit') and hasattr(estimator, 'predict'), (
            "estimator should a be an estimator implementing 'fit' and "
            "'predict' methods, %s (type %s) was passed" % (clf, type(clf))
            )
        if loss_func is None:
            assert hasattr(estimator, 'score'), ValueError(
                    "If no loss_func is specified, the estimator passed "
                    "should have a 'score' method. The estimator %s "
                    "does not." % estimator
                    )

        self.estimator = estimator
        self.param_grid = param_grid
        self.loss_func = loss_func
        self.n_jobs = n_jobs
        self.fit_params = fit_params
        self.iid = iid

    def fit(self, X, y, refit=True, cv=None, **kw):
        """Run fit with all sets of parameters
        Returns the best classifier

        Parameters
        ----------

        X: array, [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y: array, [n_samples]
            Target vector relative to X

        cv : crossvalidation generator
            see scikits.learn.cross_val module

        refit: boolean
            refit the best estimator with the entire dataset
        """
        estimator = self.estimator
        if cv is None:
            n_samples = len(X)
            if y is not None and is_classifier(estimator):
                cv = StratifiedKFold(y, k=3)
            else:
                cv = KFold(n_samples, k=3)

        grid = iter_grid(self.param_grid)
        base_clf = clone(self.estimator)
        out = Parallel(n_jobs=self.n_jobs)(
            delayed(fit_grid_point)(X, y, base_clf, clf_params,
                    cv, self.loss_func, self.iid, **self.fit_params)
                    for clf_params in grid)

        # Out is a list of pairs: score, estimator
        best_estimator = max(out)[1] # get maximum score

        if refit:
            # fit the best estimator using the entire dataset
            best_estimator.fit(X, y)

        self.best_estimator = best_estimator
        self.predict = best_estimator.predict
        if hasattr(best_estimator, 'score'):
            self.score = best_estimator.score

        # Store the computed scores
        grid = iter_grid(self.param_grid)
        self.grid_points_scores_ = dict((tuple(clf_params.items()), score)
                    for clf_params, (score, _) in zip(grid, out))

        return self


    def score(self, X, y=None):
        # This method is overridden during the fit if the best estimator
        # found has a score function.
        y_predicted = self.predict(X)
        return -self.loss_func(y_predicted, y)
