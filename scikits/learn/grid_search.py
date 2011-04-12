"""Tune the parameters of an estimator by cross-validation"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux    <gael.varoquaux@normalesup.org>
# License: BSD Style.

import copy
import time

import numpy as np
import scipy.sparse as sp

from .externals.joblib import Parallel, delayed
from .externals.joblib.logger import short_format_time
from .cross_val import KFold, StratifiedKFold
from .base import BaseEstimator, is_classifier, clone


try:
    from itertools import product
except:
    def product(*args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x + [y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)


class IterGrid(object):
    """Generators on the combination of the various parameter lists given

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
    >>> from scikits.learn.grid_search import IterGrid
    >>> param_grid = {'a':[1, 2], 'b':[True, False]}
    >>> list(IterGrid(param_grid)) #doctest: +NORMALIZE_WHITESPACE
    [{'a': 1, 'b': True}, {'a': 1, 'b': False},
     {'a': 2, 'b': True}, {'a': 2, 'b': False}]

    """
    def __init__(self, param_grid):
        self.param_grid = param_grid

    def __iter__(self):
        param_grid = self.param_grid
        if hasattr(param_grid, 'has_key'):
            param_grid = [param_grid]
        for p in param_grid:
            # Always sort the keys of a dictionary, for reproducibility
            items = sorted(p.items())
            keys, values = zip(*items)
            for v in product(*values):
                params = dict(zip(keys, v))
                yield params


def fit_grid_point(X, y, base_clf, clf_params, train, test, loss_func,
                score_func, verbose, **fit_params):
    """Run fit on one set of parameters

    Returns the score and the instance of the classifier
    """
    if verbose > 1:
        start_time = time.time()
        msg = '%s' % (', '.join('%s=%s' % (k, v)
                                     for k, v in clf_params.iteritems()))
        print "[GridSearchCV] %s %s" % (msg, (64 - len(msg)) * '.')
    # update parameters of the classifier after a copy of its base structure
    clf = copy.deepcopy(base_clf)
    clf._set_params(**clf_params)

    if isinstance(X, list) or isinstance(X, tuple):
        X_train = [X[i] for i, cond in enumerate(train) if cond]
        X_test = [X[i] for i, cond in enumerate(test) if cond]
    else:
        if sp.issparse(X):
            # For sparse matrices, slicing only works with indices
            # (no masked array). Convert to CSR format for efficiency and
            # because some sparse formats don't support row slicing.
            X = sp.csr_matrix(X)
            ind = np.arange(X.shape[0])
            train = ind[train]
            test = ind[test]
        X_train = X[train]
        X_test = X[test]
    if y is not None:
        y_test = y[test]
        y_train = y[train]
    else:
        y_test = None
        y_train = None

    clf.fit(X_train, y_train, **fit_params)

    if loss_func is not None:
        y_pred = clf.predict(X_test)
        this_score = -loss_func(y_test, y_pred)
    elif score_func is not None:
        y_pred = clf.predict(X_test)
        this_score = score_func(y_test, y_pred)
    else:
        this_score = clf.score(X_test, y_test)

    if y is not None:
        this_n_test_samples = y.shape[0]
    else:
        this_n_test_samples = X.shape[0]

    if verbose > 1:
        end_msg = "%s -%s" % (msg,
                                short_format_time(time.time() - start_time))
        print "[GridSearchCV] %s %s" % ((64 - len(end_msg)) * '.', end_msg)
    return this_score, clf, this_n_test_samples


class GridSearchCV(BaseEstimator):
    """Grid search on the parameters of a classifier

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

    score_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediciton (big is good)
        if None is passed, the score of the estimator is maximized

    fit_params : dict, optional
        parameters to pass to the fit method

    n_jobs: int, optional
        number of jobs to run in parallel (default 1)

    iid: boolean, optional
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    cv : crossvalidation generator
        see scikits.learn.cross_val module

    refit: boolean
        refit the best estimator with the entire dataset

    verbose: integer
        Controls the verbosity: the higher, the more messages.

    Examples
    --------
    >>> from scikits.learn import svm, grid_search, datasets
    >>> iris = datasets.load_iris()
    >>> parameters = {'kernel':('linear', 'rbf'), 'C':[1, 10]}
    >>> svr = svm.SVR()
    >>> clf = grid_search.GridSearchCV(svr, parameters)
    >>> clf.fit(iris.data, iris.target) # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    GridSearchCV(n_jobs=1, verbose=0, fit_params={}, loss_func=None,
                 refit=True, cv=None, iid=True,
                 estimator=SVR(kernel='rbf', C=1.0, probability=False, ...
                 ...

    Notes
    ------

    The parameters selected are those that maximize the score of the
    left out data, unless an explicit score_func is passed in which
    case it is used instead. If a loss function loss_func is passed,
    it overrides the score functions and is minimized.

    """

    def __init__(self, estimator, param_grid, loss_func=None, score_func=None,
                 fit_params={}, n_jobs=1, iid=True, refit=True, cv=None,
                 verbose=0,
                 ):
        assert hasattr(estimator, 'fit') and (hasattr(estimator, 'predict')
                        or hasattr(estimator, 'score')), (
            "estimator should a be an estimator implementing 'fit' and "
            "'predict' or 'score' methods, %s (type %s) was passed" %
                    (estimator, type(estimator)))
        if loss_func is None and score_func is None:
            assert hasattr(estimator, 'score'), ValueError(
                    "If no loss_func is specified, the estimator passed "
                    "should have a 'score' method. The estimator %s "
                    "does not." % estimator)

        self.estimator = estimator
        self.param_grid = param_grid
        self.loss_func = loss_func
        self.score_func = score_func
        self.n_jobs = n_jobs
        self.fit_params = fit_params
        self.iid = iid
        self.refit = refit
        self.cv = cv
        self.verbose = verbose

    def fit(self, X, y=None, **params):
        """Run fit with all sets of parameters

        Returns the best classifier

        Parameters
        ----------

        X: array, [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y: array, [n_samples] or None
            Target vector relative to X, None for unsupervised problems

        """
        self._set_params(**params)
        estimator = self.estimator
        cv = self.cv
        if cv is None:
            if hasattr(X, 'shape'):
                n_samples = X.shape[0]
            else:
                # support list of unstructured objects on which feature
                # extraction will be applied later in the tranformer chain
                n_samples = len(X)
            if y is not None and is_classifier(estimator):
                cv = StratifiedKFold(y, k=3)
            else:
                cv = KFold(n_samples, k=3)

        grid = IterGrid(self.param_grid)
        base_clf = clone(self.estimator)
        # XXX: Need to make use of Parallel's new pre_dispatch
        out = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
            delayed(fit_grid_point)(
                X, y, base_clf, clf_params, train, test, self.loss_func,
                self.score_func, self.verbose, **self.fit_params)
                    for clf_params in grid for train, test in cv)

        # Out is a list of triplet: score, estimator, n_test_samples
        n_grid_points = len(list(grid))
        n_fits = len(out)
        n_folds = n_fits // n_grid_points

        scores = list()
        for grid_start in range(0, n_fits, n_folds):
            n_test_samples = 0
            score = 0
            for this_score, estimator, this_n_test_samples in \
                                    out[grid_start:grid_start + n_folds]:
                if self.iid:
                    this_score *= this_n_test_samples
                score += this_score
                n_test_samples += this_n_test_samples
            if self.iid:
                score /= float(n_test_samples)
            scores.append((score, estimator))

        # Note: we do not use max(out) to make ties deterministic even if
        # comparison on estimator instances is not deterministic
        best_score = None
        for score, estimator in scores:
            if best_score is None:
                best_score = score
                best_estimator = estimator
            else:
                if score >= best_score:
                    best_score = score
                    best_estimator = estimator

        self.best_score = best_score

        if self.refit:
            # fit the best estimator using the entire dataset
            best_estimator.fit(X, y, **self.fit_params)

        self.best_estimator = best_estimator
        if hasattr(best_estimator, 'predict'):
            self.predict = best_estimator.predict
        if hasattr(best_estimator, 'score'):
            self.score = best_estimator.score

        # Store the computed scores
        # XXX: the name is too specific, it shouldn't have
        # 'grid' in it. Also, we should be retrieving/storing variance
        self.grid_points_scores_ = dict((tuple(clf_params.items()), score)
                    for clf_params, (score, _) in zip(grid, scores))

        return self

    def score(self, X, y=None):
        # This method is overridden during the fit if the best estimator
        # found has a score function.
        y_predicted = self.predict(X)
        return self.score_func(y, y_predicted)
