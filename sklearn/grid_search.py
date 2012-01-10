"""
The :mod:`sklearn.grid_search` includes utilities to fine-tune the parameters
of an estimator.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD Style.

import copy
from itertools import product
import time

import numpy as np
import scipy.sparse as sp

from .base import BaseEstimator, is_classifier, clone
from .cross_validation import check_cv
from .externals.joblib import Parallel, delayed, logger
from .utils import deprecated


class IterGrid(object):
    """Generators on the combination of the various parameter lists given

    Parameters
    -----------
    param_grid: dict of string to sequence
        The parameter grid to explore, as a dictionary mapping estimator
        parameters to sequences of allowed values.

    Returns
    -------
    params: dict of string to any
        **Yields** dictionaries mapping each estimator parameter to one of its
        allowed values.

    Examples
    ---------
    >>> from sklearn.grid_search import IterGrid
    >>> param_grid = {'a':[1, 2], 'b':[True, False]}
    >>> list(IterGrid(param_grid)) #doctest: +NORMALIZE_WHITESPACE
    [{'a': 1, 'b': True}, {'a': 1, 'b': False},
     {'a': 2, 'b': True}, {'a': 2, 'b': False}]
    """

    def __init__(self, param_grid):
        self.param_grid = param_grid

    def __iter__(self):
        param_grid = self.param_grid
        if hasattr(param_grid, 'items'):
            # wrap dictionary in a singleton list
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
    # FIXME we should be doing a clone here
    clf = copy.deepcopy(base_clf)
    clf.set_params(**clf_params)

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
        if hasattr(y, 'shape'):
            this_n_test_samples = y.shape[0]
        else:
            this_n_test_samples = len(y)
    else:
        if hasattr(X, 'shape'):
            this_n_test_samples = X.shape[0]
        else:
            this_n_test_samples = len(X)
    if verbose > 2:
        msg += ", score=%f" % this_score
    if verbose > 1:
        end_msg = "%s -%s" % (msg,
                              logger.short_format_time(time.time() -
                                                       start_time))
        print "[GridSearchCV] %s %s" % ((64 - len(end_msg)) * '.', end_msg)
    return this_score, clf, this_n_test_samples


class GridSearchCV(BaseEstimator):
    """Grid search on the parameters of a classifier

    Important members are fit, predict.

    GridSearchCV implements a "fit" method and a "predict" method like
    any classifier except that the parameters of the classifier
    used to predict is optimized by cross-validation.

    Parameters
    ----------
    estimator: object type that implements the "fit" and "predict" methods
        A object of that type is instantiated for each grid point.

    param_grid: dict
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values.

    loss_func: callable, optional
        function that takes 2 arguments and compares them in
        order to evaluate the performance of prediciton (small is good)
        if None is passed, the score of the estimator is maximized

    score_func: callable, optional
        A function that takes 2 arguments and compares them in
        order to evaluate the performance of prediction (high is good).
        If None is passed, the score of the estimator is maximized.

    fit_params : dict, optional
        parameters to pass to the fit method

    n_jobs: int, optional
        number of jobs to run in parallel (default 1)

    pre_dispatch: int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediatly
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    iid: boolean, optional
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    cv : integer or crossvalidation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific crossvalidation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    refit: boolean
        refit the best estimator with the entire dataset

    verbose: integer
        Controls the verbosity: the higher, the more messages.

    Examples
    --------
    >>> from sklearn import svm, grid_search, datasets
    >>> iris = datasets.load_iris()
    >>> parameters = {'kernel':('linear', 'rbf'), 'C':[1, 10]}
    >>> svr = svm.SVC()
    >>> clf = grid_search.GridSearchCV(svr, parameters)
    >>> clf.fit(iris.data, iris.target)
    ...                             # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    GridSearchCV(cv=None,
        estimator=SVC(C=1.0, cache_size=..., coef0=..., degree=...,
            gamma=..., kernel='rbf', probability=False,
            scale_C=False, shrinking=True, tol=...),
        fit_params={}, iid=True, loss_func=None, n_jobs=1,
            param_grid=...,
            ...)

    Attributes
    ----------
    `grid_scores_` : dict of any to float
        Contains scores for all parameter combinations in param_grid.

    `best_estimator_` : estimator
        Estimator that was choosen by grid search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data.

    `best_score_` : float
        score of best_estimator on the left out data.


    Notes
    ------
    The parameters selected are those that maximize the score of the left out
    data, unless an explicit score_func is passed in which case it is used
    instead. If a loss function loss_func is passed, it overrides the score
    functions and is minimized.

    If `n_jobs` was set to a value higher than one, the data is copied for each
    point in the grid (and not `n_jobs` times). This is done for efficiency
    reasons if individual jobs take very little time, but may raise errors if
    the dataset is large and not enough memory is available.  A workaround in
    this case is to set `pre_dispatch`. Then, the memory is copied only
    `pre_dispatch` many times. A reasonable value for `pre_dispatch` is 2 *
    `n_jobs`.

    See Also
    ---------
    IterGrid

    """

    def __init__(self, estimator, param_grid, loss_func=None, score_func=None,
                 fit_params=None, n_jobs=1, iid=True, refit=True, cv=None,
                 verbose=0, pre_dispatch='2*n_jobs',
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
        self.fit_params = fit_params if fit_params is not None else {}
        self.iid = iid
        self.refit = refit
        self.cv = cv
        self.verbose = verbose
        self.pre_dispatch = pre_dispatch

    def fit(self, X, y=None, **params):
        """Run fit with all sets of parameters

        Returns the best classifier

        Parameters
        ----------

        X: array, [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y: array-like, shape = [n_samples], optional
            Target vector relative to X for classification;
            None for unsupervised learning.

        """
        self._set_params(**params)
        estimator = self.estimator
        cv = self.cv
        if hasattr(X, 'shape'):
            n_samples = X.shape[0]
        else:
            # support list of unstructured objects on which feature
            # extraction will be applied later in the tranformer chain
            n_samples = len(X)
        if y is not None:
            if len(y) != n_samples:
                raise ValueError('Target variable (y) has a different number '
                                 'of samples (%i) than data (X: %i samples)'
                                 % (len(y), n_samples))
            y = np.asarray(y)
        cv = check_cv(cv, X, y, classifier=is_classifier(estimator))

        grid = IterGrid(self.param_grid)
        base_clf = clone(self.estimator)
        pre_dispatch = self.pre_dispatch
        out = Parallel(n_jobs=self.n_jobs, verbose=self.verbose,
                pre_dispatch=pre_dispatch)(
            delayed(fit_grid_point)(
                X, y, base_clf, clf_params, train, test, self.loss_func,
                self.score_func, self.verbose, **self.fit_params)
                    for clf_params in grid for train, test in cv)

        # Out is a list of triplet: score, estimator, n_test_samples
        n_grid_points = len(list(grid))
        n_fits = len(out)
        n_folds = n_fits // n_grid_points

        scores = list()
        cv_scores = list()
        for grid_start in range(0, n_fits, n_folds):
            n_test_samples = 0
            score = 0
            these_points = list()
            for this_score, estimator, this_n_test_samples in \
                                    out[grid_start:grid_start + n_folds]:
                these_points.append(this_score)
                if self.iid:
                    this_score *= this_n_test_samples
                score += this_score
                n_test_samples += this_n_test_samples
            if self.iid:
                score /= float(n_test_samples)
            scores.append((score, estimator))
            cv_scores.append(these_points)

        # Note: we do not use max(out) to make ties deterministic even if
        # comparison on estimator instances is not deterministic
        best_score = None
        for score, estimator in scores:
            if best_score is None:
                best_score = score
                best_estimator = estimator
            else:
                if score > best_score:
                    best_score = score
                    best_estimator = estimator

        if best_score is None:
            raise ValueError('Best score could not be found')
        self.best_score_ = best_score

        if self.refit:
            # fit the best estimator using the entire dataset
            # clone first to work around broken estimators
            best_estimator = clone(best_estimator)
            best_estimator.fit(X, y, **self.fit_params)

        self.best_estimator_ = best_estimator
        if hasattr(best_estimator, 'predict'):
            self.predict = best_estimator.predict
        if hasattr(best_estimator, 'predict_proba'):
            self.predict_proba = best_estimator.predict_proba
        if hasattr(best_estimator, 'score'):
            self.score_ = best_estimator.score

        # Store the computed scores
        # XXX: the name is too specific, it shouldn't have
        # 'grid' in it. Also, we should be retrieving/storing variance
        self.grid_scores_ = [
            (clf_params, score, all_scores)
                    for clf_params, (score, _), all_scores
                    in zip(grid, scores, cv_scores)]
        return self

    def score(self, X, y=None):
        # This method is overridden during the fit if the best estimator
        # found has a score function.
        y_predicted = self.predict(X)
        return self.score_func(y, y_predicted)

    @property
    @deprecated('GridSearchCV.best_estimator is deprecated'
                ' and will be removed in version 0.12.'
                ' Please use ``GridSearchCV.best_estimator_`` instead.')
    def best_estimator(self):
        return self.best_estimator_

    @property
    @deprecated('GridSearchCV.best_score is deprecated'
                ' and will be removed in version 0.12.'
                ' Please use ``GridSearchCV.best_score_`` instead.')
    def best_score(self):
        return self.best_score_
