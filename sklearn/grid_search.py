"""
The :mod:`sklearn.grid_search` includes utilities to fine-tune the parameters
of an estimator.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD Style.

from itertools import product
import time

import numpy as np

from .base import BaseEstimator, is_classifier, clone
from .base import MetaEstimatorMixin
from .cross_validation import check_cv
from .externals.joblib import Parallel, delayed, logger
from .utils import check_arrays, safe_mask

__all__ = ['GridSearchCV', 'IterGrid', 'fit_grid_point']


class IterGrid(object):
    """Generators on the combination of the various parameter lists given

    Parameters
    ----------
    param_grid: dict of string to sequence
        The parameter grid to explore, as a dictionary mapping estimator
        parameters to sequences of allowed values.

    Returns
    -------
    params: dict of string to any
        **Yields** dictionaries mapping each estimator parameter to one of its
        allowed values.

    Examples
    --------
    >>> from sklearn.grid_search import IterGrid
    >>> param_grid = {'a':[1, 2], 'b':[True, False]}
    >>> list(IterGrid(param_grid)) #doctest: +NORMALIZE_WHITESPACE
    [{'a': 1, 'b': True}, {'a': 1, 'b': False},
     {'a': 2, 'b': True}, {'a': 2, 'b': False}]

    See also
    --------
    :class:`GridSearchCV`:
        uses ``IterGrid`` to perform a full parallelized grid search.
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

    X, y = check_arrays(X, y, sparse_format="csr")
    # update parameters of the classifier after a copy of its base structure
    clf = clone(base_clf)
    clf.set_params(**clf_params)

    if hasattr(base_clf, 'kernel') and hasattr(base_clf.kernel, '__call__'):
        # cannot compute the kernel values with custom function
        raise ValueError(
            "Cannot use a custom kernel function. "
            "Precompute the kernel matrix instead.")

    if getattr(base_clf, "_pairwise", False):
        # X is a precomputed square kernel matrix
        if X.shape[0] != X.shape[1]:
            raise ValueError("X should be a square kernel matrix")
        X_train = X[np.ix_(train, train)]
        X_test = X[np.ix_(test, train)]
    else:
        X_train = X[safe_mask(X, train)]
        X_test = X[safe_mask(X, test)]

    if y is not None:
        y_test = y[safe_mask(y, test)]
        y_train = y[safe_mask(y, train)]
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
    return this_score, clf_params, this_n_test_samples


def _check_param_grid(param_grid):
    if hasattr(param_grid, 'items'):
        param_grid = [param_grid]

    for p in param_grid:
        for v in p.itervalues():
            if isinstance(v, np.ndarray) and v.ndim > 1:
                raise ValueError("Parameter array should be one-dimensional.")

            check = [isinstance(v, k) for k in (list, tuple, np.ndarray)]
            if not True in check:
                raise ValueError("Parameter values should be a list.")

            if len(v) == 0:
                raise ValueError("Parameter values should be a non-empty "
                        "list.")


def _has_one_grid_point(param_grid):
    if hasattr(param_grid, 'items'):
        param_grid = [param_grid]

    for p in param_grid:
        for v in p.itervalues():
            if len(v) > 1:
                return False

    return True


class GridSearchCV(BaseEstimator, MetaEstimatorMixin):
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
        refit the best estimator with the entire dataset.
        If "False", it is impossible to make predictions using
        this GridSearch instance after fitting.

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
            shrinking=True, tol=...),
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

    `best_params_` : dict
        Parameter setting that gave the best results on the hold out data.

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
    :class:`IterGrid`:
        generates all the combinations of a an hyperparameter grid.

    :func:`sklearn.cross_validation.train_test_split`:
        utility function to split the data into a development set usable
        for fitting a GridSearchCV instance and an evaluation set for
        its final evaluation.

    """

    def __init__(self, estimator, param_grid, loss_func=None, score_func=None,
                 fit_params=None, n_jobs=1, iid=True, refit=True, cv=None,
                 verbose=0, pre_dispatch='2*n_jobs',
                ):
        if not hasattr(estimator, 'fit') or \
           not (hasattr(estimator, 'predict') or hasattr(estimator, 'score')):
            raise TypeError("estimator should a be an estimator implementing"
                            " 'fit' and 'predict' or 'score' methods,"
                            " %s (type %s) was passed" %
                            (estimator, type(estimator)))
        if loss_func is None and score_func is None:
            if not hasattr(estimator, 'score'):
                raise TypeError(
                    "If no loss_func is specified, the estimator passed "
                    "should have a 'score' method. The estimator %s "
                    "does not." % estimator)

        _check_param_grid(param_grid)

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

    def _set_methods(self):
        if hasattr(self._best_estimator_, 'predict'):
            self.predict = self._best_estimator_.predict
        if hasattr(self._best_estimator_, 'predict_proba'):
            self.predict_proba = self._best_estimator_.predict_proba

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
        return self._fit(X, y)

    def _fit(self, X, y):
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

        # Return early if there is only one grid point.
        if _has_one_grid_point(self.param_grid):
            params = next(iter(grid))
            base_clf.set_params(**params)
            base_clf.fit(X, y)
            self._best_estimator_ = base_clf
            self._set_methods()
            return self

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
            for this_score, clf_params, this_n_test_samples in \
                                    out[grid_start:grid_start + n_folds]:
                these_points.append(this_score)
                if self.iid:
                    this_score *= this_n_test_samples
                score += this_score
                n_test_samples += this_n_test_samples
            if self.iid:
                score /= float(n_test_samples)
            scores.append((score, clf_params))
            cv_scores.append(these_points)

        cv_scores = np.asarray(cv_scores)

        # Note: we do not use max(out) to make ties deterministic even if
        # comparison on estimator instances is not deterministic
        best_score = -np.inf
        for score, params in scores:
            if score > best_score:
                best_score = score
                best_params = params

        if best_score is None:
            raise ValueError('Best score could not be found')
        self.best_score_ = best_score
        self.best_params_ = best_params

        if self.refit:
            # fit the best estimator using the entire dataset
            # clone first to work around broken estimators
            best_estimator = clone(base_clf).set_params(**best_params)
            best_estimator.fit(X, y, **self.fit_params)
            self._best_estimator_ = best_estimator
            self._set_methods()

        # Store the computed scores
        # XXX: the name is too specific, it shouldn't have
        # 'grid' in it. Also, we should be retrieving/storing variance
        self.grid_scores_ = [
            (clf_params, score, all_scores)
                    for clf_params, (score, _), all_scores
                    in zip(grid, scores, cv_scores)]
        return self

    def score(self, X, y=None):
        if hasattr(self.best_estimator_, 'score'):
            return self.best_estimator_.score(X, y)
        if self.score_func is None:
            raise ValueError("No score function explicitly defined, "
                             "and the estimator doesn't provide one %s"
                             % self.best_estimator_)
        y_predicted = self.predict(X)
        return self.score_func(y, y_predicted)

    # TODO around 0.13: remove this property, make it an attribute
    @property
    def best_estimator_(self):
        if hasattr(self, '_best_estimator_'):
            return self._best_estimator_
        else:
            raise RuntimeError("Grid search has to be run with 'refit=True'"
                " to make predictions or obtain an instance  of the best "
                " estimator. To obtain the best parameter settings, "
                " use ``best_params_``.")
