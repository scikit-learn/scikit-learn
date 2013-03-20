"""
The :mod:`sklearn.grid_search` includes utilities to fine-tune the parameters
of an estimator.
"""
from __future__ import print_function

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD Style.

from abc import ABCMeta, abstractmethod
from collections import Mapping, namedtuple
from functools import partial, reduce
from itertools import product
import numbers
import operator
import time
import warnings

import numpy as np

from .base import BaseEstimator, is_classifier, clone
from .base import MetaEstimatorMixin
from .cross_validation import check_cv
from .externals.joblib import Parallel, delayed, logger
from .externals.six import string_types
from .utils import safe_mask, check_random_state
from .utils.validation import _num_samples, check_arrays
from .metrics import SCORERS, Scorer

__all__ = ['GridSearchCV', 'ParameterGrid', 'fit_grid_point',
           'ParameterSampler', 'RandomizedSearchCV']


class ParameterGrid(object):
    """Grid of parameters with a discrete number of values for each.

    Can be used to iterate over parameter value combinations with the
    Python built-in function iter.

    Parameters
    ----------
    param_grid : dict of string to sequence
        The parameter grid to explore, as a dictionary mapping estimator
        parameters to sequences of allowed values.

    Examples
    --------
    >>> from sklearn.grid_search import ParameterGrid
    >>> param_grid = {'a':[1, 2], 'b':[True, False]}
    >>> list(ParameterGrid(param_grid)) #doctest: +NORMALIZE_WHITESPACE
    [{'a': 1, 'b': True}, {'a': 1, 'b': False},
     {'a': 2, 'b': True}, {'a': 2, 'b': False}]

    See also
    --------
    :class:`GridSearchCV`:
        uses ``ParameterGrid`` to perform a full parallelized parameter search.
    """

    def __init__(self, param_grid):
        if isinstance(param_grid, Mapping):
            # wrap dictionary in a singleton list
            # XXX Why? The behavior when passing a list is undocumented,
            # but not doing this breaks one of the tests.
            param_grid = [param_grid]
        self.param_grid = param_grid

    def __iter__(self):
        """Iterate over the points in the grid.

        Returns
        -------
        params : iterator over dict of string to any
            Yields dictionaries mapping each estimator parameter to one of its
            allowed values.
        """
        for p in self.param_grid:
            # Always sort the keys of a dictionary, for reproducibility
            items = sorted(p.items())
            keys, values = zip(*items)
            for v in product(*values):
                params = dict(zip(keys, v))
                yield params

    def __len__(self):
        """Number of points on the grid."""
        # Product function that can handle iterables (np.product can't).
        product = partial(reduce, operator.mul)
        return sum(product(len(v) for v in p.values())
                   for p in self.param_grid)


class IterGrid(ParameterGrid):
    """Generators on the combination of the various parameter lists given.

    This class is DEPRECATED. It was renamed to ``ParameterGrid``. The name
    ``IterGrid`` will be removed in 0.15.

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
        uses ``IterGrid`` to perform a full parallelized parameter search.
    """

    def __init__(self, param_grid):
        warnings.warn("IterGrid was renamed to ParameterGrid and will be"
                      " removed in 0.15.", DeprecationWarning)
        super(IterGrid, self).__init__(param_grid)


class ParameterSampler(object):
    """Generator on parameters sampled from given distributions.

    Parameters
    ----------
    param_distributions : dict
        Dictionary where the keys are parameters and values
        are distributions from which a parameter is to be sampled.
        Distributions either have to provide a ``rvs`` function
        to sample from them, or can be given as a list of values,
        where a uniform distribution is assumed.

    n_iter : integer
        Number of parameter settings that are produced.

    random_state : int or RandomState
            Pseudo number generator state used for random sampling.

    Returns
    -------
    params: dict of string to any
        **Yields** dictionaries mapping each estimator parameter to
        as sampled value.

    Examples
    --------
    >>> from sklearn.grid_search import ParameterSampler
    >>> from scipy.stats.distributions import expon
    >>> import numpy as np
    >>> np.random.seed(0)
    >>> param_grid = {'a':[1, 2], 'b': expon()}
    >>> list(ParameterSampler(param_grid, n_iter=4))
    ...  #doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    [{'a': 1, 'b': 0.89...}, {'a': 1, 'b': 0.92...},
     {'a': 2, 'b': 1.87...}, {'a': 2, 'b': 1.03...}]

    """
    def __init__(self, param_distributions, n_iter, random_state=None):
        self.param_distributions = param_distributions
        self.n_iter = n_iter
        self.random_state = random_state

    def __iter__(self):
        rnd = check_random_state(self.random_state)
        # Always sort the keys of a dictionary, for reproducibility
        items = sorted(self.param_distributions.items())
        for i in range(self.n_iter):
            params = dict()
            for k, v in items:
                if hasattr(v, "rvs"):
                    params[k] = v.rvs()
                else:
                    params[k] = v[rnd.randint(len(v))]
            yield params

    def __len__(self):
        """Number of points that will be sampled."""
        return self.n_iter


def fit_grid_point(X, y, base_clf, clf_params, train, test, scorer,
                   verbose, loss_func=None, **fit_params):
    """Run fit on one set of parameters.

    Parameters
    ----------
    X : array-like, sparse matrix or list
        Input data.

    y : array-like or None
        Targets for input data.

    base_clf : estimator object
        This estimator will be cloned and then fitted.

    clf_params : dict
        Parameters to be set on base_estimator clone for this grid point.

    train : ndarray, dtype int or bool
        Boolean mask or indices for training set.

    test : ndarray, dtype int or bool
        Boolean mask or indices for test set.

    scorer : callable or None.
        If provided must be a scoring object / function with signature
        ``scorer(estimator, X, y)``.

    verbose : int
        Verbosity level.

    **fit_params : kwargs
        Additional parameter passed to the fit function of the estimator.


    Returns
    -------
    score : float
        Score of this parameter setting on given training / test split.

    estimator : estimator object
        Estimator object of type base_clf that was fitted using clf_params
        and provided train / test split.

    n_samples_test : int
        Number of test samples in this split.
    """
    if verbose > 1:
        start_time = time.time()
        msg = '%s' % (', '.join('%s=%s' % (k, v)
                      for k, v in clf_params.items()))
        print("[GridSearchCV] %s %s" % (msg, (64 - len(msg)) * '.'))

    # update parameters of the classifier after a copy of its base structure
    clf = clone(base_clf)
    clf.set_params(**clf_params)

    if hasattr(base_clf, 'kernel') and callable(base_clf.kernel):
        # cannot compute the kernel values with custom function
        raise ValueError("Cannot use a custom kernel function. "
                         "Precompute the kernel matrix instead.")

    if not hasattr(X, "shape"):
        if getattr(base_clf, "_pairwise", False):
            raise ValueError("Precomputed kernels or affinity matrices have "
                             "to be passed as arrays or sparse matrices.")
        X_train = [X[idx] for idx in train]
        X_test = [X[idx] for idx in test]
    else:
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
        clf.fit(X_train, y_train, **fit_params)

        if scorer is not None:
            this_score = scorer(clf, X_test, y_test)
        else:
            this_score = clf.score(X_test, y_test)
    else:
        clf.fit(X_train, **fit_params)
        if scorer is not None:
            this_score = scorer(clf, X_test)
        else:
            this_score = clf.score(X_test)

    if not isinstance(this_score, numbers.Number):
        raise ValueError("scoring must return a number, got %s (%s)"
                         " instead." % (str(this_score), type(this_score)))

    if verbose > 2:
        msg += ", score=%f" % this_score
    if verbose > 1:
        end_msg = "%s -%s" % (msg,
                              logger.short_format_time(time.time() -
                                                       start_time))
        print("[GridSearchCV] %s %s" % ((64 - len(end_msg)) * '.', end_msg))
    return this_score, clf_params, _num_samples(X_test)


def _check_param_grid(param_grid):
    if hasattr(param_grid, 'items'):
        param_grid = [param_grid]

    for p in param_grid:
        for v in p.values():
            if isinstance(v, np.ndarray) and v.ndim > 1:
                raise ValueError("Parameter array should be one-dimensional.")

            check = [isinstance(v, k) for k in (list, tuple, np.ndarray)]
            if not True in check:
                raise ValueError("Parameter values should be a list.")

            if len(v) == 0:
                raise ValueError("Parameter values should be a non-empty "
                                 "list.")


class BaseSearchCV(BaseEstimator, MetaEstimatorMixin):
    """Base class for hyper parameter search with cross-validation.
    """
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, estimator, scoring=None, loss_func=None,
                 score_func=None, fit_params=None, n_jobs=1, iid=True,
                 refit=True, cv=None, verbose=0, pre_dispatch='2*n_jobs'):

        self.scoring = scoring
        self.estimator = estimator
        self.loss_func = loss_func
        self.score_func = score_func
        self.n_jobs = n_jobs
        self.fit_params = fit_params if fit_params is not None else {}
        self.iid = iid
        self.refit = refit
        self.cv = cv
        self.verbose = verbose
        self.pre_dispatch = pre_dispatch
        self._check_estimator()

    def score(self, X, y=None):
        """Returns the mean accuracy on the given test data and labels.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training set.

        y : array-like, shape = [n_samples], optional
            Labels for X.

        Returns
        -------
        score : float

        """
        if hasattr(self.best_estimator_, 'score'):
            return self.best_estimator_.score(X, y)
        if self.scorer_ is None:
            raise ValueError("No score function explicitly defined, "
                             "and the estimator doesn't provide one %s"
                             % self.best_estimator_)
        y_predicted = self.predict(X)
        return self.scorer(y, y_predicted)

    def _check_estimator(self):
        """Check that estimator can be fitted and score can be computed."""
        if (not hasattr(self.estimator, 'fit') or
                not (hasattr(self.estimator, 'predict')
                     or hasattr(self.estimator, 'score'))):
            raise TypeError("estimator should a be an estimator implementing"
                            " 'fit' and 'predict' or 'score' methods,"
                            " %s (type %s) was passed" %
                            (self.estimator, type(self.estimator)))
        if (self.scoring is None and self.loss_func is None and self.score_func
                is None):
            if not hasattr(self.estimator, 'score'):
                raise TypeError(
                    "If no scoring is specified, the estimator passed "
                    "should have a 'score' method. The estimator %s "
                    "does not." % self.estimator)

    def _set_methods(self):
        """Create predict and  predict_proba if present in best estimator."""
        if hasattr(self.best_estimator_, 'predict'):
            self.predict = self.best_estimator_.predict
        if hasattr(self.best_estimator_, 'predict_proba'):
            self.predict_proba = self.best_estimator_.predict_proba

    def _fit(self, X, y, parameter_iterator, **params):
        """Actual fitting,  performing the search over parameters."""
        estimator = self.estimator
        cv = self.cv

        n_samples = _num_samples(X)
        X, y = check_arrays(X, y, allow_lists=True, sparse_format='csr')

        if self.loss_func is not None:
            warnings.warn("Passing a loss function is "
                          "deprecated and will be removed in 0.15. "
                          "Either use strings or score objects."
                          "The relevant new parameter is called ''scoring''. ")
            scorer = Scorer(self.loss_func, greater_is_better=False)
        elif self.score_func is not None:
            warnings.warn("Passing function as ``score_func`` is "
                          "deprecated and will be removed in 0.15. "
                          "Either use strings or score objects."
                          "The relevant new parameter is called ''scoring''.")
            scorer = Scorer(self.score_func)
        elif isinstance(self.scoring, string_types):
            scorer = SCORERS[self.scoring]
        else:
            scorer = self.scoring

        self.scorer_ = scorer

        if y is not None:
            if len(y) != n_samples:
                raise ValueError('Target variable (y) has a different number '
                                 'of samples (%i) than data (X: %i samples)'
                                 % (len(y), n_samples))
            y = np.asarray(y)
        cv = check_cv(cv, X, y, classifier=is_classifier(estimator))

        base_clf = clone(self.estimator)

        pre_dispatch = self.pre_dispatch

        out = Parallel(
            n_jobs=self.n_jobs, verbose=self.verbose,
            pre_dispatch=pre_dispatch)(
                delayed(fit_grid_point)(
                    X, y, base_clf, clf_params, train, test, scorer,
                    self.verbose, **self.fit_params) for clf_params in
                parameter_iterator for train, test in cv)

        # Out is a list of triplet: score, estimator, n_test_samples
        n_param_points = len(list(parameter_iterator))
        n_fits = len(out)
        n_folds = n_fits // n_param_points

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
                    n_test_samples += this_n_test_samples
                score += this_score
            if self.iid:
                score /= float(n_test_samples)
            else:
                score /= float(n_folds)
            scores.append((score, clf_params))
            cv_scores.append(these_points)

        cv_scores = np.asarray(cv_scores)

        # Note: we do not use max(out) to make ties deterministic even if
        # comparison on estimator instances is not deterministic
        if scorer is not None:
            greater_is_better = scorer.greater_is_better
        else:
            greater_is_better = True

        if greater_is_better:
            best_score = -np.inf
        else:
            best_score = np.inf

        for score, params in scores:
            if ((score > best_score and greater_is_better)
                    or (score < best_score and not greater_is_better)):
                best_score = score
                best_params = params

        self.best_params_ = best_params
        self.best_score_ = best_score

        if self.refit:
            # fit the best estimator using the entire dataset
            # clone first to work around broken estimators
            best_estimator = clone(base_clf).set_params(**best_params)
            if y is not None:
                best_estimator.fit(X, y, **self.fit_params)
            else:
                best_estimator.fit(X, **self.fit_params)
            self.best_estimator_ = best_estimator
            self._set_methods()

        # Store the computed scores
        CVScoreTuple = namedtuple('CVScoreTuple', ('parameters',
                                                   'mean_validation_score',
                                                   'cv_validation_scores'))
        self.cv_scores_ = [
            CVScoreTuple(clf_params, score, all_scores)
            for clf_params, (score, _), all_scores
            in zip(parameter_iterator, scores, cv_scores)]
        return self


class GridSearchCV(BaseSearchCV):
    """Exhaustive search over specified parameter values for an estimator.

    Important members are fit, predict.

    GridSearchCV implements a "fit" method and a "predict" method like
    any classifier except that the parameters of the classifier
    used to predict is optimized by cross-validation.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        A object of that type is instantiated for each grid point.

    param_grid : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    scoring : string or callable, optional
        Either one of either a string ("zero_one", "f1", "roc_auc", ... for
        classification, "mse", "r2",... for regression) or a callable.
        See 'Scoring objects' in the model evaluation section of the user guide
        for details.

    fit_params : dict, optional
        Parameters to pass to the fit method.

    n_jobs : int, optional
        Number of jobs to run in parallel (default 1).

    pre_dispatch : int, or string, optional
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

    iid : boolean, optional
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    cv : integer or crossvalidation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific crossvalidation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    refit : boolean
        Refit the best estimator with the entire dataset.
        If "False", it is impossible to make predictions using
        this GridSearchCV instance after fitting.

    verbose : integer
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
            gamma=..., kernel='rbf', max_iter=-1, probability=False,
            shrinking=True, tol=...),
        fit_params={}, iid=True, loss_func=None, n_jobs=1,
            param_grid=...,
            ...)

    Attributes
    ----------
    `cv_scores_` : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

            * ``parameters``, a dict of parameter settings
            * ``mean_validation_score``, the mean score over the
              cross-validation folds
            * ``cv_validation_scores``, the list of scores for each fold

    `best_estimator_` : estimator
        Estimator that was choosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data.

    `best_score_` : float
        Score of best_estimator on the left out data.

    `best_params_` : dict
        Parameter setting that gave the best results on the hold out data.

    Notes
    ------
    The parameters selected are those that maximize the score of the left out
    data, unless an explicit score is passed in which case it is used instead.

    If `n_jobs` was set to a value higher than one, the data is copied for each
    point in the grid (and not `n_jobs` times). This is done for efficiency
    reasons if individual jobs take very little time, but may raise errors if
    the dataset is large and not enough memory is available.  A workaround in
    this case is to set `pre_dispatch`. Then, the memory is copied only
    `pre_dispatch` many times. A reasonable value for `pre_dispatch` is `2 *
    n_jobs`.

    See Also
    ---------
    :class:`ParameterGrid`:
        generates all the combinations of a an hyperparameter grid.

    :func:`sklearn.cross_validation.train_test_split`:
        utility function to split the data into a development set usable
        for fitting a GridSearchCV instance and an evaluation set for
        its final evaluation.

    """

    def __init__(self, estimator, param_grid, scoring=None, loss_func=None,
                 score_func=None, fit_params=None, n_jobs=1, iid=True,
                 refit=True, cv=None, verbose=0, pre_dispatch='2*n_jobs'):
        super(GridSearchCV, self).__init__(
            estimator, scoring, loss_func, score_func, fit_params, n_jobs, iid,
            refit, cv, verbose, pre_dispatch)
        self.param_grid = param_grid
        _check_param_grid(param_grid)

    @property
    def grid_scores_(self):
        warnings.warn("grid_scores_ is deprecated and will be removed in 0.15."
                      " Use cv_scores_ instead.", DeprecationWarning)
        return self.cv_scores_

    def fit(self, X, y=None, **params):
        """Run fit with all sets of parameters.

        Parameters
        ----------

        X: array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y: array-like, shape = [n_samples], optional
            Target vector relative to X for classification;
            None for unsupervised learning.

        """
        return self._fit(X, y, ParameterGrid(self.param_grid), **params)


class RandomizedSearchCV(BaseSearchCV):
    """Randomized search on hyper parameters.

    RandomizedSearchCV implements a "fit" method and a "predict" method like
    any classifier except that the parameters of the classifier
    used to predict is optimized by cross-validation.

    In constrast to GridSearchCV, not all parameter values are tried out, but
    rather a fixed number of parameter settings is sampled from the specified
    distributions. The number of parameter settings that are tried is
    given by n_iter.

    Parameters
    ----------
    estimator : object type that implements the "fit" and "predict" methods
        A object of that type is instantiated for each parameter setting.

    param_distribution : dict
        Dictionary with parameters names (string) as keys and distributions
        or lists of parameters to try. Distributions must provide a ``rvs``
        method for sampling (such as those from scipy.stats.distributions).
        If a list is given, it is sampled uniformly.

    n_iter : int, default=10
        Number of parameter settings that are sampled. n_iter trades
        off runtime vs qualitiy of the solution.

    scoring : string or callable, optional
        Either one of either a string ("zero_one", "f1", "roc_auc", ... for
        classification, "mse", "r2",... for regression) or a callable.
        See 'Scoring objects' in the model evaluation section of the user guide
        for details.

    fit_params : dict, optional
        Parameters to pass to the fit method.

    n_jobs : int, optional
        Number of jobs to run in parallel (default 1).

    pre_dispatch : int, or string, optional
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

    iid : boolean, optional
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    cv : integer or crossvalidation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific crossvalidation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    refit : boolean
        Refit the best estimator with the entire dataset.
        If "False", it is impossible to make predictions using
        this RandomizedSearchCV instance after fitting.

    verbose : integer
        Controls the verbosity: the higher, the more messages.


    Attributes
    ----------
    `cv_scores_` : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

            * ``parameters``, a dict of parameter settings
            * ``mean_validation_score``, the mean score over the
              cross-validation folds
            * ``cv_validation_scores``, the list of scores for each fold

    `best_estimator_` : estimator
        Estimator that was choosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data.

    `best_score_` : float
        Score of best_estimator on the left out data.

    `best_params_` : dict
        Parameter setting that gave the best results on the hold out data.

    Notes
    -----
    The parameters selected are those that maximize the score of the left out
    data, unless an explicit score_func is passed in which case it is used
    instead. If a loss function loss_func is passed, it overrides the score
    functions and is minimized.

    If `n_jobs` was set to a value higher than one, the data is copied for each
    parameter setting(and not `n_jobs` times). This is done for efficiency
    reasons if individual jobs take very little time, but may raise errors if
    the dataset is large and not enough memory is available.  A workaround in
    this case is to set `pre_dispatch`. Then, the memory is copied only
    `pre_dispatch` many times. A reasonable value for `pre_dispatch` is `2 *
    n_jobs`.

    See Also
    --------
    :class:`GridSearchCV`:
        Does exhaustive search over a grid of parameters.

    :class:`ParameterSampler`:
        A generator over parameter settins, constructed from
        param_distributions.

    """

    def __init__(self, estimator, param_distributions, n_iter=10, scoring=None,
                 loss_func=None, score_func=None, fit_params=None, n_jobs=1,
                 iid=True, refit=True, cv=None, verbose=0,
                 pre_dispatch='2*n_jobs'):

        self.param_distributions = param_distributions
        self.n_iter = n_iter
        super(RandomizedSearchCV, self).__init__(
            estimator, scoring, loss_func, score_func, fit_params, n_jobs, iid,
            refit, cv, verbose, pre_dispatch)

    def fit(self, X, y=None, **params):
        """Run fit on the estimator with randomly drawn parameters.

        Parameters
        ----------

        X: array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y: array-like, shape = [n_samples], optional
            Target vector relative to X for classification;
            None for unsupervised learning.

        """
        sampled_params = ParameterSampler(self.param_distributions,
                                          self.n_iter)
        return self._fit(X, y, sampled_params, **params)
