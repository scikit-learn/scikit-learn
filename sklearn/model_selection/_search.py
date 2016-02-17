"""
The :mod:`sklearn.model_selection._search` includes utilities to fine-tune the
parameters of an estimator.
"""
from __future__ import print_function

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>
#         Andreas Mueller <amueller@ais.uni-bonn.de>
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Sebastien Dubois <http://bit.ly/SebastienDubois>
#         Djalel Benbouzid <djalel.benbouzid@gmail.com>
#         Fabian Pedregosa <f@bianp.net>
# License: BSD 3 clause

from abc import ABCMeta, abstractmethod
from collections import Mapping, namedtuple, Sized
from functools import partial, reduce
from itertools import product
import operator
import warnings

import numpy as np
from scipy import stats, optimize

from ..base import BaseEstimator, is_classifier, clone
from ..base import MetaEstimatorMixin, ChangedBehaviorWarning
from ._split import check_cv
from ._validation import _fit_and_score
from ..externals.odict import OrderedDict
from ..externals.joblib import Parallel, delayed
from ..externals import six
from ..utils import check_random_state
from ..utils.fixes import sp_version
from ..utils.random import sample_without_replacement
from ..utils.validation import _num_samples, indexable
from ..utils.metaestimators import if_delegate_has_method
from ..metrics.scorer import check_scoring
from ..preprocessing import LabelEncoder
from ..gaussian_process import GaussianProcessRegressor
from ..gaussian_process.kernels import RBF, Matern


__all__ = ['GridSearchCV', 'ParameterGrid', 'fit_grid_point',
           'ParameterSampler', 'RandomizedSearchCV']


def _acquisition_func(x0, gp, prev_best, func, xi, kappa):
    x0 = np.asarray(x0)
    if x0.ndim == 1:
        x0 = np.expand_dims(x0, axis=0)
    predictions, std = gp.predict(x0, return_std=True)

    if func == 'UCB':
        acquisition_func = predictions - kappa * std
    elif func == 'EI':
        # When std is 0.0, Z is huge, safe to say the pdf at Z is 0.0
        # and cdf at Z is 1.0
        std_mask = std != 0.0
        acquisition_func = prev_best - predictions - xi
        Z = acquisition_func[std_mask] / std[std_mask]
        acquisition_func[std_mask] = std[std_mask] * (
            Z * stats.norm.cdf(Z) + stats.norm.pdf(Z))
    else:
        raise ValueError(
            'acquisition_function not implemented yet : '
            + func)

    if acquisition_func.shape == (1, 1):
        return acquisition_func[0]
    return acquisition_func


class ParameterGrid(object):
    """Grid of parameters with a discrete number of values for each.

    Can be used to iterate over parameter value combinations with the
    Python built-in function iter.

    Read more in the :ref:`User Guide <search>`.

    Parameters
    ----------
    param_grid : dict of string to sequence, or sequence of such
        The parameter grid to explore, as a dictionary mapping estimator
        parameters to sequences of allowed values.

        An empty dict signifies default parameters.

        A sequence of dicts signifies a sequence of grids to search, and is
        useful to avoid exploring parameter combinations that make no sense
        or have no effect. See the examples below.

    Examples
    --------
    >>> from sklearn.model_selection import ParameterGrid
    >>> param_grid = {'a': [1, 2], 'b': [True, False]}
    >>> list(ParameterGrid(param_grid)) == (
    ...    [{'a': 1, 'b': True}, {'a': 1, 'b': False},
    ...     {'a': 2, 'b': True}, {'a': 2, 'b': False}])
    True

    >>> grid = [{'kernel': ['linear']}, {'kernel': ['rbf'], 'gamma': [1, 10]}]
    >>> list(ParameterGrid(grid)) == [{'kernel': 'linear'},
    ...                               {'kernel': 'rbf', 'gamma': 1},
    ...                               {'kernel': 'rbf', 'gamma': 10}]
    True
    >>> ParameterGrid(grid)[1] == {'kernel': 'rbf', 'gamma': 1}
    True

    See also
    --------
    :class:`GridSearchCV`:
        Uses :class:`ParameterGrid` to perform a full parallelized parameter
        search.
    """

    def __init__(self, param_grid):
        if isinstance(param_grid, Mapping):
            # wrap dictionary in a singleton list to support either dict
            # or list of dicts
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
            if not items:
                yield {}
            else:
                keys, values = zip(*items)
                for v in product(*values):
                    params = dict(zip(keys, v))
                    yield params

    def __len__(self):
        """Number of points on the grid."""
        # Product function that can handle iterables (np.product can't).
        product = partial(reduce, operator.mul)
        return sum(product(len(v) for v in p.values()) if p else 1
                   for p in self.param_grid)

    def __getitem__(self, ind):
        """Get the parameters that would be ``ind``th in iteration

        Parameters
        ----------
        ind : int
            The iteration index

        Returns
        -------
        params : dict of string to any
            Equal to list(self)[ind]
        """
        # This is used to make discrete sampling without replacement memory
        # efficient.
        for sub_grid in self.param_grid:
            # XXX: could memoize information used here
            if not sub_grid:
                if ind == 0:
                    return {}
                else:
                    ind -= 1
                    continue

            # Reverse so most frequent cycling parameter comes first
            keys, values_lists = zip(*sorted(sub_grid.items())[::-1])
            sizes = [len(v_list) for v_list in values_lists]
            total = np.product(sizes)

            if ind >= total:
                # Try the next grid
                ind -= total
            else:
                out = {}
                for key, v_list, n in zip(keys, values_lists, sizes):
                    ind, offset = divmod(ind, n)
                    out[key] = v_list[offset]
                return out

        raise IndexError('ParameterGrid index out of range')


class ParameterSampler(object):
    """Generator on parameters sampled from given distributions.

    Non-deterministic iterable over random candidate combinations for hyper-
    parameter search. If all parameters are presented as a list,
    sampling without replacement is performed. If at least one parameter
    is given as a distribution, sampling with replacement is used.
    It is highly recommended to use continuous distributions for continuous
    parameters.

    Note that before SciPy 0.16, the ``scipy.stats.distributions`` do not
    accept a custom RNG instance and always use the singleton RNG from
    ``numpy.random``. Hence setting ``random_state`` will not guarantee a
    deterministic iteration whenever ``scipy.stats`` distributions are used to
    define the parameter search space. Deterministic behavior is however
    guaranteed from SciPy 0.16 onwards.

    Read more in the :ref:`User Guide <search>`.

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
        Pseudo random number generator state used for random uniform sampling
        from lists of possible values instead of scipy.stats distributions.

    Returns
    -------
    params : dict of string to any
        **Yields** dictionaries mapping each estimator parameter to
        as sampled value.

    Examples
    --------
    >>> from sklearn.model_selection import ParameterSampler
    >>> from scipy.stats.distributions import expon
    >>> import numpy as np
    >>> np.random.seed(0)
    >>> param_grid = {'a':[1, 2], 'b': expon()}
    >>> param_list = list(ParameterSampler(param_grid, n_iter=4))
    >>> rounded_list = [dict((k, round(v, 6)) for (k, v) in d.items())
    ...                 for d in param_list]
    >>> rounded_list == [{'b': 0.89856, 'a': 1},
    ...                  {'b': 0.923223, 'a': 1},
    ...                  {'b': 1.878964, 'a': 2},
    ...                  {'b': 1.038159, 'a': 2}]
    True
    """
    def __init__(self, param_distributions, n_iter, random_state=None):
        self.param_distributions = param_distributions
        self.n_iter = n_iter
        self.random_state = random_state

    def __iter__(self):
        # check if all distributions are given as lists
        # in this case we want to sample without replacement
        all_lists = np.all([not hasattr(v, "rvs")
                            for v in self.param_distributions.values()])
        rnd = check_random_state(self.random_state)

        if all_lists:
            # look up sampled parameter settings in parameter grid
            param_grid = ParameterGrid(self.param_distributions)
            grid_size = len(param_grid)

            if grid_size < self.n_iter:
                raise ValueError(
                    "The total space of parameters %d is smaller "
                    "than n_iter=%d." % (grid_size, self.n_iter)
                    + " For exhaustive searches, use GridSearchCV.")
            for i in sample_without_replacement(grid_size, self.n_iter,
                                                random_state=rnd):
                yield param_grid[i]

        else:
            # Always sort the keys of a dictionary, for reproducibility
            items = sorted(self.param_distributions.items())
            for _ in six.moves.range(self.n_iter):
                params = dict()
                for k, v in items:
                    if hasattr(v, "rvs"):
                        if sp_version < (0, 16):
                            params[k] = v.rvs()
                        else:
                            params[k] = v.rvs(random_state=rnd)
                    else:
                        params[k] = v[rnd.randint(len(v))]
                yield params

    def __len__(self):
        """Number of points that will be sampled."""
        return self.n_iter


def fit_grid_point(X, y, estimator, parameters, train, test, scorer,
                   verbose, error_score='raise', **fit_params):
    """Run fit on one set of parameters.

    Parameters
    ----------
    X : array-like, sparse matrix or list
        Input data.

    y : array-like or None
        Targets for input data.

    estimator : estimator object
        A object of that type is instantiated for each grid point.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``scoring`` must be passed.

    parameters : dict
        Parameters to be set on estimator for this grid point.

    train : ndarray, dtype int or bool
        Boolean mask or indices for training set.

    test : ndarray, dtype int or bool
        Boolean mask or indices for test set.

    scorer : callable or None.
        If provided must be a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    verbose : int
        Verbosity level.

    **fit_params : kwargs
        Additional parameter passed to the fit function of the estimator.

    error_score : 'raise' (default) or numeric
        Value to assign to the score if an error occurs in estimator fitting.
        If set to 'raise', the error is raised. If a numeric value is given,
        FitFailedWarning is raised. This parameter does not affect the refit
        step, which will always raise the error.

    Returns
    -------
    score : float
        Score of this parameter setting on given training / test split.

    parameters : dict
        The parameters that have been evaluated.

    n_samples_test : int
        Number of test samples in this split.
    """
    score, n_samples_test, _ = _fit_and_score(estimator, X, y, scorer, train,
                                              test, verbose, parameters,
                                              fit_params, error_score)
    return score, parameters, n_samples_test


def _check_param_grid(param_grid):
    if hasattr(param_grid, 'items'):
        param_grid = [param_grid]

    for p in param_grid:
        for v in p.values():
            if isinstance(v, np.ndarray) and v.ndim > 1:
                raise ValueError("Parameter array should be one-dimensional.")

            check = [isinstance(v, k) for k in (list, tuple, np.ndarray)]
            if True not in check:
                raise ValueError("Parameter values should be a list.")

            if len(v) == 0:
                raise ValueError("Parameter values should be a non-empty "
                                 "list.")


class _CVScoreTuple (namedtuple('_CVScoreTuple',
                                ('parameters',
                                 'mean_validation_score',
                                 'cv_validation_scores'))):
    # A raw namedtuple is very memory efficient as it packs the attributes
    # in a struct to get rid of the __dict__ of attributes in particular it
    # does not copy the string for the keys on each instance.
    # By deriving a namedtuple class just to introduce the __repr__ method we
    # would also reintroduce the __dict__ on the instance. By telling the
    # Python interpreter that this subclass uses static __slots__ instead of
    # dynamic attributes. Furthermore we don't need any additional slot in the
    # subclass so we set __slots__ to the empty tuple.
    __slots__ = ()

    def __repr__(self):
        """Simple custom repr to summarize the main info"""
        return "mean: {0:.5f}, std: {1:.5f}, params: {2}".format(
            self.mean_validation_score,
            np.std(self.cv_validation_scores),
            self.parameters)


class BaseSearchCV(six.with_metaclass(ABCMeta, BaseEstimator,
                                      MetaEstimatorMixin)):
    """Base class for hyper parameter search with cross-validation."""

    @abstractmethod
    def __init__(self, estimator, scoring=None,
                 fit_params=None, n_jobs=1, iid=True,
                 refit=True, cv=None, verbose=0, pre_dispatch='2*n_jobs',
                 error_score='raise'):

        self.scoring = scoring
        self.estimator = estimator
        self.n_jobs = n_jobs
        self.fit_params = fit_params if fit_params is not None else {}
        self.iid = iid
        self.refit = refit
        self.cv = cv
        self.verbose = verbose
        self.pre_dispatch = pre_dispatch
        self.error_score = error_score

    @property
    def _estimator_type(self):
        return self.estimator._estimator_type

    def score(self, X, y=None):
        """Returns the score on the given data, if the estimator has been refit.

        This uses the score defined by ``scoring`` where provided, and the
        ``best_estimator_.score`` method otherwise.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Input data, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.

        Returns
        -------
        score : float

        Notes
        -----
         * The long-standing behavior of this method changed in version 0.16.
         * It no longer uses the metric provided by ``estimator.score`` if the
           ``scoring`` parameter was set when fitting.

        """
        if self.scorer_ is None:
            raise ValueError("No score function explicitly defined, "
                             "and the estimator doesn't provide one %s"
                             % self.best_estimator_)
        if self.scoring is not None and hasattr(self.best_estimator_, 'score'):
            warnings.warn("The long-standing behavior to use the estimator's "
                          "score function in {0}.score has changed. The "
                          "scoring parameter is now used."
                          "".format(self.__class__.__name__),
                          ChangedBehaviorWarning)
        return self.scorer_(self.best_estimator_, X, y)

    @if_delegate_has_method(delegate='estimator')
    def predict(self, X):
        """Call predict on the estimator with the best found parameters.

        Only available if ``refit=True`` and the underlying estimator supports
        ``predict``.

        Parameters
        -----------
        X : indexable, length n_samples
            Must fulfill the input assumptions of the
            underlying estimator.

        """
        return self.best_estimator_.predict(X)

    @if_delegate_has_method(delegate='estimator')
    def predict_proba(self, X):
        """Call predict_proba on the estimator with the best found parameters.

        Only available if ``refit=True`` and the underlying estimator supports
        ``predict_proba``.

        Parameters
        -----------
        X : indexable, length n_samples
            Must fulfill the input assumptions of the
            underlying estimator.

        """
        return self.best_estimator_.predict_proba(X)

    @if_delegate_has_method(delegate='estimator')
    def predict_log_proba(self, X):
        """Call predict_log_proba on the estimator with the best found parameters.

        Only available if ``refit=True`` and the underlying estimator supports
        ``predict_log_proba``.

        Parameters
        -----------
        X : indexable, length n_samples
            Must fulfill the input assumptions of the
            underlying estimator.

        """
        return self.best_estimator_.predict_log_proba(X)

    @if_delegate_has_method(delegate='estimator')
    def decision_function(self, X):
        """Call decision_function on the estimator with the best found parameters.

        Only available if ``refit=True`` and the underlying estimator supports
        ``decision_function``.

        Parameters
        -----------
        X : indexable, length n_samples
            Must fulfill the input assumptions of the
            underlying estimator.

        """
        return self.best_estimator_.decision_function(X)

    @if_delegate_has_method(delegate='estimator')
    def transform(self, X):
        """Call transform on the estimator with the best found parameters.

        Only available if the underlying estimator supports ``transform`` and
        ``refit=True``.

        Parameters
        -----------
        X : indexable, length n_samples
            Must fulfill the input assumptions of the
            underlying estimator.

        """
        return self.best_estimator_.transform(X)

    @if_delegate_has_method(delegate='estimator')
    def inverse_transform(self, Xt):
        """Call inverse_transform on the estimator with the best found parameters.

        Only available if the underlying estimator implements
        ``inverse_transform`` and ``refit=True``.

        Parameters
        -----------
        Xt : indexable, length n_samples
            Must fulfill the input assumptions of the
            underlying estimator.

        """
        return self.best_estimator_.transform(Xt)

    def _fit(self, X, y, labels, parameter_iterable):
        """Actual fitting,  performing the search over parameters."""

        estimator = self.estimator
        cv = check_cv(self.cv, y, classifier=is_classifier(estimator))
        self.scorer_ = check_scoring(self.estimator, scoring=self.scoring)

        n_samples = _num_samples(X)
        X, y, labels = indexable(X, y, labels)

        if y is not None:
            if len(y) != n_samples:
                raise ValueError('Target variable (y) has a different number '
                                 'of samples (%i) than data (X: %i samples)'
                                 % (len(y), n_samples))
        n_splits = cv.get_n_splits(X, y, labels)

        if self.verbose > 0 and isinstance(parameter_iterable, Sized):
            n_candidates = len(parameter_iterable)
            print("Fitting {0} folds for each of {1} candidates, totalling"
                  " {2} fits".format(n_splits, n_candidates,
                                     n_candidates * n_splits))

        base_estimator = clone(self.estimator)

        pre_dispatch = self.pre_dispatch

        out = Parallel(
            n_jobs=self.n_jobs, verbose=self.verbose,
            pre_dispatch=pre_dispatch
        )(delayed(_fit_and_score)(clone(base_estimator), X, y, self.scorer_,
                                  train, test, self.verbose, parameters,
                                  self.fit_params, return_parameters=True,
                                  error_score=self.error_score)
          for parameters in parameter_iterable
          for train, test in cv.split(X, y, labels))

        # Out is a list of triplet: score, estimator, n_test_samples
        n_fits = len(out)

        scores = list()
        grid_scores = list()
        for grid_start in range(0, n_fits, n_splits):
            n_test_samples = 0
            score = 0
            all_scores = []
            for this_score, this_n_test_samples, _, parameters in \
                    out[grid_start:grid_start + n_splits]:
                all_scores.append(this_score)
                if self.iid:
                    this_score *= this_n_test_samples
                    n_test_samples += this_n_test_samples
                score += this_score
            if self.iid:
                score /= float(n_test_samples)
            else:
                score /= float(n_splits)
            scores.append((score, parameters))
            # TODO: shall we also store the test_fold_sizes?
            grid_scores.append(_CVScoreTuple(
                parameters,
                score,
                np.array(all_scores)))
        # Store the computed scores
        self.grid_scores_ = grid_scores

        # Find the best parameters by comparing on the mean validation score:
        # note that `sorted` is deterministic in the way it breaks ties
        best = sorted(grid_scores, key=lambda x: x.mean_validation_score,
                      reverse=True)[0]
        self.best_params_ = best.parameters
        self.best_score_ = best.mean_validation_score

        if self.refit:
            # fit the best estimator using the entire dataset
            # clone first to work around broken estimators
            best_estimator = clone(base_estimator).set_params(
                **best.parameters)
            if y is not None:
                best_estimator.fit(X, y, **self.fit_params)
            else:
                best_estimator.fit(X, **self.fit_params)
            self.best_estimator_ = best_estimator
        return self


class GridSearchCV(BaseSearchCV):
    """Exhaustive search over specified parameter values for an estimator.

    Important members are fit, predict.

    GridSearchCV implements a "fit" and a "score" method.
    It also implements "predict", "predict_proba", "decision_function",
    "transform" and "inverse_transform" if they are implemented in the
    estimator used.

    The parameters of the estimator used to apply these methods are optimized
    by cross-validated grid-search over a parameter grid.

    Read more in the :ref:`User Guide <grid_search>`.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``scoring`` must be passed.

    param_grid : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    fit_params : dict, optional
        Parameters to pass to the fit method.

    n_jobs : int, default=1
        Number of jobs to run in parallel.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    iid : boolean, default=True
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross validation,
          - integer, to specify the number of folds in a `(Stratified)KFold`,
          - An object to be used as a cross-validation generator.
          - An iterable yielding train, test splits.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    refit : boolean, default=True
        Refit the best estimator with the entire dataset.
        If "False", it is impossible to make predictions using
        this GridSearchCV instance after fitting.

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    error_score : 'raise' (default) or numeric
        Value to assign to the score if an error occurs in estimator fitting.
        If set to 'raise', the error is raised. If a numeric value is given,
        FitFailedWarning is raised. This parameter does not affect the refit
        step, which will always raise the error.


    Examples
    --------
    >>> from sklearn import svm, datasets
    >>> from sklearn.model_selection import GridSearchCV
    >>> iris = datasets.load_iris()
    >>> parameters = {'kernel':('linear', 'rbf'), 'C':[1, 10]}
    >>> svr = svm.SVC()
    >>> clf = GridSearchCV(svr, parameters)
    >>> clf.fit(iris.data, iris.target)
    ...                             # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    GridSearchCV(cv=None, error_score=...,
           estimator=SVC(C=1.0, cache_size=..., class_weight=..., coef0=...,
                         decision_function_shape=None, degree=..., gamma=...,
                         kernel='rbf', max_iter=-1, probability=False,
                         random_state=None, shrinking=True, tol=...,
                         verbose=False),
           fit_params={}, iid=..., n_jobs=1,
           param_grid=..., pre_dispatch=..., refit=...,
           scoring=..., verbose=...)


    Attributes
    ----------
    grid_scores_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

            * ``parameters``, a dict of parameter settings
            * ``mean_validation_score``, the mean score over the
              cross-validation folds
            * ``cv_validation_scores``, the list of scores for each fold

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.

    scorer_ : function
        Scorer function used on the held out data to choose the best
        parameters for the model.

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

    :func:`sklearn.model_selection.train_test_split`:
        utility function to split the data into a development set usable
        for fitting a GridSearchCV instance and an evaluation set for
        its final evaluation.

    :func:`sklearn.metrics.make_scorer`:
        Make a scorer from a performance metric or loss function.

    """

    def __init__(self, estimator, param_grid, scoring=None, fit_params=None,
                 n_jobs=1, iid=True, refit=True, cv=None, verbose=0,
                 pre_dispatch='2*n_jobs', error_score='raise'):

        super(GridSearchCV, self).__init__(
            estimator=estimator, scoring=scoring, fit_params=fit_params,
            n_jobs=n_jobs, iid=iid, refit=refit, cv=cv, verbose=verbose,
            pre_dispatch=pre_dispatch, error_score=error_score)
        self.param_grid = param_grid
        _check_param_grid(param_grid)

    def fit(self, X, y=None, labels=None):
        """Run fit with all sets of parameters.

        Parameters
        ----------

        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.
        """
        return self._fit(X, y, labels, ParameterGrid(self.param_grid))


class RandomizedSearchCV(BaseSearchCV):
    """Randomized search on hyper parameters.

    RandomizedSearchCV implements a "fit" and a "score" method.
    It also implements "predict", "predict_proba", "decision_function",
    "transform" and "inverse_transform" if they are implemented in the
    estimator used.

    The parameters of the estimator used to apply these methods are optimized
    by cross-validated search over parameter settings.

    In contrast to GridSearchCV, not all parameter values are tried out, but
    rather a fixed number of parameter settings is sampled from the specified
    distributions. The number of parameter settings that are tried is
    given by n_iter.

    If all parameters are presented as a list,
    sampling without replacement is performed. If at least one parameter
    is given as a distribution, sampling with replacement is used.
    It is highly recommended to use continuous distributions for continuous
    parameters.

    Read more in the :ref:`User Guide <randomized_parameter_search>`.

    Parameters
    ----------
    estimator : estimator object.
        A object of that type is instantiated for each grid point.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``scoring`` must be passed.

    param_distributions : dict
        Dictionary with parameters names (string) as keys and distributions
        or lists of parameters to try. Distributions must provide a ``rvs``
        method for sampling (such as those from scipy.stats.distributions).
        If a list is given, it is sampled uniformly.

    n_iter : int, default=10
        Number of parameter settings that are sampled. n_iter trades
        off runtime vs quality of the solution.

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    fit_params : dict, optional
        Parameters to pass to the fit method.

    n_jobs : int, default=1
        Number of jobs to run in parallel.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    iid : boolean, default=True
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross validation,
          - integer, to specify the number of folds in a `(Stratified)KFold`,
          - An object to be used as a cross-validation generator.
          - An iterable yielding train, test splits.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    refit : boolean, default=True
        Refit the best estimator with the entire dataset.
        If "False", it is impossible to make predictions using
        this RandomizedSearchCV instance after fitting.

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    random_state : int or RandomState
        Pseudo random number generator state used for random uniform sampling
        from lists of possible values instead of scipy.stats distributions.

    error_score : 'raise' (default) or numeric
        Value to assign to the score if an error occurs in estimator fitting.
        If set to 'raise', the error is raised. If a numeric value is given,
        FitFailedWarning is raised. This parameter does not affect the refit
        step, which will always raise the error.

    Attributes
    ----------
    grid_scores_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

            * ``parameters``, a dict of parameter settings
            * ``mean_validation_score``, the mean score over the
              cross-validation folds
            * ``cv_validation_scores``, the list of scores for each fold

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.

    Notes
    -----
    The parameters selected are those that maximize the score of the held-out
    data, according to the scoring parameter.

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
                 fit_params=None, n_jobs=1, iid=True, refit=True, cv=None,
                 verbose=0, pre_dispatch='2*n_jobs', random_state=None,
                 error_score='raise'):

        self.param_distributions = param_distributions
        self.n_iter = n_iter
        self.random_state = random_state
        super(RandomizedSearchCV, self).__init__(
            estimator=estimator, scoring=scoring, fit_params=fit_params,
            n_jobs=n_jobs, iid=iid, refit=refit, cv=cv, verbose=verbose,
            pre_dispatch=pre_dispatch, error_score=error_score)

    def fit(self, X, y=None, labels=None):
        """Run fit on the estimator with randomly drawn parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.
        """
        sampled_params = ParameterSampler(self.param_distributions,
                                          self.n_iter,
                                          random_state=self.random_state)
        return self._fit(X, y, labels, sampled_params)


class SequentialSearchCV(BaseSearchCV):
    """Search over hyperparameters based on a surrogate model.

    The best hyperparameters of the estimator are estimated
    by cross-validated search over the provided hyperparameter settings.

    In contrast to GridSearchCV, not all parameter values are tried out.
    Initially, a fixed number of settings are initially tried out for and
    the cross-validation scores are evaluated.
    Given these initial evaluations, the new best hyperparameter is chosen
    by maximizing a function (commonly called acquisition function).
    The number of initial random parameter settings that are tried for is given by
    ``n_init`` and the total number parameter settings that are tried for
    is given by ``n_iter``.

    Read more in the :ref:`User Guide <gp_search>`.

    Parameters
    ----------
    estimator : estimator object.
        A object of that type is instantiated for each grid point.
        This is assumed to implement the scikit-learn estimator interface.
        The estimator must implement a ``fit`` method and can either
        provide a ``score`` function, or ``scoring`` must be passed.

    param_distributions : dict of dicts
        Dictionary with parameters names (string) as keys and a dict
        containing "bounds", "scale" and "type" as values.

        - "scale": If scale is set to be "log", then sample uniformly from
           (log(upper_bound), log(lower_bounds)). This is useful in
           hyperparameters like `alpha` in `Ridge` where search has to
           be done on a log scale.

        - "bounds": Tuple containing the lower bounds and the upper bound.
           By default, the bounds are taken to be from 0.0 to 1.0
           If the scale is specified to be logscale, then the default
           bounds are set to be (10**-5, 10**5)

        - "type: int or float, specifying what type each parameter is.
           If the type is not specified, the default is assumed to
           be float.

    n_init : int, optional
        Number of random iterations to perform before the smart search.
        Default is 3.

    n_iter : int, default=10
        Number of parameter settings that are sampled. n_iter trades
        off runtime vs quality of the solution.

    acquisition_function : string, optional
        Function to maximize in order to choose the next parameter to test.
        - 'UCB : maximize the upper confidence bound
        - 'EI' : maximizes the expected improvement
        Default is 'UCB'

    scoring : string, callable or None, default=None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        If ``None``, the ``score`` method of the estimator is used.

    fit_params : dict, optional
        Parameters to pass to the fit method.

    n_jobs : int, default=1
        Number of jobs to run in parallel.

    pre_dispatch : int, or string, optional default "2*n_jobs"
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    iid : boolean, default=True
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - An object to be used as a cross-validation generator.
          - An iterable yielding train/test splits.

        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`StratifiedKFold` used. For all other cases
        :class:`KFold`` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    refit : boolean, default=True
        Refit the best estimator with the entire dataset.
        If "False", it is impossible to make predictions using
        this RandomizedSearchCV instance after fitting.

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    random_state : int or RandomState
        Pseudo random number generator state used for random uniform sampling
        from lists of possible values instead of scipy.stats distributions.

    error_score : 'raise' (default) or numeric
        Value to assign to the score if an error occurs in estimator fitting.
        If set to 'raise', the error is raised. If a numeric value is given,
        FitFailedWarning is raised. This parameter does not affect the refit
        step, which will always raise the error.

    n_candidates : int, default=500
        Number of random candidates to sample for each GP iteration
        Default is 500.

    gp_params : dict, default={}
        Parameters to pass to the GaussianProcessRegressor object.

    xi : float, default=0.01
        Exploration-exploitation trade-off parameter for the expected
        improvement acquisition function. Higher values for this parameter
        yield updates that are more conservative in which the exploitation
        criterion dominates over the exploration criterion.

    kappa: float, default=1.96
        Exploration-exploitation trade-off parameter for the upper
        confidence bound acquisition function. The default value of
        1.96 corresponds to the 95% confidence bound on the predicted
        objective function.

    Attributes
    ----------
    grid_scores_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.
        Each named tuple has the attributes:

            * ``parameters``, a dict of parameter settings
            * ``mean_validation_score``, the mean score over the
              cross-validation folds
            * ``cv_validation_scores``, the list of scores for each fold

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.

    Notes
    -----
    The parameters selected are those that maximize the score of the held-out
    data, according to the scoring parameter.

    See Also
    --------
    :class:`GridSearchCV`:
        Does exhaustive search over a grid of parameters.

    :class:`RandomizedSearchCV`:
        Randomized search on hyper parameters.

    Examples
    --------
    >>> from sklearn import linear_model, datasets
    >>> from sklearn.model_selection import SequentialSearchCV
    >>> from scipy import stats
    >>> X, y = datasets.make_regression()
    >>> clf = linear_model.Ridge()
    >>> parameters = {"alpha": {'bounds': (10**-5, 10**5), 'scale': 'log'}}
    >>> search = SequentialSearchCV(clf, parameters)
    >>> search.fit(X, y)
    ...                             # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    SequentialSearchCV(acquisition_function='UCB', cv=None, error_score='raise',
          estimator=Ridge(alpha=1.0, copy_X=True, fit_intercept=True, ...
       normalize=False, random_state=None, solver='auto', tol=0.001),
          fit_params={}, gp_params={}, iid=True, kappa=1.96, n_candidates=500,
          n_init=3, n_iter=10, n_jobs=1,
          param_distributions={'alpha': ...},
          pre_dispatch='2*n_jobs', random_state=None, refit=True, scoring=None,
          verbose=0, xi=0.01)

    """

    def __init__(self, estimator, param_distributions, n_init=3, n_iter=10,
                 scoring=None, fit_params=None, acquisition_function='UCB',
                 n_jobs=1, iid=True, refit=True,
                 cv=None, verbose=0, pre_dispatch='2*n_jobs',
                 random_state=None, error_score='raise', n_candidates=500,
                 gp_params={}, xi=0.01, kappa=1.96,
                 search='local'):
        self.acquisition_function = acquisition_function
        self.param_distributions = param_distributions
        self.n_init = n_init
        self.n_iter = n_iter
        self.random_state = random_state
        self.gp_params = gp_params
        self.xi = xi
        self.kappa = kappa
        self.n_candidates = n_candidates
        self.search = search

        super(SequentialSearchCV, self).__init__(
            estimator=estimator, scoring=scoring, fit_params=fit_params,
            n_jobs=n_jobs, iid=iid, refit=refit, cv=cv, verbose=verbose,
            pre_dispatch=pre_dispatch, error_score=error_score)

    def _evaluate_params(self, candidates, cv, X, y, labels):

        # Use the context-manager API to reuse processes for every
        # call to Parallel.
        with Parallel(
            n_jobs=self.n_jobs, verbose=self.verbose,
            pre_dispatch=self.pre_dispatch) as parallel:

            tested_parameters, cv_scores, grid_scores = [], [], []

            for i, params in enumerate(candidates):

                # We need this because we need the paramter strings
                # in the same order which params is stored.
                if i == 0:
                    params_list = list(params.keys())

                scaled_candidates = self.scale_candidates(params)
                all_scores = parallel(delayed(_fit_and_score)(
                    clone(self.estimator), X, y, self.scorer_,
                    train, test, self.verbose, scaled_candidates,
                    self.fit_params, return_parameters=False,
                    error_score=self.error_score)
                    for train, test in cv.split(X, y, labels))

                all_scores = np.asarray(all_scores)
                fold_scores = all_scores[:, 0]
                n_test_samples = all_scores[:, 1]

                if self.iid:
                    score = np.sum(
                        fold_scores * n_test_samples) / np.sum(n_test_samples)
                else:
                    score = np.mean(fold_scores)

                if self.verbose > 0:
                    print('Step ' + str(i) + ' - Hyperparameter '
                          + str(params) + ', score: ' + str(score))

                tested_parameters.append(list(params.values()))
                cv_scores.append(score)
                grid_scores.append(
                    _CVScoreTuple(scaled_candidates, score, fold_scores))

        return tested_parameters, cv_scores, grid_scores, params_list

    def scale_candidates(self, params_dict):

        scaled_candidates = dict()
        for param, val in params_dict.items():
            bounds = self._param_distributions[param]['bounds']
            scale = self._param_distributions[param]['scale']
            param_type = self._param_distributions[param]['type']

            if scale == "log":
                bounds = np.log10(bounds)
            scaled_values = val * (bounds[-1] - bounds[0]) + bounds[0]

            if scale == "log":
                scaled_values = 10**scaled_values

            if param_type == int:
                scaled_values = int(scaled_values)
            scaled_candidates[param] = scaled_values
        return scaled_candidates

    def fit(self, X, y=None, labels=None):
        """
        Run the hyper-parameter optimization process

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        """
        if self.acquisition_function not in ('UCB', 'EI', 'greedy'):
            raise ValueError(
                "acquisition_function should be one of " +
                "'UCB', 'EI' or 'greedy'. Found %s instead" %
                self.acquisition_function)

        if self.n_init > self.n_iter:
            raise ValueError('n_init cannot be larger than n_iter')

        if self.xi <= 0:
            raise ValueError('Parameter xi must be positive')

        if self.kappa <= 0:
            raise ValueError('Parameter kappa must be positive')

        # Preprocess self.param_distributions to store default values.
        param_distributions = dict()
        for param, param_describe in self.param_distributions.items():
            new_param_describe = dict()
            if not isinstance(param_describe, dict):
                raise ValueError(
                    "Expected type dict but %s provided" %
                    type(param_describe))

            scale = param_describe.get("scale", None)
            if scale is not None and scale != "log":
                raise ValueError(
                    "Only uniform and log supported. Got %s" % scale)
            new_param_describe["scale"] = scale

            bounds = param_describe.get("bounds", None)
            if bounds is None:
                if scale is None:
                    bounds = (0.0, 1.0)
                else:
                    bounds = (10**-5, 10**5)
            if len(bounds) != 2 or bounds[-1] <= bounds[0] or not np.iterable(bounds):
                raise ValueError("Provide reasonable bounds.")
            new_param_describe["bounds"] = bounds

            param_type = param_describe.get("type", float)
            if param_type not in [float, int]:
                raise ValueError(
                    "Only int and float supported. Got %s" % param_type)
            new_param_describe["type"] = param_type
            param_distributions[param] = new_param_describe

        self._param_distributions = param_distributions

        samples_uniform = dict()
        for param, _ in param_distributions.items():
            samples_uniform[param] = stats.uniform()

        cv_scores, grid_scores, tested_parameters = [], [], []
        self.scorer_ = check_scoring(self.estimator, scoring=self.scoring)
        cv = check_cv(self.cv, y, classifier=is_classifier(self.estimator))

        init_candidates = ParameterSampler(
            param_distributions=samples_uniform, n_iter=self.n_init,
            random_state=self.random_state)

        tested_parameters, cv_scores, grid_scores, params_list = \
            self._evaluate_params(
                init_candidates, cv, X, y, labels)
        n_hyperparam_dims = len(params_list)
        bounds = np.tile(
            [0.0, 1.0],
            (n_hyperparam_dims, 1))

        gp_scores = (-np.asarray(cv_scores)).tolist()
        gp_params = self.gp_params
        if not gp_params:
            length_scale = np.ones_like(tested_parameters[0])
            gp_params = {
                'kernel': Matern(length_scale=length_scale, nu=2.5),
                'normalize_y': True,
                'random_state': self.random_state}

        rng = check_random_state(self.random_state)
        for i in range(self.n_iter - self.n_init):
            # Model with a Gaussian Process
            gp = GaussianProcessRegressor(**self.gp_params)
            gp.fit(tested_parameters, gp_scores)
            best_score = np.min(gp_scores)

            if self.search == "local":
                candidates = ParameterSampler(
                    samples_uniform, self.n_candidates,
                    random_state=rng)
                candidate = list(list(candidates)[0].values())

                res = optimize.fmin_l_bfgs_b(
                    _acquisition_func,
                    np.asfortranarray(candidate),
                    args=(gp, best_score, self.acquisition_function, self.xi, self.kappa),
                    bounds=bounds, approx_grad=True, maxiter=10)
                best_candidate = dict(zip(params_list, res[0]))
            elif self.search == "sampling":
                # Sample candidates and predict their corresponding
                # acquisition values
                candidates = ParameterSampler(
                    samples_uniform, self.n_candidates,
                    random_state=rng)
                candidate_values = []
                for s in candidates:
                    candidate_values.append(list(s.values()))
                candidate_values = np.asarray(candidate_values)

                acquisition_vals = _acquisition_func(
                    candidate_values, gp, best_score,
                    self.acquisition_function,
                    self.xi, self.kappa)

                best_idx = np.argmin(acquisition_vals)
                best_candidate = dict(zip(params_list, candidate_values[best_idx]))
            else:
                raise ValueError("Search not supported.")

            best_param, best_score, best_grid_score, _ = self._evaluate_params(
                [best_candidate], cv, X, y, labels)
            tested_parameters.append(best_param[0])
            gp_scores.append(-best_score[0])
            grid_scores.append(best_grid_score[0])
            if self.verbose > 0:
                print('Step ' + str(i+self.n_init) + ' - Hyperparameter '
                      + str(best_candidate) + ', score: ' + str(best_score))

        # Store the computed scores
        self.grid_scores_ = grid_scores

        # Find the best parameters by comparing on the mean validation score:
        # note that `sorted` is deterministic in the way it breaks ties
        # that is, it picks the smallest parameter.
        # XXX: The second key should be the standard deviation of the validation
        # scores but the current grid_scores_ attribute does not allow
        # that.

        best = sorted(
            grid_scores,
            key=lambda x: (
                x.mean_validation_score
                )
            )[-1]

        self.best_params_ = best.parameters
        self.best_score_ = best.mean_validation_score

        if self.refit:
            # fit the best estimator using the entire dataset
            # clone first to work around broken estimators
            best_estimator = clone(self.estimator).set_params(
                **best.parameters)
            if y is not None:
                best_estimator.fit(X, y, **self.fit_params)
            else:
                best_estimator.fit(X, **self.fit_params)
            self.best_estimator_ = best_estimator
        return self
