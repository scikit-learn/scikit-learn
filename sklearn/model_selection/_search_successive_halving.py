from math import ceil, floor, log
from abc import abstractmethod
from collections import OrderedDict
from numbers import Integral

import numpy as np
from ._search import _check_param_grid
from ._search import BaseSearchCV
from . import ParameterGrid, ParameterSampler
from ..utils import check_random_state
from ..utils.validation import _num_samples
from ..base import is_classifier
from ._split import check_cv
from ..utils import resample


__all__ = ['HalvingGridSearchCV', 'HalvingRandomSearchCV']


def _refit_callable(results):
    # Custom refit callable to return the index of the best candidate. We want
    # the best candidate out of the last iteration. By default BaseSearchCV
    # would return the best candidate out of all iterations.

    last_iter = np.max(results['iter'])
    sorted_indices = np.argsort(results['mean_test_score'])[::-1]
    best_index = next(i for i in sorted_indices
                      if results['iter'][i] == last_iter)
    return best_index


class BaseSuccessiveHalving(BaseSearchCV):
    """Implements successive halving.

    Ref:
    Almost optimal exploration in multi-armed bandits, ICML 13
    Zohar Karnin, Tomer Koren, Oren Somekh
    """
    def __init__(self, estimator, scoring=None,
                 n_jobs=None, refit=True, cv=5, verbose=0,
                 pre_dispatch='2*n_jobs', random_state=None,
                 error_score=np.nan, return_train_score=True,
                 max_resources='auto', min_resources='auto',
                 resource='n_samples', ratio=3, aggressive_elimination=False,
                 force_exhaust_resources=False):

        refit = _refit_callable if refit else False
        super().__init__(estimator, scoring=scoring,
                         n_jobs=n_jobs, refit=refit, cv=cv,
                         verbose=verbose, pre_dispatch=pre_dispatch,
                         error_score=error_score,
                         return_train_score=return_train_score)

        self.random_state = random_state
        self.max_resources = max_resources
        self.resource = resource
        self.ratio = ratio
        self.min_resources = min_resources
        self.aggressive_elimination = aggressive_elimination
        self.force_exhaust_resources = force_exhaust_resources

    def _check_input_parameters(self, X, y, groups):

        if groups is not None:
            raise ValueError('groups are not supported.')

        if self.scoring is not None and not (isinstance(self.scoring, str)
                                             or callable(self.scoring)):
            raise ValueError('scoring parameter must be a string, '
                             'a callable or None. Multimetric scoring is not '
                             'supported.')

        if (self.resource != 'n_samples'
                and self.resource not in self.estimator.get_params()):
            raise ValueError(
                'Cannot use resource={} which is not supported '
                'by estimator {}'.format(self.resource,
                                         self.estimator.__class__.__name__))

        if (isinstance(self.max_resources, str) and
                self.max_resources != 'auto'):
            raise ValueError(
                "max_resources must be either 'auto' or a positive integer"
            )
        if self.max_resources != 'auto' and (
                not isinstance(self.max_resources, Integral) or
                self.max_resources <= 0):
            raise ValueError(
                "max_resources must be either 'auto' or a positive integer"
            )

        if (isinstance(self.min_resources, str) and
                self.min_resources != 'auto'):
            raise ValueError(
                "min_resources must be either 'auto' or a positive integer "
                "no greater than max_resources."
            )
        if self.min_resources != 'auto' and (
                not isinstance(self.min_resources, Integral) or
                self.min_resources <= 0):
            raise ValueError(
                "min_resources must be either 'auto' or a positive integer "
                "no greater than max_resources."
            )

        if self.force_exhaust_resources and self.min_resources != 'auto':
            raise ValueError(
                'min_resources must be set to auto if force_exhaust_resources'
                ' is True.'
            )

        self.min_resources_ = self.min_resources
        if self.min_resources_ == 'auto':
            if self.resource == 'n_samples':
                cv = check_cv(self.cv, y,
                              classifier=is_classifier(self.estimator))
                n_splits = cv.get_n_splits(X, y, groups)

                # please see https://gph.is/1KjihQe for a justification
                magic_factor = 2
                self.min_resources_ = n_splits * magic_factor
                if is_classifier(self.estimator):
                    n_classes = np.unique(y).shape[0]
                    self.min_resources_ *= n_classes
            else:
                self.min_resources_ = 1

        self.max_resources_ = self.max_resources
        if self.max_resources_ == 'auto':
            if not self.resource == 'n_samples':
                raise ValueError(
                    "max_resources can only be 'auto' if resource='n_samples'")
            self.max_resources_ = _num_samples(X)

        if self.min_resources_ > self.max_resources_:
            raise ValueError(
                'min_resources_={} is greater than max_resources_={}.'
                .format(self.min_resources_, self.max_resources_)
            )

    def fit(self, X, y=None, groups=None, **fit_params):
        """Run fit with all sets of parameters.

        Parameters
        ----------

        X : array-like, shape (n_samples, n_features)
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape (n_samples,) or (n_samples, n_output), optional
            Target relative to X for classification or regression;
            None for unsupervised learning.

        groups : None
            Groups are not supported

        **fit_params : dict of string -> object
            Parameters passed to the ``fit`` method of the estimator
        """
        self._check_input_parameters(
            X=X,
            y=y,
            groups=groups,
        )
        super().fit(X, y=y, groups=None, **fit_params)
        # Set best_score_: BaseSearchCV does not set it, as refit is a callable
        self.best_score_ = (
            self.cv_results_['mean_test_score'][self.best_index_])
        return self

    def _run_search(self, evaluate_candidates, X, y, **fit_params):
        rng = check_random_state(self.random_state)

        candidate_params = self._generate_candidate_params()
        # Remove duplicates (may happen with random sampling)
        candidate_params = set(tuple(d.items()) for d in candidate_params)
        candidate_params = [dict(t) for t in candidate_params]

        if self.resource != 'n_samples' and any(
                self.resource in candidate for candidate in candidate_params):
            # Can only check this now since we need the candidates list
            raise ValueError(
                "Cannot use parameter {} as the resource since it is part of "
                "the searched parameters.".format(self.resource))

        # n_required_iterations is the number of iterations needed so that the
        # last iterations evaluates less than `ratio` candidates.
        n_required_iterations = 1 + floor(log(len(candidate_params),
                                              self.ratio))

        if self.force_exhaust_resources:
            # To exhaust the resources, we want to start with the biggest
            # min_resources possible so that the last (required) iteration
            # uses as many resources as possible
            # We only force exhausting the resources if min_resources wasn't
            # specified by the user.
            last_iteration = n_required_iterations - 1
            self.min_resources_ = max(
                self.min_resources_,
                self.max_resources_ // self.ratio**last_iteration
            )

        # n_possible_iterations is the number of iterations that we can
        # actually do starting from min_resources and without exceeding the
        # max_resources. Depending on max_resources and the number of
        # candidates, this may be higher or smaller than
        # n_required_iterations.
        n_possible_iterations = 1 + floor(log(
            self.max_resources_ // self.min_resources_, self.ratio))

        if self.aggressive_elimination:
            n_iterations = n_required_iterations
        else:
            n_iterations = min(n_possible_iterations, n_required_iterations)

        if self.verbose:
            print('n_iterations: {}'.format(n_iterations))
            print('n_required_iterations: {}'.format(n_required_iterations))
            print('n_possible_iterations: {}'.format(n_possible_iterations))
            print('min_resources_: {}'.format(self.min_resources_))
            print('max_resources_: {}'.format(self.max_resources_))
            print('aggressive_elimination: {}'.format(
                self.aggressive_elimination))
            print('force_exhaust_resources: {}'.format(
                self.force_exhaust_resources))
            print('ratio: {}'.format(self.ratio))

        # list of resource_iter for each iteration, used in tests
        self._r_i_list = []
        self.n_candidates_ = []

        for iter_i in range(n_iterations):

            power = iter_i  # default
            if self.aggressive_elimination:
                # this will set resource_iter to the initial value (i.e. the
                # value of resource_iter at the first iteration) for as many
                # iterations as needed (while candidates are being
                # eliminated), and then go on as usual.
                power = max(
                    0,
                    iter_i - n_required_iterations + n_possible_iterations
                )

            resource_iter = int(self.ratio**power * self.min_resources_)
            # guard, probably not needed
            resource_iter = min(resource_iter, self.max_resources_)
            self._r_i_list.append(resource_iter)

            n_candidates = len(candidate_params)
            self.n_candidates_.append(n_candidates)

            if self.verbose:
                print('-' * 10)
                print('iter_i: {}'.format(iter_i))
                print('n_candidates: {}'.format(n_candidates))
                print('resource_iter: {}'.format(resource_iter))

            if self.resource == 'n_samples':
                # Subsample X and y as well as fit_params
                stratify = y if is_classifier(self.estimator) else None
                fit_params = OrderedDict(fit_params)
                X_iter, y_iter, *fit_params_iter_list = resample(
                    X, y, *fit_params.values(), replace=False,
                    random_state=rng, stratify=stratify,
                    n_samples=resource_iter)
                fit_params_iter = {
                    key: fit_params_iter_list[i]
                    for (i, key) in enumerate(fit_params.keys())
                }
            else:
                # Need copy so that the resource_iter of next iteration does
                # not overwrite
                candidate_params = [c.copy() for c in candidate_params]
                for candidate in candidate_params:
                    candidate[self.resource] = resource_iter
                X_iter, y_iter = X, y
                fit_params_iter = fit_params

            more_results = {'iter': [iter_i] * n_candidates,
                            'resource_iter': [resource_iter] * n_candidates}
            results = evaluate_candidates(candidate_params, X_iter, y_iter,
                                          more_results=more_results,
                                          **fit_params_iter)

            n_candidates_to_keep = ceil(n_candidates / self.ratio)
            candidate_params = self._top_k(results,
                                           n_candidates_to_keep,
                                           iter_i)

        self.n_remaining_candidates_ = len(candidate_params)
        self.n_required_iterations_ = n_required_iterations
        self.n_possible_iterations_ = n_possible_iterations
        self.n_iterations_ = n_iterations

    def _top_k(self, results, k, iter_i):
        # Return the best candidates of a given iteration
        # We need to filter out candidates from the previous iterations
        # when sorting

        best_candidates_indices = np.argsort(results['mean_test_score'])[::-1]
        best_candidates_indices = [idx for idx in best_candidates_indices
                                   if results['iter'][idx] == iter_i]
        best_candidates_indices = best_candidates_indices[:k]
        return [results['params'][idx] for idx in best_candidates_indices]

    @abstractmethod
    def _generate_candidate_params(self):
        pass


class HalvingGridSearchCV(BaseSuccessiveHalving):
    """Search over specified parameter values with successive halving.

    The search strategy starts evaluating all the candidates with a small
    amount of resources and iteratively selects the best candidates, using
    more and more resources.

    Read more in the :ref:`User guide<successive_halving_user_guide>`.

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

    scoring : string, callable, or None, default=None
        A single string (see :ref:`scoring_parameter`) or a callable
        (see :ref:`scoring`) to evaluate the predictions on the test set.
        If None, the estimator's score method is used.

    n_jobs : int or None, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

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
              as in '2*n_jobs' (default)

    cv : int, cross-validation generator or iterable, default=5
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - integer, to specify the number of folds in a `(Stratified)KFold`,
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` is used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    refit : boolean, default=True
        If True, refit an estimator using the best found parameters on the
        whole dataset.

        The refitted estimator is made available at the ``best_estimator_``
        attribute and permits using ``predict`` directly on this
        ``GridSearchCV`` instance.

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    error_score : 'raise' or numeric
        Value to assign to the score if an error occurs in estimator fitting.
        If set to 'raise', the error is raised. If a numeric value is given,
        FitFailedWarning is raised. This parameter does not affect the refit
        step, which will always raise the error. Default is ``np.nan``

    return_train_score : boolean, default=False
        If ``False``, the ``cv_results_`` attribute will not include training
        scores.
        Computing training scores is used to get insights on how different
        parameter settings impact the overfitting/underfitting trade-off.
        However computing the scores on the training set can be computationally
        expensive and is not strictly required to select the parameters that
        yield the best generalization performance.

    max_resources : int, default='auto'
        The maximum number of resources that any candidate is allowed to use
        for a given iteration. By default, this is set ``n_samples`` when
        ``resource='n_samples'`` (default), else an error is raised.

    min_resources : int, default='auto'
        The minimum amount of resource that any candidate is allowed to use for
        a given iteration. Equivalently, this defines the amount of resources
        that are allocated for each candidate at the first iteration. By
        default, this is set to:

        - ``n_splits * 2`` when ``resource='n_samples'`` for a regression
          problem
        - ``n_classes * n_splits * 2`` when ``resource='n_samples'`` for a
          regression problem
        - The highest possible value satisfying the constraint
          ``force_exhaust_resources=True``.
        - ``1`` when ``resource!='n_samples'``

        Note that the amount of resources used at each iteration is always a
        multiple of ``min_resources``.

    resource : ``'n_samples'`` or str, default='n_samples'
        Defines the resource that increases with each iteration. By default,
        the resource is the number of samples. It can also be set to any
        parameter of the base estimator that accepts positive integer
        values, e.g. 'n_iterations' or 'n_estimators' for a gradient
        boosting estimator. In this case ``max_resources`` cannot be 'auto'.

    ratio : int or float, default=3
        The 'halving' parameter, which determines the proportion of candidates
        that are selected for the next iteration. For example, ``ratio=3``
        means that only one third of the candidates are selected.

    aggressive_elimination : bool, default=False
        This is only relevant in cases where there isn't enough resources to
        eliminate enough candidates at the last iteration. If ``True``, then
        the search process will 'replay' the first iteration for as long as
        needed until the number of candidates is small enough. This is
        ``False`` by default, which means that the last iteration may evaluate
        more than ``ratio`` candidates.

    force_exhaust_resources : bool, default=False
        If True, then ``min_resources`` is set to a specific value such that
        the last iteration uses as much resources as possible. Namely, the
        last iteration uses the highest value smaller than ``max_resources``
        that is a multiple of both ``min_resources`` and ``ratio``.

    Attributes
    ----------
    n_candidates_ : list of int
        The number of candidate parameters that were evaluated at each
        iteration.

    n_remaining_candidates_ : int
        The number of candidate parameters that are left after the last
        iteration.

    max_resources_ : int
        The maximum number of resources that any candidate is allowed to use
        for a given iteration. Note that since the number of resources used
        at each iteration must be a multiple of ``min_resources_``, the
        actual number of resources used at the last iteration may be smaller
        than ``max_resources_``.

    min_resources_ : int
        The amount of resources that are allocated for each candidate at the
        first iteration.

    n_iterations_ : int
        The actual number of iterations that were run. This is equal to
        ``n_required_iterations_`` if ``aggressive_elimination`` is ``True``.
        Else, this is equal to ``min(n_possible_iterations_,
        n_required_iterations_)``.

    n_possible_iterations_ : int
        The number of iterations that are possible starting with
        ``min_resources_`` resources and without exceeding
        ``max_resources_``.

    n_required_iterations_ : int
        The number of iterations that are required to end up with less than
        ``ratio`` candidates at the last iteration, starting with
        ``min_resources_`` resources. This will be smaller than
        ``n_possible_iterations_`` when there isn't enough resources.

    cv_results_ : dict of numpy (masked) ndarrays
        A dict with keys as column headers and values as columns, that can be
        imported into a pandas ``DataFrame``. It contains many informations for
        analysing the results of a search.
        Please refer to the :ref:`User guide<successive_halving_cv_results>`
        for details.


    best_estimator_ : estimator or dict
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if ``refit=False``.

    best_score_ : float
        Mean cross-validated score of the best_estimator.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.

    best_index_ : int
        The index (of the ``cv_results_`` arrays) which corresponds to the best
        candidate parameter setting.

        The dict at ``search.cv_results_['params'][search.best_index_]`` gives
        the parameter setting for the best model, that gives the highest
        mean score (``search.best_score_``).

    scorer_ : function or a dict
        Scorer function used on the held out data to choose the best
        parameters for the model.

    n_splits_ : int
        The number of cross-validation splits (folds/iterations).

    refit_time_ : float
        Seconds used for refitting the best model on the whole dataset.

        This is present only if ``refit`` is not False.

    Examples
    --------

    >>> from sklearn.datasets import load_iris
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.model_selection import HalvingGridSearchCV
    ...
    >>> X, y = load_iris(return_X_y=True)
    >>> clf = RandomForestClassifier(random_state=0)
    ...
    >>> param_grid = {"max_depth": [3, None],
    ...               "min_samples_split": [5, 10]}
    >>> search = HalvingGridSearchCV(clf, param_grid, resource='n_estimators',
    ...                              max_resources=10,
    ...                              force_exhaust_resources=True,
    ...                              random_state=0).fit(X, y)
    >>> search.best_params_  # doctest: +SKIP
    {'max_depth': None, 'min_samples_split': 10, 'n_estimators': 9}

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
    :class:`HalvingRandomSearchCV`:
        Random search over a set of parameters using successive halving.
    """
    _required_parameters = ["estimator", "param_grid"]

    def __init__(self, estimator, param_grid, scoring=None,
                 n_jobs=None, refit=True, verbose=0, cv=5,
                 pre_dispatch='2*n_jobs', random_state=None,
                 error_score=np.nan, return_train_score=True,
                 max_resources='auto', min_resources='auto',
                 resource='n_samples', ratio=3, aggressive_elimination=False,
                 force_exhaust_resources=False):
        super().__init__(estimator, scoring=scoring,
                         n_jobs=n_jobs, refit=refit, verbose=verbose, cv=cv,
                         pre_dispatch=pre_dispatch,
                         random_state=random_state, error_score=error_score,
                         return_train_score=return_train_score,
                         max_resources=max_resources, resource=resource,
                         ratio=ratio, min_resources=min_resources,
                         aggressive_elimination=aggressive_elimination,
                         force_exhaust_resources=force_exhaust_resources)
        self.param_grid = param_grid
        _check_param_grid(self.param_grid)

    def _generate_candidate_params(self):
        return ParameterGrid(self.param_grid)


class HalvingRandomSearchCV(BaseSuccessiveHalving):
    """Randomized search on hyper parameters.

    The search strategy starts evaluating all the candidates with a small
    amount of resources and iteratively selects the best candidates, using more
    and more resources.

    Read more in the :ref:`User guide<successive_halving_user_guide>`.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``scoring`` must be passed.

    param_distributions : dict
        Dictionary with parameters names (string) as keys and distributions
        or lists of parameters to try. Distributions must provide a ``rvs``
        method for sampling (such as those from scipy.stats.distributions).
        If a list is given, it is sampled uniformly.

    n_candidates: int, default='auto'
        The number of candidate parameters to sample. By default this will
        sample enough candidates so that the last iteration uses as many
        resources as possible. Note that ``force_exhaust_resources`` has no
        effect in this case.

    scoring : string, callable, or None, default=None
        A single string (see :ref:`scoring_parameter`) or a callable
        (see :ref:`scoring`) to evaluate the predictions on the test set.
        If None, the estimator's score method is used.

    n_jobs : int or None, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

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
              as in '2*n_jobs' (default)

    cv : int, cross-validation generator or an iterable, default=5
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - integer, to specify the number of folds in a `(Stratified)KFold`,
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` is used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    refit : boolean, default=True
        If True, refit an estimator using the best found parameters on the
        whole dataset.

        The refitted estimator is made available at the ``best_estimator_``
        attribute and permits using ``predict`` directly on this
        ``GridSearchCV`` instance.

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    error_score : 'raise' or numeric
        Value to assign to the score if an error occurs in estimator fitting.
        If set to 'raise', the error is raised. If a numeric value is given,
        FitFailedWarning is raised. This parameter does not affect the refit
        step, which will always raise the error. Default is ``np.nan``

    return_train_score : boolean, default=False
        If ``False``, the ``cv_results_`` attribute will not include training
        scores.
        Computing training scores is used to get insights on how different
        parameter settings impact the overfitting/underfitting trade-off.
        However computing the scores on the training set can be computationally
        expensive and is not strictly required to select the parameters that
        yield the best generalization performance.

    max_resources : int, default='auto'
        The maximum number of resources that any candidate is allowed to use
        for a given iteration. By default, this is set ``n_samples`` when
        ``resource='n_samples'`` (default), else an error is raised.

    min_resources : int, default='auto'
        The minimum amount of resource that any candidate is allowed to use for
        a given iteration. Equivalently, this defines the amount of resources
        that are allocated for each candidate at the first iteration. By
        default, this is set to:

        - ``n_splits * 2`` when ``resource='n_samples'`` for a regression
          problem
        - ``n_classes * n_splits * 2`` when ``resource='n_samples'`` for a
          regression problem
        - The highest possible value satisfying the constraint
          ``force_exhaust_resources=True``.
        - ``1`` when ``resource!='n_samples'``

        Note that the amount of resources used at each iteration is always a
        multiple of ``min_resources``.

    resource : ``'n_samples'`` or str, default='n_samples'
        Defines the resource that increases with each iteration. By default,
        the resource is the number of samples. It can also be set to any
        parameter of the base estimator that accepts positive integer
        values, e.g. 'n_iterations' or 'n_estimators' for a gradient
        boosting estimator. In this case ``max_resources`` cannot be 'auto'.

    ratio : int or float, default=3
        The 'halving' parameter, which determines the proportion of candidates
        that are selected for the next iteration. For example, ``ratio=3``
        means that only one third of the candidates are selected.

    aggressive_elimination : bool, default=False
        This is only relevant in cases where there isn't enough resources to
        eliminate enough candidates at the last iteration. If ``True``, then
        the search process will 'replay' the first iteration for as long as
        needed until the number of candidates is small enough. This is
        ``False`` by default, which means that the last iteration may evaluate
        more than ``ratio`` candidates.

    force_exhaust_resources : bool, default=False
        If True, then ``min_resources`` is set to a specific value such that
        the last iteration uses as much resousrces as possible. Namely, the
        last iteration uses the highest value smaller than ``max_resources``
        that is a multiple of both ``min_resources`` and ``ratio``.

    Attributes
    ----------
    n_candidates_ : list of int
        The number of candidate parameters that were evaluated at each
        iteration.

    n_remaining_candidates_ : int
        The number of candidate parameters that are left after the last
        iteration.

    max_resources_ : int
        The maximum number of resources that any candidate is allowed to use
        for a given iteration. Note that since the number of resources used at
        each iteration must be a multiple of ``min_resources_``, the actual
        number of resources used at the last iteration may be smaller than
        ``max_resources_``.

    min_resources_ : int
        The amount of resources that are allocated for each candidate at the
        first iteration.

    n_iterations_ : int
        The actual number of iterations that were run. This is equal to
        ``n_required_iterations_`` if ``aggressive_elimination`` is ``True``.
        Else, this is equal to ``min(n_possible_iterations_,
        n_required_iterations_)``.

    n_possible_iterations_ : int
        The number of iterations that are possible starting with
        ``min_resources_`` resources and without exceeding
        ``max_resources_``.

    n_required_iterations_ : int
        The number of iterations that are required to end up with less than
        ``ratio`` candidates at the last iteration, starting with
        ``min_resources_`` resources. This will be smaller than
        ``n_possible_iterations_`` when there isn't enough resources.

    cv_results_ : dict of numpy (masked) ndarrays
        A dict with keys as column headers and values as columns, that can be
        imported into a pandas ``DataFrame``. It contains many informations for
        analysing the results of a search.
        Please refer to the :ref:`User guide<successive_halving_cv_results>`
        for details.

    best_estimator_ : estimator or dict
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if ``refit=False``.

    best_score_ : float
        Mean cross-validated score of the best_estimator.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.

    best_index_ : int
        The index (of the ``cv_results_`` arrays) which corresponds to the best
        candidate parameter setting.

        The dict at ``search.cv_results_['params'][search.best_index_]`` gives
        the parameter setting for the best model, that gives the highest
        mean score (``search.best_score_``).

    scorer_ : function or a dict
        Scorer function used on the held out data to choose the best
        parameters for the model.

    n_splits_ : int
        The number of cross-validation splits (folds/iterations).

    refit_time_ : float
        Seconds used for refitting the best model on the whole dataset.

        This is present only if ``refit`` is not False.

    Examples
    --------

    >>> from sklearn.datasets import load_iris
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.model_selection import HalvingRandomSearchCV
    >>> from scipy.stats import randint
    ...
    >>> X, y = load_iris(return_X_y=True)
    >>> clf = RandomForestClassifier(random_state=0)
    >>> np.random.seed(0)
    ...
    >>> param_distributions = {"max_depth": [3, None],
    ...                        "min_samples_split": randint(2, 11)}
    >>> search = HalvingRandomSearchCV(clf, param_distributions,
    ...                                resource='n_estimators',
    ...                                max_resources=10,
    ...                                force_exhaust_resources=True,
    ...                                random_state=0).fit(X, y)
    >>> search.best_params_  # doctest: +SKIP
    {'max_depth': None, 'min_samples_split': 10, 'n_estimators': 9}

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
    :class:`HalvingGridSearchCV`:
        Search over a grid of parameters using successive halving.
    """
    _required_parameters = ["estimator", "param_distributions"]

    def __init__(self, estimator, param_distributions,
                 n_candidates='auto', scoring=None, n_jobs=None, refit=True,
                 verbose=0, cv=5, pre_dispatch='2*n_jobs',
                 random_state=None, error_score=np.nan,
                 return_train_score=True, max_resources='auto',
                 min_resources='auto', resource='n_samples', ratio=3,
                 aggressive_elimination=False, force_exhaust_resources=False):
        super().__init__(estimator, scoring=scoring,
                         n_jobs=n_jobs, refit=refit, verbose=verbose, cv=cv,
                         random_state=random_state, error_score=error_score,
                         return_train_score=return_train_score,
                         max_resources=max_resources, resource=resource,
                         ratio=ratio, min_resources=min_resources,
                         aggressive_elimination=aggressive_elimination,
                         force_exhaust_resources=force_exhaust_resources)
        self.param_distributions = param_distributions
        self.n_candidates = n_candidates

    def _generate_candidate_params(self):
        n_candidates_first_iter = self.n_candidates
        if n_candidates_first_iter == 'auto':
            # This will generate enough candidate so that the last iteration
            # uses as much resources as possible
            n_candidates_first_iter = (
                self.max_resources_ // self.min_resources_)
        return ParameterSampler(self.param_distributions,
                                n_candidates_first_iter,
                                self.random_state)
