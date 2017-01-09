"""Algorithm computing the optimal value of n_clusters."""

# Authors: Thierry Guillemot <thierry.guillemot.work@gmail.com>
#          Arnaud Fouchet <foucheta@gmail.com>

import warnings
import numpy as np

from abc import ABCMeta, abstractmethod
from functools import partial
from collections import defaultdict

from ..base import clone, BaseEstimator
from ..mixture.base import BaseMixture
from ..externals import six
from ..externals.joblib import Parallel, delayed
from ..model_selection import ParameterGrid
from ..utils import check_array, check_random_state
from ..utils.fixes import rankdata
from ..utils.metaestimators import if_delegate_has_method
from ..utils.testing import ignore_warnings
from ..utils.validation import check_is_fitted


class _NClusterSearchBase(six.with_metaclass(ABCMeta, object)):
    """Base class for cluster searchers.

    This abstract class specifies an interface for all cluster searchers and
    provides some basic common methods.
    """
    def __init__(self, estimator, parameters, n_jobs=1,
                 pre_dispatch='2*n_jobs', verbose=0):
        self.estimator = estimator
        self.parameters = parameters
        self.n_jobs = n_jobs
        self.pre_dispatch = pre_dispatch
        self.verbose = verbose

    def _initialized(self, X, y):
        pass

    def _compute_values(self, n_clusters_list):
        """Define all n_clusters values needed to compute a score."""
        return np.array(n_clusters_list)

    def _parameter_grid(self, X, y):
        """Define the grid used by the cluster searcher."""
        # Define the property name controling the number of clusters
        if isinstance(self.estimator, BaseMixture):
            # Case for clustering estimator
            property_name = 'n_components'
        elif isinstance(self.estimator, BaseEstimator):
            # Case for mixture estimator
            property_name = 'n_clusters'
        else:
            raise ValueError("The estimator has to be from the type ...")

        # Extract number of cluster parameter from the parameter list
        parameters = self.parameters.copy()
        n_clusters_list = parameters.pop(property_name)
        # Avoid multiple computation of sets of same parameters
        self.n_clusters_values, self._index = np.unique(
            self._compute_values(n_clusters_list),
            return_index=True)
        self._index = self._index < len(n_clusters_list)

        # Grid definition
        for n_clusters in self.n_clusters_values:
            for params in ParameterGrid(parameters):
                params[property_name] = n_clusters
                yield params

    def _compute_score(self, estimator, X):
        """Compute the score of the estimator according to self.metric."""
        # XXX  Not good to deal with bad values to modify
        try:
            if self.metric is not None:
                # Add a warning for n_cluster = 1
                score = self.metric(X, estimator.fit_predict(X))
            elif hasattr(estimator, 'score'):
                score = estimator.fit(X).score(X)
            # else:
            #     raise ValueError(
            #         "The class %s does not have a score method. Change the "
            #         "estimator type or give an unsupervised metric to the "
            #         "``metric`` attribute.")
        except ValueError:
            # TODO Add a warning for false values
            # warnings.warn('Put a warning.')
            score = np.nan
        return score

    @abstractmethod
    def _estimator_fit(self, estimator, X, y, parameters):
        """Fit the estimator and return values used to define a score."""
        pass

    @abstractmethod
    def _compute_results(self, X, out):
        """Compute result from the values computed by _estimator_fit."""
        pass

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
        check_is_fitted(self, 'best_estimator_')
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
        check_is_fitted(self, 'best_estimator_')
        return self.best_estimator_.predict_proba(X)

    def fit(self, X, y):
        """Run fit with all sets of parameters.

        Parameters
        ----------

        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.
        """
        self._initialized(X, y)

        base_estimator = clone(self.estimator)

        # Fit the estimator for each parameter and compute output used
        # by _compute_results to define a score
        out = Parallel(
            n_jobs=self.n_jobs, verbose=self.verbose,
            pre_dispatch=self.pre_dispatch
        )(delayed(self._estimator_fit)(clone(base_estimator), X, y, parameters)
          for parameters in self._parameter_grid(X, y))

        self.results_ = self._compute_results(X, out)

        # XXX Change that ??
        ranks = np.asarray(rankdata(-self.results_['score'], method='min'),
                           dtype=np.int32)
        self.best_index_ = np.flatnonzero(ranks == 1)[0]
        self.results_['rank_score'] = ranks

        # Use one np.MaskedArray and mask all the places where the param is not
        # applicable for that candidate. Use defaultdict as each candidate may
        # not contain all the params
        param_results = defaultdict(partial(np.ma.masked_all,
                                            (len(ranks),), dtype=object))
        for cand_i, params in enumerate(self.results_['params']):
            for name, value in params.items():
                # An all masked empty array gets created for the key
                # `"param_%s" % name` at the first occurence of `name`.
                # Setting the value at an index also unmasks that index
                param_results["param_%s" % name][cand_i] = value

        self.results_.update(param_results)
        self.best_params_ = self.results_['params'][self.best_index_]
        self.best_score_ = self.results_['score'][self.best_index_]

        # fit the best estimator
        best_estimator = clone(base_estimator).set_params(
            **self.best_params_)
        if y is not None:
            best_estimator.fit(X, y)
        else:
            best_estimator.fit(X)
        self.best_estimator_ = best_estimator

        return self


# Search for unsupervised metrics
class UnsupervisedMetricSearch(_NClusterSearchBase):
    """Unsupervised metric search on hyper parameters.

    UnsupervisedMetricSearch implements an exhaustive search over specified
    parameter values for a clustering estimator using an unsupervised metric
    to define the best result.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``metric`` must be passed.

    parameters : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    metric : callable or None, defaults to None
        An unsupervised metric function with signature
        ``metric(X, y, **kwargs)``.
        If ``None``, the ``score`` method of the estimator is used.

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

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    Attributes
    ----------
    results_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.
    """
    def __init__(self, estimator, parameters, metric=None, n_jobs=-1,
                 pre_dispatch='2*n_jobs', verbose=0):
        super(UnsupervisedMetricSearch, self).__init__(
            estimator, parameters, n_jobs, pre_dispatch, verbose)
        self.metric = metric

    def _estimator_fit(self, estimator, X, y, parameters):
        """Fit the estimator and return values used to define a score."""
        estimator.set_params(**parameters)
        score = self._compute_score(estimator, X)
        return score, parameters

    def _compute_results(self, X, out):
        """Compute result from the values computed by _estimator_fit."""
        scores, parameters = zip(*out)
        return {
            'score': np.array(scores),
            'params': parameters,
        }


# Scorer for supervised metrics
class StabilitySearch(_NClusterSearchBase):
    """Unsupervised metric search on hyper parameters.

    StabilitySearch implements an exhaustive search over specified
    parameter values for a clustering estimator using a stability search
    process to define the best result.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``metric`` must be passed.

    parameters : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    metric : callable
        An unsupervised metric function with signature
        ``metric(X, y, **kwargs)``.
        If ``None``, the ``score`` method of the estimator is used.

    n_draws : int, defaults to 10
        Number of draws used to define the score.

    p_samples: float, defaults to .8
        Probability to draw point to define the sample sets.

    n_jobs : int, defaults to 1
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

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    Attributes
    ----------
    results_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.
    """
    def __init__(self, estimator, parameters, metric=None, n_draws=10,
                 p_samples=.8, random_state=None, n_jobs=-1,
                 pre_dispatch='2*n_jobs', verbose=0):
        super(StabilitySearch, self).__init__(
            estimator, parameters, n_jobs, pre_dispatch, verbose)
        self.metric = metric
        self.n_draws = n_draws
        self.p_samples = p_samples
        self.random_state = random_state

    def _initialized(self, X, y):
        """Define the subsampled sets used by _estimator_fit."""
        n_samples, _ = X.shape
        rng = check_random_state(self.random_state)
        self.data_ = (
            rng.uniform(size=(self.n_draws, 2 * n_samples)) < self.p_samples)

    def _estimator_fit(self, estimator, X, y, parameters):
        """Fit the estimator and return values used to define a score."""
        estimator.set_params(**parameters)
        # Be sure a warining is the best thing to do ?
        if parameters['n_clusters'] == 1:
            warnings.warn('StabilitySearch can not be used for n_clusers = 1.')
            return np.nan, parameters

        draw_scores = np.empty(self.n_draws)
        for l, d in enumerate(self.data_):
            p1, p2 = np.split(d, 2)

            labels1 = estimator.fit_predict(X[p1])
            labels2 = estimator.fit_predict(X[p2])

            draw_scores[l] = self.metric(labels1[p2[p1]], labels2[p1[p2]])
        score = draw_scores.mean()
        return score, parameters

    def _compute_results(self, X, out):
        """Compute result from the values computed by _estimator_fit."""
        scores, parameters = zip(*out)
        return {
            'score': np.array(scores),
            'params': parameters,
        }


# Scorer for distortion jump
class DistortionJumpSearch(_NClusterSearchBase):
    """Search on hyper parameters using the distortion jump scheme.

    DistortionJumpSearch implements an exhaustive search over specified
    parameter values for a clustering estimator using a distortion jump
    process to define the best result.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``metric`` must be passed.

    parameters : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    metric : callable
        A supervised metric function with signature
        ``metric(y_true, y_pred, **kwargs)``.

    n_draws : int, defaults to 10
        Number of draws used to define the score.

    p_samples: float, defaults to .8
        Probability to draw point to define the sample sets.

    n_jobs : int, defaults to 1
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

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    Attributes
    ----------
    results_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.
    """
    def __init__(self, estimator, parameters, metric=None, n_jobs=-1,
                 pre_dispatch='2*n_jobs', verbose=0):
        super(DistortionJumpSearch, self).__init__(
            estimator, parameters, n_jobs, pre_dispatch, verbose)
        self.metric = metric

    def _compute_values(self, n_clusters_list):
        """Define all n_clusters values needed to compute a score."""
        return np.hstack((n_clusters_list, np.array(n_clusters_list) + 1))

    def _estimator_fit(self, estimator, X, y, parameters):
        """Fit the estimator and return values used to define a score."""
        _, n_features = X.shape
        estimator.set_params(**parameters)
        score = self._compute_score(estimator, X)
        distortion = (score / n_features) ** (-n_features / 2)
        return distortion, parameters

    def _compute_results(self, X, out):
        """Compute result from the values computed by _estimator_fit."""
        distortion, parameters = zip(*out)
        distortion = np.array(distortion).reshape(
            len(self.n_clusters_values), -1)
        parameters = np.array(parameters).reshape(distortion.shape)

        return {
            'score': (distortion[self._index] -
                      distortion[1:][self._index[:-1]]),
            'params': parameters[self._index].ravel(),
        }


# Scorer for Pham
class PhamSearch(_NClusterSearchBase):
    """Unsupervised metric search on hyper parameters.

    StabilitySearch implements an exhaustive search over specified
    parameter values for a clustering estimator using a stability search
    process to define the best result.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``metric`` must be passed.

    parameters : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    metric : callable
        A supervised metric function with signature
        ``metric(y_true, y_pred, **kwargs)``.

    n_draws : int, defaults to 10
        Number of draws used to define the score.

    p_samples: float, defaults to .8
        Probability to draw point to define the sample sets.

    n_jobs : int, defaults to 1
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

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    Attributes
    ----------
    results_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.
    """
    def __init__(self, estimator, parameters, metric=None):
        super(PhamSearch, self).__init__(
            estimator, parameters)
        self.metric = metric

    def _compute_values(self, n_clusters_list):
        """Define all n_clusters values needed to compute a score."""
        return np.hstack((n_clusters_list, np.array(n_clusters_list) - 1))

    def _estimator_fit(self, estimator, X, y, parameters):
        """Fit the estimator and return values used to define a score."""
        estimator.set_params(**parameters)
        if parameters['n_clusters'] < 1:
            return np.nan, parameters
        # XXX I clearly don't like to use any unsupervised metric here.
        score = self._compute_score(estimator, X)
        return score, parameters

    def _compute_results(self, X, out):
        """Compute result from the values computed by _estimator_fit."""
        _, n_features = X.shape

        scores, parameters = zip(*out)
        scores = np.array(scores).reshape(len(self.n_clusters_values), -1)
        parameters = np.array(parameters).reshape(scores.shape)

        weights = 1. - np.exp((self.n_clusters_values[self._index] - 2) *
                              np.log(5. / 6.) + np.log(.75) -
                              np.log(n_features))

        with ignore_warnings(category=RuntimeWarning):
            scores = np.where(np.logical_and(
                np.vstack((
                    True,
                    ~np.isnan(scores[self._index])[:-1])),
                scores[self._index] > 0), scores[self._index] /
                (weights[:, np.newaxis] * scores[:-1][self._index[1:]]), 1.)

        if self.n_clusters_values[0] == 0:
            scores[0, :] = 1.

        # XXX change that to not modify the score
        with ignore_warnings(category=RuntimeWarning):
            scores[scores > 0.85] = 1.

        return {
            'score': -scores.ravel(),
            'params': parameters[self._index].ravel(),
        }


# Scorer for Gap
class GapSearch(_NClusterSearchBase):
    """Unsupervised metric search on hyper parameters.

    StabilitySearch implements an exhaustive search over specified
    parameter values for a clustering estimator using a stability search
    process to define the best result.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``metric`` must be passed.

    parameters : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    metric : callable
        A supervised metric function with signature
        ``metric(y_true, y_pred, **kwargs)``.

    n_draws : int, defaults to 10
        Number of draws used to define the score.

    p_samples: float, defaults to .8
        Probability to draw point to define the sample sets.

    n_jobs : int, defaults to 1
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

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    Attributes
    ----------
    results_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.
    """
    def __init__(self, estimator, parameters, metric=None, n_draws=10,
                 n_jobs=1, random_state=None):
        super(GapSearch, self).__init__(estimator, parameters, n_jobs)
        self.metric = metric
        self.n_draws = n_draws
        self.random_state = random_state

    def _initialized(self, X, y):
        """Define the subsampled sets used by _estimator_fit."""
        n_samples, n_features = X.shape
        rng = check_random_state(self.random_state)

        # Compute the random data once for all
        bb_min, bb_max = np.min(X, 0), np.max(X, 0)
        self._data = (rng.uniform(size=(self.n_draws, n_samples, n_features)) *
                      (bb_max - bb_min) + bb_min)

    def _compute_values(self, n_clusters_list):
        """Define all n_clusters values needed to compute a score."""
        return np.hstack((n_clusters_list, np.array(n_clusters_list) + 1))

    # def _compute_score(self, estimator, X):
    #     """Compute result from the values computed by _estimator_fit."""
    #     if self.metric is not None:
    #         # Add a warning for n_cluster = 1
    #         if estimator.n_clusters <= 1:
    #             return np.nan
    #         score = self.metric(X, estimator.fit_predict(X))
    #     elif hasattr(estimator, 'score'):
    #         score = estimator.fit(X).score(X)
    #     else:
    #         raise ValueError("The estimator does not have score_ attribute.")
    #
    #     return score

    def _estimator_fit(self, estimator, X, y, parameters):
        """Fit the estimator and return values used to define a score."""
        estimator.set_params(**parameters)

        inertia_n_draws = np.empty(self.n_draws)
        with ignore_warnings(category=RuntimeWarning):
            estimated_log_inertia = np.log(
                np.abs(self._compute_score(estimator, X)))
            for t, Xt in enumerate(self._data):
                inertia_n_draws[t] = np.abs(self._compute_score(estimator, Xt))
        inertia_n_draws = np.log(inertia_n_draws)

        expected_log_inertia = np.mean(inertia_n_draws)
        gap = expected_log_inertia - estimated_log_inertia
        safety = (np.std(inertia_n_draws) *
                  np.sqrt(1. + 1. / self.n_draws))

        return gap, safety, parameters

    def _compute_results(self, X, out):
        """Compute result from the values computed by _estimator_fit."""
        gap, safety, parameters = zip(*out)

        gap = np.array(gap).reshape(len(self.n_clusters_values), -1)
        safety = np.array(safety).reshape(gap.shape)

        with ignore_warnings(category=RuntimeWarning):
            scores = ((gap[self._index] - gap[1:][self._index[:-1]] +
                       safety[1:][self._index[:-1]]) >= 0).astype(gap.dtype)
        return {
            'gap': gap[self._index].ravel(),
            'safety': safety[self._index].ravel(),
            'score': scores.ravel(),
            'params': (
                np.array(parameters).reshape(gap.shape)[self._index].ravel()),
        }


class OptimalNClusterSearch():
    """Unsupervised metric search on hyper parameters.

    StabilitySearch implements an exhaustive search over specified
    parameter values for a clustering estimator using a stability search
    process to define the best result.

    Parameters
    ----------
    estimator : estimator object.
        This is assumed to implement the scikit-learn estimator interface.
        Either estimator needs to provide a ``score`` function,
        or ``metric`` must be passed.

    parameters : dict or list of dictionaries
        Dictionary with parameters names (string) as keys and lists of
        parameter settings to try as values, or a list of such
        dictionaries, in which case the grids spanned by each dictionary
        in the list are explored. This enables searching over any sequence
        of parameter settings.

    metric : callable
        A supervised metric function with signature
        ``metric(y_true, y_pred, **kwargs)``.

    n_draws : int, defaults to 10
        Number of draws used to define the score.

    p_samples: float, defaults to .8
        Probability to draw point to define the sample sets.

    n_jobs : int, defaults to 1
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

    verbose : integer
        Controls the verbosity: the higher, the more messages.

    Attributes
    ----------
    results_ : list of named tuples
        Contains scores for all parameter combinations in param_grid.
        Each entry corresponds to one parameter setting.

    best_estimator_ : estimator
        Estimator that was chosen by the search, i.e. estimator
        which gave highest score (or smallest loss if specified)
        on the left out data. Not available if refit=False.

    best_score_ : float
        Score of best_estimator on the left out data.

    best_params_ : dict
        Parameter setting that gave the best results on the hold out data.
    """
    def __init__(self, estimator, parameters, fitting_process='auto',
                 **kwargs):
        self.estimator = estimator
        self.parameters = parameters.copy()
        self.fitting_process = fitting_process
        self.kwargs = kwargs

    def _check_initial_parameters(self, X):
        n_samples, _ = X.shape

        X = check_array(X, dtype=[np.float64, np.float32])
        # if self.n_clusters_range.max() > n_samples:
        #     raise ValueError("Put a message")

    def _select_fitting_method(self):
        fitting_method = {
            'auto': PhamSearch,
            'pham': PhamSearch,
            'gap': GapSearch,
            'distortion_jump': DistortionJumpSearch,
            'stability': StabilitySearch,
            'unsupervised': UnsupervisedMetricSearch
        }

        self.scorer_ = fitting_method[self.fitting_process](
            estimator=self.estimator, parameters=self.parameters,
            **self.kwargs)

    def fit(self, X, y=None):
        """Run fit with all sets of parameters.

        Parameters
        ----------

        X : array-like, shape = [n_samples, n_features]
            Training vector, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples] or [n_samples, n_output], optional
            Target relative to X for classification or regression;
            None for unsupervised learning.
        """
        self._check_initial_parameters(X)

        self._select_fitting_method()
        self.scorer_.fit(X, y)

        self.results_ = self.scorer_.results_
        self.best_estimator_ = self.scorer_.best_estimator_

        return self.scorer_
