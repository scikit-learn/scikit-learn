"""
The :mod:`sklearn.model_selection.subselect` includes utilities for generating custom
refit callable objects for subselecting models from a `GridSearchCV` or
`RandomSearchCV` instances.
"""
import warnings
from functools import partial
from typing import Callable, Tuple, Dict, Optional

import numpy as np


__all__ = [
    "Refitter",
    "by_standard_error",
    "by_percentile_rank",
    "by_signed_rank",
    "by_fixed_window",
    "constrain",
    "_wrap_refit",
]


class by_standard_error:
    """A callable class that returns a window of model performance that falls within
    `sigma` standard errors of the highest performing model.

    Attributes
    ----------
    sigma : int
        Number of standard errors tolerance in the case that a standard error
        threshold is used to filter outlying scores across folds. Required if
        `rule`=='se'. Default is 1.

    """

    def __init__(self, sigma: int = 1):
        self.sigma = sigma

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        n_folds: int,
    ) -> Tuple[float, float]:
        """
        Returns a window of model performance whereby differences in model performance
        are within a specified number standard errors (`sigma`) of the best performing
        model.
        """
        # Estimate the standard error across folds for each column of the grid
        cv_se = np.array(np.nanstd(score_grid, axis=1) / np.sqrt(n_folds))

        # Determine confidence interval
        h_cutoff = cv_means[best_score_idx] + self.sigma * cv_se[best_score_idx]
        l_cutoff = cv_means[best_score_idx] - self.sigma * cv_se[best_score_idx]

        return l_cutoff, h_cutoff


class by_percentile_rank:
    """A callable class that returns a window of model performance that falls within
    `eta` percentile of the highest performing model.

    Attributes
    ----------
    eta : float
        Percentile tolerance in the case that a percentile threshold
        is used to filter outlier scores across folds. Required if
        `rule`=='percentile'. Default is 0.68.

    """

    def __init__(self, eta: float = 0.68):
        self.eta = eta

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        n_folds: int,
    ) -> Tuple[float, float]:
        """
        Returns a window of model performance whereby differences in model performance
        are within a specified percentile (`eta`) of the best performing model.
        """
        # Estimate the indicated percentile, and its inverse, across folds for
        # each column of the grid
        perc_cutoff = np.nanpercentile(
            score_grid, [100 * self.eta, 100 - 100 * self.eta], axis=1
        )

        # Determine bounds of the percentile interval
        h_cutoff = perc_cutoff[0, best_score_idx]
        l_cutoff = perc_cutoff[1, best_score_idx]

        return l_cutoff, h_cutoff


class by_signed_rank:
    """A callable class that returns a window of model performance that is not
    significantly different from the highest performing model, based on a signed
    Wilcoxon rank sum test.

    Attributes
    ----------
    alpha : float
        An alpha significance level in the case that wilcoxon rank sum
        hypothesis testing is used to filter outlying scores across folds.
        Required if `rule`=='ranksum'. Default is 0.05.
    alternative : str
        The alternative hypothesis to test against. Must be one of 'two-sided',
        'less', or 'greater'. Default is 'two-sided'. See `scipy.stats.wilcoxon` for
        more details.
    zero_method : str
        The method used to handle zero scores. Must be one of 'pratt', 'wilcox',
        'zsplit'. Default is 'zsplit'. See `scipy.stats.wilcoxon` for more details.

    """

    def __init__(
        self,
        alpha: float = 0.01,
        alternative: str = "two-sided",
        zero_method: str = "zsplit",
    ):
        self.alpha = alpha
        self.alternative = alternative
        self.zero_method = zero_method

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        n_folds: int,
    ) -> Tuple[float, float]:
        """
        Returns a window of model performance whereby the non-parametric paired
        difference in model performance is insignificant, based on a given `alpha`
        level of significance.
        """

        import itertools
        from scipy.stats import wilcoxon

        # Perform signed Wilcoxon rank sum test for each pair combination of
        # columns against the best average score column
        tests = [
            pair
            for pair in list(itertools.combinations(range(score_grid.shape[0]), 2))
            if best_score_idx in pair
        ]

        pvals = {}
        for pair in tests:
            pvals[pair] = wilcoxon(
                score_grid[pair[0]],
                score_grid[pair[1]],
                alternative=self.alternative,
                zero_method=self.zero_method,
            )[1]

        # Return the models that are insignificantly different from the best average
        # performing, else just return the best-performing model.
        surviving_ranks = [pair[0] for pair in tests if pvals[pair] > self.alpha] + [
            best_score_idx
        ]

        if len(surviving_ranks) == 0:
            surviving_ranks = [best_score_idx]
            warnings.warn(
                "The average performance of all cross-validated models is "
                "significantly different from that of the best-performing model.",
                UserWarning,
            )

        h_cutoff = np.nanmax(cv_means[surviving_ranks])
        l_cutoff = np.nanmin(cv_means[surviving_ranks])

        return l_cutoff, h_cutoff


class by_fixed_window:
    """A callable class that returns a fixed window of model performance based on
    user-specified bounds.
    """

    def __init__(
        self, lower_bound: Optional[float] = None, upper_bound: Optional[float] = None
    ):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        n_folds: int,
    ) -> Tuple[float, float]:
        l_cutoff, h_cutoff = self.lower_bound, self.upper_bound
        if l_cutoff is None:
            l_cutoff = np.nanmin(cv_means)
        if h_cutoff is None:
            h_cutoff = np.nanmax(cv_means)

        return l_cutoff, h_cutoff


class Refitter:
    """`Refitter` is an interface that generates any of a variety of tuning refit
    callables for use with `GridSearchCV` and `RandomSearchCV`. A refit callable can be
    used to select a best-performing model from among a constrained subset of models.
    This can be useful for instance the case that the user wishes to identify
    best-performing models above or below a particular threshold. It can also be useful
    for selecting alternative models whose performance is not meaningfully different
    from the best-performing model, but whose simplicity may be more preferable (e.g.
    to prevent overfitting).

    Attributes
    ----------
    cv_results_ : dict of numpy (masked) ndarrays
        A dict with keys as column headers and values as columns, that can be
        imported into a pandas ``DataFrame``.

        For instance the below given table

        +------------+-----------+------------+-----------------+---+---------+
        |param_kernel|param_gamma|param_degree|split0_test_score|...|rank_t...|
        +============+===========+============+=================+===+=========+
        |  'poly'    |     --    |      2     |       0.80      |...|    2    |
        +------------+-----------+------------+-----------------+---+---------+
        |  'poly'    |     --    |      3     |       0.70      |...|    4    |
        +------------+-----------+------------+-----------------+---+---------+
        |  'rbf'     |     0.1   |     --     |       0.80      |...|    3    |
        +------------+-----------+------------+-----------------+---+---------+
        |  'rbf'     |     0.2   |     --     |       0.93      |...|    1    |
        +------------+-----------+------------+-----------------+---+---------+

        will be represented by a ``cv_results_`` dict of::

            {
            'param_kernel': masked_array(data = ['poly', 'poly', 'rbf', 'rbf'],
                                         mask = [False False False False]...)
            'param_gamma': masked_array(data = [-- -- 0.1 0.2],
                                        mask = [ True  True False False]...),
            'param_degree': masked_array(data = [2.0 3.0 -- --],
                                         mask = [False False  True  True]...),
            'split0_test_score'  : [0.80, 0.70, 0.80, 0.93],
            'split1_test_score'  : [0.82, 0.50, 0.70, 0.78],
            'mean_test_score'    : [0.81, 0.60, 0.75, 0.85],
            'std_test_score'     : [0.01, 0.10, 0.05, 0.08],
            'rank_test_score'    : [2, 4, 3, 1],
            'split0_train_score' : [0.80, 0.92, 0.70, 0.93],
            'split1_train_score' : [0.82, 0.55, 0.70, 0.87],
            'mean_train_score'   : [0.81, 0.74, 0.70, 0.90],
            'std_train_score'    : [0.01, 0.19, 0.00, 0.03],
            'mean_fit_time'      : [0.73, 0.63, 0.43, 0.49],
            'std_fit_time'       : [0.01, 0.02, 0.01, 0.01],
            'mean_score_time'    : [0.01, 0.06, 0.04, 0.04],
            'std_score_time'     : [0.00, 0.00, 0.00, 0.01],
            'params'             : [{'kernel': 'poly', 'degree': 2}, ...],
            }
    scoring : str
        The scoring metric of interest. Must be one of the scoring metrics
        specified in `GridSearchCV` or `RandomSearchCV`.

        NOTE

        For multi-metric evaluation, the scores for all the scorers are
        available in the ``cv_results_`` dict at the keys ending with that
        scorer's name (``'_<scorer_name>'``) instead of ``'_score'`` shown
        above. ('split0_test_precision', 'mean_train_precision' etc.)

        Compatible scoring metrics must be monotonically increasing such that the best
        score is always the one with the highest value. For error metric scorers, this
        can typically be addressed by using the negated variant of the scoring metric
        of interest, i.e. ``scoring='neg_mean_squared_error'``. If this is not the case
        or a custom scorer is used, it must be wrapped with ``make_scorer`` to ensure
        that the ``greater_is_better`` parameter is set to ``True``.

    References
    ----------
    Breiman, Friedman, Olshen, and Stone. (1984) Classification and Regression
    Trees. Wadsworth.

    Notes
    -----
    A model's simplicity or complexity is here defined with respect to the
    influenced by one or more hyperparameters of interest (e.g. number of components,
    number of estimators, polynomial degree, cost, scale, number hidden units, weight
    decay, number of nearest neighbors, L1/L2 penalty, etc.).

    The `Refitter` callable API assumes that the `params` attribute of
    `cv_results_` 1) contains the indicated hyperparameter (`param`) of
    interest, and 2) contains a sequence of values (numeric, boolean, or
    categorical) that are ordered from least to most complex as defined by the user.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.svm import SVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import by_standard_error
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
            ("reduce_dim", PCA(random_state=42)),
            ("classify", SVC(random_state=42, C=0.01))
        ])
    >>> param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}
    >>> search = GridSearchCV(pipe, param_grid=param_grid, scoring="accuracy")
    >>> search.fit(X, y)
    >>> ss = Refitter(search.cv_results_)
    >>> ss.fit(by_standard_error(sigma=1))
    >>> refitted_index = ss.transform("reduce_dim__n_components")
    >>> refitted_index

    """

    def __init__(self, cv_results_: Dict, scoring: str = "score"):
        self.cv_results_ = cv_results_
        self.scoring = scoring

    def _get_splits(self):
        # Extract subgrid corresponding to the scoring metric of interest
        fitted_key_strings = "\t".join(list(self.cv_results_.keys()))
        if not all(s in fitted_key_strings for s in ["split", "params", "mean_test"]):
            raise TypeError(
                "cv_results_ must be a dict of fitted GridSearchCV or RandomSearchCV"
                " objects."
            )

        _splits = [
            i
            for i in list(self.cv_results_.keys())
            if f"test_{self.scoring}" in i and i.startswith("split")
        ]
        if len(_splits) == 0:
            raise KeyError(f"Scoring metric {self.scoring} not found in cv grid.")
        else:
            return _splits

    @property
    def _n_folds(self):
        # Extract number of folds from cv_results_. Note that we cannot get this from
        # the `n_splits_` attribute of the `cv` object because it is not exposed to the
        # refit callable.
        return len(
            list(
                set(
                    [
                        i.split("_")[0]
                        for i in list(self.cv_results_.keys())
                        if i.startswith("split")
                    ]
                )
            )
        )

    @property
    def _score_grid(self):
        # Extract subgrid corresponding to the scoring metric of interest
        return np.vstack([self.cv_results_[cv] for cv in self._get_splits()]).T

    @property
    def _cv_means(self):
        # Calculate means of subgrid corresponding to the scoring metric of interest
        return np.array(np.nanmean(self._score_grid, axis=1))

    @property
    def _lowest_score_idx(self):
        # Return index of the lowest performing model
        return np.nanargmin(self._cv_means)

    @property
    def _best_score_idx(self):
        # Return index of the highest performing model
        return np.nanargmax(self._cv_means)

    def _apply_thresh(self, param: str, l_cutoff: float, h_cutoff: float) -> int:
        if not any(param in x for x in self.cv_results_["params"][0].keys()):
            raise KeyError(f"Parameter {param} not found in cv grid.")
        else:
            hyperparam = [
                i for i in self.cv_results_["params"][0].keys() if i.endswith(param)
            ][0]

        # Mask out all grid columns that are outside the performance window
        self.cv_results_[f"param_{hyperparam}"].mask = np.where(
            (self._cv_means >= float(l_cutoff)) & (self._cv_means <= float(h_cutoff)),
            True,
            False,
        )

        if np.sum(self.cv_results_[f"param_{hyperparam}"].mask) == 0:
            print(
                f"\nLow: {l_cutoff}\nHigh: {h_cutoff}\nMeans across folds:"
                f" {self._cv_means}\nhyperparam: {hyperparam}\n"
            )
            raise ValueError(
                "No valid grid columns remain within the boundaries of the specified"
                " performance window."
            )

        # Find the highest surviving rank (i.e. the lowest performing model among those
        # remaining within the performance window)
        highest_surviving_rank = np.nanmax(
            self.cv_results_["rank_test_score"][
                self.cv_results_[f"param_{hyperparam}"].mask
            ]
        )

        # Return the index of the highest surviving rank. This is the index of the
        # simplest model that is not meaningfully different from the best-performing
        # model.
        return int(
            np.flatnonzero(
                np.isin(self.cv_results_["rank_test_score"], highest_surviving_rank)
            )[0]
        )

    def fit(self, selector: Callable) -> Tuple[float, float]:
        """Fit the refitter using the specified selector callable

        Parameters
        ----------
        selector : callable
            A callable that consumes GridSearchCV or RandomSearchCV results and
            returns a tuple of floats representing the lower and upper bounds of a
            target model performance window.

        Returns
        -------
        l_cutoff : float
            The lower bound of the target model performance window.
        h_cutoff : float
            The upper bound of the target model performance window.

        """
        if not callable(selector):
            raise TypeError(
                "selector must be a callable initialized with `score_grid`, `cv_means`,"
                " `best_score_idx`, and `n_folds` arguments."
            )

        fit_params = {
            "score_grid": self._score_grid,
            "cv_means": self._cv_means,
            "best_score_idx": self._best_score_idx,
            "n_folds": self._n_folds,
        }

        self.l_cutoff, self.h_cutoff = selector(**fit_params)
        return self.l_cutoff, self.h_cutoff

    def transform(self, param: str) -> int:
        """Re-evaluate the best-performing model under the fitted constraints with
        respect to a specified hyperparameter of interest.

        Parameters
        ----------
        param : str
            The name of the hyperparameter of interest. Note that this must be the
            name of the hyperparameter as it appears in the `cv_results_` dictionary
            of the `GridSearchCV` or `RandomSearchCV` object.

        Returns
        -------
        int
            The index of the best-performing model under the fitted constraints with
            respect to the specified hyperparameter of interest.

        Raises
        ------
        `ValueError`
            If the refitter has not been fitted before calling the `transform` method.

        """
        if not hasattr(self, "l_cutoff") or not hasattr(self, "h_cutoff"):
            raise ValueError(
                "Refitter must be fitted before calling `transform` method."
            )

        best_index_ = self._apply_thresh(param, self.l_cutoff, self.h_cutoff)
        print(f"Original best index: {self._best_score_idx}")
        print(f"Refitted best index: {best_index_}")
        print(f"Refitted best params: {self.cv_results_['params'][best_index_]}")
        print(
            f"Refitted best score: {self.cv_results_['mean_test_score'][best_index_]}"
        )
        return self._apply_thresh(param, self.l_cutoff, self.h_cutoff)


def _wrap_refit(
    cv_results_: Dict, param: str, selector: Callable, scoring: str = "score"
) -> int:
    """A wrapper for `Refitter` that returns the best index for a given parameter
    under the constraints imposed by a model selection callable. This function is used
    by `constrain` and should not be called directly.

    Parameters
    ----------
    cv_results_ : Dict
        The `cv_results_` attribute of a `GridSearchCV` or `RandomSearchCV` object.
    param : str
        A complexity-influential hyperparameter to be constrained.
    selector : Callable
        Function that returns the lower and upper bounds of an acceptable performance
        window.
    scoring : str
        The scoring metric used to select the best model. Default is "score", and only
        needs to be changed in the case of multiple scoring metrics or a custom scoring
        function.

    Returns
    -------
    int
        The index of the best model under the performance constraints imposed by the
        selector strategy.

    """
    ss = Refitter(cv_results_, scoring=scoring)
    [l_cutoff, h_cutoff] = ss.fit(selector)
    print(f"Low: {l_cutoff}\nHigh: {h_cutoff}\n")
    return ss.transform(param)


def constrain(param: str, selector: Callable) -> Callable:
    """
    Callable API for the `Refitter` class to be run as `refit` argument of
    `GridsearchCV` or `RandomSearchCV`.

    Parameters
    ----------
    param : str
        Parameter with the largest influence on model complexity.
    selector : callable
        Function that returns the lower and upper bounds of an acceptable performance
        window.

    Returns
    -------
    Callable

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.svm import SVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import constrain, by_standard_error
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
            ("reduce_dim", PCA(random_state=42)),
            ("classify", SVC(random_state=42, C=0.01))
        ])
    >>> param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}
    >>> search = GridSearchCV(pipe, param_grid=param_grid, scoring="accuracy",
    refit=constrain("reduce_dim__n_components", by_standard_error(sigma=1)))
    >>> search.fit(X, y)
    >>> search.best_params_

    """

    # avoid returning a closure in a return statement to avoid pickling issues
    best_index_callable = partial(_wrap_refit, param=param, selector=selector)
    return best_index_callable
