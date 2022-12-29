"""
The :mod:`sklearn.model_selection.refit` includes utilities for refitting a
`GridSearchCV` or `RandomSearchCV` instance based on custom strategies.
"""

from functools import partial

import numpy as np

from typing import Callable, Tuple, Dict, Optional
from ..base import BaseEstimator, TransformerMixin

__all__ = [
    "Refitter",
    "by_standard_error",
    "by_percentile_rank",
    "by_signed_rank",
    "constrain",
    "_wrap_election",
]


class by_standard_error:
    """A callable class that returns a window of model performance that falls within
    `sigma` standard errors of the highest performing model.
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
        # Estimate the standard error across folds for each column of the grid
        cv_se = np.array(np.nanstd(score_grid, axis=1) / np.sqrt(n_folds))

        # Determine confidence interval
        h_cutoff = cv_means[best_score_idx] + self.sigma * cv_se[best_score_idx]
        l_cutoff = cv_means[best_score_idx] - self.sigma * cv_se[best_score_idx]

        return l_cutoff, h_cutoff


class by_percentile_rank:
    """A callable class that returns a window of model performance that falls within
    `eta` percentile of the highest performing model.
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
            print(
                "The average performance of all cross-validated models is "
                "significantly different from that of the best-performing model."
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


class Refitter(BaseEstimator, TransformerMixin):
    def __init__(self, cv_results_: Dict, scoring: str = "score"):
        self.cv_results_ = cv_results_
        self.scoring = scoring

    def _get_splits(self):
        # Extract subgrid corresponding to the scoring metric of interest
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
                " razor"
            )

        # Find the highest surviving rank (i.e. the lowest performing model among those
        # remaining within the razor window)
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
        fit_params = {
            "score_grid": self._score_grid,
            "cv_means": self._cv_means,
            "best_score_idx": self._best_score_idx,
            "n_folds": self._n_folds,
        }

        self.l_cutoff, self.h_cutoff = selector(**fit_params)
        return self.l_cutoff, self.h_cutoff

    def transform(self, param: str) -> int:
        best_index_ = self._apply_thresh(param, self.l_cutoff, self.h_cutoff)
        print(f"Original best index: {self._best_score_idx}")
        print(f"Refitted best index: {best_index_}")
        print(f"Refitted best params: {self.cv_results_['params'][best_index_]}")
        print(
            f"Refitted best score: {self.cv_results_['mean_test_score'][best_index_]}"
        )
        return self._apply_thresh(param, self.l_cutoff, self.h_cutoff)


def _wrap_election(cv_results_: Dict, param: str, selector: Callable) -> int:
    """A wrapper for `Refitter` that returns the best index for a given parameter
    under the constraints imposed by a model selection callable.

    Parameters
    ----------
    cv_results_ : Dict
        The `cv_results_` attribute of a `GridSearchCV` or `RandomSearchCV` object.
    param : str
        A complexity-influential hyperparameter to be penalized by the razor.
    selector : Callable
        Function that returns the lower and upper bounds of the razor.

    Returns
    -------
    int
        The index of the best model under the constraints imposed by the razor.

    """
    ss = Refitter(cv_results_)
    [l_cutoff, h_cutoff] = ss.fit(selector)
    print(f"Low: {l_cutoff}\nHigh: {h_cutoff}\n")
    return ss.transform(param)


def constrain(param: str, selector: Callable) -> Callable:
    """
    Callable to be run as `refit` argument of `GridsearchCV` or `RandomSearchCV`.

    Parameters
    ----------
    param : str
        Parameter with the largest influence on model complexity.
    selector : callable
        Function that returns the lower and upper bounds of the razor.

    """
    return partial(_wrap_election, param=param, selector=selector)
