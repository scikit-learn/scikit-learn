"""
The :mod:``sklearn.model_selection.subselect`` includes refit callable factories for
subselecting models from ``GridSearchCV`` or ``RandomizedSearchCV``
"""
import warnings
from functools import partial
from typing import Callable, Tuple, Dict, Optional, Union, List

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
    """Slices a window of model performance based on standard error.

    Standard error is estimated based on a user-supplied number of standard errors,
    sigma. The resulting window of model performance represents the range of scores
    that fall within the indicated margin of error of the best performing model.

    Attributes
    ----------
    sigma : int
        Number of standard errors tolerance in the case that a standard error
        threshold is used to filter outlying scores across folds. Default is 1.

    Raises
    ------
    ValueError
        If sigma is not a positive integer.
    """

    def __init__(self, sigma: int = 1):  # noqa: D107
        self.sigma = sigma
        if not isinstance(self.sigma, int) or self.sigma < 1:
            raise ValueError("sigma must be a positive integer.")

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        lowest_score_idx: int,
        n_folds: int,
    ) -> Tuple[float, float]:
        """Returns a window of model performance based on standard error.

        Parameters
        ----------
        score_grid : np.ndarray
            A 2D array of model performance scores across folds and hyperparameter
            settings.
        cv_means : np.ndarray
            A 1D array of the average model performance across folds for each
            hyperparameter setting.
        best_score_idx : int
            The index of the highest performing hyperparameter setting.
        lowest_score_idx : int
            The index of the lowest performing hyperparameter setting.
        n_folds : int
            The number of folds used in the cross-validation.

        Returns
        -------
        min_cut : float
            The lower bound of the window of model performance.
        max_cut : float
            The upper bound of the window of model performance.
        """
        # Estimate the standard error across folds for each column of the grid
        cv_se = np.array(np.nanstd(score_grid, axis=1) / np.sqrt(n_folds))

        # Determine confidence interval
        max_cut = cv_means[best_score_idx] + self.sigma * cv_se[best_score_idx]
        min_cut = cv_means[best_score_idx] - self.sigma * cv_se[best_score_idx]
        return min_cut, max_cut

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(sigma={self.sigma})"


class by_percentile_rank:
    """Slices a window of model performance based on percentile rank.

    Percentile rank is estimated based on a user-supplied percentile threshold, eta.
    The resulting window of model performance represents the range of scores that fall
    within the indicated percentile range of the best performing model.

    Attributes
    ----------
    eta : float
        Percentile tolerance in the case that a percentile threshold is used to filter
        outlier scores across folds. Default is 0.68.

    Raises
    ------
    ValueError
        If eta is not a float between 0 and 1.
    """

    def __init__(self, eta: float = 0.68):  # noqa: D107
        self.eta = eta
        if not isinstance(self.eta, float) or self.eta < 0 or self.eta > 1:
            raise ValueError("eta must be a float between 0 and 1.")

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        lowest_score_idx: int,
        n_folds: int,
    ) -> Tuple[float, float]:
        """Returns a window of model performance based on percentile rank.

        Parameters
        ----------
        score_grid : np.ndarray
            A 2D array of model performance scores across folds and hyperparameter
            settings.
        cv_means : np.ndarray
            A 1D array of the average model performance across folds for each
            hyperparameter setting.
        best_score_idx : int
            The index of the highest performing hyperparameter setting.
        lowest_score_idx : int
            The index of the lowest performing hyperparameter setting.
        n_folds : int
            The number of folds used in the cross-validation.

        Returns
        -------
        min_cut : float
            The lower bound of the window of model performance.
        max_cut : float
            The upper bound of the window of model performance.
        """
        # Estimate the indicated percentile, and its inverse, across folds for
        # each column of the grid
        perc_cutoff = np.nanpercentile(
            score_grid, [100 * self.eta, 100 - 100 * self.eta], axis=1
        )

        # Determine bounds of the percentile interval
        max_cut = perc_cutoff[0, best_score_idx]
        min_cut = perc_cutoff[1, best_score_idx]
        return min_cut, max_cut

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(eta={self.eta})"


class by_signed_rank:
    """Slices a window of model performance based on signed rank sum.

    Signed rank sum is estimated based on a Wilcoxon rank sum test at a user-supplied
    alpha-level. The resulting window of model performance represents the range of
    scores that are not statistically different from the highest performing model.

    Attributes
    ----------
    alpha : float
        An alpha significance level in the case that wilcoxon rank sum
        hypothesis testing is used to filter outlying scores across folds.
        Default is 0.05.
    alternative : str
        The alternative hypothesis to test against. Must be one of 'two-sided',
        'less', or 'greater'. Default is 'two-sided'. See ``scipy.stats.wilcoxon`` for
        more details.
    zero_method : str
        The method used to handle zero scores. Must be one of 'pratt', 'wilcox',
        'zsplit'. Default is 'zsplit'. See ``scipy.stats.wilcoxon`` for more details.

    Raises
    ------
    ValueError
        If ``alpha`` is not a float between 0 and 1.
    """

    def __init__(
        self,
        alpha: float = 0.01,
        alternative: str = "two-sided",
        zero_method: str = "zsplit",
    ):  # noqa: D107
        self.alpha = alpha
        if not isinstance(self.alpha, float) or self.alpha < 0 or self.alpha > 1:
            raise ValueError("alpha must be a float between 0 and 1.")
        self.alternative = alternative
        self.zero_method = zero_method

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        lowest_score_idx: int,
        n_folds: int,
    ) -> Tuple[float, float]:
        """Returns a window of model performance based on signed rank sum.

        Parameters
        ----------
        score_grid : np.ndarray
            A 2D array of model performance scores across folds and hyperparameter
            settings.
        cv_means : np.ndarray
            A 1D array of the average model performance across folds for each
            hyperparameter setting.
        best_score_idx : int
            The index of the highest performing hyperparameter setting.
        lowest_score_idx : int
            The index of the lowest performing hyperparameter setting.
        n_folds : int
            The number of folds used in the cross-validation.

        Returns
        -------
        min_cut : float
            The lower bound of the window of model performance.
        max_cut : float
            The upper bound of the window of model performance.

        Raises
        ------
        ValueError
            If the number of folds is less than 3.
        """
        import itertools
        from scipy.stats import wilcoxon

        if n_folds < 3:
            raise ValueError("Number of folds must be greater than 2.")

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

        max_cut = np.nanmax(cv_means[surviving_ranks])
        min_cut = np.nanmin(cv_means[surviving_ranks])
        return min_cut, max_cut

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(alpha={self.alpha},"
            f" alternative={self.alternative}, zero_method={self.zero_method})"
        )


class by_fixed_window:
    """Slices a window of model performance based on arbitrary min/max cuts.

    Attributes
    ----------
    min_cut : float
        The lower bound of the window. Default is ``None``, which is the lowest score.
    max_cut : float
        The upper bound of the window. Default is ``None``, which is the highest score.

    Raises
    ------
    ValueError
        If ``min_cut`` is greater than ``max_cut.
    """

    def __init__(
        self, min_cut: Optional[float] = None, max_cut: Optional[float] = None
    ):  # noqa: D107
        self.min_cut = min_cut
        self.max_cut = max_cut
        if self.min_cut is not None and self.max_cut is not None:
            if self.min_cut > self.max_cut:
                raise ValueError("min_cut must be less than max_cut.")

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        lowest_score_idx: int,
        n_folds: int,
    ) -> Tuple[Union[float, None], Union[float, None]]:
        """Returns a window of performance based on min_cut and max_cut values.

        Parameters
        ----------
        score_grid : np.ndarray
            A 2D array of model performance scores across folds and hyperparameter
            settings.
        cv_means : np.ndarray
            A 1D array of the average model performance across folds for each
            hyperparameter setting.
        best_score_idx : int
            The index of the highest performing hyperparameter setting.
        lowest_score_idx : int
            The index of the lowest performing hyperparameter setting.
        n_folds : int
            The number of folds used in the cross-validation.

        Returns
        -------
        min_cut : float
            The lower bound of the window of model performance.
        max_cut : float
            The upper bound of the window of model performance.
        """
        return self.min_cut, self.max_cut

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}(min_cut={self.min_cut}, max_cut={self.max_cut})"
        )


class Refitter:
    """A refit factory for model subselection in GridSearchCV or RandomizedSearchCV.

    Model subselection can be useful for instance in the case that the user wishes to
    identify best-performing models above or below a particular threshold. It can also
    be useful for selecting alternative models whose performance is not meaningfully
    different from the best-performing model, but whose simplicity may be more
    preferable (e.g. to prevent overfitting).

    Attributes
    ----------
    cv_results_ : dict of numpy (masked) ndarrays
        A dict with keys as column headers and values as columns, as generated from
        fitting a GridSearchCV or RandomSearchCV object. See ``GridSearchCV`` or
        ``RandomSearchCV``, respectively, for more details.
    scoring : str
        The scoring metric of interest. Must be one of the scoring metrics
        specified in ``GridSearchCV`` or ``RandomSearchCV``.

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

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import Refitter, by_standard_error
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
    ...      ("reduce_dim", PCA(random_state=42)),
    ...      ("classify", LinearSVC(random_state=42, C=0.01)),
    ... ])
    >>> param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}
    >>> search = GridSearchCV(
    ...     pipe,
    ...     param_grid=param_grid,
    ...     scoring="accuracy",
    ... )
    >>> search.fit(X, y)
    GridSearchCV(estimator=Pipeline(steps=[('reduce_dim', PCA(random_state=42)),
                                           ('classify',
                                            LinearSVC(C=0.01, random_state=42))]),
                 param_grid={'reduce_dim__n_components': [6, 8, 10, 12, 14]},
                 scoring='accuracy')
    >>> ss = Refitter(search.cv_results_)
    >>> ss.fit(by_standard_error(sigma=1))
    (0.884825465639171, 0.9148526525904792)
    >>> refitted_index = ss.transform()
    Original best index: 4
    Refitted best index: 3
    Refitted best params: {'reduce_dim__n_components': 12}
    Refitted best score: 0.8926121943670691
    >>> refitted_index
    3
    """

    def __init__(self, cv_results_: Dict, scoring: str = "score"):
        self.cv_results_ = cv_results_
        self.cv_results_constrained_ = cv_results_.copy()
        self.scoring = scoring

    def _get_splits(self) -> List[str]:
        """Extracts CV splits corresponding to the specified ``scoring`` metric."""
        # Extract subgrid corresponding to the scoring metric of interest
        fitted_key_strings = "\t".join(list(self.cv_results_constrained_.keys()))
        if not all(s in fitted_key_strings for s in ["split", "params", "mean_test"]):
            raise TypeError(
                "cv_results_ must be a dict of fitted GridSearchCV or RandomSearchCV"
                " objects."
            )

        _splits = [
            i
            for i in list(self.cv_results_constrained_.keys())
            if f"test_{self.scoring}" in i and i.startswith("split")
        ]
        if len(_splits) == 0:
            raise KeyError(f"Scoring metric {self.scoring} not found in cv grid.")
        else:
            return _splits

    @property
    def _n_folds(self):
        # Extract number of folds from cv_results_. Note that we cannot get this from
        # the ``n_splits_`` attribute of the ``cv`` object because it is not exposed to
        # the refit callable.
        return len(
            list(
                set(
                    [
                        i.split("_")[0]
                        for i in list(self.cv_results_constrained_.keys())
                        if i.startswith("split")
                    ]
                )
            )
        )

    @property
    def _score_grid(self):
        # Extract subgrid corresponding to the scoring metric of interest
        return np.vstack(
            [self.cv_results_constrained_[cv] for cv in self._get_splits()]
        ).T

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

    def _apply_thresh(
        self,
        min_cut: Optional[float],
        max_cut: Optional[float],
    ) -> int:
        """Apply a performance threshold to the `_score_grid`.

        Parameters
        ----------
        min_cut : float
            The minimum performance threshold.
        max_cut : float
            The maximum performance threshold.
        """

        # Initialize a mask for the overall performance
        performance_mask = np.zeros(len(self._score_grid), dtype=bool)

        # Extract the overall performance
        if not min_cut:
            min_cut = float(np.nanmin(self._cv_means))
        if not max_cut:
            max_cut = float(np.nanmax(self._cv_means))

        # Mask all grid columns that are outside the performance window
        performance_mask = np.where(
            (self._cv_means >= float(min_cut)) & (self._cv_means <= float(max_cut)),
            True,
            False,
        )

        if np.sum(performance_mask) == 0:
            print(
                f"\nMin: {min_cut}\nMax: {max_cut}\nMeans across folds:"
                f" {self._cv_means}\n"
            )
            raise ValueError(
                "No valid grid columns remain within the boundaries of the specified"
                " performance window."
            )

        # For each hyperparameter in the grid, mask all grid columns that are outside
        # of the performance window
        for hyperparam in self.cv_results_constrained_["params"][0].keys():
            self.cv_results_constrained_[f"param_{hyperparam}"].mask = ~performance_mask

        # Among those models remaining within the performance window, find the highest
        # surviving rank (i.e. the lowest-performing model overall).
        highest_surviving_rank = np.nanmax(
            self.cv_results_constrained_["rank_test_score"][performance_mask]
        )

        # Return the index of the highest surviving rank. If the hyperparameter grid is
        # sorted sequentially from least to most complexity, the index of the
        # highest surving rank will equate to the index of the simplest model that
        # is not meaningfully different from the best-performing model.
        return int(
            np.nanargmax(
                self.cv_results_constrained_["rank_test_score"]
                == highest_surviving_rank
            )
        )

    def fit(self, selector: Callable) -> Tuple[Union[float, None], Union[float, None]]:
        """Fit the refitter using the specified selector callable.

        Parameters
        ----------
        selector : callable
            A callable that consumes GridSearchCV or RandomSearchCV results and
            returns a tuple of floats representing the lower and upper bounds of a
            target model performance window.

        Returns
        -------
        min_cut : float
            The lower bound of the target model performance window.
        max_cut : float
            The upper bound of the target model performance window.

        Raises
        ------
        ``TypeError``
            If the selector is not a callable.

        Notes
        -----
        The following keyword arguments will be automatically exposed to the selector
        by ``Refitter``:

        - best_score_idx : int
            The index of the highest performing model.
        - lowest_score_idx : int
            The index of the lowest performing model.
        - n_folds : int
            The number of cross-validation folds.
        - cv_means : array-like
            The mean performance of each model across the cross-validation folds. For
            example:

            ````
            array([0.63887341, 0.57323584, 0.50254565, 0.43688487, 0.37791086])
            ````

        - score_grid : array-like
            The performance of each model across the cross-validation folds. For
            example:

            ````
            array([[0.63888889, 0.58333333, 0.65181058, 0.66016713, 0.66016713],
                [0.53055556, 0.51111111, 0.57660167, 0.6183844 , 0.62952646],
                [0.47777778, 0.45277778, 0.46518106, 0.54874652, 0.56824513],
                [0.4       , 0.39166667, 0.41504178, 0.46518106, 0.51253482],
                [0.31666667, 0.33333333, 0.37047354, 0.40668524, 0.46239554]])
            ````
        """
        if not callable(selector):
            raise TypeError(
                f"``selector`` {selector} must be a callable but is {type(selector)}."
                " See ``Notes`` section of the"
                " :class:``~sklearn.model_selection.Refitter:fit`` API documentation"
                " for more details."
            )

        fit_params = {
            "score_grid": self._score_grid,
            "cv_means": self._cv_means,
            "best_score_idx": self._best_score_idx,
            "lowest_score_idx": self._lowest_score_idx,
            "n_folds": self._n_folds,
        }

        self.min_cut, self.max_cut = selector(**fit_params)
        return self.min_cut, self.max_cut

    def transform(self) -> int:
        """Subselects the best-performing model under the fitted constraints.

        Returns
        -------
        int
            The index of the best-performing model under the fitted constraints.

        Raises
        ------
        ``ValueError``
            If the refitter has not been fitted before calling the ``transform`` method.
        """
        if not hasattr(self, "min_cut") or not hasattr(self, "max_cut"):
            raise ValueError(
                "Refitter must be fitted before calling ``transform`` method."
            )

        best_index_ = self._apply_thresh(self.min_cut, self.max_cut)
        print(f"Original best index: {self._best_score_idx}")
        print(f"Refitted best index: {best_index_}")
        print(
            "Refitted best params:"
            f" {self.cv_results_constrained_['params'][best_index_]}"
        )
        print(
            "Refitted best score:"
            f" {self.cv_results_constrained_['mean_test_score'][best_index_]}"
        )
        return self._apply_thresh(self.min_cut, self.max_cut)


def _wrap_refit(cv_results_: Dict, selector: Callable, scoring: str = "score") -> int:
    """A wrapper function for the ``Refitter`` class.

    Should not be called directly. See the :class:``~sklearn.model_selection.Refitter``
    API documentation for more details.

    Parameters
    ----------
    cv_results_ : Dict
        The ``cv_results_`` attribute of a ``GridSearchCV`` or ``RandomSearchCV``
        object.
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
        The index of the best model under the performance constraints conferred by a
        ``selector``.
    """
    ss = Refitter(cv_results_, scoring=scoring)
    [min_cut, max_cut] = ss.fit(selector)
    print(f"Min: {min_cut}\nMax: {max_cut}")
    return ss.transform()


def constrain(selector: Callable, scoring: str = "score") -> Callable:
    """Callable returning the best index with constraints conferred by a ``selector``.

    Intended to be used as the ``refit`` parameter in ``GridSearchCV`` or
    ``RandomSearchCV``.

    Parameters
    ----------
    selector : callable
        Function that returns the lower and upper bounds of an acceptable performance
        window.
    scoring : str
        The scoring metric used to select the best model. Default is "score", and only
        needs to be changed in the case of multiple scoring metrics or a custom scoring
        function.

    Returns
    -------
    Callable
        A callable that returns the index of the best model under the performance
        constraints imposed by the selector strategy.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import constrain, by_standard_error
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
    ...      ("reduce_dim", PCA(random_state=42)),
    ...      ("classify", LinearSVC(random_state=42, C=0.01)),
    ... ])
    >>> param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14, 16, 18]}
    >>> search = GridSearchCV(
    ...     pipe,
    ...     param_grid=param_grid,
    ...     scoring="accuracy",
    ...     refit=constrain(by_standard_error(sigma=1)),
    ... )
    >>> search.fit(X, y)
    Min: 0.8898918397688278
    Max: 0.9186844524007791
    Original best index: 6
    Refitted best index: 3
    Refitted best params: {'reduce_dim__n_components': 12}
    Refitted best score: 0.8926121943670691
    ...
    >>> search.best_params_
    {'reduce_dim__n_components': 12}
    """
    # avoid returning a closure in a return statement to avoid pickling issues
    best_index_callable = partial(_wrap_refit, selector=selector, scoring=scoring)
    return best_index_callable
