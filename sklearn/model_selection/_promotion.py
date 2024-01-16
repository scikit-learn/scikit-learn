"""
The :mod:``sklearn.model_selection.promote`` includes refit callable factories for
promoteing models from ``GridSearchCV``, ``RandomizedSearchCV``, or
``HalvingRandomSearchCV``
"""

import warnings
from functools import partial
from typing import Callable, Dict, List, Optional, Tuple, Union

import numpy as np

__all__ = [
    "BaseScoreSlicer",
    "StandardErrorSlicer",
    "PercentileRankSlicer",
    "SignedRankSlicer",
    "FixedWindowSlicer",
    "FavorabilityRanker",
    "ScoreCutModelSelector",
    "promote",
    "_wrap_refit",
]


class BaseScoreSlicer:
    """A base class for classes used to define a window of proximal model performance
    based on various criteria. Should not be called directly.
    """

    def __init__(self):
        pass

    def __call__(
        self,
        score_grid: np.ndarray,
        cv_means: np.ndarray,
        best_score_idx: int,
        lowest_score_idx: int,
        n_folds: int,
    ):
        """Returns a window of model performance based on a user-defined criterion.

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
        raise NotImplementedError("Subclasses must implement this method.")

    def __repr__(self) -> str:
        slice_params = ", ".join(
            [
                f"{key}={val}"
                for key, val in self.__dict__.items()
                if key
                not in [
                    "score_grid",
                    "cv_means",
                    "best_score_idx",
                    "lowest_score_idx",
                    "n_folds",
                ]
                and not key.startswith("_")
            ]
        )
        return f"{self.__class__.__name__}({slice_params})"


class StandardErrorSlicer(BaseScoreSlicer):
    """Slices a window of model performance based on standard error.

    Standard error is estimated based on a user-supplied number of standard errors,
    sigma. The resulting window of model performance represents the range of scores
    that fall within the indicated margin of error of the best performing model.

    Parameters
    ----------
    sigma : int
        Number of standard errors tolerance in the case that a standard error
        threshold is used to filter outlying scores across folds. Default is 1.

    Raises
    ------
    ValueError
        If sigma is not a positive integer.
    """

    def __init__(self, sigma: int = 1):
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


class PercentileRankSlicer(BaseScoreSlicer):
    """Slices a window of model performance based on percentile rank.

    Percentile rank is estimated based on a user-supplied percentile threshold, eta.
    The resulting window of model performance represents the range of scores that fall
    within the indicated percentile range of the best performing model.

    Parameters
    ----------
    eta : float
        Percentile tolerance in the case that a percentile threshold is used to filter
        outlier scores across folds. Default is 0.68.

    Raises
    ------
    ValueError
        If eta is not a float between 0 and 1.
    """

    def __init__(self, eta: float = 0.68):
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


class SignedRankSlicer(BaseScoreSlicer):
    """Slices a window of model performance based on signed rank sum.

    Signed rank sum is estimated based on a Wilcoxon rank sum test at a user-supplied
    alpha-level. The resulting window of model performance represents the range of
    scores that are not statistically different from the highest performing model.

    Parameters
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
    ):
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

        if len(surviving_ranks) == 1:
            surviving_ranks = [best_score_idx]
            warnings.warn(
                (
                    "The average performance of all cross-validated models is "
                    "significantly different from that of the best-performing model."
                ),
                UserWarning,
            )

        max_cut = np.nanmax(cv_means[surviving_ranks])
        min_cut = np.nanmin(cv_means[surviving_ranks])
        return min_cut, max_cut


class FixedWindowSlicer(BaseScoreSlicer):
    """Slices a window of model performance based on arbitrary min/max cuts.

    Parameters
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
    ):
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
    ) -> Tuple[Optional[float], Optional[float]]:
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


class FavorabilityRanker:
    def __init__(
        self, favorability_rules: Dict[str, Union[Tuple[bool, float], List[str]]]
    ):
        """
        Initializes the FavorabilityRanker class with optional favorability rules.

        Parameters
        ----------
        favorability_rules: Dict[str, Union[Tuple[bool, float], List[str]]]
            A dictionary mapping hyperparameters to either a tuple or a list.
            For numeric hyperparameters, the tuple contains a boolean (indicating
            whether lower values imply higher favorability) and a float (weight).
            For string hyperparameters, the list defines the order of favorability
            from most favorable to least favorable.
        """
        self.favorability_rules = favorability_rules

    def __call__(self, params: List[Dict]) -> List[int]:
        """
        Ranks the given hyperparameter sets based on the defined favorability rules.

        Parameters
        ----------
        params: List[Dict]
            A list of dictionaries representing sets of hyperparameters.

        Returns
        -------
        List[int]
            A list of ranks corresponding to the favorability of each set of
            hyperparameters.
        """

        def calculate_favorability_score(param_set):
            favorability_score = 0
            for hyperparam, rule in self.favorability_rules.items():
                if hyperparam in param_set:
                    # Numeric hyperparameter
                    if isinstance(rule, tuple):
                        lower_is_more_complex, weight = rule
                        score_component = (
                            1 / param_set[hyperparam]
                            if lower_is_more_complex
                            else param_set[hyperparam]
                        )
                        favorability_score += weight * score_component
                    # Categorical hyperparameter
                    elif isinstance(rule, list):
                        favorability_score += rule.index(param_set[hyperparam])
            return favorability_score

        favorability_scores = [calculate_favorability_score(p) for p in params]
        ranks = [x + 1 for x in np.argsort(favorability_scores)]

        return ranks


class ScoreCutModelSelector:
    """A refit factory for promoting models in GridSearchCV, RandomizedSearchCV, or
    HalvingRandomSearchCV.

    Model promotion can be useful for instance in the case that the user wishes to
    select the most favorable alternative model whose performance is not meaningfully
    different from the globally best-performing model, but whose favorability along
    other dimensions makes it more preferable. One example of this is the case of
    selecting a simpler model that is not meaningfully different from the
    best-performing model -- the so-called "one-standard-error" rule. Other definitions
    of model favorability such as interpretability, explainability, computational
    efficiency, or other domain-specific crieria are also supported, and can be defined
    in the same way by the user using the ``FavorabilityRanker`` interface.

    Parameters
    ----------
    cv_results_ : dict of numpy (masked) ndarrays
        A dict with keys as column headers and values as columns, as generated from
        fitting a ``GridSearchCV``, ``RandomSearchCV``, or ``HalvingRandomSearchCV``
        object. For more details, see ``GridSearchCV``, ``RandomSearchCV``, or
        ``HalvingRandomSearchCV``, respectively.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import ScoreCutModelSelector,
    ... StandardErrorSlicer, FavorabilityRanker
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
    >>> ss = ScoreCutModelSelector(search.cv_results_)
    >>> ss.fit(StandardErrorSlicer(sigma=1))
    >>> favorability_rules = {
    ...     'reduce_dim__n_components': (False, 2.0)  # Lower is simpler
    ...     'classify__C': (True, 1.0) # Lower is more complex
    ... }
    >>> promoted_index = ss.transform(FavorabilityRanker(favorability_rules))
    Original best index: 4
    Promoted best index: 3
    Promoted best params: {'reduce_dim__n_components': 12}
    Promoted best score: 0.8926121943670691
    >>> promoted_index
    3
    """

    def __init__(self, cv_results_: Dict):
        self.cv_results_ = cv_results_
        self.cv_results_constrained_ = cv_results_.copy()

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
            if "test_score" in i and i.startswith("split")
        ]
        if len(_splits) == 0:
            raise KeyError("No splits found in cv grid.")
        else:
            return _splits

    @property
    def _n_folds(self) -> int:
        # Extract number of folds from cv_results_. Note that we cannot get this from
        # the ``n_splits_`` attribute of the ``cv`` object because it is not exposed to
        # the refit callable.
        return len(
            set(
                key.split("_")[0]
                for key in self.cv_results_constrained_
                if key.startswith("split")
            )
        )

    @property
    def _score_grid(self) -> np.ndarray:
        # Extract subgrid corresponding to the scoring metric of interest
        return np.vstack(
            [self.cv_results_constrained_[cv] for cv in self._get_splits()]
        ).T

    @property
    def _cv_means(self) -> np.ndarray:
        # Calculate means of subgrid corresponding to the scoring metric of interest
        return np.array(np.nanmean(self._score_grid, axis=1))

    @property
    def _lowest_score_idx(self) -> int:
        # Return index of the lowest performing model
        return np.nanargmin(self._cv_means)

    @property
    def _best_score_idx(self) -> int:
        # Return index of the highest performing model
        return np.nanargmax(self._cv_means)

    def _apply_thresh(
        self,
        min_cut: Optional[float],
        max_cut: Optional[float],
    ) -> Tuple[Optional[float], Optional[float]]:
        """Apply a performance threshold to the `_score_grid`.

        Parameters
        ----------
        min_cut : float
            The minimum performance threshold.
        max_cut : float
            The maximum performance threshold.
        """

        # Initialize a mask for the overall performance
        np.zeros(len(self._score_grid), dtype=bool)

        # Extract the overall performance
        if not min_cut:
            min_cut = float(np.nanmin(self._cv_means))
        if not max_cut:
            max_cut = float(np.nanmax(self._cv_means))

        # Mask all grid columns that are outside the performance window
        self.performance_mask = np.where(
            (self._cv_means >= float(min_cut)) & (self._cv_means <= float(max_cut)),
            True,
            False,
        )

        if np.sum(self.performance_mask) == 0:
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
            self.cv_results_constrained_[f"param_{hyperparam}"].mask = (
                ~self.performance_mask
            )
        return min_cut, max_cut

    def _select_best_favorable(self, favorability_rank_fn: Callable) -> int:
        """Selects the most favorably ranked model within the trimmed performance
        window. If multiple models are tied for highest rank, the model with the
        highest overall performance is selected.

        Parameters
        ----------
        favorability_rank_fn : callable
            A callable that consumes a hyperparameter grid in the form of a list of
            hyperparameter dictionaries and returns an list of ranked integers
            corresponding to model favorability of each hyperparameter setting, from
            least to most favorable. For example a valid ``favorability_rank_fn`` might
            consume [{"reduce_dim__n_components": [6, 8, 10, 12, 14, 16, 18]}]
            and return [0, 1, 2, 3, 4, 5, 6].
        """

        # Apply the favorability function to the hyperparameter grid and add the ensuing
        # favorability ranks to the constrained cv_results dict
        self.cv_results_constrained_["favorability_rank"] = np.array(
            favorability_rank_fn(self.cv_results_constrained_["params"])
        )

        most_favorable_surviving_rank = np.nanmin(
            self.cv_results_constrained_["favorability_rank"][self.performance_mask]
        )

        # Check if multiple models are tied for most favorable within the trimmed
        # performance
        tied_mostfavorable = (
            np.sum(
                self.cv_results_constrained_["favorability_rank"][self.performance_mask]
                == most_favorable_surviving_rank
            )
            > 1
        )

        # If only one model is definitively most favorable, return its index
        if not tied_mostfavorable:
            # Return the index of the most favorable model within the performance window
            return int(
                np.where(
                    self.cv_results_constrained_["favorability_rank"]
                    == most_favorable_surviving_rank
                )[0][0]
            )

        # Mask performance mask to only include the multiple models tied for the
        # most favorable favorability
        performance_mask_mostfavorable = np.where(
            self.cv_results_constrained_["favorability_rank"]
            == most_favorable_surviving_rank,
            True,
            False,
        )
        # Among multiple equally simple models within the trimmed performance
        # window, find the lowest surviving test performamce rank (i.e. the
        # highest-performing model overall).
        lowest_surviving_rank = np.nanmin(
            self.cv_results_constrained_["rank_test_score"][
                performance_mask_mostfavorable
            ]
        )

        # Return the index of the lowest surviving rank among equally simple models,
        # which will equate to the index of the best-performing most favorable model
        # that is not meaningfully different from the globally best-performing model.
        return int(
            np.where(
                self.cv_results_constrained_["rank_test_score"] == lowest_surviving_rank
            )[0][0]
        )

    def fit(
        self, score_slice_fn: Callable
    ) -> Tuple[Union[float, None], Union[float, None]]:
        """Fits a ScoreCutModelSelector instance with specified
        ``score_slice_fn`` callable to define a cut window of model performance.

        Parameters
        ----------
        score_slice_fn : callable
            A callable that consumes GridSearchCV, RandomSearchCV, or
            HalvingRandomSearchCV results and returns a tuple of floats representing
            the lower and upper bounds of a target model performance window.

        Returns
        -------
        min_cut : float
            The lower bound of the target model performance window.
        max_cut : float
            The upper bound of the target model performance window.

        Raises
        ------
        ``TypeError``
            If the ``score_slice_fn`` is not a callable.

        Notes
        -----
        The following keyword arguments will be automatically exposed to the
        score_slice_fn
        by ``ScoreCutModelSelector``:

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
        if not callable(score_slice_fn):
            raise TypeError(
                f"``score_slice_fn`` {score_slice_fn} must be a callable but is"
                " {type(score_slice_fn)}."
                " See ``Notes`` section of the"
                " :class:``~sklearn.model_selection.ScoreCutModelSelector:fit`` API "
                " documentation for more details."
            )

        fit_params = {
            "score_grid": self._score_grid,
            "cv_means": self._cv_means,
            "best_score_idx": self._best_score_idx,
            "lowest_score_idx": self._lowest_score_idx,
            "n_folds": self._n_folds,
        }

        self.min_cut, self.max_cut = score_slice_fn(**fit_params)
        min_cut, max_cut = self._apply_thresh(self.min_cut, self.max_cut)

        # Check that the stateful attributes have been correctly set
        assert min_cut == self.min_cut
        assert max_cut == self.max_cut

        return self.min_cut, self.max_cut

    def transform(self, favorability_rank_fn: Callable) -> int:
        """promotes the most favorable model within the constrained cut-window of
        best-performance.

        Returns
        -------
        favorability_rank_fn : callable
            A callable initialized with a dictionary of hyperparameter favorability
            ranking rules. At call, the callable consumes a hyperparameter grid in the
            form of a dictionary {"reduce_dim__n_components": [6, 8, 10, 12, 14, 16,
            18]} and returns an list of ranked integers corresponding to the
            favorability of each hyperparameter setting, from least to most complex.
            For example: [0, 1, 2, 3, 4, 5, 6].
        int
            The index of the most favorable model.

        Raises
        ------
        ValueError
            If the ScoreCutModelSelector has not been fitted before calling the
            ``transform`` method.

        """
        if not hasattr(self, "min_cut") or not hasattr(self, "max_cut"):
            raise ValueError(
                "ScoreCutModelSelector must be fitted before calling "
                "``transform`` method."
            )

        best_index_ = self._select_best_favorable(favorability_rank_fn)

        print(f"Original best index: {self._best_score_idx}")
        print(f"Promoted best index: {best_index_}")
        print(
            "Promoted best params:"
            f" {self.cv_results_constrained_['params'][best_index_]}"
        )
        print(
            "Promoted best score:"
            f" {self.cv_results_constrained_['mean_test_score'][best_index_]}"
        )
        return best_index_


def _wrap_refit(
    cv_results_: Dict, score_slice_fn: Callable, favorability_rank_fn: Callable
) -> int:
    """A wrapper function for the ``ScoreCutModelSelector`` class.

    Should not be called directly. See the :class:``~sklearn.model_selection
    ScoreCutModelSelector`` API documentation for more details.

    Parameters
    ----------
    cv_results_ : Dict
        The ``cv_results_`` attribute of a GridSearchCV, RandomSearchCV, or
            HalvingRandomSearchCV object.
    score_slice_fn : Callable
        Function that returns the lower and upper bounds of an acceptable performance
        window.
    favorability_rank_fn : Callable
        A callable that consumes a hyperparameter grid in the form of a dictionary
        {"reduce_dim__n_components": [6, 8, 10, 12, 14, 16, 18]} and returns an
        list of ranked integers corresponding to the favorability of each
        hyperparameter setting, from least to most complex. For example: [0, 1, 2,
        3, 4, 5, 6].

    Returns
    -------
    int
        The index of the best model under the performance constraints conferred by a
        ``score_slice_fn``.
    """
    ss = ScoreCutModelSelector(cv_results_)

    [min_cut, max_cut] = ss.fit(score_slice_fn)
    print(f"Min: {min_cut}\nMax: {max_cut}")
    return ss.transform(favorability_rank_fn)


def promote(score_slice_fn: Callable, favorability_rank_fn: Callable) -> Callable:
    """Callable returning the most favorable model index based on
    ``favorability_rank_fn`` with score constraints determined by a ``score_slice_fn``.

    Intended to be used as the ``refit`` parameter in ``GridSearchCV``,
    ``RandomSearchCV``, or ``HalvingRandomSearchCV``.

    Parameters
    ----------
    score_slice_fn : callable
        Function that returns the lower and upper bounds of an acceptable performance
        window.
    favorability_rank_fn : Callable
        Function that returns the favorability rank of each hyperparameter setting.

    Returns
    -------
    Callable
        A callable that returns the index of the most favorable model whose performance
        falls within the acceptable bounds imposed by the score_slice_fn rule.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import promote, StandardErrorSlicer,
    ... FavorabilityRanker
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
    ...      ("reduce_dim", PCA(random_state=42)),
    ...      ("classify", LinearSVC(random_state=42, C=0.01)),
    ... ])
    >>> param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14, 16, 18]}
    >>> favorability_rules = {
    ...     'reduce_dim__n_components': False,  # Lower is simpler
    ...     'classify__C': True # Lower is more complex
    ... }
    >>> search = GridSearchCV(
    ...     pipe,
    ...     param_grid=param_grid,
    ...     scoring="accuracy",
    ...     refit=promote(score_slice_fn=StandardErrorSlicer(sigma=1),
    ...     favorability_rank_fn=FavorabilityRanker(favorability_rules)),
    ... )
    >>> search.fit(X, y)
    Min: 0.8898918397688278
    Max: 0.9186844524007791
    Original best index: 6
    Promoted best index: 3
    Promoted best params: {'reduce_dim__n_components': 12}
    Promoted best score: 0.8926121943670691
    ...
    >>> search.best_params_
    {'reduce_dim__n_components': 12}
    """
    # avoid returning a closure in a return statement to avoid pickling issues
    best_index_callable = partial(
        _wrap_refit,
        score_slice_fn=score_slice_fn,
        favorability_rank_fn=favorability_rank_fn,
    )
    return best_index_callable
