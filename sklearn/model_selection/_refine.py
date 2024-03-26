"""
The :mod:``sklearn.model_selection._refine`` includes refit callable factories for
selecting models from ``GridSearchCV``, ``RandomizedSearchCV``, or
``HalvingRandomSearchCV`` objects based on a refined window of model performance
and favorability criteria.
"""

import itertools
import warnings
from functools import partial
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
from scipy import stats

__all__ = [
    "BaseScoreSlicer",
    "StandardErrorSlicer",
    "PercentileSlicer",
    "WilcoxonSlicer",
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
        # estimate the SE across folds for each column of the grid
        cv_se = np.array(np.nanstd(score_grid, axis=1) / np.sqrt(n_folds))

        # compute the confidence interval
        max_cut = cv_means[best_score_idx] + self.sigma * cv_se[best_score_idx]
        min_cut = cv_means[best_score_idx] - self.sigma * cv_se[best_score_idx]
        return min_cut, max_cut


class PercentileSlicer(BaseScoreSlicer):
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


class WilcoxonSlicer(BaseScoreSlicer):
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
        if n_folds < 3:
            raise ValueError("Number of folds must be greater than 2.")

        # perform signed Wilcoxon rank sum test for each pair combination of
        # columns against the best average score column
        tests = [
            pair
            for pair in list(itertools.combinations(range(score_grid.shape[0]), 2))
            if best_score_idx in pair
        ]

        pvals = {}
        for pair in tests:
            pvals[pair] = stats.wilcoxon(
                score_grid[pair[0]],
                score_grid[pair[1]],
                alternative=self.alternative,
                zero_method=self.zero_method,
            )[1]

        # return the ranks of the models that are insignificantly different from the
        # best average performing, else just return the best-performing
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
    """The FavorabilityRanker class provides a mechanism for ranking hyperparameters
    based on user-defined favorability rules. Favorability can be defined in terms
    of numerical lower-is-better, a specific order for categorical variables, or
    proximity to a statistical measure for distributions.

    Parameters
    ----------
    favorability_rules : dict
        A dictionary mapping hyperparameter names to a tuple where the first
        element is either a boolean indicating that lower numerical values are
        more favorable, a list of strings indicating the order of favorability for
        categorical variables, or a string specifying a statistical measure for
        distributions ('mean' or 'median'). The second element is a float
        indicating the relative importance of the hyperparameter.
    seed : int, optional
        An optional random seed for consistent random sampling.

    Examples
    --------
    >>> from scipy import stats
    >>> from sklearn.model_selection import FavorabilityRanker
    >>> favorability_rules = {
    ...     'reduce_dim__n_components': (True, 1.0),
    ...     'classify__degree': (True, 1.0),
    ...     'classify__kernel': (['linear', 'rbf'], 1.0),
    ...     'classify__C': ('median', 0.25),
    ... }
    >>> fr = FavorabilityRanker(favorability_rules, seed=42)
    >>> params = {
    ...     'reduce_dim__n_components': [6, 2, 8],
    ...     'classify__degree': [2, 3, 1],
    ...     'classify__kernel': ['rbf', 'linear'],
    ...     'classify__C': stats.norm(loc=0.0, scale=1.0),
    ... }
    >>> ranks = fr(params)
    >>> ranks
    [12, 11, 8, 7, 10, 9, 6, 2, 5, 1, 4, 18, 3, 14, 17, 13, 16, 15]
    """

    def __init__(
        self,
        favorability_rules: Dict[str, Tuple[Union[bool, List], float]],
        seed: Optional[int] = None,
    ):
        self.favorability_rules = favorability_rules
        self.seed = seed
        self.rng = np.random.default_rng(self.seed)
        self._validate_favorability_rules()

    def _validate_favorability_rules(self):
        for hyperparam, rule in self.favorability_rules.items():
            if not isinstance(hyperparam, str):
                raise TypeError(
                    f"Hyperparameter {hyperparam} must be a string, and correspond "
                    "to a hyperparameter in the hyperparameter grid."
                )
            if not isinstance(rule, tuple):
                raise TypeError(
                    f"Favorability rule for hyperparameter {hyperparam} must be a"
                    " tuple."
                )
            if (
                not isinstance(rule[0], bool)
                and not isinstance(rule[0], list)
                and not isinstance(rule[0], str)
            ):
                raise TypeError(
                    "First element of favorability rule for hyperparameter"
                    f" {hyperparam} must be a boolean, list, or string."
                )
            if not isinstance(rule[1], float):
                raise TypeError(
                    "Second element of favorability rule for hyperparameter"
                    f" {hyperparam} must be a float."
                )

    def _process_parameter_values(self, value: Any, rule: Tuple[Any, float]) -> Any:
        """Process a single hyperparameter value, handling distribution objects."""
        if (
            hasattr(value, "dist")
            and (
                isinstance(value.dist, stats.rv_continuous)
                or isinstance(value.dist, stats.rv_discrete)
            )
            and hasattr(value, "args")
            and hasattr(value, "kwds")
        ):
            distribution_property = rule[0]

            if distribution_property == "mean":
                return value.mean()
            elif distribution_property == "median":
                return value.median()
            elif distribution_property.startswith("percentile_"):
                percentile = float(distribution_property.split("_")[1])
                return value.ppf(percentile / 100, random_state=self.rng)
            else:
                raise ValueError(
                    f"Unsupported distribution property: {distribution_property}"
                )
        elif callable(value):
            # handle ParameterSampler or similar callable for generating parameter
            # values
            return value(self.rng)

        return value

    def __call__(self, params: Union[List[Dict], Dict]) -> List[int]:
        """
        Ranks the given hyperparameter sets based on the defined favorability rules.

        Parameters
        ----------
        params : Union[List[Dict], Dict]
            A parameter grid in the form of a list of dictionaries or a single
            dictionary.

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
                    hyperparam_value = self._process_parameter_values(
                        param_set[hyperparam], rule
                    )

                    is_numeric = isinstance(hyperparam_value, (int, float, np.number))

                    # numerical
                    if isinstance(rule[0], bool):
                        lower_is_favorable, weight = rule
                        if not is_numeric:
                            raise TypeError(
                                f"Expected numeric value for {hyperparam}, "
                                f"got {type(hyperparam_value)}"
                            )
                        if lower_is_favorable:
                            score_component = hyperparam_value
                        else:
                            score_component = (
                                0 if hyperparam_value == 0 else 1 / hyperparam_value
                            )
                        favorability_score += weight * score_component

                    # categorical
                    elif isinstance(rule[0], list):
                        rule_values, weight = rule
                        if hyperparam_value in rule_values:
                            score_component = rule_values.index(hyperparam_value)
                            favorability_score += weight * score_component
                        else:
                            raise ValueError(
                                f"Hyperparameter {hyperparam} "
                                f"must be one of {rule_values}."
                            )

            return favorability_score

        # check if 'params' is a list of dictionaries or a single dictionary
        if isinstance(params, list):
            # if it's a list, process each set of parameters separately
            favorability_scores = [calculate_favorability_score(p) for p in params]
        elif isinstance(params, dict):
            # if it's a dictionary, first check if it contains distribution objects or
            # callable parameter generators, and if so, but no seed is set, issue a
            # warning
            if self.seed is None and any(
                (
                    hasattr(v, "dist")
                    and (
                        isinstance(v.dist, stats.rv_continuous)
                        or isinstance(v.dist, stats.rv_discrete)
                    )
                    and hasattr(v, "args")
                    and hasattr(v, "kwds")
                )
                or callable(v)
                for p in params.values()
                for v in (p if isinstance(p, list) else [p])
            ):
                warnings.warn(
                    (
                        "A seed value was not set but distribution objects or callable"
                        " parameter generators are present in params. Favorability"
                        " ranks may not correspond to the actual parameter values"
                        " sampled during the SearchCV fitting."
                    ),
                    UserWarning,
                )
            # generate a grid of param combinations and
            # process them as though they were a list of dictionaries
            processed_params = {}
            for key, value in params.items():
                rule = self.favorability_rules.get(key, (None, 0))
                processed_params[key] = (
                    [self._process_parameter_values(v, rule) for v in value]
                    if isinstance(value, list)
                    else [self._process_parameter_values(value, rule)]
                )
            combinations = [
                dict(zip(processed_params.keys(), v))
                for v in itertools.product(*processed_params.values())
            ]
            favorability_scores = [
                calculate_favorability_score(p) for p in combinations
            ]
        else:
            raise ValueError(
                "params must be either a list of dictionaries or a single dictionary"
            )

        ranks = [x + 1 for x in np.argsort(favorability_scores)]
        return ranks

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.favorability_rules})"


class ScoreCutModelSelector:
    """A refit factory for refining model selection from GridSearchCV,
    RandomizedSearchCV, or HalvingRandomSearchCV.

    Model refinement can be useful in the case that the user wishes to select,
    according to one or more secondary criteria, the most favorable model from among a
    subset of similarly performant candidate models whose scores are not meaningfully
    different from the globally best-performing model, but whose favorability along the
    one or more secondary dimensions makes it nevertheless more preferable. One example
    of this is the case of selecting a simpler model that is not meaningfully different
    from the best-performing model -- the so-called "one-standard-error" rule. Other
    definitions of model favorability such as interpretability, explainability,
    computational efficiency, or other domain-specific crieria are also supported.

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
    >>> from sklearn.model_selection import ScoreCutModelSelector, \
        StandardErrorSlicer, FavorabilityRanker
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
    ...      ("reduce_dim", PCA(random_state=42)),
    ...      ("classify", LinearSVC(dual='auto', random_state=42, C=0.01)),
    ... ])
    >>> param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14]}
    >>> search = GridSearchCV(
    ...     pipe,
    ...     param_grid=param_grid,
    ...     scoring="accuracy",
    ... )
    >>> fitted = search.fit(X, y)
    >>> ss = ScoreCutModelSelector(fitted.cv_results_)
    >>> bounds = ss.fit(StandardErrorSlicer(sigma=1))
    Min: 0.884825465639171
    Max: 0.9148526525904792
    >>> favorability_rules = {
    ...     'reduce_dim__n_components': (True, 2.0), # Lower is simpler and
    ...                                              # more favorable
    ...     'classify__C': (False, 1.0) # Lower is more complex and
    ...                                 # less favorable
    ... }
    >>> favorable_index = ss.transform(FavorabilityRanker(favorability_rules))
    Original best index: 4
    Original best params: {'reduce_dim__n_components': 14}
    Original best score: 0.8998390591148251
    Promoted best index: 3
    Promoted best params: {'reduce_dim__n_components': 12}
    Promoted best score: 0.8926121943670691
    >>> favorable_index
    3
    """

    def __init__(self, cv_results_: Dict):
        self.cv_results_ = cv_results_.copy()
        self.cv_results_constrained_ = self.cv_results_.copy()

    def _get_splits(self) -> List[str]:
        """Extracts CV splits corresponding to the specified ``scoring`` metric."""
        # extract subgrid corresponding to the scoring metric of interest
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

    def _check_fitted(self, retval: Optional[Any] = None) -> Any:
        """Checks if the ``ScoreCutModelSelector`` instance has been fitted.

        Parameters
        ----------
        retval : Optional[Any], optional
            A return value, by default None.

        Returns
        -------
        Any
            The return value if the ``ScoreCutModelSelector`` instance has been
            fitted.

        Raises
        ------
        AttributeError
            If the ``ScoreCutModelSelector`` instance has not been fitted.
        """
        if hasattr(self, "min_cut_") or hasattr(self, "max_cut_"):
            return retval
        else:
            raise AttributeError(
                "The ``ScoreCutModelSelector`` instance has not been fitted. Please"
                " call the ``ScoreCutModelSelector:fit`` method first."
            )

    def _check_transformed(self, retval: Optional[Any] = None) -> Any:
        """Checks if the ``ScoreCutModelSelector`` instance has been transformed.

        Parameters
        ----------
        retval : Optional[Any], optional
            A return value, by default None.

        Returns
        -------
        Any
            The return value if the ``ScoreCutModelSelector`` instance has been
            transformed.

        Raises
        ------
        AttributeError
            If the ``ScoreCutModelSelector`` instance has not been transformed.
        AttributeError
            If the ``ScoreCutModelSelector`` instance has not been fitted.
        """
        if (hasattr(self, "min_cut_") or hasattr(self, "max_cut_")) and hasattr(
            self, "favorable_best_index_"
        ):
            return retval
        else:
            if not (hasattr(self, "min_cut_") or hasattr(self, "max_cut_")):
                raise AttributeError(
                    "The ``ScoreCutModelSelector`` instance has not been fitted. Please"
                    " call the ``ScoreCutModelSelector:fit`` method first."
                )
            else:
                raise AttributeError(
                    "The ``ScoreCutModelSelector`` instance has not been transformed."
                    " Please call the ``ScoreCutModelSelector:transform`` method first."
                )

    @property
    def _n_folds(self) -> int:
        # extract number of folds from cv_results_. Note that we cannot get this from
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
        # extract subgrid corresponding to the scoring metric of interest
        return np.vstack(
            [self.cv_results_constrained_[cv] for cv in self._get_splits()]
        ).T

    @property
    def _cv_means(self) -> np.ndarray:
        # calculate means of subgrid corresponding to the scoring metric of interest
        return np.array(np.nanmean(self._score_grid, axis=1))

    @property
    def _lowest_score_idx(self) -> int:
        # return index of the lowest performing model
        return np.nanargmin(self._cv_means)

    @property
    def _best_score_idx(self) -> int:
        # return index of the highest performing model
        return np.nanargmax(self._cv_means)

    @property
    def best_params_cut_(self) -> List[Dict]:
        return self._check_fitted(self.cv_results_constrained_["params"])

    @property
    def best_scores_cut_(self) -> np.ndarray:
        return self._check_fitted(self.cv_results_constrained_["mean_test_score"])

    @property
    def favorable_best_params_(self) -> Dict:
        return self._check_transformed(
            self.cv_results_constrained_["params"][self.favorable_best_index_]
        )

    @property
    def favorable_best_score_(self) -> float:
        return self._check_transformed(
            self.cv_results_constrained_["mean_test_score"][self.favorable_best_index_]
        )

    def _apply_thresh(
        self,
        min_cut: Optional[float],
        max_cut: Optional[float],
    ) -> Tuple[Optional[float], Optional[float]]:
        """Apply a performance threshold to the `_score_grid`.

        Parameters
        ----------
        min_cut : float
            The minimum performance threshold indicated by the ``score_slice_fn``.
        max_cut : float
            The maximum performance threshold indicated by the ``score_slice_fn``.

        Returns
        -------
        min_cut : float
            The minimum performance threshold applied to the `_score_grid`.
        max_cut : float
            The maximum performance threshold applied to the `_score_grid`.

        Raises
        ------
        ValueError
            If no valid grid columns remain within the boundaries of the specified
            performance window.
        """

        # initialize a mask for the overall performance
        np.zeros(len(self._score_grid), dtype=bool)

        # extract the overall performance
        if not min_cut:
            min_cut = float(np.nanmin(self._cv_means))
        if not max_cut:
            max_cut = float(np.nanmax(self._cv_means))

        # mask all grid columns that are outside the performance window
        self.performance_mask = np.where(
            (self._cv_means >= float(min_cut)) & (self._cv_means <= float(max_cut)),
            True,
            False,
        )

        if np.sum(self.performance_mask) == 0:
            raise ValueError(
                "No valid grid columns remain within the boundaries of the specified"
                " performance window. \nMin: {min_cut}\nMax: {max_cut}\nMeans across"
                " folds: {self._cv_means}\n"
            )

        # for each hyperparameter in the grid, mask all grid columns that are outside
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

        # apply the favorability function to the hyperparameter grid and add the ensuing
        # favorability ranks to the constrained cv_results dict
        self.cv_results_constrained_["favorability_rank"] = np.array(
            favorability_rank_fn(self.cv_results_constrained_["params"])
        )

        most_favorable_surviving_rank = np.nanmin(
            self.cv_results_constrained_["favorability_rank"][self.performance_mask]
        )

        # check if multiple models are tied for most favorable within the trimmed
        # performance
        tied_mostfavorable = (
            np.sum(
                self.cv_results_constrained_["favorability_rank"][self.performance_mask]
                == most_favorable_surviving_rank
            )
            > 1
        )

        # if only one model is definitively most favorable, return its index
        if not tied_mostfavorable:
            # Return the index of the most favorable model within the performance window
            return int(
                np.where(
                    self.cv_results_constrained_["favorability_rank"]
                    == most_favorable_surviving_rank
                )[0][0]
            )

        # mask performance mask to only include the multiple models tied for the
        # most favorable favorability
        performance_mask_mostfavorable = np.where(
            self.cv_results_constrained_["favorability_rank"]
            == most_favorable_surviving_rank,
            True,
            False,
        )
        # among multiple equally simple models within the trimmed performance
        # window, find the lowest surviving test performamce rank (i.e. the
        # highest-performing model overall).
        lowest_surviving_rank = np.nanmin(
            self.cv_results_constrained_["rank_test_score"][
                performance_mask_mostfavorable
            ]
        )

        # return the index of the lowest surviving rank among equally simple models,
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
        TypeError
            If the ``score_slice_fn`` is not a callable.

        Notes
        -----
        The following keyword arguments will be automatically exposed to the
        ``score_slice_fn`` by ``ScoreCutModelSelector``:

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

        self.min_cut_, self.max_cut_ = score_slice_fn(**fit_params)
        print(f"Min: {self.min_cut_}\nMax: {self.max_cut_}")

        return self._apply_thresh(self.min_cut_, self.max_cut_)

    def transform(self, favorability_rank_fn: Callable) -> int:
        """promotes the most favorable model within the constrained cut-window of
        best-performance.

        Parameters
        ----------
        favorability_rank_fn : callable
            A callable initialized with a dictionary of hyperparameter favorability
            ranking rules. At call, the callable consumes a hyperparameter grid in the
            form of a dictionary {"reduce_dim__n_components": [6, 8, 10, 12, 14, 16,
            18]} and returns an list of ranked integers corresponding to the
            favorability of each hyperparameter setting, from least to most complex.
            For example: [0, 1, 2, 3, 4, 5, 6].

        Returns
        -------
        int
            The index of the most favorable model.

        Raises
        ------
        TypeError
            If the ``favorability_rank_fn`` is not a callable.

        """
        if not callable(favorability_rank_fn):
            raise TypeError(
                f"``favorability_rank_fn`` {favorability_rank_fn} must be a callable"
                f" but is {type(favorability_rank_fn)}. See ``Notes`` section of the"
                " :class:``~sklearn.model_selection.ScoreCutModelSelector:fit`` API "
                " documentation for more details."
            )

        self._check_fitted()

        self.favorable_best_index_ = self._select_best_favorable(favorability_rank_fn)

        print(
            f"Original best index: {self._best_score_idx}\nOriginal best "
            f"params: {self.cv_results_['params'][self._best_score_idx]}\nOriginal "
            "best score: "
            f"{self.cv_results_['mean_test_score'][self._best_score_idx]}\nPromoted "
            f"best index: {self.favorable_best_index_}\nPromoted best params: "
            f"{self.favorable_best_params_}\nPromoted best "
            f"score: {self.favorable_best_score_}"
        )

        return self.favorable_best_index_


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
    ss.fit(score_slice_fn)
    return ss.transform(favorability_rank_fn)


def promote(score_slice_fn: Callable, favorability_rank_fn: Callable) -> Callable:
    """Return a callable to select the most favorable model within performance bounds.

    This function is designed to be used as the `refit` parameter in `GridSearchCV`,
    `RandomSearchCV`, or `HalvingRandomSearchCV`, allowing for custom model selection
    based on `favorability_rank_fn` within score constraints defined by
    `score_slice_fn`.

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

    Raises
    ------
    TypeError
        If ``score_slice_fn`` or ``favorability_rank_fn`` are not callable.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.decomposition import PCA
    >>> from sklearn.svm import LinearSVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import promote, StandardErrorSlicer, \
        FavorabilityRanker
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
    ...      ("reduce_dim", PCA(random_state=42)),
    ...      ("classify", LinearSVC(dual='auto', random_state=42, C=0.01)),
    ... ])
    >>> param_grid = {"reduce_dim__n_components": [6, 8, 10, 12, 14, 16, 18]}
    >>> favorability_rules = {
    ...     'reduce_dim__n_components': (True, 1.0),  # Lower is simpler and
    ...                                        # more favorable
    ...     'classify__C': (False, 1.0) # Lower is more complex and
    ...                          # less favorable
    ... }
    >>> search = GridSearchCV(
    ...     pipe,
    ...     param_grid=param_grid,
    ...     scoring="accuracy",
    ...     refit=promote(score_slice_fn=StandardErrorSlicer(sigma=1),
    ...     favorability_rank_fn=FavorabilityRanker(favorability_rules)),
    ... )
    >>> fitted = search.fit(X, y)
    Min: 0.8898918397688278
    Max: 0.9186844524007791
    Original best index: 6
    Original best params: {'reduce_dim__n_components': 18}
    Original best score: 0.9042881460848035
    Promoted best index: 3
    Promoted best params: {'reduce_dim__n_components': 12}
    Promoted best score: 0.8926121943670691
    >>> fitted.best_params_
    {'reduce_dim__n_components': 12}
    """
    if not callable(score_slice_fn) or not callable(favorability_rank_fn):
        raise TypeError(
            "``score_slice_fn`` and ``favorability_rank_fn`` must be callables"
        )  # this error gets raised initially here to avoid
        # delaying an equivalent error internal to
        # `ScoreCutModelSelector` until after the SearchCV
        # object has been fit.

    # avoid returning a closure in a return statement to avoid pickling issues
    best_index_callable = partial(
        _wrap_refit,
        score_slice_fn=score_slice_fn,
        favorability_rank_fn=favorability_rank_fn,
    )
    return best_index_callable
