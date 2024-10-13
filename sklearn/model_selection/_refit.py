"""
The :mod:``sklearn.model_selection._refit`` includes refit callable factories for
selecting models from ``GridSearchCV``, ``RandomizedSearchCV``, or
``HalvingRandomSearchCV`` objects based on a constrained window of model performance
and favorability criteria.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

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
    TypeError
        If ``sigma`` is not an integer.
    ValueError
        If ``sigma`` is less than 1.
    """

    def __init__(self, sigma: int = 1):
        self.sigma = sigma
        if not isinstance(self.sigma, int):
            raise TypeError("sigma must be an integer.")
        if self.sigma < 1:
            raise ValueError("sigma must be positive.")

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
    TypeError
        If ``eta`` is not a float.
    ValueError
        If ``eta`` is not between 0 and 1.
    """

    def __init__(self, eta: float = 0.68):
        self.eta = eta
        if not isinstance(self.eta, float):
            raise TypeError("eta must be a float.")
        if self.eta < 0 or self.eta > 1:
            raise ValueError("eta must be between 0 and 1.")

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
    TypeError
        If ``alpha`` is not a float.
    ValueError
        If ``alpha`` is not between 0 and 1.
    """

    def __init__(
        self,
        alpha: float = 0.01,
        alternative: str = "two-sided",
        zero_method: str = "zsplit",
    ):
        self.alpha = alpha
        if not isinstance(self.alpha, float):
            raise TypeError("alpha must be a float.")
        if self.alpha < 0 or self.alpha > 1:
            raise ValueError("alpha must be between 0 and 1.")
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
                "The average performance of all cross-validated models is "
                "significantly different from that of the best-performing model.",
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
    based on user-specified favorability rules. Favorability can be defined in terms
    of numerical lower-is-better, or a specific order for categorical variables.

    Parameters
    ----------
    favorability_rules : dict
        A dictionary mapping hyperparameter names to a tuple where the first
        element is either a boolean indicating that lower numerical values are
        more favorable, or a list of strings indicating the order of favorability for
        categorical variables. The second element is a float indicating the relative
        importance of the hyperparameter among all hyperparameters in the grid.

    Examples
    --------
    >>> from scipy import stats
    >>> from sklearn.model_selection import FavorabilityRanker
    >>> favorability_rules = {
    ...     'reduce_dim__n_components': (True, 1.0),
    ...     'classify__degree': (True, 1.0),
    ...     'classify__kernel': (['linear', 'rbf'], 1.0),
    ... }
    >>> fr = FavorabilityRanker(favorability_rules)
    >>> params = {
    ...     'reduce_dim__n_components': [6, 2, 8],
    ...     'classify__degree': [2, 3, 1],
    ...     'classify__kernel': ['rbf', 'linear'],
    ... }
    >>> ranks = fr(params)
    >>> ranks
    [17, 5, 14, 2, 11, 8, 16, 4, 13, 1, 10, 18, 6, 7, 15, 3, 12, 9]
    """

    def __init__(
        self,
        favorability_rules: Dict[str, Tuple[Union[bool, List], float]],
    ):
        self.favorability_rules = {
            k: favorability_rules[k] for k in sorted(favorability_rules)
        }
        self._validate_favorability_rules()

    def _validate_favorability_rules(self):
        """Validates the favorability rules provided by the user."""
        for hyperparam, rule in self.favorability_rules.items():
            if not isinstance(hyperparam, str):
                raise TypeError(f"Hyperparameter {hyperparam} must be a string.")
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
                    "First element of favorability rule tuple for hyperparameter"
                    f" {hyperparam} must be a boolean, list, or string."
                )
            if not isinstance(rule[1], float):
                raise TypeError(
                    "Second element of favorability rule tuple for hyperparameter"
                    f" {hyperparam} must be a float."
                )

            if rule[1] < 0:
                raise ValueError(
                    f"Weight for hyperparameter {hyperparam} must be non-negative."
                )

    def _parse_param_values(self, value: Any) -> Any:
        """Parses a single hyperparameter value, for a variety of data types."""
        if isinstance(value, (int, float, str, np.number)):
            return value
        else:
            raise ValueError(
                "FavorabilityRanker only supports numeric or string values for "
                f"hyperparameters. The provided value {value} is not supported."
            )

    @staticmethod
    def get_favorability_score(
        param_set: Dict, favorability_rules: Dict, parse_param_values: Callable
    ) -> float:
        """
        Calculates the favorability score for a given hyperparameter set.

        Parameters
        ----------
        param_set : Dict
            A dictionary of hyperparameters and their values.
        favorability_rules : Dict
            A dictionary mapping hyperparameter names to their favorability rules. Each
            entry in the dictionary corresponds to a hyperparameter and its rule, which
            can either be a boolean indicating whether a lower numerical value is more
            favorable or a list of categorical values sorted by their favorability.
        parse_param_values : Callable
            A function that parses hyperparameter values, ensuring they are in a format
            suitable for favorability calculation. This can involve converting data
            types or handling special cases.

        Returns
        -------
        float
            The calculated favorability score for the given set of hyperparameters.
            This score is used to rank the hyperparameter sets in terms of their
            favorability according to the specified rules.
        """
        favorability_score = 0
        for hyperparam, rule in favorability_rules.items():
            if hyperparam in param_set:
                hyperparam_value = parse_param_values(param_set[hyperparam])

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
                            f"Hyperparameter {hyperparam} must be one of {rule_values}."
                        )

        return favorability_score

    def __call__(self, params: Union[List[Dict], Dict]) -> List[int]:
        """
        Ranks the provided hyperparameter set based on the specified favorability rules.

        Parameters
        ----------
        params : Union[List[Dict], Dict]
            A parameter grid in the form of a list of dictionaries or a single
            dictionary.

        Returns
        -------
        List[int]
            A list of favorability ranks corresponding to each unique combination of
            hyperparameter values in the parameter grid.
        """

        # check if 'params' is a list of dictionaries or a single dictionary
        if isinstance(params, list):
            # if it's a list, process each dictionary separately in the order they
            # appear in the list
            favorability_scores = [
                self.get_favorability_score(
                    p,
                    self.favorability_rules,
                    lambda x: self._parse_param_values(x),
                )
                for p in params
            ]
        elif isinstance(params, dict):
            # if it's a single dictionary with hyperparameters as keys and values as
            # lists of hyperparameter values, generate all unique combinations of
            # hyperparameter values and calculate the favorability score for each
            from sklearn.model_selection import ParameterGrid

            # convert the params dictionary to a cleaned dictionary of lists, where each
            # key is a hyperparameter and each value is a list of hyperparameter
            # values, all of which are numeric or categorical, to be ranked.
            sorted_params = {k: params[k] for k in sorted(params)}
            processed_params = {}
            for key, value in sorted_params.items():
                processed_params[key] = (
                    [self._parse_param_values(v) for v in value]
                    if isinstance(value, list)
                    else [self._parse_param_values(value)]
                )

            # generate a ParameterGrid object from the processed_params dictionary
            # we opt to use ParameterGrid here to ensure that the order of the
            # hyperparameters is preserved in the output and is consistent with
            # the order as it is handled in the searchCV objects
            favorability_scores = [
                self.get_favorability_score(
                    p,
                    self.favorability_rules,
                    lambda x: self._parse_param_values(x),
                )
                for p in list(ParameterGrid(processed_params))
            ]
        else:
            raise ValueError(
                "`params` must be either a list of dictionaries or a single dictionary"
            )

        return [int(x + 1) for x in np.argsort(favorability_scores, kind="stable")]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.favorability_rules})"


class ScoreCutModelSelector:
    """A refit factory for refining model selection from GridSearchCV,
    RandomizedSearchCV, or HalvingRandomSearchCV.

    Constrained model selection can be useful in the case that a user wishes to select,
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
    >>> from sklearn.svm import SVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import ScoreCutModelSelector, \
        StandardErrorSlicer, FavorabilityRanker
    >>> import numpy as np
    >>> np.random.seed(42)  # set seed
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
    ...      ("reduce_dim", PCA(random_state=42)),
    ...      ("classify", SVC(kernel='linear', random_state=42, C=0.01)),
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
    Min: 0.9225
    Max: 0.9462
    >>> favorability_rules = {
    ...     'reduce_dim__n_components': (True, 2.0),  # Lower is simpler and
    ...                                               # more favorable
    ...     'classify__C': (False, 1.0)  # Lower is more complex and
    ...                                  # less favorable
    ... }
    >>> favorable_index = ss.transform(FavorabilityRanker(favorability_rules))
    Original best index: 4
    Original best params: {'reduce_dim__n_components': 14}
    Original best score: 0.9344
    Promoted best index: 2
    Promoted best params: {'reduce_dim__n_components': 10}
    Promoted best score: 0.9255
    >>> favorable_index
    2
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
        self._check_transformed()
        return self._check_transformed(
            self.cv_results_constrained_["params"][self.favorable_best_index_]
        )

    @property
    def favorable_best_score_(self) -> float:
        self._check_transformed()
        return self._check_transformed(
            self.cv_results_constrained_["mean_test_score"][self.favorable_best_index_]
        )

    def _apply_thresh(
        self,
        min_cut: Optional[float],
        max_cut: Optional[float],
        cv_results_constrained: Dict,
    ) -> Tuple[Dict, np.ndarray, Optional[float], Optional[float]]:
        """Apply a performance threshold to the `_score_grid`."""
        # initialize a mask for the overall performance
        _performance_mask = np.zeros(len(self._score_grid), dtype=bool)

        # extract the overall performance
        if not min_cut:
            min_cut = float(np.nanmin(self._cv_means))
        if not max_cut:
            max_cut = float(np.nanmax(self._cv_means))

        if min_cut > max_cut:
            raise ValueError(
                f"min_cut ({min_cut}) must be less than or equal to max_cut"
                f" ({max_cut})."
            )

        # mask all grid columns that are outside the performance window
        _performance_mask = np.where(
            (self._cv_means >= float(min_cut)) & (self._cv_means <= float(max_cut)),
            True,
            False,
        )

        if np.sum(_performance_mask) == 0:
            raise ValueError(
                "No valid grid columns remain within the boundaries of the specified"
                " performance window. \nMin: {min_cut}\nMax: {max_cut}\nMeans across"
                " folds: {self._cv_means}\n"
            )

        # for each hyperparameter in the grid, mask all grid columns that are outside
        # of the performance window
        for hyperparam in cv_results_constrained["params"][0].keys():
            cv_results_constrained[f"param_{hyperparam}"] = np.ma.masked_array(
                cv_results_constrained[f"param_{hyperparam}"],
                mask=~_performance_mask,
            )

        return cv_results_constrained, _performance_mask, min_cut, max_cut

    def _apply_favorability_ranks(
        self, favorability_rank_fn: Callable, cv_results_constrained: Dict
    ) -> Dict:
        """Apply the favorability function to the hyperparameter grid."""
        cv_results_constrained = cv_results_constrained.copy()
        cv_results_constrained["favorability_rank"] = np.array(
            favorability_rank_fn(cv_results_constrained["params"])
        )
        return cv_results_constrained

    def _get_most_favorable_surviving_rank(
        self, cv_results_constrained: Dict, _performance_mask: np.ndarray
    ) -> int:
        """Get the most favorable surviving rank within trimmed performance window."""
        return np.nanmin(cv_results_constrained["favorability_rank"][_performance_mask])

    def _check_tied_most_favorable(
        self,
        cv_results_constrained: Dict,
        _performance_mask: np.ndarray,
        most_favorable_surviving_rank: int,
    ) -> bool:
        """Check if multiple models are tied for most favorable within the trimmed
        performance window."""
        return (
            np.sum(
                cv_results_constrained["favorability_rank"][_performance_mask]
                == most_favorable_surviving_rank
            )
            > 1
        )

    def _get_most_favorable_index(
        self, cv_results_constrained: Dict, most_favorable_surviving_rank: int
    ) -> int:
        """Return the index of the most favorable model within the trimmed performance
        window."""
        return int(
            np.where(
                cv_results_constrained["favorability_rank"]
                == most_favorable_surviving_rank
            )[0][0]
        )

    def _get_mostfavorable_mask(
        self, cv_results_constrained: Dict, most_favorable_surviving_rank: int
    ) -> np.ndarray:
        """Get the performance mask for the most favorable models."""
        return np.where(
            cv_results_constrained["favorability_rank"]
            == most_favorable_surviving_rank,
            True,
            False,
        )

    def _get_lowest_surviving_rank(
        self, cv_results_constrained: Dict, mostfavorable_mask: np.ndarray
    ) -> int:
        """Get the lowest surviving test performance rank among the most favorable
        models."""
        return np.nanmin(cv_results_constrained["rank_test_score"][mostfavorable_mask])

    def _get_best_most_favorable_index(
        self, cv_results_constrained: Dict, lowest_surviving_rank: int
    ) -> int:
        """Return the index of the best-performing most favorable model."""
        return int(
            np.where(
                cv_results_constrained["rank_test_score"] == lowest_surviving_rank
            )[0][0]
        )

    def _select_best_favorable(
        self,
        favorability_rank_fn: Callable,
        cv_results_constrained: Dict,
        _performance_mask: np.ndarray,
    ) -> int:
        """Selects the most favorably ranked model within the trimmed performance
        window."""
        cv_results_constrained = self._apply_favorability_ranks(
            favorability_rank_fn, cv_results_constrained
        )
        most_favorable_surviving_rank = self._get_most_favorable_surviving_rank(
            cv_results_constrained, _performance_mask
        )
        tied_mostfavorable = self._check_tied_most_favorable(
            cv_results_constrained, _performance_mask, most_favorable_surviving_rank
        )

        if not tied_mostfavorable:
            return self._get_most_favorable_index(
                cv_results_constrained, most_favorable_surviving_rank
            )

        _mostfavorable_mask = self._get_mostfavorable_mask(
            cv_results_constrained, most_favorable_surviving_rank
        )
        lowest_surviving_rank = self._get_lowest_surviving_rank(
            cv_results_constrained, _mostfavorable_mask
        )
        return self._get_best_most_favorable_index(
            cv_results_constrained, lowest_surviving_rank
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
        print(f"Min: {self.min_cut_:.4f}\nMax: {self.max_cut_:.4f}")

        cv_results_constrained = self.cv_results_.copy()
        (
            cv_results_constrained,
            _performance_mask,
            self.min_cut_,
            self.max_cut_,
        ) = self._apply_thresh(self.min_cut_, self.max_cut_, cv_results_constrained)
        self.cv_results_constrained_ = cv_results_constrained
        self._performance_mask = _performance_mask

        return self.min_cut_, self.max_cut_

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

        cv_results_constrained = self.cv_results_constrained_.copy()
        _performance_mask = self._performance_mask.copy()
        self.favorable_best_index_ = self._select_best_favorable(
            favorability_rank_fn, cv_results_constrained, _performance_mask
        )

        print(
            f"Original best index: {self._best_score_idx}\nOriginal best "
            f"params: {self.cv_results_['params'][self._best_score_idx]}\nOriginal "
            "best score: "
            f"{self.cv_results_['mean_test_score'][self._best_score_idx]:.4f}"
            "\nPromoted best index: "
            f"{self.favorable_best_index_}\nPromoted best params: "
            f"{self.favorable_best_params_}\nPromoted best score: "
            f"{self.favorable_best_score_:.4f}"
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
    >>> from sklearn.svm import SVC
    >>> from sklearn.pipeline import Pipeline
    >>> from sklearn.model_selection import promote, StandardErrorSlicer, \
        FavorabilityRanker
    >>> import numpy as np
    >>> np.random.seed(42)
    >>> X, y = load_digits(return_X_y=True)
    >>> pipe = Pipeline([
    ...      ("reduce_dim", PCA(random_state=42)),
    ...      ("classify", SVC(kernel='linear', random_state=42)),
    ... ])
    >>> param_grid = {"reduce_dim__n_components": [18, 24, 30, 36],
    ... "classify__C": [0.0001, 0.001, 0.01, 1, 10]}
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
    Min: 0.9401
    Max: 0.9620
    Original best index: 7
    Original best params: {'classify__C': 0.001, 'reduce_dim__n_components': 36}
    Original best score: 0.9510
    Promoted best index: 17
    Promoted best params: {'classify__C': 10, 'reduce_dim__n_components': 24}
    Promoted best score: 0.9455
    >>> fitted.best_params_
    {'classify__C': 10, 'reduce_dim__n_components': 24}
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
