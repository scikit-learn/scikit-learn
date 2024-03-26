import numpy as np
import pytest
import scipy.stats
from _pytest.mark.structures import ParameterSet

from sklearn.datasets import make_classification
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.experimental import enable_halving_search_cv  # noqa: F401
from sklearn.model_selection import (
    FavorabilityRanker,
    FixedWindowSlicer,
    GridSearchCV,
    HalvingRandomSearchCV,
    PercentileSlicer,
    RandomizedSearchCV,
    ScoreCutModelSelector,
    StandardErrorSlicer,
    WilcoxonSlicer,
    promote,
)
from sklearn.pipeline import Pipeline
from sklearn.svm import SVC, LinearSVC
from sklearn.utils._testing import ignore_warnings


@pytest.fixture(scope="function")
def grid_search_simulated():
    X, y = make_classification(n_samples=50, n_features=4, random_state=42)

    n_splits = 3
    params = [
        dict(
            kernel=[
                "rbf",
            ],
            C=[1, 10],
            gamma=[0.1, 1],
        ),
        dict(
            kernel=[
                "poly",
            ],
            degree=[1, 2],
        ),
    ]

    search = GridSearchCV(
        SVC(), cv=n_splits, param_grid=params, return_train_score=True
    )
    search.fit(X, y)

    cv_results = search.cv_results_

    yield {"cv_results": cv_results, "n_splits": n_splits}


@pytest.fixture(scope="function")
def generate_fit_params(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]
    n_splits = grid_search_simulated["n_splits"]
    ss = ScoreCutModelSelector(cv_results)

    yield {
        "score_grid": ss._score_grid,
        "n_folds": n_splits,
        "cv_means": ss._cv_means,
        "best_score_idx": ss._best_score_idx,
        "lowest_score_idx": ss._lowest_score_idx,
    }


def mock_favorability_ranker(hyperparams):
    print("Received hyperparameters:", hyperparams)
    ranks = []
    for params in hyperparams:
        rank = 0
        if params["kernel"] == "rbf":
            rank += 1  # Assuming 'rbf' is more favorable
        if params.get("C", 1) < 5:
            rank += 1  # Assuming lower 'C' is more favorable
        ranks.append(rank)
    # Debug: print ranks
    for idx, rank in enumerate(ranks):
        print(f"Hyperparameters at index {idx} have rank {rank}")
    return [
        rank for rank, _ in sorted(enumerate(ranks), key=lambda x: x[1], reverse=True)
    ]


def test_ScoreCutModelSelector_methods(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]
    n_splits = grid_search_simulated["n_splits"]

    ss = ScoreCutModelSelector(cv_results)

    # Test that the _get_splits method extracts the correct subgrid
    assert len(ss._get_splits()) == n_splits

    # Test that the _n_folds property returns the correct number of folds
    assert ss._n_folds == n_splits

    # Test that the _score_grid property returns the correct subgrid of scores
    assert ss._score_grid.shape == (6, n_splits)

    # Test that the _cv_means property returns the correct array of mean scores
    assert ss._cv_means.shape == (6,)

    # Test that the _lowest_score_idx property returns the correct index
    assert ss._lowest_score_idx == 5

    # Test that the _best_score_idx property returns the correct index
    assert ss._best_score_idx == 0

    assert ss._apply_thresh(0.93, 0.96) == (0.93, 0.96)

    # Omit min_thresh
    assert ss._apply_thresh(None, 0.99) == (0.4803921568627451, 0.99)

    # Omit max_thresh
    assert ss._apply_thresh(0.80, None) == (0.8, 0.9583333333333334)

    # Test that the fit method returns the correct score cuts
    assert ss.fit(StandardErrorSlicer(sigma=1)) == (
        0.9243126424613448,
        0.9923540242053219,
    )

    assert ss.transform(mock_favorability_ranker) == 0

    assert ss.best_params_cut_ == [
        {"C": 1, "gamma": 0.1, "kernel": "rbf"},
        {"C": 1, "gamma": 1, "kernel": "rbf"},
        {"C": 10, "gamma": 0.1, "kernel": "rbf"},
        {"C": 10, "gamma": 1, "kernel": "rbf"},
        {"degree": 1, "kernel": "poly"},
        {"degree": 2, "kernel": "poly"},
    ]

    assert (
        np.testing.assert_almost_equal(
            ss.best_scores_cut_,
            np.array(
                [0.95833333, 0.93872549, 0.93872549, 0.93995098, 0.95833333, 0.48039216]
            ),
        )
        is None
    )


def test_ScoreCutModelSelector_errors(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]
    n_splits = grid_search_simulated["n_splits"]

    with pytest.raises(ValueError):
        ss = ScoreCutModelSelector(cv_results)
        assert ss._apply_thresh(0.98, 0.99) == 1

    with pytest.raises(TypeError):
        ss = ScoreCutModelSelector(cv_results)
        assert ss.fit("Not_a_rule") == (0.9243126424613448, 0.9923540242053219)

    with pytest.raises(TypeError):
        ss = ScoreCutModelSelector(cv_results)
        ss.fit(StandardErrorSlicer(sigma=1))
        assert ss.transform("Not_a_rule") == 1

    del cv_results["params"]
    ss = ScoreCutModelSelector(cv_results)
    with pytest.raises(TypeError):
        assert len(ss._get_splits()) == n_splits


@ignore_warnings
@pytest.mark.parametrize(
    "param",
    [
        "reduce_dim__n_components",
        None,
    ],
)
@pytest.mark.parametrize(
    "scoring,score_slice_rule",
    [
        ("roc_auc", StandardErrorSlicer(sigma=1)),
        ("roc_auc", WilcoxonSlicer(alpha=0.01)),
        ("roc_auc", PercentileSlicer(eta=0.68)),
        ("roc_auc", FixedWindowSlicer(min_cut=0.96, max_cut=0.97)),
        ("roc_auc", "Not_a_rule"),
        ("neg_log_loss", StandardErrorSlicer(sigma=1)),
        ("neg_log_loss", WilcoxonSlicer(alpha=0.01)),
        ("neg_log_loss", PercentileSlicer(eta=0.68)),
        (
            "neg_log_loss",
            pytest.param(
                FixedWindowSlicer(min_cut=0.96, max_cut=0.97), marks=pytest.mark.xfail
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "favorability_rank_rule",
    [
        FavorabilityRanker(
            {
                "reduce_dim__n_components": (True, 2.0),  # Lower is more favorable
                "classify__C": (False, 1.0),  # Lower is less favorable
            }
        ),
        FavorabilityRanker(
            {
                "reduce_dim__n_components": ([4, 8, 12], 2.0),  # String-based rule
                "classify__C": (False, 1.0),  # Lower is less favorable
            }
        ),
        "Not_a_rule",
    ],
)
@pytest.mark.parametrize(
    "search_cv",
    ["GridSearchCV", "RandomizedSearchCV"],
)
def test_promote(param, scoring, score_slice_rule, favorability_rank_rule, search_cv):
    """
    A function that tests the promote function by comparing the results of a
    refitted grid and random search object to those of a non-refitted grid and random
    search object, respectively.
    """

    if search_cv == "GridSearchCV":
        search_cv = GridSearchCV
    else:
        search_cv = RandomizedSearchCV

    X, y = make_classification(n_samples=350, n_features=16, random_state=42)

    # Instantiate a pipeline with parameter grid representing different levels of
    # favorability
    clf = LinearSVC(dual="auto", random_state=42)
    if param == "reduce_dim__n_components":
        param_grid = {"reduce_dim__n_components": [4, 8, 12]}
        pipe = Pipeline([("reduce_dim", PCA(random_state=42)), ("classify", clf)])
    else:
        param_grid = {"classify__C": [0.1, 1], "reduce_dim__n_components": [4, 8, 12]}
        pipe = Pipeline(
            [("reduce_dim", PCA(random_state=42)), ("classify", SVC(random_state=42))]
        )

    # Instantiate a non-refitted grid search object for comparison
    grid = search_cv(pipe, param_grid, scoring=scoring, n_jobs=1)
    grid.fit(X, y)

    score_slice_rule = (
        score_slice_rule.values[0]
        if isinstance(score_slice_rule, ParameterSet)
        else score_slice_rule
    )
    favorability_rank_rule = (
        favorability_rank_rule.values[0]
        if isinstance(favorability_rank_rule, ParameterSet)
        else favorability_rank_rule
    )

    if score_slice_rule == "Not_a_rule" or favorability_rank_rule == "Not_a_rule":
        with pytest.raises(TypeError):
            promote(score_slice_rule, favorability_rank_rule)
        return

    # Instantiate a refitted grid search object
    grid_refitted = search_cv(
        pipe,
        param_grid,
        scoring=scoring,
        refit=promote(score_slice_rule, favorability_rank_rule),
    )

    # If the cv results were not all NaN, then we can test the refit callable
    if not np.isnan(grid.fit(X, y).cv_results_["mean_test_score"]).all():
        if param and param not in favorability_rank_rule.favorability_rules:
            with pytest.raises(ValueError):
                grid_refitted.fit(X, y)
        else:
            grid_refitted.fit(X, y)
        simplified_best_score_ = grid_refitted.cv_results_["mean_test_score"][
            grid_refitted.best_index_
        ]
        # Ensure that if the refit callable promoted a lower scoring model,
        # it was because it was only because it was a more favorable model.
        if abs(grid.best_score_) > abs(simplified_best_score_):
            assert grid.best_index_ != grid_refitted.best_index_
            if param:
                assert grid.best_params_[param] > grid_refitted.best_params_[param]
        elif grid.best_score_ == simplified_best_score_:
            assert grid.best_index_ == grid_refitted.best_index_
            assert grid.best_params_ == grid_refitted.best_params_
        else:
            assert grid.best_index_ != grid_refitted.best_index_
            assert grid.best_params_ != grid_refitted.best_params_
            assert grid.best_score_ > simplified_best_score_


@ignore_warnings
@pytest.mark.parametrize(
    "param",
    [
        "max_depth",
        None,
    ],
)
@pytest.mark.parametrize(
    "favorability_rank_rule",
    [
        FavorabilityRanker(
            {
                "min_samples_split": (True, 1.0),  # Lower is more favorable
                "max_depth": (True, 1.0),  # Lower is more favorable
            }
        ),
        FavorabilityRanker(
            {
                "min_samples_split": ([2, 4, 6, 8], 1.0),  # String-based rule
                "max_depth": (True, 1.0),  # Lower is more favorable
            }
        ),
        "Not_a_rule",
    ],
)
@pytest.mark.parametrize(
    "scoring,score_slice_rule",
    [
        ("roc_auc", StandardErrorSlicer(sigma=1)),
        ("roc_auc", WilcoxonSlicer(alpha=0.01)),
        ("roc_auc", PercentileSlicer(eta=0.68)),
        ("roc_auc", FixedWindowSlicer(min_cut=0.96, max_cut=0.97)),
        ("roc_auc", "Not_a_rule"),
        ("neg_log_loss", StandardErrorSlicer(sigma=1)),
        ("neg_log_loss", WilcoxonSlicer(alpha=0.01)),
        ("neg_log_loss", PercentileSlicer(eta=0.68)),
        (
            "neg_log_loss",
            pytest.param(
                FixedWindowSlicer(min_cut=0.96, max_cut=0.97), marks=pytest.mark.xfail
            ),
        ),
    ],
)
def test_promote_successive_halving(
    param, scoring, score_slice_rule, favorability_rank_rule
):
    """
    A function that tests the promote function using HalvingRandomSearchCV by
    comparing the results of a refitted search object to those of a non-refitted
    search object.
    """

    X, y = make_classification(n_samples=350, n_features=16, random_state=42)

    # Instantiate the classifier
    clf = RandomForestClassifier(random_state=0)

    # Define parameter distributions
    if param == "max_depth":
        param_distributions = {
            "max_depth": [3, None],
            "min_samples_split": [2, 4, 6, 8],
        }
    else:
        param_distributions = {"min_samples_split": [2, 4, 6, 8]}

    # Instantiate a non-refitted HalvingRandomSearchCV object for comparison
    search = HalvingRandomSearchCV(
        clf,
        param_distributions,
        resource="n_estimators",
        max_resources=10,
        scoring=scoring,
        random_state=0,
    )
    search.fit(X, y)

    score_slice_rule = (
        score_slice_rule.values[0]
        if isinstance(score_slice_rule, ParameterSet)
        else score_slice_rule
    )
    favorability_rank_rule = (
        favorability_rank_rule.values[0]
        if isinstance(favorability_rank_rule, ParameterSet)
        else favorability_rank_rule
    )

    if score_slice_rule == "Not_a_rule" or favorability_rank_rule == "Not_a_rule":
        with pytest.raises(TypeError):
            promote(score_slice_rule, favorability_rank_rule)
        return

    # Instantiate a refitted HalvingRandomSearchCV object
    search_simplified = HalvingRandomSearchCV(
        clf,
        param_distributions,
        resource="n_estimators",
        max_resources=10,
        scoring=scoring,
        refit=promote(score_slice_rule, favorability_rank_rule),
        random_state=0,
    )

    # If the cv results were not all NaN, then we can test the refit callable
    if not np.isnan(search.fit(X, y).cv_results_["mean_test_score"]).all():
        if param and param not in favorability_rank_rule.favorability_rules:
            with pytest.raises(ValueError):
                search_simplified.fit(X, y)
        else:
            search_simplified.fit(X, y)
        simplified_best_score_ = search_simplified.cv_results_["mean_test_score"][
            search_simplified.best_index_
        ]
        # Ensure that if the refit callable promoted a lower scoring model,
        # it was only because it was a more favorable model.
        if abs(search.best_score_) > abs(simplified_best_score_):
            assert search.best_index_ != search_simplified.best_index_
            if param:
                assert (
                    search.best_params_[param] > search_simplified.best_params_[param]
                )
        elif search.best_score_ == simplified_best_score_:
            assert search.best_index_ == search_simplified.best_index_
            assert search.best_params_ == search_simplified.best_params_
        else:
            assert search.best_index_ != search_simplified.best_index_
            assert search.best_params_ != search_simplified.best_params_
            assert search.best_score_ > simplified_best_score_


def test_score_cut_model_selector_tied_ranks(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]

    def mock_favorability_ranker(params):
        favorability_ranks = {
            "rbf": 1,  # assuming 'rbf' kernel is more favorable
            "poly": 2,  # 'poly' kernel is less favorable
        }
        default_rank = len(favorability_ranks) + 1
        ranks = []
        for p in params:
            if "kernel" in p:
                rank = favorability_ranks.get(p["kernel"], default_rank)
            else:
                rank = default_rank
            ranks.append(rank)
        return ranks

    ss = ScoreCutModelSelector(cv_results)
    ss.fit(FixedWindowSlicer(min_cut=0.90, max_cut=0.94))
    best_index = ss.transform(mock_favorability_ranker)

    assert best_index in [0, 2, 3]


def test_standard_error_slicer(generate_fit_params):
    # Test that the StandardErrorSlicer function returns the correct score_slice_rule
    assert pytest.approx(
        StandardErrorSlicer(sigma=1).__call__(**generate_fit_params), rel=1e-2
    ) == (
        0.9243126424613448,
        0.9923540242053219,
    )

    assert StandardErrorSlicer(sigma=1).__repr__() == "StandardErrorSlicer(sigma=1)"

    # Test that the StandardErrorSlicer function raises a ValueError
    with pytest.raises(ValueError):
        StandardErrorSlicer(sigma=-1)


def test_signed_rank_slicer(generate_fit_params):
    # Test that the WilcoxonSlicer function returns the correct score_slice_rule
    assert pytest.approx(
        WilcoxonSlicer(alpha=0.01).__call__(**generate_fit_params), rel=1e-2
    ) == (
        0.9583333333333334,
        0.9583333333333334,
    )

    assert (
        WilcoxonSlicer(alpha=0.01).__repr__()
        == "WilcoxonSlicer(alpha=0.01, alternative=two-sided, zero_method=zsplit)"
    )

    # Test that the WilcoxonSlicer function raises a ValueError if alpha is not
    # between 0 and 1
    with pytest.raises(ValueError):
        WilcoxonSlicer(alpha=-1)

    # Test that the WilcoxonSlicer function raises a ValueError if the number of
    # folds is less than 3
    with pytest.raises(ValueError):
        generate_mod_fit_params = generate_fit_params.copy()
        generate_mod_fit_params.update({"n_folds": 2})
        WilcoxonSlicer(alpha=0.01)(**generate_mod_fit_params)

    # Select rows 0, 1, and 5 from the score grid
    score_grid = generate_fit_params["score_grid"][[0, 1, 5]]
    cv_means = np.mean(score_grid, axis=1)
    best_score_idx = 0
    lowest_score_idx = 2
    n_folds = 3

    with pytest.warns(UserWarning):
        assert WilcoxonSlicer(alpha=0.5)(
            score_grid=score_grid,
            cv_means=cv_means,
            best_score_idx=best_score_idx,
            lowest_score_idx=lowest_score_idx,
            n_folds=n_folds,
        )


def test_percentile_rank_slicer(generate_fit_params):
    # Test that the PercentileSlicer function returns the correct score_slice_rule
    assert pytest.approx(
        PercentileSlicer(eta=0.68).__call__(**generate_fit_params), rel=1e-2
    ) == (0.955, 1.0)

    assert PercentileSlicer(eta=0.68).__repr__() == "PercentileSlicer(eta=0.68)"

    # Test that the PercentileSlicer function raises a ValueError
    with pytest.raises(ValueError):
        PercentileSlicer(eta=-1)


def test_fixed_window_slicer(generate_fit_params):
    # Test that the FixedWindowSlicer function returns the correct score_slice_rule
    assert FixedWindowSlicer(min_cut=0.80, max_cut=0.91).__call__(
        **generate_fit_params
    ) == (
        0.8,
        0.91,
    )

    # No min_cut
    assert FixedWindowSlicer(max_cut=0.91).__call__(**generate_fit_params) == (
        None,
        0.91,
    )

    # No max_cut
    assert FixedWindowSlicer(min_cut=0.80).__call__(**generate_fit_params) == (
        0.8,
        None,
    )

    assert (
        FixedWindowSlicer(min_cut=0.80, max_cut=0.91).__repr__()
        == "FixedWindowSlicer(min_cut=0.8, max_cut=0.91)"
    )

    # Test that the FixedWindowSlicer function raises a ValueError
    with pytest.raises(ValueError):
        FixedWindowSlicer(min_cut=0.99, max_cut=0.92)


def test_favorability_ranker():
    rng = 42

    ranker = FavorabilityRanker(
        {
            "param1": (True, 1.0),
            "param2": (["low", "medium", "high"], 1.0),
            "param3": ("mean", 1.0),
        },
        seed=rng,
    )

    params = [
        {"param1": 10, "param2": "low", "param3": 0.5},
        {"param1": 5, "param2": "medium", "param3": 1.5},
        {"param1": 1, "param2": "high", "param3": 2.5},
    ]

    assert ranker(params) == [3, 2, 1]

    def param_generator(seed_or_rng):
        if isinstance(seed_or_rng, int):
            rng = np.random.default_rng(seed_or_rng)
        else:
            rng = seed_or_rng
        return rng.choice([10, 5, 1])

    params_with_callable = [
        {
            "param1": lambda rng: param_generator(rng),
            "param2": "high",
            "param3": scipy.stats.norm(loc=0, scale=1),
        },
        {
            "param1": lambda rng: param_generator(rng),
            "param2": "medium",
            "param3": scipy.stats.norm(loc=1, scale=1),
        },
        {
            "param1": lambda rng: param_generator(rng),
            "param2": "low",
            "param3": scipy.stats.norm(loc=2, scale=1),
        },
    ]

    expected_ranks = [3, 2, 1]
    assert ranker(params_with_callable) == expected_ranks

    assert (
        repr(ranker)
        == "FavorabilityRanker({'param1': (True, 1.0), 'param2': (['low', 'medium',"
        " 'high'], 1.0), 'param3': ('mean', 1.0)})"
    )


def test_favorability_ranker_validation():
    # invalid hyperparameter type (not a string)
    with pytest.raises(TypeError):
        FavorabilityRanker({1: (True, 1.0)})

    # invalid rule type (not a tuple)
    with pytest.raises(TypeError):
        FavorabilityRanker({"param": "not_a_tuple"})

    # invalid first element in rule (not a bool, list, or string)
    with pytest.raises(TypeError):
        FavorabilityRanker({"param": (1, 1.0)})

    # invalid second element in rule (not a float)
    with pytest.raises(TypeError):
        FavorabilityRanker({"param": (True, "not_a_float")})


def test_favorability_ranker_process_parameter_values():
    ranker = FavorabilityRanker(
        {
            "param_continuous_mean": ("mean", 1.0),
            "param_continuous_median": ("median", 1.0),
            "param_continuous_percentile": ("percentile_75", 1.0),
        },
        seed=42,
    )

    params = {
        "param_continuous_mean": scipy.stats.norm(loc=0, scale=1),
        "param_continuous_median": scipy.stats.norm(loc=0, scale=1),
        "param_continuous_percentile": scipy.stats.norm(loc=0, scale=1),
    }

    # expects: mean, median, and 75th percentile of a normal dist
    expected_values = {
        "param_continuous_mean": 0.0,
        "param_continuous_median": 0.0,
        "param_continuous_percentile": scipy.stats.norm(loc=0, scale=1).ppf(0.75),
    }

    processed_values = {}
    for key, value in params.items():
        processed_res = ranker._process_parameter_values(
            value, ranker.favorability_rules[key]
        )
        processed_values[key] = processed_res

    for key in processed_values:
        assert np.isclose(processed_values[key], expected_values[key]), (
            f"Processed value {processed_values[key]} does not match expected value"
            f" {expected_values[key]}"
        )

    # test case where `distribution_property` is not a supported one
    with pytest.raises(ValueError):
        ranker._process_parameter_values(scipy.stats.norm(loc=0, scale=1), "foo")


def test_favorability_ranker_with_distribution_handling_corrected():
    favorability_rules = {
        "param1": ("mean", 1.0),
        "param2": (["low", "medium", "high"], 1.0),
    }

    params = {
        "param1": scipy.stats.norm(loc=0, scale=1),
        "param2": "medium",
    }

    ranker = FavorabilityRanker(favorability_rules, seed=42)

    ranks = ranker([params])
    assert isinstance(ranks, list), "Expected output to be a list"
    assert len(ranks) == 1, "Expected a single rank output for a single parameter set"


def test_favorability_ranker_with_categorical_combinations():
    favorability_rules = {
        "param2": (["a", "b", "c"], 1.0),
    }

    # Correcting the test to reflect how FavorabilityRanker expects its input
    params = [
        {"param2": value}  # Each value is treated as a separate parameter set
        for value in ["a", "b", "c"]
    ]

    ranker = FavorabilityRanker(favorability_rules, seed=42)

    # Directly use the list of parameter dictionaries without additional wrapping
    ranks = ranker(params)

    assert isinstance(ranks, list), "Ranks should be a list"
    assert len(ranks) == len(params), "Should produce a rank for each parameter set"


def test_favorability_ranker_simple_warning_emission():
    favorability_rules = {
        "param_simple": ("mean", 1.0),
    }

    params_simple = {
        "param_simple": scipy.stats.norm(loc=0, scale=1),
    }

    ranker = FavorabilityRanker(favorability_rules, seed=None)

    with pytest.warns(UserWarning, match="A seed value was not set"):
        ranker(params_simple)


def is_user_warning(warning_record):
    """Check if the warning record is a UserWarning."""
    return issubclass(warning_record.category, UserWarning)
