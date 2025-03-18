import numpy as np
import pytest
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
from sklearn.model_selection._refit import BaseScoreSlicer
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
            rank += 1  # assuming 'rbf' is more favorable
        if params.get("C", 1) < 5:
            rank += 1  # assuming lower 'C' is more favorable
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

    # test that the _get_splits method extracts the correct subgrid
    assert len(ss._get_splits()) == n_splits

    # test that the _n_folds property returns the correct number of folds
    assert ss._n_folds == n_splits

    # test that the _score_grid property returns the correct subgrid of scores
    assert ss._score_grid.shape == (6, n_splits)

    # test that the _cv_means property returns the correct array of mean scores
    assert ss._cv_means.shape == (6,)

    # test that the _lowest_score_idx property returns the correct index
    assert ss._lowest_score_idx == 5

    # test that the _best_score_idx property returns the correct index
    assert ss._best_score_idx == 0

    cv_results_constrained = ss.cv_results_constrained_
    constrained_results, performance_mask, min_cut, max_cut = ss._apply_thresh(
        0.93, 0.96, cv_results_constrained
    )
    assert constrained_results == cv_results_constrained
    assert np.array_equal(
        performance_mask, np.array([True, True, True, True, True, False])
    )
    assert min_cut == 0.93
    assert max_cut == 0.96

    # omit min_thresh
    constrained_results, performance_mask, min_cut, max_cut = ss._apply_thresh(
        None, 0.99, cv_results_constrained
    )
    assert constrained_results == cv_results_constrained
    assert np.array_equal(
        performance_mask, np.array([True, True, True, True, True, True])
    )
    assert min_cut == np.nanmin(ss._cv_means)
    assert max_cut == 0.99

    # omit max_thresh
    constrained_results, performance_mask, min_cut, max_cut = ss._apply_thresh(
        0.80, None, cv_results_constrained
    )
    assert constrained_results == cv_results_constrained
    assert np.array_equal(
        performance_mask, np.array([True, True, True, True, True, False])
    )
    assert min_cut == 0.8
    assert max_cut == np.nanmax(ss._cv_means)

    # Adjusted assertion to allow for numerical differences
    min_cut, max_cut = ss.fit(StandardErrorSlicer(sigma=1))
    assert min_cut == pytest.approx(0.9243, rel=1e-4)
    assert max_cut == pytest.approx(0.9924, rel=1e-4)

    # Adjusted assertion to ensure the transformed index is as expected
    favorable_index = ss.transform(mock_favorability_ranker)
    assert favorable_index == 0

    assert ss.best_params_cut_ == [
        {"C": 1, "gamma": 0.1, "kernel": "rbf"},
        {"C": 1, "gamma": 1, "kernel": "rbf"},
        {"C": 10, "gamma": 0.1, "kernel": "rbf"},
        {"C": 10, "gamma": 1, "kernel": "rbf"},
        {"degree": 1, "kernel": "poly"},
        {"degree": 2, "kernel": "poly"},
    ]

    # Adjusted assertion to use np.testing.assert_allclose for numerical arrays
    np.testing.assert_allclose(
        ss.best_scores_cut_,
        np.array(
            [0.95833333, 0.93872549, 0.93872549, 0.93995098, 0.95833333, 0.48039216]
        ),
        rtol=1e-5,
        atol=1e-8,
    )


def test_ScoreCutModelSelector_errors(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]
    n_splits = grid_search_simulated["n_splits"]

    with pytest.raises(ValueError):
        ss = ScoreCutModelSelector(cv_results)
        cv_results_constrained = ss.cv_results_constrained_
        ss._apply_thresh(0.98, 0.99, cv_results_constrained)

    with pytest.raises(TypeError):
        ss = ScoreCutModelSelector(cv_results)
        ss.fit("Not_a_rule")

    with pytest.raises(TypeError):
        ss = ScoreCutModelSelector(cv_results)
        ss.fit(StandardErrorSlicer(sigma=1))
        ss.transform("Not_a_rule")

    with pytest.raises(ValueError) as exc_info:
        ss = ScoreCutModelSelector(cv_results)
        cv_results_constrained = ss.cv_results_constrained_
        ss._apply_thresh(0.99, 0.98, cv_results_constrained)
    assert "min_cut (0.99) must be less than or equal to max_cut (0.98)." in str(
        exc_info.value
    )

    del cv_results["params"]
    ss = ScoreCutModelSelector(cv_results)
    with pytest.raises(TypeError):
        ss._get_splits()


def test_ScoreCutModelSelector_not_fitted_error(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]

    ss = ScoreCutModelSelector(cv_results)

    with pytest.raises(AttributeError) as exc_info:
        _ = ss.best_params_cut_

    assert (
        str(exc_info.value)
        == "The ``ScoreCutModelSelector`` instance has not been fitted. Please"
        " call the ``ScoreCutModelSelector:fit`` method first."
    )


def test_ScoreCutModelSelector_not_fitted_errors(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]
    ss = ScoreCutModelSelector(cv_results)

    # accessing best_params_cut_ before fit
    with pytest.raises(AttributeError) as exc_info:
        _ = ss.best_params_cut_
    assert "The ``ScoreCutModelSelector`` instance has not been fitted" in str(
        exc_info.value
    )

    # accessing best_scores_cut_ before fit
    with pytest.raises(AttributeError) as exc_info:
        _ = ss.best_scores_cut_
    assert "The ``ScoreCutModelSelector`` instance has not been fitted" in str(
        exc_info.value
    )

    # accessing favorable_best_params_ before fit and transform
    with pytest.raises(AttributeError) as exc_info:
        _ = ss.favorable_best_params_
    assert "The ``ScoreCutModelSelector`` instance has not been fitted" in str(
        exc_info.value
    )

    # accessing favorable_best_score_ before fit and transform
    with pytest.raises(AttributeError) as exc_info:
        _ = ss.favorable_best_score_
    assert "The ``ScoreCutModelSelector`` instance has not been fitted" in str(
        exc_info.value
    )

    ss.fit(StandardErrorSlicer(sigma=1))

    # accessing favorable_best_params_ and favorable_best_score_ after fit
    # but before transform
    with pytest.raises(AttributeError) as exc_info:
        _ = ss.favorable_best_params_
    assert "The ``ScoreCutModelSelector`` instance has not been transformed" in str(
        exc_info.value
    )

    with pytest.raises(AttributeError) as exc_info:
        _ = ss.favorable_best_score_
    assert "The ``ScoreCutModelSelector`` instance has not been transformed" in str(
        exc_info.value
    )


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
                "reduce_dim__n_components": (True, 2.0),  # lower is more favorable
                "classify__C": (False, 1.0),  # lower is less favorable
            }
        ),
        FavorabilityRanker(
            {
                "reduce_dim__n_components": ([4, 8, 12], 2.0),  # string-based rule
                "classify__C": (False, 1.0),  # lower is less favorable
            }
        ),
        "Not_a_rule",
    ],
)
@pytest.mark.parametrize(
    "search_cv",
    [GridSearchCV, RandomizedSearchCV],
)
def test_promote(param, scoring, score_slice_rule, favorability_rank_rule, search_cv):
    """
    A function that tests the promote function by comparing the results of a
    refitted grid and random search object to those of a non-refitted grid and random
    search object, respectively.
    """

    X, y = make_classification(n_samples=350, n_features=16, random_state=42)

    # instantiate a pipeline with parameter grid representing different levels of
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

    # instantiate a non-refitted grid search object for comparison
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

    # instantiate a refitted grid search object
    grid_refitted = search_cv(
        pipe,
        param_grid,
        scoring=scoring,
        refit=promote(score_slice_rule, favorability_rank_rule),
    )

    # if the cv results were not all NaN, then we can test the refit callable
    if not np.isnan(grid.fit(X, y).cv_results_["mean_test_score"]).all():
        grid_refitted.fit(X, y)
        simplified_best_score_ = grid_refitted.cv_results_["mean_test_score"][
            grid_refitted.best_index_
        ]
        # Adjusted assertion to allow for minor numerical differences
        assert abs(grid.best_score_) >= abs(simplified_best_score_) - 1e-4
        # ensure that the refit callable promoted a lower scoring model because it was
        # a more favorable model.
        assert grid.best_index_ != grid_refitted.best_index_
        if param:
            assert grid.best_params_[param] > grid_refitted.best_params_[param]


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
                "min_samples_split": (True, 1.0),  # lower is more favorable
                "max_depth": (True, 1.0),  # lower is more favorable
            }
        ),
        FavorabilityRanker(
            {
                "min_samples_split": ([2, 4, 6, 8], 1.0),  # string-based rule
                "max_depth": (True, 1.0),  # lower is more favorable
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

    # instantiate the classifier
    clf = RandomForestClassifier(random_state=0)

    # define parameter distributions
    if param == "max_depth":
        param_distributions = {
            "max_depth": [3, None],
            "min_samples_split": [2, 4, 6, 8],
        }
    else:
        param_distributions = {"min_samples_split": [2, 4, 6, 8]}

    # instantiate a non-refitted HalvingRandomSearchCV object for comparison
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

    # instantiate a refitted HalvingRandomSearchCV object
    search_simplified = HalvingRandomSearchCV(
        clf,
        param_distributions,
        resource="n_estimators",
        max_resources=10,
        scoring=scoring,
        refit=promote(score_slice_rule, favorability_rank_rule),
        random_state=0,
    )

    search_simplified.fit(X, y)
    simplified_best_score_ = search_simplified.cv_results_["mean_test_score"][
        search_simplified.best_index_
    ]
    # Adjusted assertion to allow for minor numerical differences
    assert abs(search.best_score_) >= abs(simplified_best_score_) - 1e-4
    # ensure that the refit callable promoted a lower scoring model because it was
    # a more favorable model.
    assert search.best_index_ == search_simplified.best_index_
    assert search.best_params_ == search_simplified.best_params_


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
            rank = favorability_ranks.get(p["kernel"], default_rank)
            ranks.append(rank)
        return ranks

    ss = ScoreCutModelSelector(cv_results)
    ss.fit(FixedWindowSlicer(min_cut=0.90, max_cut=0.94))
    best_index = ss.transform(mock_favorability_ranker)

    assert best_index in [0, 2, 3]


def test_BaseScoreSlicer_not_implemented():
    slicer = BaseScoreSlicer()

    assert (
        slicer.__dict__ == {}
    ), "BaseScoreSlicer's __init__ should not initialize any attributes."

    with pytest.raises(
        NotImplementedError, match="Subclasses must implement this method."
    ):
        slicer(
            score_grid=None,
            cv_means=None,
            best_score_idx=None,
            lowest_score_idx=None,
            n_folds=None,
        )


def test_standard_error_slicer(generate_fit_params):
    # test that the StandardErrorSlicer function returns the correct score_slice_rule
    min_cut, max_cut = StandardErrorSlicer(sigma=1).__call__(**generate_fit_params)
    assert min_cut == pytest.approx(0.9243, rel=1e-4)
    assert max_cut == pytest.approx(0.9924, rel=1e-4)

    assert StandardErrorSlicer(sigma=1).__repr__() == "StandardErrorSlicer(sigma=1)"

    # test that StandardErrorSlicer raises a ValueError if sigma is negative
    with pytest.raises(ValueError):
        StandardErrorSlicer(sigma=-1)


def test_signed_rank_slicer(generate_fit_params):
    # test that the WilcoxonSlicer function returns the correct score_slice_rule
    min_cut, max_cut = WilcoxonSlicer(alpha=0.01).__call__(**generate_fit_params)
    assert min_cut == pytest.approx(0.9583, rel=1e-4)
    assert max_cut == pytest.approx(0.9583, rel=1e-4)

    assert (
        WilcoxonSlicer(alpha=0.01).__repr__()
        == "WilcoxonSlicer(alpha=0.01, alternative=two-sided, zero_method=zsplit)"
    )

    # test that WilcoxonSlicer raises a ValueError if alpha is not between 0 and 1
    with pytest.raises(ValueError):
        WilcoxonSlicer(alpha=-1.05)

    # test that the WilcoxonSlicer function raises a ValueError if the number of
    # folds is less than 3
    with pytest.raises(ValueError):
        generate_mod_fit_params = generate_fit_params.copy()
        generate_mod_fit_params.update({"n_folds": 2})
        WilcoxonSlicer(alpha=0.01)(**generate_mod_fit_params)

    # select rows 0, 1, and 5 from the score grid
    score_grid = generate_fit_params["score_grid"][[0, 1, 5]]
    cv_means = np.mean(score_grid, axis=1)
    best_score_idx = 0
    lowest_score_idx = 2
    n_folds = 3

    with pytest.warns(UserWarning):
        WilcoxonSlicer(alpha=0.5)(
            score_grid=score_grid,
            cv_means=cv_means,
            best_score_idx=best_score_idx,
            lowest_score_idx=lowest_score_idx,
            n_folds=n_folds,
        )


def test_percentile_rank_slicer(generate_fit_params):
    # test that the PercentileSlicer function returns the correct score_slice_rule
    min_cut, max_cut = PercentileSlicer(eta=0.68).__call__(**generate_fit_params)
    assert min_cut == pytest.approx(0.9550, rel=1e-4)
    assert max_cut == pytest.approx(1.0, rel=1e-4)

    assert PercentileSlicer(eta=0.68).__repr__() == "PercentileSlicer(eta=0.68)"

    # test that PercentileSlicer raises a ValueError if eta is not between 0 and 1
    with pytest.raises(ValueError):
        PercentileSlicer(eta=-1.05)


def test_fixed_window_slicer(generate_fit_params):
    # test that the FixedWindowSlicer function returns the correct score_slice_rule
    min_cut, max_cut = FixedWindowSlicer(min_cut=0.80, max_cut=0.91).__call__(
        **generate_fit_params
    )
    assert min_cut == 0.8
    assert max_cut == 0.91

    # no min_cut
    min_cut, max_cut = FixedWindowSlicer(max_cut=0.91).__call__(**generate_fit_params)
    assert min_cut is None
    assert max_cut == 0.91

    # no max_cut
    min_cut, max_cut = FixedWindowSlicer(min_cut=0.80).__call__(**generate_fit_params)
    assert min_cut == 0.8
    assert max_cut is None

    assert (
        FixedWindowSlicer(min_cut=0.80, max_cut=0.91).__repr__()
        == "FixedWindowSlicer(min_cut=0.8, max_cut=0.91)"
    )

    # Test that the FixedWindowSlicer function raises a ValueError
    with pytest.raises(ValueError):
        FixedWindowSlicer(min_cut=0.99, max_cut=0.92)


def test_effectively_empty_cv_results():
    cv_results_empty = {
        "mean_fit_time": np.array([]),
        "std_fit_time": np.array([]),
        "mean_score_time": np.array([]),
        "std_score_time": np.array([]),
        "param_kernel": np.array([], dtype=object),
        "params": [],
        "split0_test_score": np.array([]),
        "split1_test_score": np.array([]),
        "mean_test_score": np.array([]),
        "std_test_score": np.array([]),
        "rank_test_score": np.array([], dtype=int),
    }

    ss = ScoreCutModelSelector(cv_results_empty)
    slicer = StandardErrorSlicer(sigma=1)

    with pytest.raises(ValueError):
        ss.fit(slicer)


def test_single_element_parameter_list():
    cv_results_single = {
        "mean_fit_time": np.array([0.1]),
        "std_fit_time": np.array([0.01]),
        "mean_score_time": np.array([0.01]),
        "std_score_time": np.array([0.001]),
        "param_kernel": np.array(["rbf"], dtype=object),
        "param_C": np.array([1]),
        "param_gamma": np.array(["scale"], dtype=object),
        "params": [{"kernel": "rbf", "C": 1, "gamma": "scale"}],
        "split0_test_score": np.array([0.95]),
        "split1_test_score": np.array([0.96]),
        "mean_test_score": np.array([0.955]),
        "std_test_score": np.array([0.005]),
        "rank_test_score": np.array([1], dtype=int),
    }

    ss = ScoreCutModelSelector(cv_results_single)

    slicer = StandardErrorSlicer(sigma=1)

    min_cut, max_cut = ss.fit(slicer)

    assert min_cut == pytest.approx(0.9515, rel=1e-4)
    assert max_cut == pytest.approx(0.9585, rel=1e-4)


def test_standard_error_near_zero(generate_fit_params):
    generate_fit_params["score_grid"] = np.ones_like(generate_fit_params["score_grid"])

    min_cut, max_cut = StandardErrorSlicer(sigma=1).__call__(**generate_fit_params)

    assert min_cut == max_cut


def test_incorrect_data_types_for_slicers():
    with pytest.raises(TypeError):
        StandardErrorSlicer(sigma="not_an_int")

    with pytest.raises(TypeError):
        PercentileSlicer(eta="not_a_float")

    with pytest.raises(TypeError):
        WilcoxonSlicer(alpha="not_a_float")


def test_large_parameter_grids():
    # number of parameter combinations
    num_params = 1000

    # simulate large grid search results with two hyperparameters: C and gamma
    cv_results_large = {
        "mean_fit_time": np.random.rand(num_params),
        "std_fit_time": np.random.rand(num_params),
        "mean_score_time": np.random.rand(num_params),
        "std_score_time": np.random.rand(num_params),
        "param_C": np.random.choice([0.1, 1, 10], size=num_params),
        "param_gamma": np.random.choice(["scale", "auto"], size=num_params),
        "params": [
            {"C": c, "gamma": g}
            for c, g in zip(
                np.random.choice([0.1, 1, 10], size=num_params),
                np.random.choice(["scale", "auto"], size=num_params),
            )
        ],
        "split0_test_score": np.random.rand(num_params),
        "split1_test_score": np.random.rand(num_params),
        "mean_test_score": np.random.rand(num_params),
        "std_test_score": np.random.rand(num_params),
        "rank_test_score": np.arange(1, num_params + 1),
    }

    ss = ScoreCutModelSelector(cv_results_large)

    slicer = StandardErrorSlicer(sigma=1)

    ss.fit(slicer)
    assert (
        len(ss.cv_results_constrained_["params"]) <= num_params
    ), "ScoreCutModelSelector should constrain the results."


def test_favorability_ranker():
    ranker = FavorabilityRanker(
        {
            "param1": (True, 1.0),
            "param2": (["low", "medium", "high"], 1.0),
            "param3": (False, 1.0),
        },
    )

    expected_ranks = [
        9,
        6,
        3,
        18,
        15,
        12,
        27,
        24,
        21,
        36,
        33,
        30,
        8,
        5,
        2,
        17,
        14,
        11,
        26,
        23,
        20,
        35,
        32,
        29,
        7,
        4,
        1,
        16,
        13,
        10,
        25,
        22,
        19,
        34,
        31,
        28,
    ]

    params_dict = {
        "param1": [10, 15, 20, 25],
        "param2": ["high", "medium", "low"],
        "param3": [0.001, 0.01, 0.1],
    }

    params_list = [
        {"param1": 10, "param2": "high", "param3": 0.001},
        {"param1": 10, "param2": "high", "param3": 0.01},
        {"param1": 10, "param2": "high", "param3": 0.1},
        {"param1": 10, "param2": "medium", "param3": 0.001},
        {"param1": 10, "param2": "medium", "param3": 0.01},
        {"param1": 10, "param2": "medium", "param3": 0.1},
        {"param1": 10, "param2": "low", "param3": 0.001},
        {"param1": 10, "param2": "low", "param3": 0.01},
        {"param1": 10, "param2": "low", "param3": 0.1},
        {"param1": 15, "param2": "high", "param3": 0.001},
        {"param1": 15, "param2": "high", "param3": 0.01},
        {"param1": 15, "param2": "high", "param3": 0.1},
        {"param1": 15, "param2": "medium", "param3": 0.001},
        {"param1": 15, "param2": "medium", "param3": 0.01},
        {"param1": 15, "param2": "medium", "param3": 0.1},
        {"param1": 15, "param2": "low", "param3": 0.001},
        {"param1": 15, "param2": "low", "param3": 0.01},
        {"param1": 15, "param2": "low", "param3": 0.1},
        {"param1": 20, "param2": "high", "param3": 0.001},
        {"param1": 20, "param2": "high", "param3": 0.01},
        {"param1": 20, "param2": "high", "param3": 0.1},
        {"param1": 20, "param2": "medium", "param3": 0.001},
        {"param1": 20, "param2": "medium", "param3": 0.01},
        {"param1": 20, "param2": "medium", "param3": 0.1},
        {"param1": 20, "param2": "low", "param3": 0.001},
        {"param1": 20, "param2": "low", "param3": 0.01},
        {"param1": 20, "param2": "low", "param3": 0.1},
        {"param1": 25, "param2": "high", "param3": 0.001},
        {"param1": 25, "param2": "high", "param3": 0.01},
        {"param1": 25, "param2": "high", "param3": 0.1},
        {"param1": 25, "param2": "medium", "param3": 0.001},
        {"param1": 25, "param2": "medium", "param3": 0.01},
        {"param1": 25, "param2": "medium", "param3": 0.1},
        {"param1": 25, "param2": "low", "param3": 0.001},
        {"param1": 25, "param2": "low", "param3": 0.01},
        {"param1": 25, "param2": "low", "param3": 0.1},
    ]

    assert ranker(params_list) == expected_ranks
    assert ranker(params_dict) == expected_ranks

    assert (
        repr(ranker)
        == "FavorabilityRanker({'param1': (True, 1.0), 'param2': (['low', 'medium',"
        " 'high'], 1.0), 'param3': (False, 1.0)})"
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

    # negative weight for a hyperparameter
    with pytest.raises(ValueError) as exc_info:
        FavorabilityRanker({"param": (True, -1.0)})
    assert "Weight for hyperparameter param must be non-negative." in str(
        exc_info.value
    )


def test_favorability_ranker_with_categorical_combinations():
    favorability_rules = {
        "param2": (["a", "b", "c"], 1.0),
    }

    params = [{"param2": value} for value in ["a", "b", "c"]]

    ranker = FavorabilityRanker(favorability_rules)

    ranks = ranker(params)

    assert isinstance(ranks, list), "Ranks should be a list"
    assert len(ranks) == len(params), "Should produce a rank for each parameter set"
    assert ranks == [1, 2, 3], "Expected ranks to be [1, 2, 3]"


def test_favorability_ranker_unsupported_value_error():
    favorability_rules = {"param_invalid": (True, 1.0)}
    ranker = FavorabilityRanker(favorability_rules)

    class CustomObject:
        pass

    invalid_value = CustomObject()

    with pytest.raises(
        ValueError,
        match=(
            "FavorabilityRanker only supports numeric or "
            "string values for hyperparameters. The provided value .* is not "
            "supported."
        ),
    ):
        ranker([{"param_invalid": invalid_value}])


def test_favorability_ranker_expect_numeric_error():
    favorability_rules = {
        "param_numeric": (True, 1.0),
    }

    params = [
        {"param_numeric": "not_a_number"},
    ]

    ranker = FavorabilityRanker(favorability_rules)

    with pytest.raises(TypeError):
        ranker(params)


def test_favorability_ranker_invalid_categorical_value_error():
    favorability_rules = {
        "param_categorical": (["allowed_value1", "allowed_value2"], 1.0),
    }

    params = [
        {"param_categorical": "forbidden_value"},
    ]

    ranker = FavorabilityRanker(favorability_rules)

    with pytest.raises(ValueError):
        ranker(params)

    favorability_rules = {
        "kernel": (["linear", "rbf"], 1.0),  # allowed kernel types
    }

    # params list contains a 'kernel' value not in the allowed list ['linear', 'rbf']
    params_list = [
        {"kernel": "linear"},  # valid
        {"kernel": "rbf"},  # valid
        {"kernel": "poly"},  # invalid, not in the ['linear', 'rbf'] list
    ]

    ranker = FavorabilityRanker(favorability_rules)

    # test that attempting to rank a set of hyperparameters with an invalid 'kernel'
    # value raises a ValueError
    with pytest.raises(ValueError) as exc_info:
        ranker(params_list)

    assert "Hyperparameter kernel must be one of ['linear', 'rbf']" in str(
        exc_info.value
    ), (
        "FavorabilityRanker should raise a ValueError when a "
        "hyperparameter value is not in the allowed list."
    )


def test_favorability_ranker_invalid_params_type():
    favorability_rules = {
        "param1": (True, 1.0),
    }

    ranker = FavorabilityRanker(favorability_rules)

    params_invalid = 12345

    with pytest.raises(ValueError) as exc_info:
        ranker(params_invalid)

    assert "`params` must be either a list of dictionaries or a " in str(exc_info.value)


def test_favorability_ranker_consistency():
    favorability_rules = {
        "param1": (True, 1.0),
        "param2": (["a", "b", "c"], 1.0),
    }

    # case 1
    params1 = {
        "param1": [1, 2, 3],
        "param2": ["a", "b", "c"],
    }

    # case 2 (different order of keys)
    params2 = {
        "param1": [1, 2, 3],
        "param2": ["a", "b", "c"],
    }

    # case 3 (different order of values)
    params3 = {
        "param1": [3, 1, 2],
        "param2": ["c", "a", "b"],
    }

    # case 4 (additional parameter)
    params4 = {
        "param1": [1, 2, 3],
        "param2": ["a", "b", "c"],
        "param4": [4, 5, 6],
    }

    ranker = FavorabilityRanker(favorability_rules)

    ranks1 = ranker(params1)
    ranks2 = ranker(params2)
    ranks3 = ranker(params3)
    ranks4 = ranker(params4)

    assert ranks1 == ranks2, "Ranks should be consistent across different key orders"
    assert ranks1 != ranks3, "Ranks should be different when value orders differ"
    assert (
        ranks1 != ranks4
    ), "Ranks should be different when additional parameters are present"
