import numpy as np
import pytest

from sklearn.utils._testing import ignore_warnings

from sklearn.datasets import make_classification
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import (
    Refitter,
    by_standard_error,
    by_percentile_rank,
    by_signed_rank,
    by_fixed_window,
    constrain,
)

from sklearn.svm import LinearSVC, SVC
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline


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
    ss = Refitter(cv_results)

    yield {
        "score_grid": ss._score_grid,
        "n_folds": n_splits,
        "cv_means": ss._cv_means,
        "best_score_idx": ss._best_score_idx,
        "lowest_score_idx": ss._lowest_score_idx,
    }


def test_refitter_methods(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]
    n_splits = grid_search_simulated["n_splits"]

    ss = Refitter(cv_results)

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

    assert ss._apply_thresh(0.93, 0.96) == 1

    # Omit min_thresh
    assert ss._apply_thresh(None, 0.99) == 5

    # Omit max_thresh
    assert ss._apply_thresh(0.80, None) == 1

    # Test that the fit method returns the correct scores
    assert ss.fit(by_standard_error(sigma=1)) == (
        0.9243126424613448,
        0.9923540242053219,
    )

    # Test that the transform method returns the correct model
    assert ss.transform() == 1


def test_refitter_errors(grid_search_simulated):
    cv_results = grid_search_simulated["cv_results"]
    n_splits = grid_search_simulated["n_splits"]

    with pytest.raises(ValueError):
        ss = Refitter(cv_results)
        assert ss._apply_thresh(0.98, 0.99) == 1

    with pytest.raises(TypeError):
        ss = Refitter(cv_results)
        assert ss.fit("Not_a_rule") == (0.9243126424613448, 0.9923540242053219)

    with pytest.raises(ValueError):
        ss = Refitter(cv_results)
        assert ss.transform() == 1

    del cv_results["params"]
    ss = Refitter(cv_results)
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
    "scoring,rule",
    [
        ("roc_auc", by_standard_error(sigma=1)),
        ("roc_auc", by_signed_rank(alpha=0.01)),
        ("roc_auc", by_percentile_rank(eta=0.68)),
        ("roc_auc", by_fixed_window(min_cut=0.96, max_cut=0.97)),
        ("roc_auc", "Not_a_rule"),
        ("neg_log_loss", by_standard_error(sigma=1)),
        ("neg_log_loss", by_signed_rank(alpha=0.01)),
        ("neg_log_loss", by_percentile_rank(eta=0.68)),
        (
            "neg_log_loss",
            pytest.param(
                by_fixed_window(min_cut=0.96, max_cut=0.97), marks=pytest.mark.xfail
            ),
        ),
    ],
)
@pytest.mark.parametrize(
    "search_cv",
    [GridSearchCV, RandomizedSearchCV],
)
def test_constrain(param, scoring, rule, search_cv):
    """
    A function tests `refit=callable` interface where the callable is the `simplify`
    method of the `Refitter` refit class that returnsthe most parsimonious,
    highest-performing model.
    """

    X, y = make_classification(n_samples=350, n_features=16, random_state=42)

    # Instantiate a pipeline with parameter grid representing different levels of
    # complexity
    clf = LinearSVC(random_state=42)
    if param == "reduce_dim__n_components":
        param_grid = {"reduce_dim__n_components": [4, 8, 12]}
        pipe = Pipeline([("reduce_dim", PCA(random_state=42)), ("classify", clf)])
    else:
        param_grid = {"classify__C": [0.1, 1], "reduce_dim__n_components": [4, 8, 12]}
        pipe = Pipeline(
            [("reduce_dim", PCA(random_state=42)), ("classify", SVC(random_state=42))]
        )

    # Instantiate a non-refitted grid search object for comparison
    grid = search_cv(pipe, param_grid, scoring=scoring, n_jobs=-1)
    grid.fit(X, y)

    # Instantiate a refitted grid search object
    grid_simplified = search_cv(
        pipe,
        param_grid,
        scoring=scoring,
        refit=constrain(rule),
    )

    # If the cv results were not all NaN, then we can test the refit callable
    if not np.isnan(grid.fit(X, y).cv_results_["mean_test_score"]).all():
        if rule == "Not_a_rule":
            with pytest.raises(TypeError):
                grid_simplified.fit(X, y)  # pragma: no cover
        else:
            grid_simplified.fit(X, y)  # pragma: no cover
            simplified_best_score_ = grid_simplified.cv_results_["mean_test_score"][
                grid_simplified.best_index_
            ]  # pragma: no cover
            # Ensure that if the refit callable subselected a lower scoring model,
            # it was because it was only because it was a simpler model.
            if abs(grid.best_score_) > abs(simplified_best_score_):  # pragma: no cover
                assert (
                    grid.best_index_ != grid_simplified.best_index_
                )  # pragma: no cover
                if param:
                    assert (
                        grid.best_params_[param] > grid_simplified.best_params_[param]
                    )  # pragma: no cover
            elif grid.best_score_ == simplified_best_score_:  # pragma: no cover
                assert grid.best_index_ == grid_simplified.best_index_
                assert grid.best_params_ == grid_simplified.best_params_
            else:  # pragma: no cover
                assert (
                    grid.best_index_ != grid_simplified.best_index_
                )  # pragma: no cover
                assert (
                    grid.best_params_ != grid_simplified.best_params_
                )  # pragma: no cover
                assert grid.best_score_ > simplified_best_score_  # pragma: no cover


def test_by_standard_error(generate_fit_params):
    # Test that the by_standard_error function returns the correct rule
    assert pytest.approx(
        by_standard_error(sigma=1).__call__(**generate_fit_params), rel=1e-2
    ) == (
        0.9243126424613448,
        0.9923540242053219,
    )

    assert by_standard_error(sigma=1).__repr__() == "by_standard_error(sigma=1)"

    # Test that the by_standard_error function raises a ValueError
    with pytest.raises(ValueError):
        by_standard_error(sigma=-1)


def test_by_signed_rank(generate_fit_params):
    # Test that the by_signed_rank function returns the correct rule
    assert pytest.approx(
        by_signed_rank(alpha=0.01).__call__(**generate_fit_params), rel=1e-2
    ) == (
        0.9583333333333334,
        0.9583333333333334,
    )

    assert (
        by_signed_rank(alpha=0.01).__repr__()
        == "by_signed_rank(alpha=0.01, alternative=two-sided, zero_method=zsplit)"
    )

    # Test that the by_signed_rank function raises a ValueError if alpha is not
    # between 0 and 1
    with pytest.raises(ValueError):
        by_signed_rank(alpha=-1)

    # Test that the by_signed_rank function raises a ValueError if the number of
    # folds is less than 3
    with pytest.raises(ValueError):
        generate_mod_fit_params = generate_fit_params.copy()
        generate_mod_fit_params.update({"n_folds": 2})
        by_signed_rank(alpha=0.01)(**generate_mod_fit_params)

    # The average performance of all cross-validated models is significantly different
    # from that of the best-performing model

    # Select rows 0, 1, and 5 from the score grid
    score_grid = generate_fit_params["score_grid"][[0, 1, 5]]
    cv_means = np.mean(score_grid, axis=1)
    best_score_idx = 0
    lowest_score_idx = 2
    n_folds = 3

    with pytest.warns(UserWarning):
        assert by_signed_rank(alpha=0.5)(
            score_grid=score_grid,
            cv_means=cv_means,
            best_score_idx=best_score_idx,
            lowest_score_idx=lowest_score_idx,
            n_folds=n_folds,
        )


def test_by_percentile_rank(generate_fit_params):
    # Test that the by_percentile_rank function returns the correct rule
    assert pytest.approx(
        by_percentile_rank(eta=0.68).__call__(**generate_fit_params), rel=1e-2
    ) == (0.955, 1.0)

    assert by_percentile_rank(eta=0.68).__repr__() == "by_percentile_rank(eta=0.68)"

    # Test that the by_percentile_rank function raises a ValueError
    with pytest.raises(ValueError):
        by_percentile_rank(eta=-1)


def test_by_fixed_window(generate_fit_params):
    # Test that the by_fixed_window function returns the correct rule
    assert by_fixed_window(min_cut=0.80, max_cut=0.91).__call__(
        **generate_fit_params
    ) == (
        0.8,
        0.91,
    )

    # No min_cut
    assert by_fixed_window(max_cut=0.91).__call__(**generate_fit_params) == (None, 0.91)

    # No max_cut
    assert by_fixed_window(min_cut=0.80).__call__(**generate_fit_params) == (0.8, None)

    assert (
        by_fixed_window(min_cut=0.80, max_cut=0.91).__repr__()
        == "by_fixed_window(min_cut=0.8, max_cut=0.91)"
    )

    # Test that the by_fixed_window function raises a ValueError
    with pytest.raises(ValueError):
        by_fixed_window(min_cut=0.99, max_cut=0.92)
