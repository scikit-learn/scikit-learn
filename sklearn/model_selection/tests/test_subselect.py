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
from sklearn.metrics import make_scorer
from sklearn.pipeline import Pipeline


def test_subselector():
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

    # Test that the _apply_thresh method returns the correct index
    assert ss._apply_thresh("C", 0.93, 0.96) == 1

    # Test that the fit method returns the correct scores
    assert ss.fit(by_standard_error(sigma=1)) == (
        0.9243126424613448,
        0.9923540242053219,
    )

    # Test that the transform method returns the correct model
    assert ss.transform("C") == 1


@ignore_warnings
@pytest.mark.parametrize(
    "param",
    [
        "reduce_dim__n_components",
        "poly__degree",
    ],
)
@pytest.mark.parametrize(
    "scoring",
    ["roc_auc", "neg_log_loss"],
)
@pytest.mark.parametrize(  # iterate over extra liberal thresholds
    "rule",
    [
        by_standard_error(sigma=1),
        by_signed_rank(alpha=0.01),
        by_percentile_rank(eta=0.68),
        by_fixed_window(lower_bound=0.9, upper_bound=0.95),
        pytest.mark.xfail("Not_a_rule"),
    ],
)
@pytest.mark.parametrize(
    "search_cv",
    [GridSearchCV, RandomizedSearchCV],
)
def test_refit_callable_constrain(param, scoring, rule, search_cv):
    """
    A function tests `refit=callable` interface where the callable is the `simplify`
    method of the `Refitter` refit class that returnsthe most parsimonious,
    highest-performing model.
    """

    X, y = make_classification(n_samples=350, n_features=16, random_state=42)

    # Instantiate a pipeline with parameter grid representing different levels of
    # complexity
    clf = LinearSVC(random_state=42)
    if param == "poly__degree":
        from sklearn.preprocessing import PolynomialFeatures

        param_grid = {"poly__degree": [2, 4, 6]}
        pipe = Pipeline([("poly", PolynomialFeatures()), ("classify", clf)])
    elif param == "reduce_dim__n_components":
        from sklearn.decomposition import PCA

        param_grid = {"reduce_dim__n_components": [4, 8, 12]}
        pipe = Pipeline([("reduce_dim", PCA(random_state=42)), ("classify", clf)])
    else:
        raise NotImplementedError(f"{param} not recognized.")

    scoring = make_scorer(scoring, greater_is_better=True)

    # Instantiate a refitted grid search object using `ModelElection`
    grid_simplified = search_cv(
        pipe,
        param_grid,
        scoring=scoring,
        refit=constrain(param, rule),
    )

    # Instantiate a non-refitted grid search object for comparison
    grid = search_cv(pipe, param_grid, scoring=scoring, n_jobs=-1)
    grid.fit(X, y)

    # If the cv results were not all NaN, then we can test the refit callable
    if not np.isnan(grid.fit(X, y).cv_results_["split0_test_score"]).all():
        grid_simplified.fit(X, y)
        best_score_ = grid_simplified.cv_results_["mean_test_score"][
            grid_simplified.best_index_
        ]

        # Ensure that if the `razors` refit callable chose a lower scoring model, that
        # it was because it was a simpler model.
        if abs(grid.best_score_) > abs(best_score_):
            assert (
                getattr(grid, "best_params_")[param]
                >= getattr(grid_simplified, "best_params_")[param]
            )
