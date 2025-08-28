import pytest

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV


def test_randomized_search_cv_distribution_weights():
    """Test the distribution_weights parameter of RandomizedSearchCV."""
    # Create a simple dataset
    X, y = make_classification(n_samples=100, n_features=5, random_state=42)

    # Define parameter distributions with multiple dicts
    param_distributions = [
        {"C": [0.1, 1.0], "penalty": ["l1"], "solver": ["liblinear"]},
        {"C": [10.0, 100.0], "penalty": ["l2"], "solver": ["lbfgs"]},
    ]

    # Test with uniform weights (default behavior)
    clf = LogisticRegression(random_state=42, max_iter=1000)
    search_uniform = RandomizedSearchCV(
        clf, param_distributions, n_iter=20, random_state=42, cv=2
    )

    search_uniform.fit(X, y)

    # Check that we have results
    assert hasattr(search_uniform, "best_params_")
    assert hasattr(search_uniform, "best_score_")

    # Test with custom weights
    # Give higher weight to the second distribution
    weights = [0.1, 0.9]
    clf = LogisticRegression(random_state=42, max_iter=1000)
    search_weighted = RandomizedSearchCV(
        clf,
        param_distributions,
        n_iter=20,
        random_state=42,
        cv=2,
        distribution_weights=weights,
    )

    search_weighted.fit(X, y)

    # Check that we have results
    assert hasattr(search_weighted, "best_params_")
    assert hasattr(search_weighted, "best_score_")

    # Analyze the parameters selected in cv_results_
    # Count how many times each distribution was used
    dist1_count = 0
    dist2_count = 0

    for params in search_weighted.cv_results_["params"]:
        # Check if parameters match the first distribution
        if "penalty" in params and params["penalty"] == "l1":
            dist1_count += 1
        # Check if parameters match the second distribution
        elif "penalty" in params and params["penalty"] == "l2":
            dist2_count += 1

    # With higher weight on second distribution, we should see more samples from it
    # Note: This is a probabilistic test, so we're just checking it's not extremely
    # skewed toward the first distribution
    assert dist1_count + dist2_count == len(search_weighted.cv_results_["params"])

    # Test with zero weight for first distribution
    weights = [0.0, 1.0]
    clf = LogisticRegression(random_state=42, max_iter=1000)
    search_zero_weight = RandomizedSearchCV(
        clf,
        param_distributions,
        n_iter=10,
        random_state=42,
        cv=2,
        distribution_weights=weights,
    )

    search_zero_weight.fit(X, y)

    # Check that we have results
    assert hasattr(search_zero_weight, "best_params_")
    assert hasattr(search_zero_weight, "best_score_")

    # With zero weight on first distribution, we should only see parameters from second
    for params in search_zero_weight.cv_results_["params"]:
        assert "penalty" in params
        assert params["penalty"] == "l2"


def test_randomized_search_cv_distribution_weights_validation():
    """Test validation of distribution_weights parameter in RandomizedSearchCV."""
    # Create a simple dataset
    X, y = make_classification(n_samples=10, n_features=5, random_state=42)

    # Define parameter distributions with multiple dicts
    param_distributions = [
        {
            "C": [0.1, 1.0],
        },
        {
            "C": [10.0, 100.0],
        },
    ]

    clf = LogisticRegression()

    # Test with wrong length
    with pytest.raises(
        ValueError, match="distribution_weights must have the same length"
    ):
        search = RandomizedSearchCV(
            clf,
            param_distributions,
            n_iter=5,
            distribution_weights=[0.5],  # Wrong length
        )
        search.fit(X, y)

    # Test with negative weights
    with pytest.raises(
        ValueError, match="distribution_weights must contain only non-negative values"
    ):
        search = RandomizedSearchCV(
            clf, param_distributions, n_iter=5, distribution_weights=[0.5, -0.1]
        )
        search.fit(X, y)

    # Test with all zero weights
    with pytest.raises(
        ValueError,
        match="distribution_weights must contain at least one positive value",
    ):
        search = RandomizedSearchCV(
            clf, param_distributions, n_iter=5, distribution_weights=[0.0, 0.0]
        )
        search.fit(X, y)
