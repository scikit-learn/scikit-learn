import matplotlib.pyplot as plt
import pytest

from sklearn.datasets import make_classification
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.linear_model import LogisticRegression


class TestDecisionBoundaryDisplayLimits:
    """Test xlim and ylim params for DecisionBoundaryDisplay.from_estimator"""
    def setup_method(self):
        """Set up test data and classifier"""
        X, y = make_classification(
            n_samples=100, n_features=2, n_redundant=0, n_informative=2,
            random_state=42, n_clusters_per_class=1
        )
        self.X = X
        self.y = y

        self.classifier = LogisticRegression(random_state=42)
        self.classifier.fit(X, y)

    def test_xlim_ylim_basic_functionality(self):
        """Test that custom xlim and ylim are properly applied"""
        custom_xlim = (-2, 2)
        custom_ylim = (-1.5, 1.5)

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            response_method="predict",
        )

        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)
        assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
        assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)

        plt.close("all")

    def test_xlim_only(self):
        """Test providing only xlim parameter"""
        custom_xlim = (-3, 3)

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            response_method="predict",
        )

        # xlim should be custom
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)

        # ylim should be auto-calculated (data range + eps)
        expected_y_min = self.X[:, 1].min() - 1.0
        expected_y_max = self.X[:, 1].max() + 1.0
        assert display.xx1.min() == pytest.approx(expected_y_min, abs=1e-10)
        assert display.xx1.max() == pytest.approx(expected_y_max, abs=1e-10)

        plt.close("all")

    def test_ylim_only(self):
        """Test providing only ylim parameter"""
        custom_ylim = (-2, 2)

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            ylim=custom_ylim,
            response_method="predict",
        )

        # ylim should be custom
        assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
        assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)

        # xlim should be auto-calculated (data range + eps)
        expected_x_min = self.X[:, 0].min() - 1.0
        expected_x_max = self.X[:, 0].max() + 1.0
        assert display.xx0.min() == pytest.approx(expected_x_min, abs=1e-10)
        assert display.xx0.max() == pytest.approx(expected_x_max, abs=1e-10)

        plt.close("all")

    def test_xlim_validation_not_tuple(self):
        """Test xlim validation: must be tuple"""
        with pytest.raises(
            ValueError,
            match="xlim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                xlim=[0, 1],
                response_method="predict",
            )

    def test_ylim_validation_not_tuple(self):
        """Test ylim validation: must be tuple"""
        with pytest.raises(
            ValueError,
            match="ylim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                ylim=[0, 1],
                response_method="predict",
            )

    def test_xlim_validation_wrong_length(self):
        """Test xlim validation: must have exactly 2 values"""
        with pytest.raises(
            ValueError,
            match="xlim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                xlim=(0, 1, 2),
                response_method="predict",
            )

    def test_ylim_validation_wrong_length(self):
        """Test ylim validation: must have exactly 2 values"""
        with pytest.raises(
            ValueError,
            match="ylim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                ylim=(0,),
                response_method="predict",
            )

    def test_xlim_validation_min_greater_than_max(self):
        """Test xlim validation: min must be less than max"""
        with pytest.raises(
            ValueError,
            match="xlim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                xlim=(2, 1),
                response_method="predict",
            )

    def test_ylim_validation_min_greater_than_max(self):
        """Test ylim validation: min must be less than max"""
        with pytest.raises(
            ValueError,
            match="ylim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                ylim=(3, 1),
                response_method="predict",
            )

    def test_xlim_validation_min_equal_to_max(self):
        """Test xlim validation: min must be strictly less than max"""
        with pytest.raises(
            ValueError,
            match="xlim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                xlim=(1.5, 1.5),
                response_method="predict",
            )

    def test_ylim_validation_min_equal_to_max(self):
        """Test ylim validation: min must be strictly less than max"""
        with pytest.raises(
            ValueError,
            match="ylim must be a tuple of \\(min, max\\) with min < max"
        ):
            DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                ylim=(1.5, 1.5),
                response_method="predict",
            )

    def test_limits_with_different_response_methods(self):
        """Test that limits work with different response methods"""
        custom_xlim = (-2, 2)
        custom_ylim = (-1, 1)

        # Test with predict_proba
        display_proba = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            response_method="predict_proba",
        )

        # Test with decision_function
        display_decision = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            response_method="decision_function",
        )

        # Both should have the same grid bounds
        for display in [display_proba, display_decision]:
            assert display.xx0.min() == pytest.approx(
                custom_xlim[0], abs=1e-10
            )
            assert display.xx0.max() == pytest.approx(
                custom_xlim[1], abs=1e-10
            )
            assert display.xx1.min() == pytest.approx(
                custom_ylim[0], abs=1e-10
            )
            assert display.xx1.max() == pytest.approx(
                custom_ylim[1], abs=1e-10
            )

        plt.close("all")

    def test_limits_with_eps_ignored(self):
        """Test that eps is ignored when both xlim and ylim are provided"""
        custom_xlim = (-1, 1)
        custom_ylim = (-1, 1)
        large_eps = 10.0

        # With custom limits, eps should be ignored
        display_with_limits = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            eps=large_eps,
            response_method="predict",
        )

        # Bounds should not be affected by eps
        assert display_with_limits.xx0.min() == pytest.approx(
            custom_xlim[0], abs=1e-10
        )
        assert display_with_limits.xx0.max() == pytest.approx(
            custom_xlim[1], abs=1e-10
        )
        assert display_with_limits.xx1.min() == pytest.approx(
            custom_ylim[0], abs=1e-10
        )
        assert display_with_limits.xx1.max() == pytest.approx(
            custom_ylim[1], abs=1e-10
        )

        plt.close("all")

    def test_limits_with_different_grid_resolution(self):
        """Test that limits work with different grid resolutions"""
        custom_xlim = (-2, 2)
        custom_ylim = (-1, 1)

        for grid_res in [50, 100, 200]:
            display = DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                xlim=custom_xlim,
                ylim=custom_ylim,
                grid_resolution=grid_res,
                response_method="predict",
            )

            # Grid bounds should be the same regardless of resolution
            assert display.xx0.min() == pytest.approx(
                custom_xlim[0], abs=1e-10
            )
            assert display.xx0.max() == pytest.approx(
                custom_xlim[1], abs=1e-10
            )
            assert display.xx1.min() == pytest.approx(
                custom_ylim[0], abs=1e-10
            )
            assert display.xx1.max() == pytest.approx(
                custom_ylim[1], abs=1e-10
            )

            # But grid shape should change
            assert display.xx0.shape == (grid_res, grid_res)
            assert display.xx1.shape == (grid_res, grid_res)

        plt.close("all")

    def test_negative_limits(self):
        """Test that negative limits work correctly"""
        custom_xlim = (-5, -1)
        custom_ylim = (-3, -0.5)

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            response_method="predict",
        )
        # Grid bounds should be the same
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)
        assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
        assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)

        plt.close("all")

    def teardown_method(self):
        """Clean up after each test"""
        plt.close("all")


def test_limits_integration_with_real_plot():
    """Integration test that actually creates and verifies the plot"""
    from sklearn.datasets import make_classification
    from sklearn.svm import SVC

    # Create synthetic data
    X, y = make_classification(
        n_samples=50,
        n_features=2,
        n_redundant=0,
        n_informative=2,
        random_state=42
    )

    classifier = SVC(kernel="rbf", random_state=42)
    classifier.fit(X, y)

    fig, ax = plt.subplots(figsize=(8, 6))

    custom_xlim = (-3, 3)
    custom_ylim = (-2, 2)

    display = DecisionBoundaryDisplay.from_estimator(
        classifier,
        X,
        xlim=custom_xlim,
        ylim=custom_ylim,
        response_method="predict",
        ax=ax,
        alpha=0.8,
    )

    ax.scatter(X[:, 0], X[:, 1], c=y, cmap="viridis", edgecolors="black")

    # Verify the decision boundary grid limits (not the matplotlib axis limits)
    # Decision boundary grid should exactly match our custom limits
    assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
    assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)
    assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
    assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)

    # Verify that the display object was created successfully
    assert display.response is not None
    assert display.xx0.shape == display.xx1.shape
    assert display.response.shape == display.xx0.shape

    plt.close("all")
