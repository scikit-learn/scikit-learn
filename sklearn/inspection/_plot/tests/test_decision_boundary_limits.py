import pytest

from sklearn.datasets import make_classification
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.linear_model import LogisticRegression

# Test data
X, y = make_classification(
    n_informative=1,
    n_redundant=1,
    n_clusters_per_class=1,
    n_features=2,
    random_state=42,
)


@pytest.fixture(scope="module")
def fitted_clf():
    return LogisticRegression().fit(X, y)


class TestXlimYlimValidation:
    """Test input validation for xlim and ylim parameters."""

    @pytest.mark.parametrize(
        "xlim, ylim, expected_error",
        [
            # xlim validation errors
            ([1, 2], None, "xlim must be a tuple of \\(min, max\\) with min < max"),
            ((1,), None, "xlim must be a tuple of \\(min, max\\) with min < max"),
            ((1, 2, 3), None, "xlim must be a tuple of \\(min, max\\) with min < max"),
            ((2, 1), None, "xlim must be a tuple of \\(min, max\\) with min < max"),
            ((1, 1), None, "xlim must be a tuple of \\(min, max\\) with min < max"),
            # ylim validation errors
            (None, [1, 2], "ylim must be a tuple of \\(min, max\\) with min < max"),
            (None, (1,), "ylim must be a tuple of \\(min, max\\) with min < max"),
            (None, (1, 2, 3), "ylim must be a tuple of \\(min, max\\) with min < max"),
            (None, (2, 1), "ylim must be a tuple of \\(min, max\\) with min < max"),
            (None, (1, 1), "ylim must be a tuple of \\(min, max\\) with min < max"),
            # Both invalid
            ((2, 1), (2, 1), "xlim must be a tuple of \\(min, max\\) with min < max"),
        ],
    )
    def test_invalid_xlim_ylim(self, pyplot, fitted_clf, xlim, ylim, expected_error):
        """Test that invalid xlim/ylim values raise ValueError."""
        with pytest.raises(ValueError, match=expected_error):
            DecisionBoundaryDisplay.from_estimator(fitted_clf, X, xlim=xlim, ylim=ylim)

    def test_valid_xlim_ylim(self, pyplot, fitted_clf):
        """Test that valid xlim/ylim values don't raise errors."""
        # These should not raise
        DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=(-5, 5), ylim=(-3, 3)
        )
        DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=(-1.5, 1.5), ylim=None
        )
        DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=None, ylim=(-2.5, 2.5)
        )


class TestXlimYlimFunctionality:
    """Test the functional behavior of xlim and ylim parameters."""

    def test_xlim_ylim_override_eps(self, pyplot, fitted_clf):
        """Test that xlim/ylim override eps-based calculations."""
        xlim = (-10, 10)
        ylim = (-8, 8)
        eps = 2.0  # Should be ignored when xlim/ylim are provided

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim, eps=eps, grid_resolution=5
        )

        # Check that the mesh grid uses exact xlim/ylim values
        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    def test_xlim_only(self, pyplot, fitted_clf):
        """Test behavior when only xlim is provided."""
        xlim = (-5, 5)
        eps = 1.5

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, eps=eps, grid_resolution=5
        )

        # xlim should be used exactly
        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])

        # ylim should be calculated from data + eps
        x1 = X[:, 1]
        expected_y_min = x1.min() - eps
        expected_y_max = x1.max() + eps
        assert disp.xx1.min() == pytest.approx(expected_y_min)
        assert disp.xx1.max() == pytest.approx(expected_y_max)

    def test_ylim_only(self, pyplot, fitted_clf):
        """Test behavior when only ylim is provided."""
        ylim = (-3, 3)
        eps = 1.5

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, ylim=ylim, eps=eps, grid_resolution=5
        )

        # ylim should be used exactly
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

        # xlim should be calculated from data + eps
        x0 = X[:, 0]
        expected_x_min = x0.min() - eps
        expected_x_max = x0.max() + eps
        assert disp.xx0.min() == pytest.approx(expected_x_min)
        assert disp.xx0.max() == pytest.approx(expected_x_max)

    def test_neither_xlim_nor_ylim(self, pyplot, fitted_clf):
        """Test default behavior when neither xlim nor ylim are provided."""
        eps = 2.0

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, eps=eps, grid_resolution=5
        )

        # Both should be calculated from data + eps
        x0, x1 = X[:, 0], X[:, 1]
        expected_x_min = x0.min() - eps
        expected_x_max = x0.max() + eps
        expected_y_min = x1.min() - eps
        expected_y_max = x1.max() + eps

        assert disp.xx0.min() == pytest.approx(expected_x_min)
        assert disp.xx0.max() == pytest.approx(expected_x_max)
        assert disp.xx1.min() == pytest.approx(expected_y_min)
        assert disp.xx1.max() == pytest.approx(expected_y_max)

    @pytest.mark.parametrize("grid_resolution", [5, 10, 20])
    def test_xlim_ylim_with_different_resolutions(
        self, pyplot, fitted_clf, grid_resolution
    ):
        """Test that xlim/ylim work correctly with different grid resolutions."""
        xlim = (-4, 4)
        ylim = (-3, 3)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim, grid_resolution=grid_resolution
        )

        # Check grid shape
        assert disp.xx0.shape == (grid_resolution, grid_resolution)
        assert disp.xx1.shape == (grid_resolution, grid_resolution)

        # Check bounds
        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    @pytest.mark.parametrize(
        "response_method", ["predict", "predict_proba", "decision_function"]
    )
    def test_xlim_ylim_with_different_response_methods(
        self, pyplot, fitted_clf, response_method
    ):
        """Test that xlim/ylim work with different response methods."""
        xlim = (-6, 6)
        ylim = (-4, 4)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim, response_method=response_method
        )

        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    @pytest.mark.parametrize("plot_method", ["contourf", "contour", "pcolormesh"])
    def test_xlim_ylim_with_different_plot_methods(
        self, pyplot, fitted_clf, plot_method
    ):
        """Test that xlim/ylim work with different plot methods."""
        xlim = (-2, 2)
        ylim = (-1.5, 1.5)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim, plot_method=plot_method
        )

        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    def test_xlim_ylim_with_negative_values(self, pyplot, fitted_clf):
        """Test xlim/ylim with negative values."""
        xlim = (-10, -5)
        ylim = (-8, -2)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim
        )

        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    def test_xlim_ylim_with_large_values(self, pyplot, fitted_clf):
        """Test xlim/ylim with large values."""
        xlim = (100, 200)
        ylim = (500, 1000)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim
        )

        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    def test_xlim_ylim_response_shape_consistency(self, pyplot, fitted_clf):
        """Test that response array shape is consistent with xlim/ylim bounds."""
        xlim = (-3, 3)
        ylim = (-2, 2)
        grid_resolution = 10

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim, grid_resolution=grid_resolution
        )

        # Response should have the same shape as the meshgrid
        assert disp.response.shape == (grid_resolution, grid_resolution)
        assert disp.response.shape == disp.xx0.shape
        assert disp.response.shape == disp.xx1.shape

    def test_xlim_ylim_preserves_other_parameters(self, pyplot, fitted_clf):
        """Test that xlim/ylim don't interfere with other parameters."""
        xlim = (-4, 4)
        ylim = (-3, 3)

        # Test with various other parameters
        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf,
            X,
            xlim=xlim,
            ylim=ylim,
            grid_resolution=15,
            plot_method="contour",
            response_method="predict_proba",
        )

        # Verify xlim/ylim are respected
        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

        # Verify other parameters work
        assert disp.response.shape == (15, 15)


class TestXlimYlimEdgeCases:
    """Test edge cases and special scenarios."""

    def test_xlim_ylim_with_zero_values(self, pyplot, fitted_clf):
        """Test xlim/ylim with zero values."""
        xlim = (-1, 0)
        ylim = (0, 1)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim
        )

        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    def test_xlim_ylim_with_very_small_range(self, pyplot, fitted_clf):
        """Test xlim/ylim with very small ranges."""
        xlim = (0.0, 0.001)
        ylim = (-0.0005, 0.0005)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim
        )

        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])

    def test_xlim_ylim_with_floating_point_precision(self, pyplot, fitted_clf):
        """Test xlim/ylim with floating point precision values."""
        xlim = (-1.23456789, 1.23456789)
        ylim = (-0.98765432, 0.98765432)

        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim
        )

        assert disp.xx0.min() == pytest.approx(xlim[0])
        assert disp.xx0.max() == pytest.approx(xlim[1])
        assert disp.xx1.min() == pytest.approx(ylim[0])
        assert disp.xx1.max() == pytest.approx(ylim[1])


class TestXlimYlimDocstringExamples:
    """Test examples that would be in documentation."""

    def test_basic_usage_example(self, pyplot, fitted_clf):
        """Test basic usage example for documentation."""
        # Example: Setting custom plot boundaries
        disp = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=(-3, 3), ylim=(-2, 2), grid_resolution=50
        )

        # Verify the plot uses the specified boundaries
        assert disp.xx0.min() == pytest.approx(-3)
        assert disp.xx0.max() == pytest.approx(3)
        assert disp.xx1.min() == pytest.approx(-2)
        assert disp.xx1.max() == pytest.approx(2)

    def test_comparison_with_without_limits(self, pyplot, fitted_clf):
        """Test comparison between using limits vs not using them."""
        eps = 1.0

        # Without limits (default behavior)
        disp_default = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, eps=eps, grid_resolution=10
        )

        # With custom limits
        xlim = (-5, 5)
        ylim = (-4, 4)
        disp_custom = DecisionBoundaryDisplay.from_estimator(
            fitted_clf, X, xlim=xlim, ylim=ylim, grid_resolution=10
        )

        # They should have different bounds
        assert disp_default.xx0.min() != disp_custom.xx0.min()
        assert disp_default.xx0.max() != disp_custom.xx0.max()
        assert disp_default.xx1.min() != disp_custom.xx1.min()
        assert disp_default.xx1.max() != disp_custom.xx1.max()

        # Custom limits should match exactly
        assert disp_custom.xx0.min() == pytest.approx(xlim[0])
        assert disp_custom.xx0.max() == pytest.approx(xlim[1])
        assert disp_custom.xx1.min() == pytest.approx(ylim[0])
        assert disp_custom.xx1.max() == pytest.approx(ylim[1])
