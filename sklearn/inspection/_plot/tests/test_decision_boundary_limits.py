import numpy as np
import pytest

from sklearn.datasets import make_classification
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.utils._testing import assert_allclose


class TestDecisionBoundaryDisplayLimits:
    """Test xlim and ylim params for DecisionBoundaryDisplay.from_estimator"""

    def setup_method(self):
        """Set up test data and classifier"""
        X, y = make_classification(
            n_samples=100,
            n_features=2,
            n_redundant=0,
            n_informative=2,
            random_state=42,
            n_clusters_per_class=1,
        )
        self.X = X
        self.y = y
        self.classifier = LogisticRegression(random_state=42)
        self.classifier.fit(X, y)

    def test_xlim_ylim_grid_bounds(self):
        """Test that custom xlim and ylim create correct grid bounds"""
        custom_xlim = (-2, 2)
        custom_ylim = (-1.5, 1.5)
        grid_resolution = 50

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            grid_resolution=grid_resolution,
            response_method="predict",
        )

        # Verify grid bounds match exactly
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)
        assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
        assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)

        # Verify grid shape
        assert display.xx0.shape == (grid_resolution, grid_resolution)
        assert display.response.shape == (grid_resolution, grid_resolution)

    def test_xlim_only_auto_ylim(self):
        """Test xlim with auto-calculated ylim"""
        custom_xlim = (-3, 3)
        eps = 1.0

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            eps=eps,
            response_method="predict",
        )

        # xlim should be custom
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)

        # ylim should be auto-calculated from data
        expected_y_min = self.X[:, 1].min() - eps
        expected_y_max = self.X[:, 1].max() + eps
        assert display.xx1.min() == pytest.approx(expected_y_min, abs=1e-10)
        assert display.xx1.max() == pytest.approx(expected_y_max, abs=1e-10)

    def test_ylim_only_auto_xlim(self):
        """Test ylim with auto-calculated xlim"""
        custom_ylim = (-2, 2)
        eps = 1.0

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            ylim=custom_ylim,
            eps=eps,
            response_method="predict",
        )

        # ylim should be custom
        assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
        assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)

        # xlim should be auto-calculated from data
        expected_x_min = self.X[:, 0].min() - eps
        expected_x_max = self.X[:, 0].max() + eps
        assert display.xx0.min() == pytest.approx(expected_x_min, abs=1e-10)
        assert display.xx0.max() == pytest.approx(expected_x_max, abs=1e-10)

    def test_response_values_match_manual_prediction(self):
        """Test that response values match manual prediction on same grid"""
        custom_xlim = (-1, 1)
        custom_ylim = (-1, 1)
        grid_resolution = 20

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            grid_resolution=grid_resolution,
            response_method="predict",
        )

        # Create expected grid manually
        xx0, xx1 = np.meshgrid(
            np.linspace(custom_xlim[0], custom_xlim[1], grid_resolution),
            np.linspace(custom_ylim[0], custom_ylim[1], grid_resolution),
        )

        # Manual prediction
        grid_points = np.c_[xx0.ravel(), xx1.ravel()]
        expected_response = self.classifier.predict(grid_points)
        expected_response = expected_response.reshape(xx0.shape)

        # Compare with display results
        assert_allclose(display.xx0, xx0)
        assert_allclose(display.xx1, xx1)
        assert_allclose(display.response, expected_response)


    def test_multiclass_with_limits(self):
        """Test multiclass classification with custom limits"""
        X, y = make_classification(
            n_samples=150,
            n_features=2,
            n_classes=3,
            n_informative=2,
            n_redundant=0,
            n_clusters_per_class=1,
            random_state=42,
        )

        clf = LogisticRegression(random_state=42)
        clf.fit(X, y)

        custom_xlim = (-2, 2)
        custom_ylim = (-2, 2)
        grid_resolution = 30

        display = DecisionBoundaryDisplay.from_estimator(
            clf,
            X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            grid_resolution=grid_resolution,
            response_method="predict",
        )

        # Check grid shape
        assert len(np.unique(display.response)) <= 3
        assert display.response.shape == (grid_resolution, grid_resolution)

        # Check grid bounds
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)

    def test_predict_proba_with_limits(self):
        """Test predict_proba response method with custom limits"""
        custom_xlim = (-1, 1)
        custom_ylim = (-1, 1)
        grid_resolution = 25

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            grid_resolution=grid_resolution,
            response_method="predict_proba",
        )

        # Resulting probability values should be between 0 and 1
        assert np.all(display.response >= 0)
        assert np.all(display.response <= 1)
        assert display.response.shape == (grid_resolution, grid_resolution)

    def test_svm_decision_function_with_limits(self):
        """Test SVM decision function with custom limits"""
        svm_clf = SVC(kernel="rbf", random_state=42)
        svm_clf.fit(self.X, self.y)

        custom_xlim = (-2, 2)
        custom_ylim = (-1, 1)
        grid_resolution = 40

        display = DecisionBoundaryDisplay.from_estimator(
            svm_clf,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            grid_resolution=grid_resolution,
            response_method="decision_function",
        )

        # Decision function returns continuous values
        assert display.response.shape == (grid_resolution, grid_resolution)
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)

    def test_different_grid_resolutions(self):
        """Test that different grid resolutions work with custom limits"""
        custom_xlim = (-1, 1)
        custom_ylim = (-1, 1)

        for grid_res in [10, 50, 100]:
            display = DecisionBoundaryDisplay.from_estimator(
                self.classifier,
                self.X,
                xlim=custom_xlim,
                ylim=custom_ylim,
                grid_resolution=grid_res,
                response_method="predict",
            )

            # Grid shape should match resolution
            assert display.xx0.shape == (grid_res, grid_res)
            assert display.xx1.shape == (grid_res, grid_res)
            assert display.response.shape == (grid_res, grid_res)

            # Bounds should remain the same
            assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
            assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)

    def test_edge_case_small_limits(self):
        """Test very small custom limits"""
        custom_xlim = (-0.1, 0.1)
        custom_ylim = (-0.05, 0.05)

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            response_method="predict",
        )

        # Should still work with small ranges
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)
        assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
        assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)

    def test_large_negative_limits(self):
        """Test large negative custom limits"""
        custom_xlim = (-100, -50)
        custom_ylim = (-200, -100)

        display = DecisionBoundaryDisplay.from_estimator(
            self.classifier,
            self.X,
            xlim=custom_xlim,
            ylim=custom_ylim,
            response_method="predict",
        )

        # Should handle large negative ranges
        assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
        assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)
        assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
        assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)


def test_limits_integration():
    """Integration test for limits functionality"""
    # Create synthetic dataset
    X, y = make_classification(
        n_samples=100,
        n_features=2,
        n_classes=2,
        n_informative=2,
        n_redundant=0,
        random_state=123,
    )

    # Train classifier
    clf = LogisticRegression(random_state=123)
    clf.fit(X, y)

    # Custom limits
    custom_xlim = (-3, 3)
    custom_ylim = (-2, 2)
    grid_resolution = 50

    # Create display
    display = DecisionBoundaryDisplay.from_estimator(
        clf,
        X,
        xlim=custom_xlim,
        ylim=custom_ylim,
        grid_resolution=grid_resolution,
        response_method="predict",
    )

    # Manual grid creation for comparison
    xx0, xx1 = np.meshgrid(
        np.linspace(custom_xlim[0], custom_xlim[1], grid_resolution),
        np.linspace(custom_ylim[0], custom_ylim[1], grid_resolution),
    )

    # Manual prediction
    grid_points = np.c_[xx0.ravel(), xx1.ravel()]
    expected_response = clf.predict(grid_points).reshape(xx0.shape)

    # Verify everything matches
    assert_allclose(display.xx0, xx0)
    assert_allclose(display.xx1, xx1)
    assert_allclose(display.response, expected_response)

    # Verify bounds
    assert display.xx0.min() == pytest.approx(custom_xlim[0], abs=1e-10)
    assert display.xx0.max() == pytest.approx(custom_xlim[1], abs=1e-10)
    assert display.xx1.min() == pytest.approx(custom_ylim[0], abs=1e-10)
    assert display.xx1.max() == pytest.approx(custom_ylim[1], abs=1e-10)
