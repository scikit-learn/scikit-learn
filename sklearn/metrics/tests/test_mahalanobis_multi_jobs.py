import numpy as np
import pytest

from sklearn.metrics import DistanceMetric
from sklearn.neighbors import LocalOutlierFactor


def make_dataset(rng, n_inliers, n_outliers, n_features, outlier_shift):
    X_in = rng.normal(0.0, 1.0, size=(n_inliers, n_features))
    X_out = rng.normal(outlier_shift, 0.5, size=(n_outliers, n_features))
    return np.vstack([X_in, X_out])


def mahalanobis_metric(X):
    V = np.cov(X.T)
    return DistanceMetric.get_metric("mahalanobis", V=V)


def fit_predict_lof(X, metric, n_neighbors, n_jobs):
    lof = LocalOutlierFactor(
        n_neighbors=n_neighbors,
        metric=metric,
        n_jobs=n_jobs,
    )
    return lof.fit_predict(X)


@pytest.mark.parametrize("n_jobs", [2, 4, 8])
@pytest.mark.parametrize(
    "seed,n_inliers,n_outliers,n_features,outlier_shift,n_neighbors",
    [
        (0, 200, 10, 2, 10.0, 20),
        (1, 500, 25, 2, 8.0, 35),
        (2, 300, 30, 5, 12.0, 25),
        (3, 1000, 50, 10, 15.0, 40),
        (4, 150, 15, 3, 6.0, 15),
    ],
)
def test_lof_mahalanobis_parallel_matches_serial(
    n_jobs, seed, n_inliers, n_outliers, n_features, outlier_shift, n_neighbors
):
    rng = np.random.RandomState(seed)
    X = make_dataset(rng, n_inliers, n_outliers, n_features, outlier_shift)
    metric = mahalanobis_metric(X)

    y_serial = fit_predict_lof(X, metric, n_neighbors=n_neighbors, n_jobs=1)
    y_parallel = fit_predict_lof(X, metric, n_neighbors=n_neighbors, n_jobs=n_jobs)

    np.testing.assert_array_equal(y_parallel, y_serial)


@pytest.mark.parametrize("n_jobs", [2, 4, 8])
def test_lof_mahalanobis_parallel_stable_across_runs(n_jobs):
    rng = np.random.RandomState(123)
    X = make_dataset(
        rng, n_inliers=400, n_outliers=20, n_features=5, outlier_shift=11.0
    )
    metric = mahalanobis_metric(X)

    y1 = fit_predict_lof(X, metric, n_neighbors=30, n_jobs=n_jobs)
    y2 = fit_predict_lof(X, metric, n_neighbors=30, n_jobs=n_jobs)

    np.testing.assert_array_equal(y1, y2)
