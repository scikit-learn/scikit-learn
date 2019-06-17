import numpy as np
import pytest

from sklearn.impute import KNNImputer
from sklearn.metrics.pairwise import nan_euclidean_distances
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.neighbors import KNeighborsRegressor
from sklearn.utils.mask import _get_missing_mask
from sklearn.utils.testing import assert_array_almost_equal


def test_knn_imputer_shape():
    # Verify the shapes of the imputed matrix for different weights and
    # number of neighbors.
    n_rows = 10
    n_cols = 2
    X = np.random.rand(n_rows, n_cols)
    X[0, 0] = np.nan

    for weights in ['uniform', 'distance']:
        for n_neighbors in range(1, 6):
            imputer = KNNImputer(n_neighbors=n_neighbors, weights=weights)
            X_imputed = imputer.fit_transform(X)
            assert X_imputed.shape == (n_rows, n_cols)


def test_knn_imputer_errors():
    missing_values = 0
    n_neighbors = 2
    imputer = KNNImputer(missing_values=missing_values,
                         n_neighbors=n_neighbors,
                         weights="uniform")

    # Test with missing_values=0 when NaN present
    X = np.array([
        [np.nan, 0, 0, 0, 5],
        [np.nan, 1, 0, np.nan, 3],
        [np.nan, 2, 0, 0, 0],
        [np.nan, 6, 0, 5, 13],
    ])
    msg = (r"Input contains NaN, infinity or a value too large for "
           r"dtype\('float64'\)")
    with pytest.raises(ValueError, match=msg):
        imputer.fit(X)

    # Test with % zeros in samples > sample_max_missing
    X = np.array([
        [1, 0, 0, 0, 5],
        [2, 1, 0, 2, 3],
        [3, 2, 0, 0, 0],
        [4, 6, 0, 5, 13],
    ])
    msg = "Some columns have more than {}% missing values".format(
        imputer.sample_max_missing * 100)
    with pytest.raises(ValueError, match=msg):
        imputer.fit(X)


def _missing_mean(X, missing_value):
    masked_X = np.ma.array(X, mask=_get_missing_mask(X, missing_value))
    return np.ma.average(masked_X, axis=0).data


@pytest.mark.parametrize("na", [np.nan, -1])
def test_knn_imputer_zero_nan_imputes_the_same(na):
    # Test with an imputable matrix and also compare with missing_values=np.NaN
    X_zero = np.array([
        [1, 0, 1, 1, 1.],
        [2, 2, 2, 2, 2],
        [3, 3, 3, 3, 0],
        [6, 6, 0, 6, 6],
    ])

    X_nan = np.array([
        [1, na, 1, 1, 1.],
        [2, 2, 2, 2, 2],
        [3, 3, 3, 3, na],
        [6, 6, na, 6, 6],
    ])
    statistics_mean = _missing_mean(X_nan, na)

    X_imputed = np.array([
        [1, 2.5, 1, 1, 1.],
        [2, 2, 2, 2, 2],
        [3, 3, 3, 3, 1.5],
        [6, 6, 2.5, 6, 6],
    ])

    imputer_zero = KNNImputer(missing_values=0, n_neighbors=2,
                              weights="uniform")

    imputer_nan = KNNImputer(missing_values=na,
                             n_neighbors=2,
                             weights="uniform")

    assert_array_almost_equal(imputer_zero.fit_transform(X_zero), X_imputed)
    assert_array_almost_equal(imputer_zero.statistics_, statistics_mean)
    assert_array_almost_equal(imputer_zero.fit_transform(X_zero),
                              imputer_nan.fit_transform(X_nan))


@pytest.mark.parametrize("na", [np.nan, -1])
def test_knn_imputer_verify(na):
    # Test imputation with default parameter values

    # Test with an imputable matrix
    X = np.array([
        [1, 0, 0, 1],
        [2, 1, 2, na],
        [3, 2, 3, na],
        [na, 4, 5, 5],
        [6, na, 6, 7],
        [8, 8, 8, 8],
        [16, 15, 18, 19],
    ])
    statistics_mean = _missing_mean(X, na)

    X_imputed = np.array([
        [1, 0, 0, 1],
        [2, 1, 2, 8],
        [3, 2, 3, 8],
        [4, 4, 5, 5],
        [6, 3, 6, 7],
        [8, 8, 8, 8],
        [16, 15, 18, 19],
    ])

    imputer = KNNImputer(missing_values=na)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed)
    assert_array_almost_equal(imputer.statistics_, statistics_mean)

    # Test with % missing in features > feature_max_missing
    X = np.array([
        [1, 0, 0, 1],
        [2, 1, 2, na],
        [3, 2, 3, na],
        [na, 4, 5, 5],
        [6, na, 6, 7],
        [8, 8, 8, 8],
        [19, 19, 19, 19],
        [na, na, na, 19],
    ])
    statistics_mean = _missing_mean(X, na)
    r7c0, r7c1, r7c2, _ = statistics_mean

    X_imputed = np.array([
        [1, 0, 0, 1],
        [2, 1, 2, 8],
        [3, 2, 3, 8],
        [4, 4, 5, 5],
        [6, 3, 6, 7],
        [8, 8, 8, 8],
        [19, 19, 19, 19],
        [r7c0, r7c1, r7c2, 19],
    ])

    imputer = KNNImputer(missing_values=na)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed, decimal=6)
    assert_array_almost_equal(imputer.statistics_, statistics_mean, decimal=6)

    # Test with all neighboring donors also having missing feature values
    X = np.array([
        [1, 0, 0, na],
        [2, 1, 2, na],
        [3, 2, 3, na],
        [4, 4, 5, na],
        [6, 7, 6, na],
        [8, 8, 8, na],
        [20, 20, 20, 20],
        [22, 22, 22, 22]
    ])
    statistics_mean = _missing_mean(X, na)

    X_imputed = np.array([
        [1, 0, 0, 21],
        [2, 1, 2, 21],
        [3, 2, 3, 21],
        [4, 4, 5, 21],
        [6, 7, 6, 21],
        [8, 8, 8, 21],
        [20, 20, 20, 20],
        [22, 22, 22, 22]
    ])

    imputer = KNNImputer(missing_values=na)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed)
    assert_array_almost_equal(imputer.statistics_, statistics_mean)

    # Test when data in fit() and transform() are different
    X = np.array([
        [0, 0],
        [na, 2],
        [4, 3],
        [5, 6],
        [7, 7],
        [9, 8],
        [11, 16]
    ])
    statistics_mean = _missing_mean(X, na)

    Y = np.array([
        [1, 0],
        [3, 2],
        [4, na]
    ])

    Y_imputed = np.array([
        [1, 0],
        [3, 2],
        [4, 4.8]
    ])

    imputer = KNNImputer(missing_values=na)
    assert_array_almost_equal(imputer.fit(X).transform(Y), Y_imputed)
    assert_array_almost_equal(imputer.statistics_, statistics_mean)


@pytest.mark.parametrize("na", [np.nan, -1])
def test_knn_imputer_default_with_invalid_input(na):
    # Test imputation with default values and invalid input
    # Test with % missing in a samples > sample_max_missing
    X = np.array([
        [na, 0, 0, 0, 5],
        [na, 1, 0, na, 3],
        [na, 2, 0, 0, 0],
        [na, 6, 0, 5, 13],
        [na, 7, 0, 7, 8],
        [na, 8, 0, 8, 9],
    ])
    imputer = KNNImputer(missing_values=na)
    msg = "Some columns have more than {}% missing values".format(
        imputer.sample_max_missing * 100)
    with pytest.raises(ValueError, match=msg):
        imputer.fit(X)

    # Test with insufficient number of neighbors
    X = np.array([
        [1, 1, 1, 2, na],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [6, 6, 2, 5, 13],
    ])
    msg = "There are only %d samples, but n_neighbors=%d" % \
          (X.shape[0], imputer.n_neighbors)
    with pytest.raises(ValueError, match=msg):
        imputer.fit(X)

    # Test with inf present
    X = np.array([
        [np.inf, 1, 1, 2, na],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [na, 6, 0, 5, 13],
        [na, 7, 0, 7, 8],
        [6, 6, 2, 5, 7],
    ])
    with pytest.raises(ValueError, match="Input contains (infinity|NaN)"):
        KNNImputer(missing_values=na).fit(X)

    # Test with inf present in matrix passed in transform()
    X = np.array([
        [np.inf, 1, 1, 2, na],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [na, 6, 0, 5, 13],
        [na, 7, 0, 7, 8],
        [6, 6, 2, 5, 7],
    ])

    X_fit = np.array([
        [0, 1, 1, 2, na],
        [2, 1, 2, 2, 3],
        [3, 2, 3, 3, 8],
        [na, 6, 0, 5, 13],
        [na, 7, 0, 7, 8],
        [6, 6, 2, 5, 7],
    ])
    with pytest.raises(ValueError, match="Input contains (infinity|NaN)"):
        KNNImputer(missing_values=na).fit(X_fit).transform(X)


@pytest.mark.parametrize("na", [np.nan, -1])
def test_knn_imputer_n_neighbors(na):

    X = np.array([
        [0, 0],
        [na, 2],
        [4, 3],
        [5, na],
        [7, 7],
        [na, 8],
        [14, 13]
    ])
    statistics_mean = _missing_mean(X, na)

    # Test with 1 neighbor
    X_imputed_1NN = np.array([
        [0, 0],
        [4, 2],
        [4, 3],
        [5, 3],
        [7, 7],
        [7, 8],
        [14, 13]
    ])

    n_neighbors = 1
    imputer = KNNImputer(n_neighbors=n_neighbors, missing_values=na)

    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_1NN)
    assert_array_almost_equal(imputer.statistics_, statistics_mean)

    # Test with 6 neighbors
    X = np.array([
        [0, 0],
        [na, 2],
        [4, 3],
        [5, na],
        [7, 7],
        [na, 8],
        [14, 13]
    ])

    X_imputed_6NN = np.array([
        [0, 0],
        [6, 2],
        [4, 3],
        [5, 5.5],
        [7, 7],
        [6, 8],
        [14, 13]
    ])

    n_neighbors = 6
    imputer = KNNImputer(n_neighbors=6, missing_values=na)
    imputer_plus1 = KNNImputer(n_neighbors=n_neighbors + 1,
                               missing_values=na)

    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_6NN)
    assert_array_almost_equal(imputer.statistics_, statistics_mean)
    assert_array_almost_equal(imputer.fit_transform(X), imputer_plus1.fit(
        X).transform(X))


@pytest.mark.parametrize("na", [np.nan, -1])
def test_knn_imputer_weight_uniform(na):

    X = np.array([
        [0, 0],
        [na, 2],
        [4, 3],
        [5, 6],
        [7, 7],
        [9, 8],
        [11, 10]
    ])

    # Test with "uniform" weight (or unweighted)
    X_imputed_uniform = np.array([
        [0, 0],
        [5, 2],
        [4, 3],
        [5, 6],
        [7, 7],
        [9, 8],
        [11, 10]
    ])

    imputer = KNNImputer(weights="uniform", missing_values=na)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_uniform)

    # Test with "callable" weight
    def no_weight(dist=None):
        return None

    imputer = KNNImputer(weights=no_weight, missing_values=na)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_uniform)


def test_knn_imputer_weight_distance():
    na = np.nan
    X = np.array([
        [0, 0],
        [na, 2],
        [4, 3],
        [5, 6],
        [7, 7],
        [9, 8],
        [11, 10]
    ])

    # Test with "distance" weight
    nn = KNeighborsRegressor(metric="euclidean", weights="distance")
    nn.fit(np.delete(X, 1, axis=0)[:, 1:], np.delete(X, 1, axis=0)[:, 0])
    imputed = nn.predict(X[1:2, 1:])

    # Manual calculation
    X_imputed_distance1 = np.array([
        [0, 0],
        [3.850394, 2],
        [4, 3],
        [5, 6],
        [7, 7],
        [9, 8],
        [11, 10]
    ])

    # NearestNeighbor calculation
    X_imputed_distance2 = np.array([
        [0, 0],
        [imputed, 2],
        [4, 3],
        [5, 6],
        [7, 7],
        [9, 8],
        [11, 10]
    ])

    imputer = KNNImputer(weights="distance")
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_distance1,
                              decimal=6)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed_distance2,
                              decimal=6)

    # Test with weights = "distance" and n_neighbors=2
    X = np.array([
        [na, 0, 0],
        [2, 1, 2],
        [3, 2, 3],
        [4, 5, 5],
    ])
    statistics_mean = np.nanmean(X, axis=0)

    X_imputed = np.array([
        [2.3828, 0, 0],
        [2, 1, 2],
        [3, 2, 3],
        [4, 5, 5],
    ])

    imputer = KNNImputer(n_neighbors=2, weights="distance")
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed,
                              decimal=4)
    assert_array_almost_equal(imputer.statistics_, statistics_mean)

    # Test with varying missingness patterns
    X = np.array([
        [1, 0, 0, 1],
        [0, na, 1, na],
        [1, 1, 1, na],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [1, 0, 1, 1],
        [10, 10, 10, 10],
    ])
    statistics_mean = np.nanmean(X, axis=0)

    # Get weights of donor neighbors
    dist = nan_euclidean_distances(X)
    r1c1_nbor_dists = dist[1, [0, 2, 3, 4, 5]]
    r1c3_nbor_dists = dist[1, [0, 3, 4, 5, 6]]
    r1c1_nbor_wt = (1/r1c1_nbor_dists)
    r1c3_nbor_wt = (1 / r1c3_nbor_dists)

    r2c3_nbor_dists = dist[2, [0, 3, 4, 5, 6]]
    r2c3_nbor_wt = 1/r2c3_nbor_dists

    # Collect donor values
    col1_donor_values = np.ma.masked_invalid(X[[0, 2, 3, 4, 5], 1]).copy()
    col3_donor_values = np.ma.masked_invalid(X[[0, 3, 4, 5, 6], 3]).copy()

    # Final imputed values
    r1c1_imp = np.ma.average(col1_donor_values, weights=r1c1_nbor_wt)
    r1c3_imp = np.ma.average(col3_donor_values, weights=r1c3_nbor_wt)
    r2c3_imp = np.ma.average(col3_donor_values, weights=r2c3_nbor_wt)

    X_imputed = np.array([
        [1, 0, 0, 1],
        [0, r1c1_imp, 1, r1c3_imp],
        [1, 1, 1, r2c3_imp],
        [0, 1, 0, 0],
        [0, 0, 0, 0],
        [1, 0, 1, 1],
        [10, 10, 10, 10],
    ])

    imputer = KNNImputer(weights="distance")
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed, decimal=6)
    assert_array_almost_equal(imputer.statistics_, statistics_mean)


def test_knn_imputer_metric_type():
    X = np.array([
        [0, 0],
        [np.nan, 2],
    ])

    # Test with a metric type without NaN support
    imputer = KNNImputer(metric="euclidean")
    bad_metric_msg = "The selected metric does not support NaN values"
    with pytest.raises(ValueError, match=bad_metric_msg):
        imputer.fit(X)


def test_knn_imputer_callable_metric():

    # Define callable metric that returns the l1 norm:
    def custom_callable(x, y, missing_values=np.nan, squared=False):
        x = np.ma.array(x, mask=np.isnan(x))
        y = np.ma.array(y, mask=np.isnan(y))
        dist = np.nansum(np.abs(x-y))
        return dist

    X = np.array([
        [4, 3, 3, np.nan],
        [6, 9, 6, 9],
        [4, 8, 6, 9],
        [np.nan, 9, 11, 10.]
    ])

    X_imputed = np.array([
        [4, 3, 3, 9],
        [6, 9, 6, 9],
        [4, 8, 6, 9],
        [5, 9, 11, 10.]
    ])

    imputer = KNNImputer(n_neighbors=2, metric=custom_callable)
    assert_array_almost_equal(imputer.fit_transform(X), X_imputed)


def test_knn_imputer_with_simple_example():

    na = np.nan

    X = np.array([
        [0, na, 0, na],
        [1, 1, 1, na],
        [2, 2, na,  2],
        [3, 3, 3, 3],
        [4, 4, 4, 4],
        [5, 5, 5, 5],
        [6, 6, 6, 6],
        [na, 7, 7, 7]
    ])

    r0c1 = np.mean(X[1:6, 1])
    r0c3 = np.mean(X[2:-1, -1])
    r1c3 = np.mean(X[2:-1, -1])
    r2c2 = np.nanmean(X[:6, 2])
    r7c0 = np.mean(X[2:-1, 0])

    X_imputed = np.array([
        [0, r0c1, 0, r0c3],
        [1, 1, 1, r1c3],
        [2, 2, r2c2, 2],
        [3, 3, 3, 3],
        [4, 4, 4, 4],
        [5, 5, 5, 5],
        [6, 6, 6, 6],
        [r7c0, 7, 7, 7]
    ])

    imputer_comp = KNNImputer(missing_values=na)
    assert_array_almost_equal(imputer_comp.fit_transform(X), X_imputed)


@pytest.mark.parametrize("na", [-1, np.nan])
def test_knn_imputer_with_weighted_features(na):

    X = np.array([
        [0, 0, 0, na],
        [1, 1, 1, na],
        [2, 2, na, 2],
        [3, 3, 3, 3],
        [4, 4, 4, 4],
        [5, 5, 5, 5],
        [6, 6, 6, 6],
        [na, 7, 7, 7]
    ])

    dist = pairwise_distances(X, metric="nan_euclidean", squared=False,
                              missing_values=na)

    # Calculate weights
    r0c3_w = 1.0 / dist[0, 2:-1]
    r1c3_w = 1.0 / dist[1, 2:-1]
    r2c2_w = 1.0 / dist[2, (0, 1, 3, 4, 5)]
    r7c0_w = 1.0 / dist[7, 2:7]

    # Calculate weighted averages
    r0c3 = np.average(X[2:-1, -1], weights=r0c3_w)
    r1c3 = np.average(X[2:-1, -1], weights=r1c3_w)
    r2c2 = np.average(X[(0, 1, 3, 4, 5), 2], weights=r2c2_w)
    r7c0 = np.average(X[2:7, 0], weights=r7c0_w)

    X_imputed = np.array([
        [0, 0, 0, r0c3],
        [1, 1, 1, r1c3],
        [2, 2, r2c2, 2],
        [3, 3, 3, 3],
        [4, 4, 4, 4],
        [5, 5, 5, 5],
        [6, 6, 6, 6],
        [r7c0, 7, 7, 7]
    ])

    imputer_comp_wt = KNNImputer(missing_values=na, weights="distance")
    assert_array_almost_equal(imputer_comp_wt.fit_transform(X), X_imputed)
