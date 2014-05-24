import numpy as np
from numpy.testing import assert_equal, assert_raises
from numpy.testing import assert_array_almost_equal
from scipy import sparse

from sklearn.utils.testing import assert_less
from sklearn.linear_model import LinearRegression, RANSACRegressor
from sklearn.linear_model.ransac import _dynamic_max_trials


# Generate coordinates of line
X = np.arange(-200, 200)
y = 0.2 * X + 20
data = np.column_stack([X, y])

# Add some faulty data
outliers = np.array((10, 30, 200))
data[outliers[0], :] = (1000, 1000)
data[outliers[1], :] = (-1000, -1000)
data[outliers[2], :] = (-100, -50)

X = data[:, 0][:, np.newaxis]
y = data[:, 1]


def test_ransac_inliers_outliers():

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, random_state=0)

    # Estimate parameters of corrupted data
    ransac_estimator.fit(X, y)

    # Ground truth / reference inlier mask
    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_
                                   ).astype(np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_is_data_valid():
    def is_data_valid(X, y):
        assert_equal(X.shape[0], 2)
        assert_equal(y.shape[0], 2)
        return False

    X = np.random.rand(10, 2)
    y = np.random.rand(10, 1)

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5,
                                       is_data_valid=is_data_valid,
                                       random_state=0)

    assert_raises(ValueError, ransac_estimator.fit, X, y)


def test_ransac_is_model_valid():
    def is_model_valid(estimator, X, y):
        assert_equal(X.shape[0], 2)
        assert_equal(y.shape[0], 2)
        return False

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5,
                                       is_model_valid=is_model_valid,
                                       random_state=0)

    assert_raises(ValueError, ransac_estimator.fit, X, y)


def test_ransac_max_trials():
    base_estimator = LinearRegression()

    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, max_trials=0,
                                       random_state=0)
    assert_raises(ValueError, ransac_estimator.fit, X, y)

    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, max_trials=11,
                                       random_state=0)
    assert getattr(ransac_estimator, 'n_trials_', None) is None
    ransac_estimator.fit(X, y)
    assert_equal(ransac_estimator.n_trials_, 2)


def test_ransac_stop_n_inliers():
    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, stop_n_inliers=2,
                                       random_state=0)
    ransac_estimator.fit(X, y)

    assert_equal(ransac_estimator.n_trials_, 1)


def test_ransac_stop_score():
    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, stop_score=0,
                                       random_state=0)
    ransac_estimator.fit(X, y)

    assert_equal(ransac_estimator.n_trials_, 1)


def test_ransac_score():
    X = np.arange(100)[:, None]
    y = np.zeros((100, ))
    y[0] = 1
    y[1] = 100

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=0.5, random_state=0)
    ransac_estimator.fit(X, y)

    assert_equal(ransac_estimator.score(X[2:], y[2:]), 1)
    assert_less(ransac_estimator.score(X[:2], y[:2]), 1)


def test_ransac_predict():
    X = np.arange(100)[:, None]
    y = np.zeros((100, ))
    y[0] = 1
    y[1] = 100

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=0.5, random_state=0)
    ransac_estimator.fit(X, y)

    assert_equal(ransac_estimator.predict(X), np.zeros((100, 1)))


def test_ransac_sparse_coo():
    X_sparse = sparse.coo_matrix(X)

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, random_state=0)
    ransac_estimator.fit(X_sparse, y)

    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_
                                   ).astype(np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_sparse_csr():
    X_sparse = sparse.csr_matrix(X)

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, random_state=0)
    ransac_estimator.fit(X_sparse, y)

    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_
                                   ).astype(np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_sparse_csc():
    X_sparse = sparse.csc_matrix(X)

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, random_state=0)
    ransac_estimator.fit(X_sparse, y)

    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_
                                   ).astype(np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_none_estimator():

    base_estimator = LinearRegression()

    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, random_state=0)
    ransac_none_estimator = RANSACRegressor(None, 2, 5, random_state=0)

    ransac_estimator.fit(X, y)
    ransac_none_estimator.fit(X, y)

    assert_array_almost_equal(ransac_estimator.predict(X),
                              ransac_none_estimator.predict(X))


def test_ransac_min_n_samples():
    base_estimator = LinearRegression()
    ransac_estimator1 = RANSACRegressor(base_estimator, min_samples=2,
                                        residual_threshold=5,  random_state=0)
    ransac_estimator2 = RANSACRegressor(base_estimator,
                                        min_samples=2. / X.shape[0],
                                        residual_threshold=5, random_state=0)
    ransac_estimator3 = RANSACRegressor(base_estimator, min_samples=-1,
                                        residual_threshold=5, random_state=0)
    ransac_estimator4 = RANSACRegressor(base_estimator, min_samples=5.2,
                                        residual_threshold=5, random_state=0)
    ransac_estimator5 = RANSACRegressor(base_estimator, min_samples=2.0,
                                        residual_threshold=5, random_state=0)
    ransac_estimator6 = RANSACRegressor(base_estimator,
                                        residual_threshold=5, random_state=0)
    ransac_estimator7 = RANSACRegressor(base_estimator,
                                        min_samples=X.shape[0] + 1,
                                        residual_threshold=5, random_state=0)

    ransac_estimator1.fit(X, y)
    ransac_estimator2.fit(X, y)
    ransac_estimator5.fit(X, y)
    ransac_estimator6.fit(X, y)

    assert_array_almost_equal(ransac_estimator1.predict(X),
                              ransac_estimator2.predict(X))
    assert_array_almost_equal(ransac_estimator1.predict(X),
                              ransac_estimator5.predict(X))
    assert_array_almost_equal(ransac_estimator1.predict(X),
                              ransac_estimator6.predict(X))

    assert_raises(ValueError, ransac_estimator3.fit, X, y)
    assert_raises(ValueError, ransac_estimator4.fit, X, y)
    assert_raises(ValueError, ransac_estimator7.fit, X, y)


def test_ransac_multi_dimensional_targets():

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       residual_threshold=5, random_state=0)

    # 3-D target values
    yyy = np.column_stack([y, y, y])

    # Estimate parameters of corrupted data
    ransac_estimator.fit(X, yyy)

    # Ground truth / reference inlier mask
    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_
                                   ).astype(np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_residual_metric():
    residual_metric1 = lambda dy: np.sum(np.abs(dy), axis=1)
    residual_metric2 = lambda dy: np.sum(dy ** 2, axis=1)

    yyy = np.column_stack([y, y, y])

    base_estimator = LinearRegression()
    ransac_estimator0 = RANSACRegressor(base_estimator, min_samples=2,
                                        residual_threshold=5, random_state=0)
    ransac_estimator1 = RANSACRegressor(base_estimator, min_samples=2,
                                        residual_threshold=5, random_state=0,
                                        residual_metric=residual_metric1)
    ransac_estimator2 = RANSACRegressor(base_estimator, min_samples=2,
                                        residual_threshold=5, random_state=0,
                                        residual_metric=residual_metric2)

    # multi-dimensional
    ransac_estimator0.fit(X, yyy)
    ransac_estimator1.fit(X, yyy)
    ransac_estimator2.fit(X, yyy)
    assert_array_almost_equal(ransac_estimator0.predict(X),
                              ransac_estimator1.predict(X))
    assert_array_almost_equal(ransac_estimator0.predict(X),
                              ransac_estimator2.predict(X))

    # one-dimensional
    ransac_estimator0.fit(X, y)
    ransac_estimator2.fit(X, y)
    assert_array_almost_equal(ransac_estimator0.predict(X),
                              ransac_estimator2.predict(X))


def test_ransac_default_residual_threshold():
    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       random_state=0)

    # Estimate parameters of corrupted data
    ransac_estimator.fit(X, y)

    # Ground truth / reference inlier mask
    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_
                                   ).astype(np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_dynamic_max_trials():
    # Numbers hand-calculated and confirmed on page 119 (Table 4.3) in
    #   Hartley, R.~I. and Zisserman, A., 2004,
    #   Multiple View Geometry in Computer Vision, Second Edition,
    #   Cambridge University Press, ISBN: 0521540518

    # e = 0%, min_samples = X
    assert_equal(_dynamic_max_trials(100, 100, 2, 0.99), 1)

    # e = 5%, min_samples = 2
    assert_equal(_dynamic_max_trials(95, 100, 2, 0.99), 2)
    # e = 10%, min_samples = 2
    assert_equal(_dynamic_max_trials(90, 100, 2, 0.99), 3)
    # e = 30%, min_samples = 2
    assert_equal(_dynamic_max_trials(70, 100, 2, 0.99), 7)
    # e = 50%, min_samples = 2
    assert_equal(_dynamic_max_trials(50, 100, 2, 0.99), 17)

    # e = 5%, min_samples = 8
    assert_equal(_dynamic_max_trials(95, 100, 8, 0.99), 5)
    # e = 10%, min_samples = 8
    assert_equal(_dynamic_max_trials(90, 100, 8, 0.99), 9)
    # e = 30%, min_samples = 8
    assert_equal(_dynamic_max_trials(70, 100, 8, 0.99), 78)
    # e = 50%, min_samples = 8
    assert_equal(_dynamic_max_trials(50, 100, 8, 0.99), 1177)

    # e = 0%, min_samples = 10
    assert_equal(_dynamic_max_trials(1, 100, 10, 0), 0)
    assert_equal(_dynamic_max_trials(1, 100, 10, 1), float('inf'))

    base_estimator = LinearRegression()
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       stop_probability=-0.1)
    assert_raises(ValueError, ransac_estimator.fit, X, y)
    ransac_estimator = RANSACRegressor(base_estimator, min_samples=2,
                                       stop_probability=1.1)
    assert_raises(ValueError, ransac_estimator.fit, X, y)


if __name__ == "__main__":
    np.testing.run_module_suite()
