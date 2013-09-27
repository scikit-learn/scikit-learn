import numpy as np
from numpy.testing import assert_equal, assert_raises
from scipy import sparse

from sklearn import linear_model


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

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           random_state=0)

    # Estimate parameters of corrupted data
    ransac_estimator.fit(X, y)

    # Ground truth / reference inlier mask
    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_,
                                   dtype=np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_is_data_valid():
    def is_data_valid(X, y):
        assert X.shape[0] == y.shape[0] == 2
        return False

    X = np.random.rand(10, 2)
    y = np.random.rand(10, 1)

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           is_data_valid=is_data_valid,
                                           random_state=0)

    assert_raises(ValueError, ransac_estimator.fit, X, y)


def test_ransac_is_model_valid():
    def is_model_valid(estimator, X, y):
        assert X.shape[0] == y.shape[0] == 2
        return False

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           is_model_valid=is_model_valid,
                                           random_state=0)

    assert_raises(ValueError, ransac_estimator.fit, X, y)


def test_ransac_max_trials():
    base_estimator = linear_model.LinearRegression()

    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5, max_trials=0,
                                           random_state=0)
    assert_raises(ValueError, ransac_estimator.fit, X, y)

    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5, max_trials=11,
                                           random_state=0)
    assert getattr(ransac_estimator, 'n_trials_', None) is None
    ransac_estimator.fit(X, y)
    assert ransac_estimator.n_trials_ == 11


def test_ransac_stop_n_inliers():
    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           stop_n_inliers=2, random_state=0)
    ransac_estimator.fit(X, y)

    assert ransac_estimator.n_trials_ == 1


def test_ransac_stop_score():
    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           stop_score=0, random_state=0)
    ransac_estimator.fit(X, y)

    assert ransac_estimator.n_trials_ == 1


def test_ransac_score():
    X = np.arange(100)[:, None]
    y = np.zeros((100, ))
    y[0] = 1
    y[1] = 100

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 0.5,
                                           random_state=0)
    ransac_estimator.fit(X, y)

    assert ransac_estimator.score(X[2:], y[2:]) == 1
    assert ransac_estimator.score(X[:2], y[:2]) < 1


def test_ransac_predict():
    X = np.arange(100)[:, None]
    y = np.zeros((100, ))
    y[0] = 1
    y[1] = 100

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 0.5,
                                           random_state=0)
    ransac_estimator.fit(X, y)

    assert_equal(ransac_estimator.predict(X), np.zeros((100, 1)))


def test_ransac_sparse_coo():
    X_sparse = sparse.coo_matrix(X)

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           random_state=0)
    ransac_estimator.fit(X_sparse, y)

    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_,
                                   dtype=np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_sparse_csr():
    X_sparse = sparse.csr_matrix(X)

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           random_state=0)
    ransac_estimator.fit(X_sparse, y)

    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_,
                                   dtype=np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_sparse_csc():
    X_sparse = sparse.csc_matrix(X)

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           random_state=0)
    ransac_estimator.fit(X_sparse, y)

    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_,
                                   dtype=np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


def test_ransac_none_estimator():

    base_estimator = linear_model.LinearRegression()

    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           random_state=0)
    ransac_none_estimator = linear_model.RANSAC(None, 2, 5, random_state=0)

    ransac_estimator.fit(X, y)
    ransac_none_estimator.fit(X, y)

    assert_equal(ransac_estimator.predict(X), ransac_none_estimator.predict(X))
    assert_raises(ValueError, ransac_none_estimator.fit, X, y.astype(int))


def test_ransac_min_n_samples():
    base_estimator = linear_model.LinearRegression()
    ransac_estimator1 = linear_model.RANSAC(base_estimator, 2, 5,
                                            random_state=0)
    ransac_estimator2 = linear_model.RANSAC(base_estimator, 2. / X.shape[0],
                                            5, random_state=0)
    ransac_estimator3 = linear_model.RANSAC(base_estimator, -1,
                                            5, random_state=0)

    ransac_estimator1.fit(X, y)
    ransac_estimator2.fit(X, y)

    assert_equal(ransac_estimator1.predict(X), ransac_estimator2.predict(X))
    assert_raises(ValueError, ransac_estimator3.fit, X, y)


def test_ransac_multi_dimensional_targets():

    base_estimator = linear_model.LinearRegression()
    ransac_estimator = linear_model.RANSAC(base_estimator, 2, 5,
                                           random_state=0)

    # 3-D target values
    yyy = np.column_stack([y, y, y])

    # Estimate parameters of corrupted data
    ransac_estimator.fit(X, yyy)

    # Ground truth / reference inlier mask
    ref_inlier_mask = np.ones_like(ransac_estimator.inlier_mask_,
                                   dtype=np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(ransac_estimator.inlier_mask_, ref_inlier_mask)


if __name__ == "__main__":
    np.testing.run_module_suite()
