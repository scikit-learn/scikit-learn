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

    assert_equal(ransac_estimator.predict(X), np.zeros((100, )))


def test_ransac_sparse():
    X_sparse = sparse.coo_matrix(X)

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


if __name__ == "__main__":
    np.testing.run_module_suite()
