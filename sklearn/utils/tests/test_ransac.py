import numpy as np
from numpy.testing import assert_equal, assert_raises

from sklearn import linear_model
from sklearn.utils import ransac


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
    np.random.seed(1)

    # Estimate parameters of corrupted data
    _, inlier_mask = ransac(X, y, linear_model.LinearRegression(), 2, 5)

    # Ground truth / reference inlier mask
    ref_inlier_mask = np.ones_like(inlier_mask, dtype=np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(inlier_mask, ref_inlier_mask)


def test_ransac_is_data_valid():
    def is_data_valid(X, y):
        return False

    X = np.random.rand(10, 2)
    y = np.random.rand(10, 1)
    estimator = linear_model.LinearRegression()
    assert_raises(ValueError, ransac, X, y, estimator, 2, 5,
                  is_data_valid=is_data_valid)


def test_ransac_is_model_valid():
    def is_model_valid(estimator, X, y):
        return False
    estimator = linear_model.LinearRegression()
    assert_raises(ValueError, ransac, X, y, estimator, 2, 5,
                  is_model_valid=is_model_valid)


def test_ransac_max_trials():
    estimator = linear_model.LinearRegression()
    assert_raises(ValueError, ransac, X, y, estimator, 2, 5, max_trials=0)
    assert ransac(X, y, estimator, 2, 5, max_trials=11)[0] == 11


if __name__ == "__main__":
    np.testing.run_module_suite()
