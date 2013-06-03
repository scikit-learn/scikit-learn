import numpy as np
from numpy.testing import assert_equal

from sklearn import linear_model
from sklearn.utils import ransac


def test_ransac_inliers_outliers():
    np.random.seed(1)

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

    # Estimate parameters of corrupted data
    inlier_mask = ransac(X, y, linear_model.LinearRegression(), 3, 5)

    # Ground truth / reference inlier mask
    ref_inlier_mask = np.ones_like(inlier_mask, dtype=np.bool_)
    ref_inlier_mask[outliers] = False

    assert_equal(inlier_mask, ref_inlier_mask)


if __name__ == "__main__":
    np.testing.run_module_suite()
