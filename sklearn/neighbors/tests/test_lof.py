# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

from math import sqrt
import numpy as np
from sklearn import neighbors

from numpy.testing import assert_array_equal

from sklearn import metrics
from sklearn.metrics import roc_auc_score

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_warns_message

from sklearn.datasets import load_iris


# load the iris dataset
# and randomly permute it
rng = check_random_state(0)
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_lof():
    # Toy sample (the last two samples are outliers):
    X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1], [5, 3], [-4, 2]]

    # Test LocalOutlierFactor:
    clf = neighbors.LocalOutlierFactor(n_neighbors=5)
    score = clf.fit(X).negative_outlier_factor_
    assert_array_equal(clf._fit_X, X)

    # Assert largest outlier score is smaller than smallest inlier score:
    assert_greater(np.min(score[:-2]), np.max(score[-2:]))

    # Assert predict() works:
    clf = neighbors.LocalOutlierFactor(contamination=0.25,
                                       n_neighbors=5).fit(X)
    assert_array_equal(clf._predict(), 6 * [1] + 2 * [-1])


def test_lof_performance():
    # Generate train/test data
    rng = check_random_state(2)
    X = 0.3 * rng.randn(120, 2)
    X_train = np.r_[X + 2, X - 2]
    X_train = X[:100]

    # Generate some abnormal novel observations
    X_outliers = rng.uniform(low=-4, high=4, size=(20, 2))
    X_test = np.r_[X[100:], X_outliers]
    y_test = np.array([0] * 20 + [1] * 20)

    # fit the model
    clf = neighbors.LocalOutlierFactor().fit(X_train)

    # predict scores (the lower, the more normal)
    y_pred = -clf._decision_function(X_test)

    # check that roc_auc is good
    assert_greater(roc_auc_score(y_test, y_pred), .99)


def test_lof_values():
    # toy samples:
    X_train = [[1, 1], [1, 2], [2, 1]]
    clf = neighbors.LocalOutlierFactor(n_neighbors=2).fit(X_train)
    s_0 = 2. * sqrt(2.) / (1. + sqrt(2.))
    s_1 = (1. + sqrt(2)) * (1. / (4. * sqrt(2.)) + 1. / (2. + 2. * sqrt(2)))
    # check predict()
    assert_array_almost_equal(-clf.negative_outlier_factor_, [s_0, s_1, s_1])
    # check predict(one sample not in train)
    assert_array_almost_equal(-clf._decision_function([[2., 2.]]), [s_0])
    # # check predict(one sample already in train)
    assert_array_almost_equal(-clf._decision_function([[1., 1.]]), [s_1])


def test_lof_precomputed(random_state=42):
    """Tests LOF with a distance matrix."""
    # Note: smaller samples may result in spurious test success
    rng = np.random.RandomState(random_state)
    X = rng.random_sample((10, 4))
    Y = rng.random_sample((3, 4))
    DXX = metrics.pairwise_distances(X, metric='euclidean')
    DYX = metrics.pairwise_distances(Y, X, metric='euclidean')
    # As a feature matrix (n_samples by n_features)
    lof_X = neighbors.LocalOutlierFactor(n_neighbors=3)
    lof_X.fit(X)
    pred_X_X = lof_X._predict()
    pred_X_Y = lof_X._predict(Y)

    # As a dense distance matrix (n_samples by n_samples)
    lof_D = neighbors.LocalOutlierFactor(n_neighbors=3, algorithm='brute',
                                         metric='precomputed')
    lof_D.fit(DXX)
    pred_D_X = lof_D._predict()
    pred_D_Y = lof_D._predict(DYX)

    assert_array_almost_equal(pred_X_X, pred_D_X)
    assert_array_almost_equal(pred_X_Y, pred_D_Y)


def test_n_neighbors_attribute():
    X = iris.data
    clf = neighbors.LocalOutlierFactor(n_neighbors=500).fit(X)
    assert_equal(clf.n_neighbors_, X.shape[0] - 1)

    clf = neighbors.LocalOutlierFactor(n_neighbors=500)
    assert_warns_message(UserWarning,
                         "n_neighbors will be set to (n_samples - 1)",
                         clf.fit, X)
    assert_equal(clf.n_neighbors_, X.shape[0] - 1)
