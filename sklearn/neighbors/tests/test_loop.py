# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Valentino Constantinou <vc@valentino.io>
# License: BSD 3 clause

from math import log10
import numpy as np
from sklearn import neighbors

from numpy.testing import assert_array_equal

from sklearn import metrics
from sklearn.metrics import roc_auc_score

from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_warns_message, assert_raises

from sklearn.datasets import load_iris

# load the iris dataset
# and randomly permute it
rng = check_random_state(0)
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_loop():
    # Toy sample (the last two samples are outliers):
    X = [[-2, -1], [-1, -1], [-1, -2], [1, 2], [1, 2], [2, 1], [5, 3], [-4, 2]]

    # Test LocalOutlierFactor:
    clf = neighbors.LocalOutlierProbability(n_neighbors=5)
    score = -1. * clf.fit(X).negative_local_outlier_probability_
    assert_array_equal(clf._fit_X, X)

    # Assert smallest outlier score is greater than largest inlier score:
    assert_greater(np.min(score[-2:]), np.max(score[:-2]))

    # Assert predict() works:
    clf = neighbors.LocalOutlierProbability(n_neighbors=5, norm_factor=0.8).fit(X)
    assert_array_equal(clf._predict(), 6 * [1] + 2 * [-1])


def test_loop_performance():
    # Generate train/test data
    rng = check_random_state(2)
    X = 0.3 * rng.randn(120, 2)
    X_train = X[:100]

    # Generate some abnormal novel observations
    X_outliers = rng.uniform(low=-4, high=4, size=(20, 2))
    X_test = np.r_[X[100:], X_outliers]
    y_test = np.array([0] * 20 + [1] * 20)

    # fit the model
    clf = neighbors.LocalOutlierProbability().fit(X_train)

    # predict scores (the lower, the more normal)
    y_pred = -clf._decision_function(X_test)

    # check that roc_auc is good
    assert_greater(roc_auc_score(y_test, y_pred), .99)


def test_loop_values():
    # toy samples:
    X_train = [[1, 1], [1, 2], [2, 1]]
    clf1 = neighbors.LocalOutlierProbability(n_neighbors=2, norm_factor=0.95).fit(X_train)
    clf2 = neighbors.LocalOutlierProbability(n_neighbors=2).fit(X_train)
    # define test values and labels
    s_0 = 0.
    s_1 = log10(3.493965)
    s_0_label = 1
    s_1_label = -1
    # check predict()
    assert_array_almost_equal(-clf1.negative_local_outlier_probability_, [s_0, s_1, s_1])
    assert_array_almost_equal(-clf2.negative_local_outlier_probability_, [s_0, s_1, s_1])
    # check predict(one sample not in train)
    assert_array_almost_equal(-clf1._score_samples([[2., 2.]], mode='loop'), [s_0_label])
    assert_array_almost_equal(-clf2._score_samples([[2., 2.]], mode='loop'), [s_0_label])
    # check predict(one sample already in train)
    assert_array_almost_equal(-clf1._score_samples([[1., 1.]], mode='loop'), [s_1_label])
    assert_array_almost_equal(-clf2._score_samples([[1., 1.]], mode='loop'), [s_1_label])


def test_loop_precomputed(random_state=42):
    """Tests LoOP with a distance matrix."""
    # Note: smaller samples may result in spurious test success
    rng = np.random.RandomState(random_state)
    X = rng.random_sample((10, 4))
    Y = rng.random_sample((3, 4))
    DXX = metrics.pairwise_distances(X, metric='euclidean')
    DYX = metrics.pairwise_distances(Y, X, metric='euclidean')
    # As a feature matrix (n_samples by n_features)
    loop_X = neighbors.LocalOutlierProbability(n_neighbors=3)
    loop_X.fit(X)
    pred_X_X = loop_X._predict()
    pred_X_Y = loop_X._predict(Y)

    # As a dense distance matrix (n_samples by n_samples)
    loop_D = neighbors.LocalOutlierProbability(n_neighbors=3, algorithm='brute',
                                               metric='precomputed')
    loop_D.fit(DXX)
    pred_D_X = loop_D._predict()
    pred_D_Y = loop_D._predict(DYX)

    assert_array_almost_equal(pred_X_X, pred_D_X)
    assert_array_almost_equal(pred_X_Y, pred_D_Y)


def test_n_neighbors_attribute():
    X = iris.data
    clf = neighbors.LocalOutlierProbability(n_neighbors=500).fit(X)
    assert_equal(clf.n_neighbors_, X.shape[0] - 1)

    clf = neighbors.LocalOutlierProbability(n_neighbors=500)
    assert_warns_message(UserWarning,
                         "n_neighbors will be set to (n_samples - 1)",
                         clf.fit, X)
    assert_equal(clf.n_neighbors_, X.shape[0] - 1)


def test_score_samples():
    X_train = [[1, 1], [1, 2], [2, 1]]
    clf1 = neighbors.LocalOutlierProbability(n_neighbors=2,
                                             norm_factor=0.9).fit(X_train)
    clf2 = neighbors.LocalOutlierProbability(n_neighbors=2).fit(X_train)
    assert_array_equal(clf1._score_samples([[2., 2.]], mode='loop'),
                       clf1._decision_function([[2., 2.]]))

    assert_array_equal(clf2._score_samples([[2., 2.]], mode='loop'),
                       clf2._decision_function([[2., 2.]]))
    assert_array_equal(clf1._score_samples([[2., 2.]], mode='loop'),
                       clf2._score_samples([[2., 2.]], mode='loop'))


def test_norm_factor():
    X = [[1, 1], [1, 0]]
    clf = neighbors.LocalOutlierProbability(n_neighbors=2, norm_factor=1.2)
    assert_raises(ValueError, clf.fit, X)
