# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

from math import sqrt

import pytest
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
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.estimator_checks import check_estimator

from sklearn.datasets import load_iris


# load the iris dataset
# and randomly permute it
rng = check_random_state(0)
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
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
    assert_array_equal(clf.fit_predict(X), 6 * [1] + 2 * [-1])


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_lof_performance():
    # Generate train/test data
    rng = check_random_state(2)
    X = 0.3 * rng.randn(120, 2)
    X_train = X[:100]

    # Generate some abnormal novel observations
    X_outliers = rng.uniform(low=-4, high=4, size=(20, 2))
    X_test = np.r_[X[100:], X_outliers]
    y_test = np.array([0] * 20 + [1] * 20)

    # fit the model for novelty detection
    clf = neighbors.LocalOutlierFactor(novelty=True).fit(X_train)

    # predict scores (the lower, the more normal)
    y_pred = -clf.decision_function(X_test)

    # check that roc_auc is good
    assert_greater(roc_auc_score(y_test, y_pred), .99)


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_lof_values():
    # toy samples:
    X_train = [[1, 1], [1, 2], [2, 1]]
    clf1 = neighbors.LocalOutlierFactor(n_neighbors=2,
                                        contamination=0.1,
                                        novelty=True).fit(X_train)
    clf2 = neighbors.LocalOutlierFactor(n_neighbors=2,
                                        novelty=True).fit(X_train)
    s_0 = 2. * sqrt(2.) / (1. + sqrt(2.))
    s_1 = (1. + sqrt(2)) * (1. / (4. * sqrt(2.)) + 1. / (2. + 2. * sqrt(2)))
    # check predict()
    assert_array_almost_equal(-clf1.negative_outlier_factor_, [s_0, s_1, s_1])
    assert_array_almost_equal(-clf2.negative_outlier_factor_, [s_0, s_1, s_1])
    # check predict(one sample not in train)
    assert_array_almost_equal(-clf1.score_samples([[2., 2.]]), [s_0])
    assert_array_almost_equal(-clf2.score_samples([[2., 2.]]), [s_0])
    # check predict(one sample already in train)
    assert_array_almost_equal(-clf1.score_samples([[1., 1.]]), [s_1])
    assert_array_almost_equal(-clf2.score_samples([[1., 1.]]), [s_1])


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_lof_precomputed(random_state=42):
    """Tests LOF with a distance matrix."""
    # Note: smaller samples may result in spurious test success
    rng = np.random.RandomState(random_state)
    X = rng.random_sample((10, 4))
    Y = rng.random_sample((3, 4))
    DXX = metrics.pairwise_distances(X, metric='euclidean')
    DYX = metrics.pairwise_distances(Y, X, metric='euclidean')
    # As a feature matrix (n_samples by n_features)
    lof_X = neighbors.LocalOutlierFactor(n_neighbors=3, novelty=True)
    lof_X.fit(X)
    pred_X_X = lof_X._predict()
    pred_X_Y = lof_X.predict(Y)

    # As a dense distance matrix (n_samples by n_samples)
    lof_D = neighbors.LocalOutlierFactor(n_neighbors=3, algorithm='brute',
                                         metric='precomputed', novelty=True)
    lof_D.fit(DXX)
    pred_D_X = lof_D._predict()
    pred_D_Y = lof_D.predict(DYX)

    assert_array_almost_equal(pred_X_X, pred_D_X)
    assert_array_almost_equal(pred_X_Y, pred_D_Y)


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_n_neighbors_attribute():
    X = iris.data
    clf = neighbors.LocalOutlierFactor(n_neighbors=500).fit(X)
    assert_equal(clf.n_neighbors_, X.shape[0] - 1)

    clf = neighbors.LocalOutlierFactor(n_neighbors=500)
    assert_warns_message(UserWarning,
                         "n_neighbors will be set to (n_samples - 1)",
                         clf.fit, X)
    assert_equal(clf.n_neighbors_, X.shape[0] - 1)


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_score_samples():
    X_train = [[1, 1], [1, 2], [2, 1]]
    clf1 = neighbors.LocalOutlierFactor(n_neighbors=2,
                                        contamination=0.1,
                                        novelty=True).fit(X_train)
    clf2 = neighbors.LocalOutlierFactor(n_neighbors=2,
                                        novelty=True).fit(X_train)
    assert_array_equal(clf1.score_samples([[2., 2.]]),
                       clf1.decision_function([[2., 2.]]) + clf1.offset_)
    assert_array_equal(clf2.score_samples([[2., 2.]]),
                       clf2.decision_function([[2., 2.]]) + clf2.offset_)
    assert_array_equal(clf1.score_samples([[2., 2.]]),
                       clf2.score_samples([[2., 2.]]))


def test_contamination():
    X = [[1, 1], [1, 0]]
    clf = neighbors.LocalOutlierFactor(contamination=0.6)
    assert_raises(ValueError, clf.fit, X)


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_novelty_errors():
    X = iris.data

    # check errors for novelty=False
    clf = neighbors.LocalOutlierFactor()
    clf.fit(X)
    # predict, decision_function and score_samples raise ValueError
    for method in ['predict', 'decision_function', 'score_samples']:
        msg = ('{} is not available when novelty=False'.format(method))
        assert_raises_regex(AttributeError, msg, getattr, clf, method)

    # check errors for novelty=True
    clf = neighbors.LocalOutlierFactor(novelty=True)
    msg = 'fit_predict is not available when novelty=True'
    assert_raises_regex(AttributeError, msg, getattr, clf, 'fit_predict')


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_novelty_training_scores():
    # check that the scores of the training samples are still accessible
    # when novelty=True through the negative_outlier_factor_ attribute
    X = iris.data

    # fit with novelty=False
    clf_1 = neighbors.LocalOutlierFactor()
    clf_1.fit(X)
    scores_1 = clf_1.negative_outlier_factor_

    # fit with novelty=True
    clf_2 = neighbors.LocalOutlierFactor(novelty=True)
    clf_2.fit(X)
    scores_2 = clf_2.negative_outlier_factor_

    assert_array_almost_equal(scores_1, scores_2)


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_hasattr_prediction():
    # check availability of prediction methods depending on novelty value.
    X = [[1, 1], [1, 2], [2, 1]]

    # when novelty=True
    clf = neighbors.LocalOutlierFactor(novelty=True)
    clf.fit(X)
    assert hasattr(clf, 'predict')
    assert hasattr(clf, 'decision_function')
    assert hasattr(clf, 'score_samples')
    assert not hasattr(clf, 'fit_predict')

    # when novelty=False
    clf = neighbors.LocalOutlierFactor(novelty=False)
    clf.fit(X)
    assert hasattr(clf, 'fit_predict')
    assert not hasattr(clf, 'predict')
    assert not hasattr(clf, 'decision_function')
    assert not hasattr(clf, 'score_samples')


@pytest.mark.filterwarnings(
    'ignore:default contamination parameter 0.1:FutureWarning')
# XXX: Remove in 0.22
def test_novelty_true_common_tests():

    # the common tests are run for the default LOF (novelty=False).
    # here we run these common tests for LOF when novelty=True
    check_estimator(neighbors.LocalOutlierFactor(novelty=True))


def test_contamination_future_warning():
    X = [[1, 1], [1, 2], [2, 1]]
    assert_warns_message(FutureWarning,
                         'default contamination parameter 0.1 will change '
                         'in version 0.22 to "auto"',
                         neighbors.LocalOutlierFactor().fit, X)
