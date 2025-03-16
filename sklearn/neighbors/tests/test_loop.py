# Authors: Erich Schubert <erich.schubert@tu-dortmund.de>
#          Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import re

import numpy as np
import pytest

from sklearn import metrics, neighbors
from sklearn.datasets import load_iris
from sklearn.metrics import roc_auc_score
from sklearn.utils import check_random_state
from sklearn.utils._testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_equal,
)
from sklearn.utils.estimator_checks import (
    check_outlier_corruption,
    parametrize_with_checks,
)
from sklearn.utils.fixes import CSR_CONTAINERS

# load the iris dataset
# and randomly permute it
rng = check_random_state(0)
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_loop(global_dtype):
    # Toy sample (the last two samples are outliers):
    X = np.asarray(
        [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1], [5, 3], [-4, 2]],
        dtype=global_dtype,
    )
    elki_scores = -1 * np.asarray(
        [
            0.0,
            0.0,
            0.12387384545791083,
            0.0,
            0.0,
            0.0,
            0.7713640481622056,
            0.5318645769367495,
        ]
    )
    elki_nplof = 0.7595192477675781

    # Test LocalOutlierProbabilities:
    clf = neighbors.LocalOutlierProbabilities(n_neighbors=5)
    score = clf.fit(X).negative_outlier_probabilities_
    assert_array_equal(clf._fit_X, X)

    # compare to ELKI results
    assert_almost_equal(clf._nplof, elki_nplof)
    assert_allclose(score, elki_scores)

    # Assert largest outlier score is smaller than smallest inlier score:
    assert np.min(score[:-2]) > np.max(score[-2:])

    # Assert predict() works:
    clf = neighbors.LocalOutlierProbabilities(contamination=0.25, n_neighbors=5).fit(X)
    expected_predictions = 6 * [1] + 2 * [-1]
    assert_array_equal(clf._predict(), expected_predictions)
    assert_array_equal(clf.fit_predict(X), expected_predictions)


def test_loop_performance(global_dtype):
    # Generate train/test data
    rng = check_random_state(2)
    X = 0.3 * rng.randn(120, 2).astype(global_dtype, copy=False)
    X_train = X[:100]

    # Generate some abnormal novel observations
    X_outliers = rng.uniform(low=-4, high=4, size=(20, 2)).astype(
        global_dtype, copy=False
    )
    X_test = np.r_[X[100:], X_outliers]
    y_test = np.array([0] * 20 + [1] * 20)

    # fit the model for novelty detection
    clf = neighbors.LocalOutlierProbabilities(novelty=True).fit(X_train)

    # predict scores (the lower, the more normal)
    y_pred = -clf.decision_function(X_test)

    # check that roc_auc is good
    assert roc_auc_score(y_test, y_pred) > 0.99


def test_loop_precomputed(global_dtype, random_state=42):
    """Tests LOF with a distance matrix."""
    # Note: smaller samples may result in spurious test success
    rng = np.random.RandomState(random_state)
    X = rng.random_sample((10, 4)).astype(global_dtype, copy=False)
    Y = rng.random_sample((3, 4)).astype(global_dtype, copy=False)
    DXX = metrics.pairwise_distances(X, metric="euclidean")
    DYX = metrics.pairwise_distances(Y, X, metric="euclidean")
    # As a feature matrix (n_samples by n_features)
    loop_X = neighbors.LocalOutlierProbabilities(n_neighbors=3, novelty=True)
    loop_X.fit(X)
    pred_X_X = loop_X._predict()
    pred_X_Y = loop_X.predict(Y)

    # As a dense distance matrix (n_samples by n_samples)
    loop_D = neighbors.LocalOutlierProbabilities(
        n_neighbors=3, algorithm="brute", metric="precomputed", novelty=True
    )
    loop_D.fit(DXX)
    pred_D_X = loop_D._predict()
    pred_D_Y = loop_D.predict(DYX)

    assert_allclose(pred_X_X, pred_D_X)
    assert_allclose(pred_X_Y, pred_D_Y)


def test_n_neighbors_attribute():
    X = iris.data
    clf = neighbors.LocalOutlierProbabilities(n_neighbors=500).fit(X)
    assert clf.n_neighbors_ == X.shape[0] - 1

    clf = neighbors.LocalOutlierProbabilities(n_neighbors=500)
    msg = "n_neighbors will be set to (n_samples - 1)"
    with pytest.warns(UserWarning, match=re.escape(msg)):
        clf.fit(X)
    assert clf.n_neighbors_ == X.shape[0] - 1


def test_score_samples(global_dtype):
    X_train = np.asarray([[1, 1], [1, 2], [2, 1]], dtype=global_dtype)
    X_test = np.asarray([[2.0, 2.0]], dtype=global_dtype)
    clf1 = neighbors.LocalOutlierProbabilities(
        n_neighbors=2, contamination=0.1, novelty=True
    ).fit(X_train)
    clf2 = neighbors.LocalOutlierProbabilities(n_neighbors=2, novelty=True).fit(X_train)

    clf1_scores = clf1.score_samples(X_test)
    clf1_decisions = clf1.decision_function(X_test)

    clf2_scores = clf2.score_samples(X_test)
    clf2_decisions = clf2.decision_function(X_test)

    assert_allclose(
        clf1_scores,
        clf1_decisions + clf1.offset_,
    )
    assert_allclose(
        clf2_scores,
        clf2_decisions + clf2.offset_,
    )
    assert_allclose(clf1_scores, clf2_scores)


def test_novelty_errors():
    X = iris.data

    # check errors for novelty=False
    clf = neighbors.LocalOutlierProbabilities()
    clf.fit(X)
    # predict, decision_function and score_samples raise ValueError
    for method in ["predict", "decision_function", "score_samples"]:
        msg = "{} is not available when novelty=False".format(method)
        with pytest.raises(AttributeError, match=msg):
            getattr(clf, method)

    # check errors for novelty=True
    clf = neighbors.LocalOutlierProbabilities(novelty=True)
    msg = "fit_predict is not available when novelty=True"
    with pytest.raises(AttributeError, match=msg):
        getattr(clf, "fit_predict")


def test_novelty_training_scores(global_dtype):
    # check that the scores of the training samples are still accessible
    # when novelty=True through the negative_outlier_probabilities_ attribute
    X = iris.data.astype(global_dtype)

    # fit with novelty=False
    clf_1 = neighbors.LocalOutlierProbabilities()
    clf_1.fit(X)
    scores_1 = clf_1.negative_outlier_probabilities_

    # fit with novelty=True
    clf_2 = neighbors.LocalOutlierProbabilities(novelty=True)
    clf_2.fit(X)
    scores_2 = clf_2.negative_outlier_probabilities_

    assert_allclose(scores_1, scores_2)


def test_hasattr_prediction():
    # check availability of prediction methods depending on novelty value.
    X = [[1, 1], [1, 2], [2, 1]]

    # when novelty=True
    clf = neighbors.LocalOutlierProbabilities(novelty=True)
    clf.fit(X)
    assert hasattr(clf, "predict")
    assert hasattr(clf, "decision_function")
    assert hasattr(clf, "score_samples")
    assert not hasattr(clf, "fit_predict")

    # when novelty=False
    clf = neighbors.LocalOutlierProbabilities(novelty=False)
    clf.fit(X)
    assert hasattr(clf, "fit_predict")
    assert not hasattr(clf, "predict")
    assert not hasattr(clf, "decision_function")
    assert not hasattr(clf, "score_samples")


@parametrize_with_checks([neighbors.LocalOutlierProbabilities(novelty=True)])
def test_novelty_true_common_tests(estimator, check):
    # the common tests are run for the default LOF (novelty=False).
    # here we run these common tests for LOF when novelty=True
    check(estimator)


@pytest.mark.parametrize("expected_outliers", [30, 53])
def test_predicted_outlier_number(expected_outliers):
    # the number of predicted outliers should be equal to the number of
    # expected outliers unless there are ties in the abnormality scores.
    X = iris.data
    n_samples = X.shape[0]
    contamination = float(expected_outliers) / n_samples

    clf = neighbors.LocalOutlierProbabilities(contamination=contamination)
    y_pred = clf.fit_predict(X)

    num_outliers = np.sum(y_pred != 1)
    if num_outliers != expected_outliers:
        y_dec = clf.negative_outlier_probabilities_
        check_outlier_corruption(num_outliers, expected_outliers, y_dec)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_sparse(csr_container):
    # LocalOutlierProbabilities must support CSR inputs
    # TODO: compare results on dense and sparse data as proposed in:
    # https://github.com/scikit-learn/scikit-learn/pull/23585#discussion_r968388186
    X = csr_container(iris.data)

    lof = neighbors.LocalOutlierProbabilities(novelty=True)
    lof.fit(X)
    lof.predict(X)
    lof.score_samples(X)
    lof.decision_function(X)

    lof = neighbors.LocalOutlierProbabilities(novelty=False)
    lof.fit_predict(X)


@pytest.mark.parametrize("algorithm", ["auto", "ball_tree", "kd_tree", "brute"])
@pytest.mark.parametrize("novelty", [True, False])
@pytest.mark.parametrize("contamination", [0.5, "auto"])
def test_loop_input_dtype_preservation(global_dtype, algorithm, contamination, novelty):
    """Check that the fitted attributes are stored using the data type of X."""
    X = iris.data.astype(global_dtype, copy=False)

    iso = neighbors.LocalOutlierProbabilities(
        n_neighbors=5, algorithm=algorithm, contamination=contamination, novelty=novelty
    )
    iso.fit(X)

    assert iso.negative_outlier_probabilities_.dtype == global_dtype

    for method in ("score_samples", "decision_function"):
        if hasattr(iso, method):
            y_pred = getattr(iso, method)(X)
            assert y_pred.dtype == global_dtype


@pytest.mark.parametrize("algorithm", ["auto", "ball_tree", "kd_tree", "brute"])
@pytest.mark.parametrize("novelty", [True, False])
@pytest.mark.parametrize("contamination", [0.5, "auto"])
def test_loop_dtype_equivalence(algorithm, novelty, contamination):
    """Check the equivalence of the results with 32 and 64 bits input."""

    inliers = iris.data[:50]  # setosa iris are really distinct from others
    outliers = iris.data[-5:]  # virginica will be considered as outliers
    # lower the precision of the input data to check that we have an equivalence when
    # making the computation in 32 and 64 bits.
    X = np.concatenate([inliers, outliers], axis=0).astype(np.float32)

    loop_32 = neighbors.LocalOutlierProbabilities(
        algorithm=algorithm, novelty=novelty, contamination=contamination
    )
    X_32 = X.astype(np.float32, copy=True)
    loop_32.fit(X_32)

    loop_64 = neighbors.LocalOutlierProbabilities(
        algorithm=algorithm, novelty=novelty, contamination=contamination
    )
    X_64 = X.astype(np.float64, copy=True)
    loop_64.fit(X_64)

    assert_allclose(
        loop_32.negative_outlier_probabilities_, loop_64.negative_outlier_probabilities_
    )

    for method in ("score_samples", "decision_function", "predict", "fit_predict"):
        if hasattr(loop_32, method):
            y_pred_32 = getattr(loop_32, method)(X_32)
            y_pred_64 = getattr(loop_64, method)(X_64)
            assert_allclose(y_pred_32, y_pred_64, atol=0.0002)
