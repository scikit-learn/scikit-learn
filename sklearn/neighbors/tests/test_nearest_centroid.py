"""
Testing for the nearest centroid module.
"""

import numpy as np
import pytest

from sklearn import datasets
from sklearn.neighbors import NearestCentroid
from sklearn.utils._testing import (
    assert_allclose,
    assert_array_almost_equal,
    assert_array_equal,
)
from sklearn.utils.fixes import CSR_CONTAINERS

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]
true_result_prior1 = [-1, 1, 1]

true_discriminant_scores = [-32, 64, 80]
true_proba = [[1, 1.26642e-14], [1.60381e-28, 1], [1.80485e-35, 1]]


# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = np.random.RandomState(1)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_classification_toy(csr_container):
    # Check classification on a toy dataset, including sparse versions.
    X_csr = csr_container(X)
    T_csr = csr_container(T)

    # Check classification on a toy dataset, including sparse versions.
    clf = NearestCentroid()
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.decision_function(T), true_discriminant_scores)
    assert_array_almost_equal(clf.predict_proba(T), true_proba)

    # Test uniform priors
    clf = NearestCentroid(priors="uniform")
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.decision_function(T), true_discriminant_scores)
    assert_array_almost_equal(clf.predict_proba(T), true_proba)

    clf = NearestCentroid(priors="empirical")
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.decision_function(T), true_discriminant_scores)
    assert_array_almost_equal(clf.predict_proba(T), true_proba)

    # Test custom priors
    clf = NearestCentroid(priors=[0.25, 0.75])
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result_prior1)

    # Same test, but with a sparse matrix to fit and test.
    clf = NearestCentroid()
    clf.fit(X_csr, y)
    assert_array_equal(clf.predict(T_csr), true_result)

    # Fit with sparse, test with non-sparse
    clf = NearestCentroid()
    clf.fit(X_csr, y)
    assert_array_equal(clf.predict(T), true_result)

    # Fit with non-sparse, test with sparse
    clf = NearestCentroid()
    clf.fit(X, y)
    assert_array_equal(clf.predict(T_csr), true_result)

    # Fit and predict with non-CSR sparse matrices
    clf = NearestCentroid()
    clf.fit(X_csr.tocoo(), y)
    assert_array_equal(clf.predict(T_csr.tolil()), true_result)


def test_iris():
    # Check consistency on dataset iris.
    for metric in ("euclidean", "manhattan"):
        clf = NearestCentroid(metric=metric).fit(iris.data, iris.target)
        score = np.mean(clf.predict(iris.data) == iris.target)
        assert score > 0.9, "Failed with score = " + str(score)


def test_iris_shrinkage():
    # Check consistency on dataset iris, when using shrinkage.
    for metric in ("euclidean", "manhattan"):
        for shrink_threshold in [None, 0.1, 0.5]:
            clf = NearestCentroid(metric=metric, shrink_threshold=shrink_threshold)
            clf = clf.fit(iris.data, iris.target)
            score = np.mean(clf.predict(iris.data) == iris.target)
            assert score > 0.8, "Failed with score = " + str(score)


def test_pickle():
    import pickle

    # classification
    obj = NearestCentroid()
    obj.fit(iris.data, iris.target)
    score = obj.score(iris.data, iris.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert type(obj2) == obj.__class__
    score2 = obj2.score(iris.data, iris.target)
    assert_array_equal(
        score,
        score2,
        "Failed to generate same score after pickling (classification).",
    )


def test_shrinkage_correct():
    # Ensure that the shrinking is correct.
    # The expected result is calculated by R (pamr),
    # which is implemented by the author of the original paper.
    # (One need to modify the code to output the new centroid in pamr.predict)

    X = np.array([[0, 1], [1, 0], [1, 1], [2, 0], [6, 8]])
    y = np.array([1, 1, 2, 2, 2])
    clf = NearestCentroid(shrink_threshold=0.1)
    clf.fit(X, y)
    expected_result = np.array([[0.7787310, 0.8545292], [2.814179, 2.763647]])
    np.testing.assert_array_almost_equal(clf.centroids_, expected_result)


def test_shrinkage_threshold_decoded_y():
    clf = NearestCentroid(shrink_threshold=0.01)
    y_ind = np.asarray(y)
    y_ind[y_ind == -1] = 0
    clf.fit(X, y_ind)
    centroid_encoded = clf.centroids_
    clf.fit(X, y)
    assert_array_equal(centroid_encoded, clf.centroids_)


def test_predict_translated_data():
    # Test that NearestCentroid gives same results on translated data

    rng = np.random.RandomState(0)
    X = rng.rand(50, 50)
    y = rng.randint(0, 3, 50)
    noise = rng.rand(50)
    clf = NearestCentroid(shrink_threshold=0.1)
    clf.fit(X, y)
    y_init = clf.predict(X)
    clf = NearestCentroid(shrink_threshold=0.1)
    X_noise = X + noise
    clf.fit(X_noise, y)
    y_translate = clf.predict(X_noise)
    assert_array_equal(y_init, y_translate)


@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_manhattan_metric(csr_container):
    # Test the manhattan metric.
    X_csr = csr_container(X)

    clf = NearestCentroid(metric="manhattan")
    clf.fit(X, y)
    dense_centroid = clf.centroids_
    clf.fit(X_csr, y)
    assert_array_equal(clf.centroids_, dense_centroid)
    assert_array_equal(dense_centroid, [[-1, -1], [1, 1]])


def test_features_zero_var():
    # Test that features with 0 variance throw error

    X = np.empty((10, 2))
    X[:, 0] = -0.13725701
    X[:, 1] = -0.9853293
    y = np.zeros((10))
    y[0] = 1

    clf = NearestCentroid(shrink_threshold=0.1)
    with pytest.raises(ValueError):
        clf.fit(X, y)


def test_negative_priors_error():
    """Check that we raise an error when the user-defined priors are negative."""
    clf = NearestCentroid(priors=[-2, 4])
    with pytest.raises(ValueError, match="priors must be non-negative"):
        clf.fit(X, y)


def test_warn_non_normalized_priors():
    """Check that we raise a warning and normalize the user-defined priors when they
    don't sum to 1.
    """
    priors = [2, 4]
    clf = NearestCentroid(priors=priors)
    with pytest.warns(
        UserWarning,
        match="The priors do not sum to 1. Normalizing such that it sums to one.",
    ):
        clf.fit(X, y)

    assert_allclose(clf.class_prior_, np.asarray(priors) / np.asarray(priors).sum())


@pytest.mark.parametrize(
    "response_method", ["decision_function", "predict_proba", "predict_log_proba"]
)
def test_method_not_available_with_manhattan(response_method):
    """Check that we raise an AttributeError with Manhattan metric when trying
    to call a non-thresholded response method.
    """
    clf = NearestCentroid(metric="manhattan").fit(X, y)
    with pytest.raises(AttributeError):
        getattr(clf, response_method)(T)


@pytest.mark.parametrize("array_constructor", [np.array] + CSR_CONTAINERS)
def test_error_zero_variances(array_constructor):
    """Check that we raise an error when the variance for all features is zero."""
    X = np.ones((len(y), 2))
    X[:, 1] *= 2
    X = array_constructor(X)

    clf = NearestCentroid()
    with pytest.raises(ValueError, match="All features have zero variance"):
        clf.fit(X, y)


def test_nearest_centroid_sample_weight():
    """Test that NearestCentroid correctly applies sample_weight."""
    # Toy dataset: 3 samples for class -1, 3 for class 1
    X_data = np.array(
        [[-2.0, -1.0], [-1.0, -1.0], [-1.0, -2.0], [1.0, 1.0], [1.0, 2.0], [2.0, 1.0]]
    )
    y_data = np.array([-1, -1, -1, 1, 1, 1])

    # -- Test 1: Uniform weights == unweighted --
    clf_no_weight = NearestCentroid()
    clf_no_weight.fit(X_data, y_data)

    clf_uniform = NearestCentroid()
    clf_uniform.fit(X_data, y_data, sample_weight=np.ones(6))

    np.testing.assert_array_almost_equal(
        clf_no_weight.centroids_, clf_uniform.centroids_
    )

    # -- Test 2: Heavily weight one extreme sample in class 1 --
    # Extreme sample: [2, 1] with weight 10. Other samples in class 1: weight 1.
    # Weighted mean for class 1 should shift toward [2, 1].
    weights_biased = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 10.0])
    clf_biased = NearestCentroid()
    clf_biased.fit(X_data, y_data, sample_weight=weights_biased)

    # The class 1 centroid should shift toward [2, 1] on both coordinates.
    # x-coordinate increases (toward 2), y-coordinate decreases (toward 1).
    assert (
        clf_biased.centroids_[1, 0] > clf_no_weight.centroids_[1, 0]
    ), "Centroid x for class 1 should shift right with high weight on [2,1]"
    assert (
        clf_biased.centroids_[1, 1] < clf_no_weight.centroids_[1, 1]
    ), "Centroid y for class 1 should shift down with high weight on [2,1]"

    # -- Test 3: Zero-weight class raises ValueError --
    weights_zero_class = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    clf_zero = NearestCentroid()
    with pytest.raises(ValueError, match="zero total sample weight"):
        clf_zero.fit(X_data, y_data, sample_weight=weights_zero_class)

    # -- Test 4: Sparse input with sample_weight gives same result as dense --
    for csr_container in CSR_CONTAINERS:
        X_sparse = csr_container(X_data)
        clf_dense = NearestCentroid()
        clf_dense.fit(X_data, y_data, sample_weight=weights_biased)

        clf_sparse = NearestCentroid()
        clf_sparse.fit(X_sparse, y_data, sample_weight=weights_biased)

        np.testing.assert_array_almost_equal(
            clf_dense.centroids_, clf_sparse.centroids_, decimal=5
        )


def test_nearest_centroid_weighted_shrinkage():
    """Test that shrinkage works correctly with sample_weight.

    A heavily weighted outlier should pull the shrunken centroid noticeably
    compared to the uniform-weight case, confirming the nk/total_weight
    m-parameter works end-to-end under shrinkage.
    """
    # Dataset matches test_shrinkage_correct for easy cross-referencing
    X_shrink = np.array(
        [[0.0, 1.0], [1.0, 0.0], [1.0, 1.0], [2.0, 0.0], [6.0, 8.0]]
    )
    y_shrink = np.array([1, 1, 2, 2, 2])

    # Uniform weights must reproduce unweighted result
    clf_no_weight = NearestCentroid(shrink_threshold=0.1)
    clf_no_weight.fit(X_shrink, y_shrink)

    clf_uniform = NearestCentroid(shrink_threshold=0.1)
    clf_uniform.fit(X_shrink, y_shrink, sample_weight=np.ones(5))

    np.testing.assert_allclose(
        clf_uniform.centroids_, clf_no_weight.centroids_, rtol=1e-5
    )

    # Heavily weight the outlier [6, 8] in class 2 (last sample, index 4).
    # The shrunken centroid for class 2 should shift toward [6, 8].
    weights_outlier = np.array([1.0, 1.0, 1.0, 1.0, 20.0])
    clf_weighted = NearestCentroid(shrink_threshold=0.1)
    clf_weighted.fit(X_shrink, y_shrink, sample_weight=weights_outlier)

    # Class index for label 2 is 1 (LabelEncoder sorts: 1 < 2)
    assert clf_weighted.centroids_[1, 0] > clf_uniform.centroids_[1, 0], (
        "Weighted shrunken centroid x for class 2 should move toward outlier [6,8]"
    )
    assert clf_weighted.centroids_[1, 1] > clf_uniform.centroids_[1, 1], (
        "Weighted shrunken centroid y for class 2 should move toward outlier [6,8]"
    )


def test_nearest_centroid_denominator_guard():
    """Test the within-class variance denominator guard.

    The pooled variance estimate (total_weight - n_classes) must be positive.
    This guard fires for degenerate cases like one unit-weight sample per class.
    """
    # Case 1: Exactly one unit-weight sample per class (total_weight = 2, n_classes = 2)
    # -> denominator = 0 -> should raise.
    X = np.array([[0.0, 0.0], [1.0, 1.0]])
    y = np.array([0, 1])

    clf = NearestCentroid(shrink_threshold=0.1)
    with pytest.raises(ValueError, match="effective degrees of freedom"):
        clf.fit(X, y, sample_weight=np.ones(2))

    # Case 2: Fractional weights summing to <= n_classes (total_weight = 1.0, n_classes = 2)
    # -> denominator = -1.0 -> should raise.
    with pytest.raises(ValueError, match="effective degrees of freedom"):
        clf.fit(X, y, sample_weight=np.array([0.5, 0.5]))

    # Case 3: Positive path. total_weight > n_classes (total_weight = 2.1, n_classes = 2)
    # -> denominator = 0.1 -> should pass.
    clf.fit(X, y, sample_weight=np.array([1.05, 1.05]))
    assert hasattr(clf, "centroids_")
