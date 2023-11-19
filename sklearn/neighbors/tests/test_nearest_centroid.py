"""
Testing for the nearest centroid module.
"""
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal
from scipy import sparse as sp

from sklearn import datasets
from sklearn.neighbors import NearestCentroid
from sklearn.utils import check_random_state

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
X_csr = sp.csr_matrix(X)  # Sparse matrix
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
T_nan = np.asarray([[-1, -1], [2, 2], [3, 2]], dtype=np.float64)
T_nan[0][0] = float("nan")
T_csr = sp.csr_matrix(T)
T_nan_csr = sp.csr_matrix(T_nan)
T_zero_var = [[1, 2], [1, 2], [1, 2], [1, 2], [1, 2], [1, 2]]
true_result = [-1, 1, 1]
true_result_prior1 = [-1, 1, 1]

true_discriminant_scores = [-32, 64, 80]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = np.random.RandomState(1)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_classification_toy():
    # Check classification on a toy dataset, including sparse versions.
    clf = NearestCentroid()
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.decision_function(T), true_discriminant_scores)

    # Test uniform priors
    clf = NearestCentroid(priors="uniform")
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.decision_function(T), true_discriminant_scores)

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


@pytest.mark.parametrize("n_classes", [2, 3])
def test_predict_proba(n_classes):
    # Fit and predict probability estimates
    # compare with results from pamr package
    def generate_dataset(n_samples, centers, covariances, random_state=None):
        """Generate a multivariate normal data given some centers and
        covariances"""
        rng = check_random_state(random_state)
        X = np.vstack(
            [
                rng.multivariate_normal(mean, cov, size=n_samples // len(centers))
                for mean, cov in zip(centers, covariances)
            ]
        )
        y = np.hstack(
            [[clazz] * (n_samples // len(centers)) for clazz in range(len(centers))]
        )
        return X, y

    blob_centers = np.array([[0, 0], [-10, 40], [-30, 30]])[:n_classes]
    blob_stds = np.array([[[10, 10], [10, 100]]] * len(blob_centers))
    X, y = generate_dataset(
        n_samples=90000, centers=blob_centers, covariances=blob_stds, random_state=42
    )
    clf = NearestCentroid().fit(X, y)
    probabilities = clf.predict_proba(X)
    assert probabilities.shape == (X.shape[0], n_classes)
    assert_array_almost_equal(probabilities.sum(axis=1), np.ones(X.shape[0]))
    assert 0 <= probabilities.all() <= 1


# TODO(1.5): Remove filterwarnings when support for some metrics is removed
@pytest.mark.filterwarnings("ignore:Support for distance metrics:FutureWarning:sklearn")
def test_iris():
    # Check consistency on dataset iris.
    for metric in ("euclidean", "cosine"):
        clf = NearestCentroid(metric=metric).fit(iris.data, iris.target)
        score = np.mean(clf.predict(iris.data) == iris.target)
        assert score > 0.9, "Failed with score = " + str(score)


# TODO(1.5): Remove filterwarnings when support for some metrics is removed
@pytest.mark.filterwarnings("ignore:Support for distance metrics:FutureWarning:sklearn")
def test_iris_shrinkage():
    # Check consistency on dataset iris, when using shrinkage.
    for metric in ("euclidean", "cosine"):
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


def test_manhattan_metric():
    # Test the manhattan metric.

    clf = NearestCentroid(metric="manhattan")
    clf.fit(X, y)
    dense_centroid = clf.centroids_
    clf.fit(X_csr, y)
    assert_array_equal(clf.centroids_, dense_centroid)
    assert_array_equal(dense_centroid, [[-1, -1], [1, 1]])


# TODO(1.5): remove this test
@pytest.mark.parametrize(
    "metric", sorted(list(NearestCentroid._valid_metrics - {"manhattan", "euclidean"}))
)
def test_deprecated_distance_metric_supports(metric):
    # Check that a warning is raised for all deprecated distance metric supports
    clf = NearestCentroid(metric=metric)
    with pytest.warns(
        FutureWarning,
        match="Support for distance metrics other than euclidean and manhattan",
    ):
        clf.fit(X, y)


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


def test_neg_priors():
    # Test negative priors.

    clf = NearestCentroid(priors=[-2, 4])
    with pytest.raises(ValueError):
        clf.fit(X, y)


def test_wrong_priors():
    # Test normalizing priors that don't sum to 1.
    clf = NearestCentroid(priors=[2, 4])
    with pytest.warns(
        UserWarning,
        match="The priors do not sum to 1. Normalizing such that it sums to one.",
    ):
        clf.fit(X, y)


def test_manhattan_decision_func_error():
    # Make sure error is raised when calling decision_function with
    # manhattan metric.

    clf = NearestCentroid(metric="manhattan", priors=[0.2, 0.8])
    clf.fit(X, y)
    with pytest.raises(AttributeError):
        clf.decision_function(T)


def test_manhattan_pred_proba_error():
    # Make sure error is raised when calling predict_proba with
    # manhattan metric.

    clf = NearestCentroid(metric="manhattan", priors=[0.2, 0.8])
    clf.fit(X, y)
    with pytest.raises(AttributeError):
        clf.predict_proba(T)


def test_manhattan_pred_log_proba_error():
    # Make sure error is raised when calling predict_log_proba with
    # manhattan metric.

    clf = NearestCentroid(metric="manhattan", priors=[0.2, 0.8])
    clf.fit(X, y)
    with pytest.raises(AttributeError):
        clf.predict_log_proba(T)


def test_zero_var():
    clf = NearestCentroid(priors=[0.2, 0.8])
    with pytest.raises(ValueError):
        clf.fit(T_zero_var, y)


def test_nan():
    clf = NearestCentroid(priors=[0.2, 0.8], shrink_threshold=0.5)
    clf.fit(X, y)
    with pytest.raises(ValueError):
        clf.decision_function(T_nan)


def test_nan_sparse():
    clf = NearestCentroid(priors=[0.2, 0.8], shrink_threshold=0.5)
    clf.fit(X, y)
    with pytest.raises(ValueError):
        clf.decision_function(T_nan_csr)
