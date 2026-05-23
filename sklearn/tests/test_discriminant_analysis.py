import warnings

import numpy as np
import pytest
from scipy import linalg

from sklearn.cluster import KMeans
from sklearn.covariance import LedoitWolf, ShrunkCovariance, ledoit_wolf
from sklearn.datasets import make_blobs
from sklearn.discriminant_analysis import (
    LinearDiscriminantAnalysis,
    QuadraticDiscriminantAnalysis,
    _cov,
)
from sklearn.model_selection import ShuffleSplit, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.utils import check_random_state
from sklearn.utils._testing import (
    _convert_container,
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_equal,
)

# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]], dtype="f")
y = np.array([1, 1, 1, 2, 2, 2])
y3 = np.array([1, 1, 2, 2, 3, 3])

# Degenerate data with only one feature (still should be separable)
X1 = np.array(
    [[-2], [-1], [-1], [1], [1], [2]],
    dtype="f",
)

# Data is just 9 separable points in the plane
X6 = np.array(
    [[0, 0], [-2, -2], [-2, -1], [-1, -1], [-1, -2], [1, 3], [1, 2], [2, 1], [2, 2]]
)
y6 = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2])
y7 = np.array([1, 2, 3, 2, 3, 1, 2, 3, 1])

# Degenerate data with 1 feature (still should be separable)
X7 = np.array([[-3], [-2], [-1], [-1], [0], [1], [1], [2], [3]])

# Data that has zero variance in one dimension and needs regularization
X2 = np.array(
    [[-3, 0], [-2, 0], [-1, 0], [-1, 0], [0, 0], [1, 0], [1, 0], [2, 0], [3, 0]]
)

# One element class
y4 = np.array([1, 1, 1, 1, 1, 1, 1, 1, 2])

solver_shrinkage = [
    ("svd", None),
    ("lsqr", None),
    ("eigen", None),
    ("lsqr", "auto"),
    ("lsqr", 0),
    ("lsqr", 0.43),
    ("eigen", "auto"),
    ("eigen", 0),
    ("eigen", 0.43),
]


def test_lda_predict():
    # Test LDA classification.
    # This checks that LDA implements fit and predict and returns correct
    # values for simple toy data.
    for test_case in solver_shrinkage:
        solver, shrinkage = test_case
        clf = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
        y_pred = clf.fit(X, y).predict(X)
        assert_array_equal(y_pred, y, "solver %s" % solver)

        # Assert that it works with 1D data
        y_pred1 = clf.fit(X1, y).predict(X1)
        assert_array_equal(y_pred1, y, "solver %s" % solver)

        # Test probability estimates
        y_proba_pred1 = clf.predict_proba(X1)
        assert_array_equal((y_proba_pred1[:, 1] > 0.5) + 1, y, "solver %s" % solver)
        y_log_proba_pred1 = clf.predict_log_proba(X1)
        assert_allclose(
            np.exp(y_log_proba_pred1),
            y_proba_pred1,
            rtol=1e-6,
            atol=1e-6,
            err_msg="solver %s" % solver,
        )

        # Primarily test for commit 2f34950 -- "reuse" of priors
        y_pred3 = clf.fit(X, y3).predict(X)
        # LDA shouldn't be able to separate those
        assert np.any(y_pred3 != y3), "solver %s" % solver

    clf = LinearDiscriminantAnalysis(solver="svd", shrinkage="auto")
    with pytest.raises(NotImplementedError):
        clf.fit(X, y)

    clf = LinearDiscriminantAnalysis(
        solver="lsqr", shrinkage=0.1, covariance_estimator=ShrunkCovariance()
    )
    with pytest.raises(
        ValueError,
        match=(
            "covariance_estimator and shrinkage "
            "parameters are not None. "
            "Only one of the two can be set."
        ),
    ):
        clf.fit(X, y)

    # test bad solver with covariance_estimator
    clf = LinearDiscriminantAnalysis(solver="svd", covariance_estimator=LedoitWolf())
    with pytest.raises(
        ValueError, match="covariance estimator is not supported with svd"
    ):
        clf.fit(X, y)

    # test bad covariance estimator
    clf = LinearDiscriminantAnalysis(
        solver="lsqr", covariance_estimator=KMeans(n_clusters=2, n_init="auto")
    )
    with pytest.raises(ValueError):
        clf.fit(X, y)


@pytest.mark.parametrize("n_classes", [2, 3])
@pytest.mark.parametrize("solver", ["svd", "lsqr", "eigen"])
def test_lda_predict_proba(solver, n_classes):
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
    lda = LinearDiscriminantAnalysis(
        solver=solver, store_covariance=True, shrinkage=None
    ).fit(X, y)
    # check that the empirical means and covariances are close enough to the
    # one used to generate the data
    assert_allclose(lda.means_, blob_centers, atol=1e-1)
    assert_allclose(lda.covariance_, blob_stds[0], atol=1)

    # implement the method to compute the probability given in The Elements
    # of Statistical Learning (cf. p.127, Sect. 4.4.5 "Logistic Regression
    # or LDA?")
    precision = linalg.inv(blob_stds[0])
    alpha_k = []
    alpha_k_0 = []
    for clazz in range(len(blob_centers) - 1):
        alpha_k.append(
            np.dot(precision, (blob_centers[clazz] - blob_centers[-1])[:, np.newaxis])
        )
        alpha_k_0.append(
            np.dot(
                -0.5 * (blob_centers[clazz] + blob_centers[-1])[np.newaxis, :],
                alpha_k[-1],
            )
        )

    sample = np.array([[-22, 22]])

    def discriminant_func(sample, coef, intercept, clazz):
        return np.exp(intercept[clazz] + np.dot(sample, coef[clazz])).item()

    prob = np.array(
        [
            float(
                discriminant_func(sample, alpha_k, alpha_k_0, clazz)
                / (
                    1
                    + sum(
                        [
                            discriminant_func(sample, alpha_k, alpha_k_0, clazz)
                            for clazz in range(n_classes - 1)
                        ]
                    )
                )
            )
            for clazz in range(n_classes - 1)
        ]
    )

    prob_ref = 1 - np.sum(prob)

    # check the consistency of the computed probability
    # all probabilities should sum to one
    prob_ref_2 = float(
        1
        / (
            1
            + sum(
                [
                    discriminant_func(sample, alpha_k, alpha_k_0, clazz)
                    for clazz in range(n_classes - 1)
                ]
            )
        )
    )

    assert prob_ref == pytest.approx(prob_ref_2)
    # check that the probability of LDA are close to the theoretical
    # probabilities
    assert_allclose(
        lda.predict_proba(sample), np.hstack([prob, prob_ref])[np.newaxis], atol=1e-2
    )


def test_lda_priors():
    # Test priors (negative priors)
    priors = np.array([0.5, -0.5])
    clf = LinearDiscriminantAnalysis(priors=priors)
    msg = "priors must be non-negative"

    with pytest.raises(ValueError, match=msg):
        clf.fit(X, y)

    # Test that priors passed as a list are correctly handled (run to see if
    # failure)
    clf = LinearDiscriminantAnalysis(priors=[0.5, 0.5])
    clf.fit(X, y)

    # Test that priors always sum to 1
    priors = np.array([0.5, 0.6])
    prior_norm = np.array([0.45, 0.55])
    clf = LinearDiscriminantAnalysis(priors=priors)

    with pytest.warns(UserWarning):
        clf.fit(X, y)

    assert_array_almost_equal(clf.priors_, prior_norm, 2)


def test_lda_coefs():
    # Test if the coefficients of the solvers are approximately the same.
    n_features = 2
    n_classes = 2
    n_samples = 1000
    X, y = make_blobs(
        n_samples=n_samples, n_features=n_features, centers=n_classes, random_state=11
    )

    clf_lda_svd = LinearDiscriminantAnalysis(solver="svd")
    clf_lda_lsqr = LinearDiscriminantAnalysis(solver="lsqr")
    clf_lda_eigen = LinearDiscriminantAnalysis(solver="eigen")

    clf_lda_svd.fit(X, y)
    clf_lda_lsqr.fit(X, y)
    clf_lda_eigen.fit(X, y)

    assert_array_almost_equal(clf_lda_svd.coef_, clf_lda_lsqr.coef_, 1)
    assert_array_almost_equal(clf_lda_svd.coef_, clf_lda_eigen.coef_, 1)
    assert_array_almost_equal(clf_lda_eigen.coef_, clf_lda_lsqr.coef_, 1)


def test_lda_transform():
    # Test LDA transform.
    clf = LinearDiscriminantAnalysis(solver="svd", n_components=1)
    X_transformed = clf.fit(X, y).transform(X)
    assert X_transformed.shape[1] == 1
    clf = LinearDiscriminantAnalysis(solver="eigen", n_components=1)
    X_transformed = clf.fit(X, y).transform(X)
    assert X_transformed.shape[1] == 1

    clf = LinearDiscriminantAnalysis(solver="lsqr", n_components=1)
    clf.fit(X, y)
    msg = "transform not implemented for 'lsqr'"

    with pytest.raises(NotImplementedError, match=msg):
        clf.transform(X)


def test_lda_explained_variance_ratio():
    # Test if the sum of the normalized eigen vectors values equals 1,
    # Also tests whether the explained_variance_ratio_ formed by the
    # eigen solver is the same as the explained_variance_ratio_ formed
    # by the svd solver

    state = np.random.RandomState(0)
    X = state.normal(loc=0, scale=100, size=(40, 20))
    y = state.randint(0, 3, size=(40,))

    clf_lda_eigen = LinearDiscriminantAnalysis(solver="eigen")
    clf_lda_eigen.fit(X, y)
    assert_almost_equal(clf_lda_eigen.explained_variance_ratio_.sum(), 1.0, 3)
    assert clf_lda_eigen.explained_variance_ratio_.shape == (2,), (
        "Unexpected length for explained_variance_ratio_"
    )

    clf_lda_svd = LinearDiscriminantAnalysis(solver="svd")
    clf_lda_svd.fit(X, y)
    assert_almost_equal(clf_lda_svd.explained_variance_ratio_.sum(), 1.0, 3)
    assert clf_lda_svd.explained_variance_ratio_.shape == (2,), (
        "Unexpected length for explained_variance_ratio_"
    )

    assert_array_almost_equal(
        clf_lda_svd.explained_variance_ratio_, clf_lda_eigen.explained_variance_ratio_
    )


def test_lda_orthogonality():
    # arrange four classes with their means in a kite-shaped pattern
    # the longer distance should be transformed to the first component, and
    # the shorter distance to the second component.
    means = np.array([[0, 0, -1], [0, 2, 0], [0, -2, 0], [0, 0, 5]])

    # We construct perfectly symmetric distributions, so the LDA can estimate
    # precise means.
    scatter = np.array(
        [
            [0.1, 0, 0],
            [-0.1, 0, 0],
            [0, 0.1, 0],
            [0, -0.1, 0],
            [0, 0, 0.1],
            [0, 0, -0.1],
        ]
    )

    X = (means[:, np.newaxis, :] + scatter[np.newaxis, :, :]).reshape((-1, 3))
    y = np.repeat(np.arange(means.shape[0]), scatter.shape[0])

    # Fit LDA and transform the means
    clf = LinearDiscriminantAnalysis(solver="svd").fit(X, y)
    means_transformed = clf.transform(means)

    d1 = means_transformed[3] - means_transformed[0]
    d2 = means_transformed[2] - means_transformed[1]
    d1 /= np.sqrt(np.sum(d1**2))
    d2 /= np.sqrt(np.sum(d2**2))

    # the transformed within-class covariance should be the identity matrix
    assert_almost_equal(np.cov(clf.transform(scatter).T), np.eye(2))

    # the means of classes 0 and 3 should lie on the first component
    assert_almost_equal(np.abs(np.dot(d1[:2], [1, 0])), 1.0)

    # the means of classes 1 and 2 should lie on the second component
    assert_almost_equal(np.abs(np.dot(d2[:2], [0, 1])), 1.0)


def test_lda_scaling():
    # Test if classification works correctly with differently scaled features.
    n = 100
    rng = np.random.RandomState(1234)
    # use uniform distribution of features to make sure there is absolutely no
    # overlap between classes.
    x1 = rng.uniform(-1, 1, (n, 3)) + [-10, 0, 0]
    x2 = rng.uniform(-1, 1, (n, 3)) + [10, 0, 0]
    x = np.vstack((x1, x2)) * [1, 100, 10000]
    y = [-1] * n + [1] * n

    for solver in ("svd", "lsqr", "eigen"):
        clf = LinearDiscriminantAnalysis(solver=solver)
        # should be able to separate the data perfectly
        assert clf.fit(x, y).score(x, y) == 1.0, "using covariance: %s" % solver


def test_lda_store_covariance():
    # Test for solver 'lsqr' and 'eigen'
    # 'store_covariance' has no effect on 'lsqr' and 'eigen' solvers
    for solver in ("lsqr", "eigen"):
        clf = LinearDiscriminantAnalysis(solver=solver).fit(X6, y6)
        assert hasattr(clf, "covariance_")

        # Test the actual attribute:
        clf = LinearDiscriminantAnalysis(solver=solver, store_covariance=True).fit(
            X6, y6
        )
        assert hasattr(clf, "covariance_")

        assert_array_almost_equal(
            clf.covariance_, np.array([[0.422222, 0.088889], [0.088889, 0.533333]])
        )

    # Test for SVD solver, the default is to not set the covariances_ attribute
    clf = LinearDiscriminantAnalysis(solver="svd").fit(X6, y6)
    assert not hasattr(clf, "covariance_")

    # Test the actual attribute:
    clf = LinearDiscriminantAnalysis(solver=solver, store_covariance=True).fit(X6, y6)
    assert hasattr(clf, "covariance_")

    assert_array_almost_equal(
        clf.covariance_, np.array([[0.422222, 0.088889], [0.088889, 0.533333]])
    )


@pytest.mark.parametrize("seed", range(10))
def test_lda_shrinkage(seed):
    # Test that shrunk covariance estimator and shrinkage parameter behave the
    # same
    rng = np.random.RandomState(seed)
    X = rng.rand(100, 10)
    y = rng.randint(3, size=(100))
    c1 = LinearDiscriminantAnalysis(store_covariance=True, shrinkage=0.5, solver="lsqr")
    c2 = LinearDiscriminantAnalysis(
        store_covariance=True,
        covariance_estimator=ShrunkCovariance(shrinkage=0.5),
        solver="lsqr",
    )
    c1.fit(X, y)
    c2.fit(X, y)
    assert_allclose(c1.means_, c2.means_)
    assert_allclose(c1.covariance_, c2.covariance_)


def test_lda_ledoitwolf():
    # When shrinkage="auto" current implementation uses ledoitwolf estimation
    # of covariance after standardizing the data. This checks that it is indeed
    # the case
    class StandardizedLedoitWolf:
        def fit(self, X):
            sc = StandardScaler()  # standardize features
            X_sc = sc.fit_transform(X)
            s = ledoit_wolf(X_sc)[0]
            # rescale
            s = sc.scale_[:, np.newaxis] * s * sc.scale_[np.newaxis, :]
            self.covariance_ = s

    rng = np.random.RandomState(0)
    X = rng.rand(100, 10)
    y = rng.randint(3, size=(100,))
    c1 = LinearDiscriminantAnalysis(
        store_covariance=True, shrinkage="auto", solver="lsqr"
    )
    c2 = LinearDiscriminantAnalysis(
        store_covariance=True,
        covariance_estimator=StandardizedLedoitWolf(),
        solver="lsqr",
    )
    c1.fit(X, y)
    c2.fit(X, y)
    assert_allclose(c1.means_, c2.means_)
    assert_allclose(c1.covariance_, c2.covariance_)


@pytest.mark.parametrize("n_features", [3, 5])
@pytest.mark.parametrize("n_classes", [5, 3])
def test_lda_dimension_warning(n_classes, n_features):
    rng = check_random_state(0)
    n_samples = 10
    X = rng.randn(n_samples, n_features)
    # we create n_classes labels by repeating and truncating a
    # range(n_classes) until n_samples
    y = np.tile(range(n_classes), n_samples // n_classes + 1)[:n_samples]
    max_components = min(n_features, n_classes - 1)

    for n_components in [max_components - 1, None, max_components]:
        # if n_components <= min(n_classes - 1, n_features), no warning
        lda = LinearDiscriminantAnalysis(n_components=n_components)
        lda.fit(X, y)

    for n_components in [max_components + 1, max(n_features, n_classes - 1) + 1]:
        # if n_components > min(n_classes - 1, n_features), raise error.
        # We test one unit higher than max_components, and then something
        # larger than both n_features and n_classes - 1 to ensure the test
        # works for any value of n_component
        lda = LinearDiscriminantAnalysis(n_components=n_components)
        msg = "n_components cannot be larger than "
        with pytest.raises(ValueError, match=msg):
            lda.fit(X, y)


@pytest.mark.parametrize(
    "data_type, expected_type",
    [
        (np.float32, np.float32),
        (np.float64, np.float64),
        (np.int32, np.float64),
        (np.int64, np.float64),
    ],
)
def test_lda_dtype_match(data_type, expected_type):
    for solver, shrinkage in solver_shrinkage:
        clf = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
        clf.fit(X.astype(data_type), y.astype(data_type))
        assert clf.coef_.dtype == expected_type


def test_lda_numeric_consistency_float32_float64():
    for solver, shrinkage in solver_shrinkage:
        clf_32 = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
        clf_32.fit(X.astype(np.float32), y.astype(np.float32))
        clf_64 = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
        clf_64.fit(X.astype(np.float64), y.astype(np.float64))

        # Check value consistency between types
        rtol = 1e-6
        assert_allclose(clf_32.coef_, clf_64.coef_, rtol=rtol)


@pytest.mark.parametrize("solver", ["svd", "eigen"])
def test_qda(solver):
    # QDA classification.
    # This checks that QDA implements fit and predict and returns
    # correct values for a simple toy dataset.
    clf = QuadraticDiscriminantAnalysis(solver=solver)
    y_pred = clf.fit(X6, y6).predict(X6)
    assert_array_equal(y_pred, y6)

    # Assure that it works with 1D data
    y_pred1 = clf.fit(X7, y6).predict(X7)
    assert_array_equal(y_pred1, y6)

    # Test probas estimates
    y_proba_pred1 = clf.predict_proba(X7)
    assert_array_equal((y_proba_pred1[:, 1] > 0.5) + 1, y6)
    y_log_proba_pred1 = clf.predict_log_proba(X7)
    assert_array_almost_equal(np.exp(y_log_proba_pred1), y_proba_pred1, 8)

    y_pred3 = clf.fit(X6, y7).predict(X6)
    # QDA shouldn't be able to separate those
    assert np.any(y_pred3 != y7)

    # Classes should have at least 2 elements
    with pytest.raises(ValueError):
        clf.fit(X6, y4)


def test_qda_covariance_estimator():
    # Test that the correct errors are raised when using inappropriate
    # covariance estimators or shrinkage parameters with QDA.
    clf = QuadraticDiscriminantAnalysis(solver="svd", shrinkage="auto")
    with pytest.raises(NotImplementedError):
        clf.fit(X, y)

    clf = QuadraticDiscriminantAnalysis(
        solver="eigen", shrinkage=0.1, covariance_estimator=ShrunkCovariance()
    )
    with pytest.raises(
        ValueError,
        match=(
            "covariance_estimator and shrinkage parameters are not None. "
            "Only one of the two can be set."
        ),
    ):
        clf.fit(X, y)

    # test bad solver with covariance_estimator
    clf = QuadraticDiscriminantAnalysis(solver="svd", covariance_estimator=LedoitWolf())
    with pytest.raises(
        ValueError, match="covariance_estimator is not supported with solver='svd'"
    ):
        clf.fit(X, y)

    # test bad covariance estimator
    clf = QuadraticDiscriminantAnalysis(
        solver="eigen", covariance_estimator=KMeans(n_clusters=2, n_init="auto")
    )
    with pytest.raises(ValueError):
        clf.fit(X, y)


def test_qda_ledoitwolf(global_random_seed):
    # When shrinkage="auto" current implementation uses ledoitwolf estimation
    # of covariance after standardizing the data. This checks that it is indeed
    # the case
    class StandardizedLedoitWolf:
        def fit(self, X):
            sc = StandardScaler()  # standardize features
            X_sc = sc.fit_transform(X)
            s = ledoit_wolf(X_sc)[0]
            # rescale
            s = sc.scale_[:, np.newaxis] * s * sc.scale_[np.newaxis, :]
            self.covariance_ = s

    rng = np.random.RandomState(global_random_seed)
    X = rng.rand(100, 10)
    y = rng.randint(3, size=(100,))
    c1 = QuadraticDiscriminantAnalysis(
        store_covariance=True, shrinkage="auto", solver="eigen"
    )
    c2 = QuadraticDiscriminantAnalysis(
        store_covariance=True,
        covariance_estimator=StandardizedLedoitWolf(),
        solver="eigen",
    )
    c1.fit(X, y)
    c2.fit(X, y)
    assert_allclose(c1.means_, c2.means_)
    assert_allclose(c1.covariance_, c2.covariance_)


def test_qda_coefs(global_random_seed):
    # Test if the coefficients of the solvers are approximately the same.
    n_features = 2
    n_classes = 2
    n_samples = 3000
    X, y = make_blobs(
        n_samples=n_samples,
        n_features=n_features,
        centers=n_classes,
        cluster_std=[1.0, 3.0],
        random_state=global_random_seed,
    )

    clf_svd = QuadraticDiscriminantAnalysis(solver="svd")
    clf_eigen = QuadraticDiscriminantAnalysis(solver="eigen")

    clf_svd.fit(X, y)
    clf_eigen.fit(X, y)

    for class_idx in range(n_classes):
        assert_allclose(
            np.abs(clf_svd.rotations_[class_idx]),
            np.abs(clf_eigen.rotations_[class_idx]),
            rtol=1e-3,
            err_msg=f"SVD and Eigen rotations differ for class {class_idx}",
        )
        assert_allclose(
            clf_svd.scalings_[class_idx],
            clf_eigen.scalings_[class_idx],
            rtol=1e-3,
            err_msg=f"SVD and Eigen scalings differ for class {class_idx}",
        )


def test_qda_priors():
    clf = QuadraticDiscriminantAnalysis()
    y_pred = clf.fit(X6, y6).predict(X6)
    n_pos = np.sum(y_pred == 2)

    neg = 1e-10
    clf = QuadraticDiscriminantAnalysis(priors=np.array([neg, 1 - neg]))
    y_pred = clf.fit(X6, y6).predict(X6)
    n_pos2 = np.sum(y_pred == 2)

    assert n_pos2 > n_pos


@pytest.mark.parametrize("priors_type", ["list", "tuple", "array"])
def test_qda_prior_type(priors_type):
    """Check that priors accept array-like."""
    priors = [0.5, 0.5]
    clf = QuadraticDiscriminantAnalysis(
        priors=_convert_container([0.5, 0.5], priors_type)
    ).fit(X6, y6)
    assert isinstance(clf.priors_, np.ndarray)
    assert_array_equal(clf.priors_, priors)


def test_qda_prior_copy():
    """Check that altering `priors` without `fit` doesn't change `priors_`"""
    priors = np.array([0.5, 0.5])
    qda = QuadraticDiscriminantAnalysis(priors=priors).fit(X, y)

    # we expect the following
    assert_array_equal(qda.priors_, qda.priors)

    # altering `priors` without `fit` should not change `priors_`
    priors[0] = 0.2
    assert qda.priors_[0] != qda.priors[0]


def test_qda_store_covariance():
    # The default is to not set the covariances_ attribute
    clf = QuadraticDiscriminantAnalysis().fit(X6, y6)
    assert not hasattr(clf, "covariance_")

    # Test the actual attribute:
    clf = QuadraticDiscriminantAnalysis(store_covariance=True).fit(X6, y6)
    assert hasattr(clf, "covariance_")

    assert_array_almost_equal(clf.covariance_[0], np.array([[0.7, 0.45], [0.45, 0.7]]))

    assert_array_almost_equal(
        clf.covariance_[1],
        np.array([[0.33333333, -0.33333333], [-0.33333333, 0.66666667]]),
    )


@pytest.mark.parametrize("solver", ["svd", "eigen"])
def test_qda_regularization(global_random_seed, solver):
    # The default is reg_param=0. and will cause issues when there is a
    # constant variable.
    rng = np.random.default_rng(global_random_seed)

    # Fitting on data with constant variable without regularization
    # triggers a LinAlgError.
    msg = r"The covariance matrix of class .+ is not full rank."
    clf = QuadraticDiscriminantAnalysis(solver=solver)
    with pytest.raises(linalg.LinAlgError, match=msg):
        clf.fit(X2, y6)

    with pytest.raises(AttributeError):
        y_pred = clf.predict(X2)

    # Adding a little regularization fixes the fit time error.
    if solver == "svd":
        clf = QuadraticDiscriminantAnalysis(solver=solver, reg_param=0.01)
    elif solver == "eigen":
        clf = QuadraticDiscriminantAnalysis(solver=solver, shrinkage=0.01)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
    clf.fit(X2, y6)
    y_pred = clf.predict(X2)
    assert_array_equal(y_pred, y6)

    # LinAlgError should also be there for the n_samples_in_a_class <
    # n_features case.
    X = rng.normal(size=(9, 4))
    y = np.array([1, 1, 1, 1, 1, 1, 2, 2, 2])

    clf = QuadraticDiscriminantAnalysis(solver=solver)
    if solver == "svd":
        msg2 = msg + " When using `solver='svd'`"
    elif solver == "eigen":
        msg2 = msg

    with pytest.raises(linalg.LinAlgError, match=msg2):
        clf.fit(X, y)

    # The error will persist even with regularization for SVD
    # because the number of singular values is limited by n_samples_in_a_class.
    if solver == "svd":
        clf = QuadraticDiscriminantAnalysis(solver=solver, reg_param=0.3)
        with pytest.raises(linalg.LinAlgError, match=msg2):
            clf.fit(X, y)
    # The warning will be gone for Eigen with regularization, because
    # the covariance matrix will be full-rank.
    elif solver == "eigen":
        clf = QuadraticDiscriminantAnalysis(solver=solver, shrinkage=0.3)
        clf.fit(X, y)


def test_covariance():
    x, y = make_blobs(n_samples=100, n_features=5, centers=1, random_state=42)

    # make features correlated
    x = np.dot(x, np.arange(x.shape[1] ** 2).reshape(x.shape[1], x.shape[1]))

    c_e = _cov(x, "empirical")
    assert_almost_equal(c_e, c_e.T)

    c_s = _cov(x, "auto")
    assert_almost_equal(c_s, c_s.T)


@pytest.mark.parametrize("solver", ["svd", "lsqr", "eigen"])
def test_raises_value_error_on_same_number_of_classes_and_samples(solver):
    """
    Tests that if the number of samples equals the number
    of classes, a ValueError is raised.
    """
    X = np.array([[0.5, 0.6], [0.6, 0.5]])
    y = np.array(["a", "b"])
    clf = LinearDiscriminantAnalysis(solver=solver)
    with pytest.raises(ValueError, match="The number of samples must be more"):
        clf.fit(X, y)


@pytest.mark.parametrize("solver", ["svd", "eigen"])
def test_raises_value_error_on_one_sample_per_class(solver):
    """
    Tests that if a class has one sample, a ValueError is raised.
    """
    X = np.array([[0.5, 0.6], [0.6, 0.5], [0.4, 0.4], [0.6, 0.5]])
    y = np.array(["a", "a", "a", "b"])
    clf = QuadraticDiscriminantAnalysis(solver=solver)
    with pytest.raises(ValueError, match="y has only 1 sample in class"):
        clf.fit(X, y)


def test_get_feature_names_out():
    """Check get_feature_names_out uses class name as prefix."""

    est = LinearDiscriminantAnalysis().fit(X, y)
    names_out = est.get_feature_names_out()

    class_name_lower = "LinearDiscriminantAnalysis".lower()
    expected_names_out = np.array(
        [
            f"{class_name_lower}{i}"
            for i in range(est.explained_variance_ratio_.shape[0])
        ],
        dtype=object,
    )
    assert_array_equal(names_out, expected_names_out)


@pytest.mark.parametrize("n_features", [25])
@pytest.mark.parametrize("train_size", [100])
@pytest.mark.parametrize("solver_no_shrinkage", ["svd", "eigen"])
def test_qda_shrinkage_performance(
    global_random_seed, n_features, train_size, solver_no_shrinkage
):
    # Test that QDA with shrinkage performs better than without shrinkage on
    # a case where there's a small number of samples per class relative to
    # the number of features.
    n_samples = 1000
    n_features = n_features

    rng = np.random.default_rng(global_random_seed)

    # Sample from two Gaussians with different variances and same null means.
    vars1 = rng.uniform(2.0, 3.0, size=n_features)
    vars2 = rng.uniform(0.2, 1.0, size=n_features)

    X = np.concatenate(
        [
            np.random.randn(n_samples // 2, n_features) * np.sqrt(vars1),
            np.random.randn(n_samples // 2, n_features) * np.sqrt(vars2),
        ],
        axis=0,
    )
    y = np.array([0] * (n_samples // 2) + [1] * (n_samples // 2))

    # Use small training sets to illustrate the regularization effect of
    # covariance shrinkage.
    cv = ShuffleSplit(n_splits=5, train_size=train_size, random_state=0)
    qda_shrinkage = QuadraticDiscriminantAnalysis(solver="eigen", shrinkage="auto")
    qda_no_shrinkage = QuadraticDiscriminantAnalysis(
        solver=solver_no_shrinkage, shrinkage=None
    )

    scores_no_shrinkage = cross_val_score(
        qda_no_shrinkage, X, y, cv=cv, scoring="d2_brier_score"
    )
    scores_shrinkage = cross_val_score(
        qda_shrinkage, X, y, cv=cv, scoring="d2_brier_score"
    )

    assert scores_shrinkage.mean() > 0.9
    assert scores_no_shrinkage.mean() < 0.6
