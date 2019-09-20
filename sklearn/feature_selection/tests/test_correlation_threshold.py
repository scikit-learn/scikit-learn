import pytest
import numpy as np
from scipy.sparse import csr_matrix
from numpy.testing import assert_allclose

from sklearn.feature_selection import CorrelationThreshold
from sklearn.datasets import load_boston
from sklearn.datasets import load_iris
from sklearn.datasets import load_diabetes
from sklearn.datasets import load_digits
from sklearn.datasets import load_wine
from sklearn.datasets import load_breast_cancer
from sklearn.utils.testing import assert_allclose_dense_sparse


def test_sparse_matrix_non_pearson():
    X = np.array([[0, 0.2], [1, 1.8], [2, 2.2], [3, 2.8]])
    X = csr_matrix(X)

    msg = "only pearson correlation is supported with sparse matrices"
    with pytest.raises(ValueError, match=msg):
        CorrelationThreshold(kind='spearmanr').fit(X)


@pytest.mark.parametrize("threshold", [-1, 2])
def test_threshold_out_of_bounds(threshold):
    msg = r"threshold must be in \[0.0, 1.0\], got {}".format(threshold)
    with pytest.raises(ValueError, match=msg):
        CorrelationThreshold(threshold=threshold).fit([[0, 1]])


def test_incorrect_kind_parameter():
    X = np.array([[0, 0.2], [1, 1.8], [2, 2.2], [3, 2.8]])
    msg = "kind must be 'pearson' or 'spearmanr'"
    with pytest.raises(ValueError, match=msg):
        CorrelationThreshold(kind='bad').fit(X)


@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
def test_correlated_features_are_removed(toarray):
    rng = np.random.RandomState(0)
    X_shape = (1000, 3)

    X_uncorr = rng.normal(size=X_shape)
    X = np.c_[X_uncorr,
              X_uncorr + 2 * rng.normal(scale=0.05, size=X_shape),
              X_uncorr + 3 * rng.normal(scale=0.05, size=X_shape)]
    X = toarray(X)

    cor_thres = CorrelationThreshold()
    X_tran = cor_thres.fit_transform(X)

    assert X_tran.shape[1] == X_shape[1]


@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
@pytest.mark.parametrize("flip_cols", [True, False])
def test_keeps_feature_with_the_most_variance(toarray, flip_cols):
    X = np.array([[0, 0.2], [1, 1.8], [2, 2.2], [3, 2.8]])
    if flip_cols:
        X = np.c_[X[:, 1], X[:, 0]]

    higher_variance_idx = np.argmax(np.std(X, axis=0))
    X = toarray(X)

    cor_thres = CorrelationThreshold(threshold=0.5)
    X_tran = cor_thres.fit_transform(X)

    assert_allclose_dense_sparse(X[:, [higher_variance_idx]], X_tran)


@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
def test_uncorrelated_features_are_kept(toarray):
    rng = np.random.RandomState(0)
    X_uncorr = toarray(rng.normal(size=(1000, 3)))

    cor_thres = CorrelationThreshold()
    X_tran = cor_thres.fit_transform(X_uncorr)

    assert_allclose_dense_sparse(X_tran, X_uncorr)


@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
def test_all_correlated_features_are_removed(toarray):
    X = np.linspace(0, 10, 100)
    X = toarray(np.c_[X, 2 * X, 0.5 * X])

    cor_thres = CorrelationThreshold()
    X_tran = cor_thres.fit_transform(X)

    assert X_tran.shape[1] == 1


@pytest.mark.parametrize("load_dataset", [
    load_boston, load_iris, load_diabetes, load_digits, load_wine,
    load_breast_cancer
])
@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
def test_increasing_threshold_removes_features_consistently(load_dataset,
                                                            toarray):
    X, _ = load_dataset(return_X_y=True)
    X = toarray(X)

    cor_thresholds = []
    for threshold in np.linspace(0, 1, 20):
        cor_thres = CorrelationThreshold(threshold=threshold)
        X_tran = cor_thres.fit_transform(X)
        assert_allclose_dense_sparse(X_tran, X[:, cor_thres.support_mask_])

        cor_thresholds.append(cor_thres)

    # lower threshold produces a mask that is a subset of a higher threshold
    for lower_cor_thres, higher_cor_thres in zip(cor_thresholds,
                                                 cor_thresholds[1:]):
        assert np.all(lower_cor_thres.support_mask_ <=
                      higher_cor_thres.support_mask_)


@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
def test_threshold_one_keeps_all_features(toarray):
    X = np.linspace(0, 10, 100)
    X = toarray(np.c_[X, 2 * X, 0.5 * X])

    cor_thres = CorrelationThreshold(threshold=1.0)
    X_tran = cor_thres.fit_transform(X)

    assert_allclose_dense_sparse(X_tran, X)


@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
def test_threshold_zero_keeps_one_feature(toarray):
    X = np.linspace(0, 10, 100)
    X = toarray(np.c_[X, 2 * X, 0.5 * X])

    cor_thres = CorrelationThreshold(threshold=0.0)
    X_tran = cor_thres.fit_transform(X)

    # column 1 has the most variance
    assert_allclose_dense_sparse(X_tran, X[:, [1]])


@pytest.mark.parametrize("toarray", [np.asarray, csr_matrix])
def test_constant_features(toarray):
    X = np.array([[0, 0, 0], [1, 1, 1], [1, 2, 3], [2, 3, 4.1]]).T
    X = toarray(X)

    cor_thres = CorrelationThreshold()
    X_tran = cor_thres.fit_transform(X)
    assert X_tran.shape == (3, 1)

    # Feature with the most variance is kept
    assert_allclose_dense_sparse(X_tran, X[:, [3]])


def test_spearmanr_invariant_only_scaling():
    rng = np.random.RandomState(0)
    X_shape = (1000, 3)

    X_uncorr = rng.normal(size=X_shape)
    X = np.c_[X_uncorr,
              X_uncorr + 2 * rng.normal(scale=0.05, size=X_shape),
              X_uncorr + 3 * rng.normal(scale=0.05, size=X_shape)]

    scales = [0.5, 1, 2, 10]

    cor_thresholds = []
    for scale in scales:
        cor_thres = CorrelationThreshold(kind='spearmanr')
        cor_thresholds.append(cor_thres.fit(scale * X))

    for cor_thres_l, cor_thres_r in zip(cor_thresholds[1:],
                                        cor_thresholds[:-1]):
        assert np.all(cor_thres_l.support_mask_ == cor_thres_r.support_mask_)


@pytest.mark.parametrize("flip_cols", [True, False])
def test_spearmanr_keeps_feature_with_the_most_variance(flip_cols):
    X = np.array([[0, 0.2], [1, 1.8], [2, 2.2], [3, 2.8]])
    if flip_cols:
        X = np.c_[X[:, 1], X[:, 0]]

    higher_variance_idx = np.argmax(np.std(X, axis=0))

    cor_thres = CorrelationThreshold(threshold=0.5, kind='spearmanr')
    X_tran = cor_thres.fit_transform(X)

    assert_allclose(X[:, [higher_variance_idx]], X_tran)


def test_spearmanr_uncorrelated_features_are_kept():
    rng = np.random.RandomState(0)
    X_uncorr = rng.normal(size=(1000, 3))

    cor_thres = CorrelationThreshold(kind='spearmanr')
    X_tran = cor_thres.fit_transform(X_uncorr)

    assert_allclose(X_tran, X_uncorr)


@pytest.mark.parametrize("load_dataset", [
    load_boston, load_iris, load_diabetes, load_digits, load_wine,
    load_breast_cancer
])
def test_spearmanrincreasing_threshold_removes_features_consistently(
        load_dataset):
    X, _ = load_dataset(return_X_y=True)

    cor_thresholds = []
    for threshold in np.linspace(0, 1, 20):
        cor_thres = CorrelationThreshold(threshold=threshold, kind='spearmanr')
        X_tran = cor_thres.fit_transform(X)
        assert_allclose_dense_sparse(X_tran, X[:, cor_thres.support_mask_])

        cor_thresholds.append(cor_thres)

    # lower threshold produces a mask that is a subset of a higher threshold
    for lower_cor_thres, higher_cor_thres in zip(cor_thresholds,
                                                 cor_thresholds[1:]):
        assert np.all(lower_cor_thres.support_mask_ <=
                      higher_cor_thres.support_mask_)


def test_spearmanr_threshold_one_keeps_all_features():
    X = np.linspace(0, 10, 100)
    X = np.c_[X, 2 * X, 0.5 * X]

    cor_thres = CorrelationThreshold(threshold=1.0, kind='spearmanr')
    X_tran = cor_thres.fit_transform(X)

    assert_allclose_dense_sparse(X_tran, X)


def test_spearmanr_threshold_zero_keeps_one_feature():
    X = np.linspace(0, 10, 100)
    X = np.c_[X, 2 * X, 0.5 * X]

    cor_thres = CorrelationThreshold(threshold=0.0, kind='spearmanr')
    X_tran = cor_thres.fit_transform(X)

    # column 1 has the most variance
    assert_allclose_dense_sparse(X_tran, X[:, [1]])
