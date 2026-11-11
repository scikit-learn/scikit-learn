import pickle

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from scipy import sparse

from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper
from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower

SPARSE_CONSTRUCTORS = [
    sparse.csr_matrix,
    sparse.csc_matrix,
    sparse.csr_array,
    sparse.csc_array,
]


@pytest.mark.parametrize("problem", ["format", "canonical", "zero_bin"])
def test_sparse_binned_input_validation(problem):
    X = sparse.csr_matrix(np.array([[1], [2]], dtype=np.uint8))
    zero_bin_idx = np.array([0], dtype=np.uint8)
    if problem == "format":
        X = X.tocsc()
        message = "CSR format"
    elif problem == "canonical":
        X = sparse.csr_matrix(
            (np.array([1, 2], dtype=np.uint8), [0, 0], [0, 2, 2]), shape=(2, 1)
        )
        message = "canonical CSR"
    else:
        zero_bin_idx = None
        message = "zero_bin_idx"

    with pytest.raises(ValueError, match=message):
        TreeGrower(
            X,
            gradients=np.ones(2, dtype=np.float32),
            hessians=np.ones(2, dtype=np.float32),
            zero_bin_idx=zero_bin_idx,
        )

    model = HistGradientBoostingRegressor(
        max_iter=1, min_samples_leaf=1, early_stopping=False
    ).fit(np.array([[0.0], [1.0]]), [0.0, 1.0])
    predictor = model._predictors[0][0]
    with pytest.raises(ValueError, match=message):
        predictor.predict_binned(X, 255, 1, zero_bin_idx=zero_bin_idx)


def test_sparse_prediction_rejected_after_categorical_fit():
    X = np.array([[0.0], [1.0], [2.0], [1.0]])
    model = HistGradientBoostingRegressor(
        categorical_features=[0], max_iter=1, early_stopping=False
    ).fit(X, np.arange(4))
    with pytest.raises(ValueError, match="Sparse input is not supported"):
        model.predict(sparse.csr_matrix(X))


def _decode_sparse_bins(X_binned, zero_bin_idx, missing_values_bin_idx):
    decoded = X_binned.toarray()
    for feature_idx, zero_bin in enumerate(zero_bin_idx):
        encoded = decoded[:, feature_idx].copy()
        is_missing = encoded == missing_values_bin_idx
        decoded[:, feature_idx] = np.where(
            encoded == 0,
            zero_bin,
            np.where(encoded <= zero_bin, encoded - 1, encoded),
        )
        decoded[is_missing, feature_idx] = missing_values_bin_idx
    return decoded.astype(np.uint8)


@pytest.mark.parametrize("sparse_constructor", SPARSE_CONSTRUCTORS)
@pytest.mark.parametrize("n_bins", [3, 8])
def test_sparse_bin_mapping_preserves_zero_and_missing(sparse_constructor, n_bins):
    X = np.array(
        [
            [0, 0, 0, 2, -3],
            [0, 1, 0, 1, -2],
            [0, 2, 0, 0, -1],
            [0, 0, 0, -1, 0],
            [0, -1, 0, 0, 1],
            [np.nan, 0, 0, 3, 2],
        ],
        dtype=np.float64,
    )
    X_sparse = sparse_constructor(X)
    mapper = _BinMapper(n_bins=n_bins, random_state=0, n_threads=1).fit(X_sparse)

    X_binned_sparse = mapper.transform(X_sparse)
    X_binned_dense = mapper.transform(X)
    assert sparse.issparse(X_binned_sparse)
    assert X_binned_sparse.format == "csr"
    assert_array_equal(
        _decode_sparse_bins(
            X_binned_sparse,
            mapper.zero_bin_idx_,
            mapper.missing_values_bin_idx_,
        ),
        X_binned_dense,
    )
    assert mapper.zero_bin_idx_[0] == 0
    if n_bins > 3:
        assert mapper.zero_bin_idx_[1] > 0
        assert mapper.zero_bin_idx_[4] > 0
    assert X_binned_sparse.data.dtype == np.uint8


def test_sparse_bin_mapping_coalesces_duplicates_and_ignores_explicit_zeros():
    X_dense = np.array(
        [[0.0, 0.0], [1.0, -2.0], [-1.0, 0.0], [0.5, 3.0], [np.nan, 0.0]]
    )
    clean = sparse.csr_matrix(X_dense)
    # The first coordinate is an explicit zero and the second has duplicate
    # entries whose sum is the corresponding value in X_dense.
    data = np.array([0.0, 0.25, 0.75, -2.0, -1.0, 0.5, 3.0, np.nan])
    noncanonical = sparse.csr_matrix(
        (
            data,
            np.array([0, 0, 0, 1, 0, 0, 1, 0]),
            np.array([0, 1, 4, 5, 7, 8]),
        ),
        shape=X_dense.shape,
    )
    assert not noncanonical.has_canonical_format

    clean_mapper = _BinMapper(n_bins=8, random_state=0, n_threads=1).fit(clean)
    duplicate_mapper = _BinMapper(n_bins=8, random_state=0, n_threads=1).fit(
        noncanonical
    )
    assert_array_equal(
        clean_mapper.bin_thresholds_[0], duplicate_mapper.bin_thresholds_[0]
    )
    assert_array_equal(
        clean_mapper.bin_thresholds_[1], duplicate_mapper.bin_thresholds_[1]
    )
    assert_array_equal(
        clean_mapper.transform(clean).toarray(),
        duplicate_mapper.transform(noncanonical).toarray(),
    )
    assert not noncanonical.has_canonical_format


@pytest.mark.parametrize("sparse_constructor", SPARSE_CONSTRUCTORS)
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize(
    "Estimator, y_factory",
    [
        (
            HistGradientBoostingClassifier,
            lambda X: (np.nan_to_num(X[:, 0]) - X[:, 1] > 0).astype(int),
        ),
        (
            HistGradientBoostingRegressor,
            lambda X: (
                np.nan_to_num(X[:, 0])
                + 2 * np.nan_to_num(X[:, 1])
                - np.nan_to_num(X[:, 2])
            ),
        ),
    ],
)
def test_sparse_fit_and_sparse_dense_prediction_agree(
    sparse_constructor, dtype, Estimator, y_factory
):
    rng = np.random.default_rng(42)
    X = rng.normal(size=(96, 7)).astype(dtype)
    X[rng.random(X.shape) < 0.75] = 0
    X[3, 4] = np.nan
    X[26, 1] = np.nan
    y = y_factory(X)
    X_sparse = sparse_constructor(X)
    sample_weight = np.linspace(0.5, 1.5, X.shape[0])

    model = Estimator(
        max_iter=4,
        max_leaf_nodes=5,
        min_samples_leaf=3,
        max_bins=16,
        early_stopping=False,
        random_state=0,
    ).fit(X_sparse, y, sample_weight=sample_weight)

    assert_array_equal(model.predict(X_sparse), model.predict(X))
    assert model.__sklearn_tags__().input_tags.sparse
    assert model._bin_mapper.transform(X_sparse).nnz <= X_sparse.nnz
    assert all(
        np.allclose(sparse_stage, dense_stage)
        for sparse_stage, dense_stage in zip(
            model.staged_predict(X_sparse), model.staged_predict(X)
        )
    )
    if isinstance(model, HistGradientBoostingClassifier):
        np.testing.assert_allclose(
            model.predict_proba(X_sparse), model.predict_proba(X)
        )
        np.testing.assert_allclose(model._raw_predict(X_sparse), model._raw_predict(X))

    restored = pickle.loads(pickle.dumps(model))
    assert_array_equal(restored.predict(X_sparse), model.predict(X_sparse))


def test_sparse_multiclass_with_explicit_validation_and_interaction_constraints():
    rng = np.random.default_rng(0)
    X_train = rng.normal(size=(80, 6))
    X_train[rng.random(X_train.shape) < 0.7] = 0
    X_val = rng.normal(size=(24, 6))
    X_val[rng.random(X_val.shape) < 0.7] = 0
    y_train = np.argmax(
        np.column_stack((X_train[:, 0], X_train[:, 2], -X_train[:, 4])), axis=1
    )
    y_val = np.argmax(np.column_stack((X_val[:, 0], X_val[:, 2], -X_val[:, 4])), axis=1)
    X_train_sparse = sparse.csr_array(X_train)
    X_val_sparse = sparse.csc_array(X_val)
    sample_weight = np.linspace(0.5, 1.5, X_train.shape[0])

    model = HistGradientBoostingClassifier(
        max_iter=4,
        max_leaf_nodes=5,
        min_samples_leaf=3,
        max_bins=16,
        early_stopping=True,
        n_iter_no_change=2,
        interaction_cst=[[0, 1, 2], [3, 4, 5]],
        max_features=0.8,
        l2_regularization=0.5,
        random_state=0,
    ).fit(
        X_train_sparse,
        y_train,
        sample_weight=sample_weight,
        X_val=X_val_sparse,
        y_val=y_val,
    )

    np.testing.assert_allclose(
        model.predict_proba(X_val_sparse), model.predict_proba(X_val)
    )
    assert model.classes_.shape == (3,)
    assert model._raw_predict(X_train_sparse).shape == (80, 3)


def test_sparse_implicit_early_stopping_split_with_imbalanced_classes():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(100, 5))
    X[rng.random(X.shape) < 0.8] = 0
    y = np.zeros(100, dtype=int)
    y[:4] = 1
    X[:4, 0] = 3
    sample_weight = np.linspace(0.5, 1.5, X.shape[0])

    model = HistGradientBoostingClassifier(
        max_iter=3,
        max_leaf_nodes=4,
        min_samples_leaf=2,
        early_stopping=True,
        validation_fraction=0.2,
        n_iter_no_change=2,
        random_state=0,
    ).fit(sparse.csr_matrix(X), y, sample_weight=sample_weight)

    assert model.validation_score_.size > 0
    assert_array_equal(model.predict(sparse.csr_matrix(X)), model.predict(X))


def test_sparse_warm_start_and_monotonic_constraint():
    X = np.arange(80, dtype=np.float64).reshape(-1, 1)
    y = (X[:, 0] > 36).astype(int)
    X_sparse = sparse.csr_matrix(X)
    model = HistGradientBoostingClassifier(
        max_iter=2,
        max_leaf_nodes=3,
        min_samples_leaf=2,
        early_stopping=False,
        warm_start=True,
        monotonic_cst=[1],
        random_state=0,
    ).fit(X_sparse, y)
    model.set_params(max_iter=4).fit(X_sparse, y)

    probabilities = model.predict_proba(X_sparse)[:, 1]
    assert np.all(np.diff(probabilities) >= -1e-14)
    assert model.n_iter_ == 4
    np.testing.assert_allclose(model.predict_proba(X_sparse), model.predict_proba(X))


def test_sparse_categorical_features_are_rejected_explicitly():
    X = sparse.csr_matrix([[0.0, 1.0], [1.0, 0.0], [0.0, 2.0]])
    y = np.array([0, 1, 1])
    model = HistGradientBoostingClassifier(
        categorical_features=[0], max_iter=2, early_stopping=False
    )
    with pytest.raises(ValueError, match="Categorical features are not supported"):
        model.fit(X, y)


def test_all_zero_sparse_training_data():
    X = sparse.csr_matrix((12, 1), dtype=np.float64)
    y = np.array([0, 1] * 6)
    model = HistGradientBoostingClassifier(
        max_iter=2,
        max_leaf_nodes=3,
        min_samples_leaf=2,
        early_stopping=False,
        random_state=0,
    ).fit(X, y)

    assert_array_equal(model.predict(X), model.predict(X.toarray()))
    assert model._bin_mapper.bin_thresholds_[0].size == 0


def test_sparse_fit_and_prediction_do_not_call_toarray(monkeypatch):
    X = sparse.csr_matrix([[0.0, 1.0], [1.0, 0.0], [0.0, -1.0], [-1.0, 0.0]] * 4)
    y = np.array([0, 1, 1, 0] * 4)

    def fail_toarray(*args, **kwargs):
        raise AssertionError("Sparse HGB input must not be densified.")

    monkeypatch.setattr(sparse.csr_matrix, "toarray", fail_toarray)
    model = HistGradientBoostingClassifier(
        max_iter=2, min_samples_leaf=1, early_stopping=False, random_state=0
    ).fit(X, y)
    model.predict(X)


def test_sparse_splitter_and_histogram_accept_int64_indices():
    X = sparse.csr_matrix(
        [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, -1.0], [-1.0, 0.0, 0.0]]
    )
    X.indices = X.indices.astype(np.int64)
    X.indptr = X.indptr.astype(np.int64)
    y = np.array([0, 1, 1, 0])
    model = HistGradientBoostingClassifier(
        max_iter=2,
        max_leaf_nodes=3,
        min_samples_leaf=1,
        early_stopping=False,
        random_state=0,
    ).fit(X, y)
    np.testing.assert_allclose(model.predict_proba(X), model.predict_proba(X.toarray()))
