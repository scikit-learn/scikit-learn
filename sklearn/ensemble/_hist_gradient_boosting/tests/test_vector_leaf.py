"""Tests for trees with vector leaves (vector_leaf=True)."""

import itertools
import pickle

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from sklearn.datasets import (
    make_classification,
    make_multilabel_classification,
    make_regression,
)
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    HistGradientBoostingRegressor,
)
from sklearn.ensemble._hist_gradient_boosting import gradient_boosting as gb
from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper
from sklearn.ensemble._hist_gradient_boosting.common import G_H_DTYPE, X_BINNED_DTYPE
from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower
from sklearn.ensemble._hist_gradient_boosting.histogram import (
    HistogramBuilder,
    VectorHistogramBuilder,
)
from sklearn.metrics import log_loss
from sklearn.utils._openmp_helpers import _openmp_effective_n_threads

n_threads = _openmp_effective_n_threads()

PREDICTOR_STRUCTURE_FIELDS = [
    "feature_idx",
    "bin_threshold",
    "missing_go_to_left",
    "is_leaf",
    "left",
    "right",
    "count",
    "is_categorical",
    "bitset_idx",
]


def _make_binned_data(n_samples=2000, n_features=5, seed=0, with_missing=True):
    rng = np.random.RandomState(seed)
    X = rng.randn(n_samples, n_features)
    if with_missing:
        X[rng.rand(n_samples, n_features) < 0.05] = np.nan
    bin_mapper = _BinMapper(n_bins=64)
    X_binned = bin_mapper.fit_transform(X)
    return X, X_binned, bin_mapper


def _make_gradients(n_samples, n_outputs, constant_hessian, seed=0):
    rng = np.random.RandomState(seed)
    gradients = rng.randn(n_samples, n_outputs).astype(G_H_DTYPE)
    if constant_hessian:
        hessians = np.ones((1, n_outputs), dtype=G_H_DTYPE)
    else:
        hessians = (rng.rand(n_samples, n_outputs) + 0.1).astype(G_H_DTYPE)
    return np.asfortranarray(gradients), np.asfortranarray(hessians)


@pytest.mark.parametrize("constant_hessian", [True, False])
def test_vector_histograms_match_scalar_histograms(constant_hessian):
    n_samples, n_outputs, n_bins = 1000, 3, 32
    rng = np.random.RandomState(0)
    X_binned = np.asfortranarray(
        rng.randint(0, n_bins, size=(n_samples, 4), dtype=X_BINNED_DTYPE)
    )
    gradients, hessians = _make_gradients(n_samples, n_outputs, constant_hessian)
    sample_indices = np.sort(rng.choice(n_samples, 300, replace=False)).astype(
        np.uint32
    )
    sibling_indices = np.setdiff1d(np.arange(n_samples), sample_indices).astype(
        np.uint32
    )

    vector_builder = VectorHistogramBuilder(
        X_binned,
        n_bins,
        np.ascontiguousarray(gradients),
        np.ascontiguousarray(hessians),
        constant_hessian,
        n_threads,
    )
    vector_root = vector_builder.compute_histograms_brute(
        np.arange(n_samples, dtype=np.uint32)
    ).copy()
    vector_node = vector_builder.compute_histograms_brute(sample_indices)
    vector_sibling = vector_builder.compute_histograms_subtraction(
        vector_root.copy(), vector_node
    )
    n_stats = n_outputs + 1 if constant_hessian else 2 * n_outputs + 1
    assert vector_node.shape == (4, n_bins, n_stats)

    for k in range(n_outputs):
        hessians_k = hessians[:, k] if not constant_hessian else hessians[:1, k]
        scalar_builder = HistogramBuilder(
            X_binned,
            n_bins,
            np.ascontiguousarray(gradients[:, k]),
            np.ascontiguousarray(hessians_k),
            constant_hessian,
            n_threads,
        )
        for indices, vector_hist in [
            (sample_indices, vector_node),
            (sibling_indices, vector_sibling),
        ]:
            scalar_hist = np.asarray(scalar_builder.compute_histograms_brute(indices))
            assert_allclose(vector_hist[:, :, k], scalar_hist["sum_gradients"])
            assert_array_equal(vector_hist[:, :, -1], scalar_hist["count"])
            if not constant_hessian:
                assert_allclose(
                    vector_hist[:, :, n_outputs + k], scalar_hist["sum_hessians"]
                )


@pytest.mark.parametrize("constant_hessian", [True, False])
@pytest.mark.parametrize(
    "grower_params",
    [
        dict(),
        dict(l2_regularization=1.0, min_samples_leaf=50),
        dict(interaction_cst=[{0, 1}, {2, 3, 4}]),
        dict(feature_fraction_per_split=0.5),
        dict(categorical=True),
    ],
)
def test_single_output_vector_tree_matches_scalar_tree(constant_hessian, grower_params):
    """A vector tree with a single output is exactly the scalar tree."""
    grower_params = dict(grower_params)
    X, X_binned, bin_mapper = _make_binned_data()
    n_samples = X.shape[0]
    is_categorical = None
    if grower_params.pop("categorical", False):
        X_binned = X_binned.copy(order="F")
        X_binned[:, 0] = X_binned[:, 0] % 20
        is_categorical = np.array([1, 0, 0, 0, 0], dtype=np.uint8)
    gradients, hessians = _make_gradients(n_samples, 1, constant_hessian)
    common_params = dict(
        n_bins=65,
        n_bins_non_missing=bin_mapper.n_bins_non_missing_,
        has_missing_values=np.isnan(X).any(axis=0),
        is_categorical=is_categorical,
        max_leaf_nodes=20,
        min_samples_leaf=10,
        shrinkage=0.1,
        n_threads=n_threads,
    )
    common_params.update(grower_params)

    scalar_grower = TreeGrower(
        X_binned,
        np.ascontiguousarray(gradients[:, 0]),
        np.ascontiguousarray(hessians[:, 0]),
        rng=np.random.default_rng(0),
        **common_params,
    )
    scalar_grower.grow()
    vector_grower = TreeGrower(
        X_binned, gradients, hessians, rng=np.random.default_rng(0), **common_params
    )
    vector_grower.grow()
    assert vector_grower.n_outputs == 1

    scalar_predictor = scalar_grower.make_predictor(bin_mapper.bin_thresholds_)
    vector_predictor = vector_grower.make_predictor(bin_mapper.bin_thresholds_)
    assert vector_predictor.n_outputs == 1
    assert_array_equal(
        scalar_predictor.nodes[PREDICTOR_STRUCTURE_FIELDS],
        vector_predictor.nodes[PREDICTOR_STRUCTURE_FIELDS],
    )
    assert_array_equal(scalar_predictor.nodes["gain"], vector_predictor.nodes["gain"])
    assert_array_equal(
        scalar_predictor.nodes["value"], vector_predictor.node_values[:, 0]
    )
    assert_array_equal(
        scalar_predictor.binned_left_cat_bitsets,
        vector_predictor.binned_left_cat_bitsets,
    )

    scalar_pred = scalar_predictor.predict_binned(
        X_binned, bin_mapper.missing_values_bin_idx_, n_threads
    )
    vector_pred = vector_predictor.predict_binned(
        X_binned, bin_mapper.missing_values_bin_idx_, n_threads
    )
    assert vector_pred.shape == (n_samples, 1)
    assert_array_equal(scalar_pred, vector_pred[:, 0])


@pytest.mark.parametrize("constant_hessian", [True, False])
def test_vector_tree_with_identical_outputs(constant_hessian):
    """Replicating the gradients of one output gives the scalar tree structure."""
    X, X_binned, bin_mapper = _make_binned_data()
    n_samples, n_outputs = X.shape[0], 3
    gradients, hessians = _make_gradients(n_samples, 1, constant_hessian)
    common_params = dict(
        n_bins=65,
        n_bins_non_missing=bin_mapper.n_bins_non_missing_,
        has_missing_values=np.isnan(X).any(axis=0),
        max_leaf_nodes=20,
        min_samples_leaf=10,
        l2_regularization=0.5,
        shrinkage=0.1,
        n_threads=n_threads,
    )
    scalar_grower = TreeGrower(
        X_binned,
        np.ascontiguousarray(gradients[:, 0]),
        np.ascontiguousarray(hessians[:, 0]),
        **common_params,
    )
    scalar_grower.grow()
    vector_grower = TreeGrower(
        X_binned,
        np.asfortranarray(np.tile(gradients, n_outputs)),
        np.asfortranarray(np.tile(hessians, n_outputs)),
        **common_params,
    )
    vector_grower.grow()

    scalar_predictor = scalar_grower.make_predictor(bin_mapper.bin_thresholds_)
    vector_predictor = vector_grower.make_predictor(bin_mapper.bin_thresholds_)
    assert_array_equal(
        scalar_predictor.nodes[PREDICTOR_STRUCTURE_FIELDS],
        vector_predictor.nodes[PREDICTOR_STRUCTURE_FIELDS],
    )
    is_split = ~scalar_predictor.nodes["is_leaf"].astype(bool)
    assert_allclose(
        scalar_predictor.nodes["gain"][is_split] * n_outputs,
        vector_predictor.nodes["gain"][is_split],
    )
    for k in range(n_outputs):
        assert_allclose(
            scalar_predictor.nodes["value"], vector_predictor.node_values[:, k]
        )

    X_test = np.random.RandomState(1).randn(100, X.shape[1])
    known_cat_bitsets, f_idx_map = bin_mapper.make_known_categories_bitsets()
    scalar_pred = scalar_predictor.predict(X_test, known_cat_bitsets, f_idx_map, 1)
    vector_pred = vector_predictor.predict(X_test, known_cat_bitsets, f_idx_map, 1)
    assert_allclose(np.tile(scalar_pred[:, None], n_outputs), vector_pred)

    # predictions can be accumulated in place
    out = np.ones((X_test.shape[0], n_outputs), order="F")
    vector_predictor.predict(X_test, known_cat_bitsets, f_idx_map, 1, out=out)
    assert_allclose(out, vector_pred + 1)


def test_vector_tree_categorical_split():
    """Categories are grouped by their joint effect on all outputs."""
    rng = np.random.RandomState(0)
    n_samples, n_categories, n_outputs = 3000, 8, 3
    X_binned = np.asfortranarray(
        rng.randint(0, n_categories, size=(n_samples, 1), dtype=X_BINNED_DTYPE)
    )
    # Each category shifts the gradients in a direction along the vector
    # (1, -1, 0): categories {0, 2, 4, 6} and {1, 3, 5, 7} have opposite
    # effects, so the optimal split must separate odd and even categories.
    sign = np.where(X_binned[:, 0] % 2 == 0, 1.0, -1.0)
    gradients = sign[:, None] * np.array([1.0, -1.0, 0.0]) + 0.1 * rng.randn(
        n_samples, n_outputs
    )
    hessians = np.ones((n_samples, n_outputs))
    grower = TreeGrower(
        X_binned,
        np.asfortranarray(gradients, dtype=G_H_DTYPE),
        np.asfortranarray(hessians, dtype=G_H_DTYPE),
        n_bins=n_categories + 1,
        n_bins_non_missing=np.array([n_categories], dtype=np.uint32),
        is_categorical=np.array([1], dtype=np.uint8),
        max_leaf_nodes=2,
        shrinkage=1.0,
    )
    grower.grow()
    split_info = grower.root.split_info
    assert split_info.is_categorical
    left_categories = {
        c for c in range(n_categories) if split_info.left_cat_bitset[0] & (1 << c)
    }
    assert left_categories in ({0, 2, 4, 6}, {1, 3, 5, 7})
    for leaf in grower.finalized_leaves:
        # leaf values are close to -mean(gradients) of their samples
        assert abs(abs(leaf.value[0]) - 1) < 0.05
        assert abs(leaf.value[0] + leaf.value[1]) < 0.05
        assert abs(leaf.value[2]) < 0.05


def test_vector_tree_monotonic_cst_not_supported():
    X, X_binned, _ = _make_binned_data(with_missing=False)
    gradients, hessians = _make_gradients(X.shape[0], 2, False)
    with pytest.raises(ValueError, match="Monotonic constraints are not supported"):
        TreeGrower(X_binned, gradients, hessians, monotonic_cst=[1, 0, 0, 0, 0])


def _multiclass_data(n_samples=2000, n_classes=4, seed=0):
    return make_classification(
        n_samples=n_samples,
        n_features=10,
        n_informative=6,
        n_classes=n_classes,
        random_state=seed,
    )


def _multilabel_data(n_samples=2000, n_labels=4, seed=0):
    return make_multilabel_classification(
        n_samples=n_samples, n_features=10, n_classes=n_labels, random_state=seed
    )


def _regression_data(n_samples=2000, seed=0):
    X, y = make_regression(
        n_samples=n_samples, n_features=5, noise=1.0, random_state=seed
    )
    # heteroscedastic noise
    rng = np.random.RandomState(seed)
    y = y + (1 + np.abs(X[:, 0])) * 20 * rng.standard_normal(n_samples)
    return X, y


def test_multiclass_vector_leaf():
    X, y = _multiclass_data()
    X_train, X_test, y_train, y_test = X[:1500], X[1500:], y[:1500], y[1500:]
    losses = {}
    for vector_leaf in [False, True]:
        clf = HistGradientBoostingClassifier(
            vector_leaf=vector_leaf, max_iter=50, early_stopping=False, random_state=0
        ).fit(X_train, y_train)
        n_trees = 1 if vector_leaf else 4
        assert clf.n_trees_per_iteration_ == n_trees
        assert all(len(predictors) == n_trees for predictors in clf._predictors)
        proba = clf.predict_proba(X_test)
        assert proba.shape == (X_test.shape[0], 4)
        assert_allclose(proba.sum(axis=1), 1)
        assert_array_equal(clf.predict(X_test), clf.classes_[proba.argmax(axis=1)])
        losses[vector_leaf] = log_loss(y_test, proba)
    baseline = log_loss(y_test, np.tile(np.bincount(y_train) / len(y_train), (500, 1)))
    assert losses[True] < 0.5 * baseline
    assert losses[True] == pytest.approx(losses[False], rel=0.2)


@pytest.mark.parametrize("vector_leaf", [False, True])
def test_multilabel(vector_leaf):
    X, Y = _multilabel_data()
    X_train, X_test, Y_train, Y_test = X[:1500], X[1500:], Y[:1500], Y[1500:]
    clf = HistGradientBoostingClassifier(
        vector_leaf=vector_leaf, max_iter=50, early_stopping=False, random_state=0
    ).fit(X_train, Y_train)
    assert_array_equal(clf.classes_, np.arange(4))
    assert clf.n_trees_per_iteration_ == (1 if vector_leaf else 4)
    proba = clf.predict_proba(X_test)
    assert proba.shape == Y_test.shape
    assert np.all((proba >= 0) & (proba <= 1))
    pred = clf.predict(X_test)
    assert pred.shape == Y_test.shape
    assert_array_equal(pred, (proba > 0.5).astype(int))
    assert_allclose(clf.decision_function(X_test), np.log(proba / (1 - proba)))
    for k in range(4):
        baseline = log_loss(Y_test[:, k], np.full(500, Y_train[:, k].mean()))
        assert log_loss(Y_test[:, k], proba[:, k]) < 0.8 * baseline


def test_multilabel_one_tree_per_label_matches_independent_classifiers():
    X, Y = _multilabel_data(n_samples=500)
    params = dict(max_iter=10, early_stopping=False, random_state=0)
    clf = HistGradientBoostingClassifier(vector_leaf=False, **params).fit(X, Y)
    for k in range(Y.shape[1]):
        clf_k = HistGradientBoostingClassifier(**params).fit(X, Y[:, k])
        assert_allclose(clf.decision_function(X)[:, k], clf_k.decision_function(X))


def test_multiple_quantiles_one_tree_per_quantile_matches_single_quantiles():
    X, y = _regression_data(n_samples=500)
    quantiles = [0.1, 0.5, 0.9]
    params = dict(loss="quantile", max_iter=10, early_stopping=False, random_state=0)
    reg = HistGradientBoostingRegressor(
        quantile=quantiles, vector_leaf=False, **params
    ).fit(X, y)
    assert reg.n_trees_per_iteration_ == 3
    pred = reg.predict(X)
    assert pred.shape == (X.shape[0], 3)
    for k, quantile in enumerate(quantiles):
        reg_k = HistGradientBoostingRegressor(quantile=quantile, **params).fit(X, y)
        assert_allclose(pred[:, k], reg_k.predict(X))


@pytest.mark.parametrize("vector_leaf", [False, True])
def test_multiple_quantiles_coverage(vector_leaf):
    X, y = _regression_data(n_samples=4000)
    X_train, X_test, y_train, y_test = X[:3000], X[3000:], y[:3000], y[3000:]
    quantiles = (0.1, 0.5, 0.9)
    reg = HistGradientBoostingRegressor(
        loss="quantile",
        quantile=quantiles,
        vector_leaf=vector_leaf,
        max_iter=100,
        early_stopping=False,
        random_state=0,
    ).fit(X_train, y_train)
    assert reg.n_trees_per_iteration_ == (1 if vector_leaf else 3)
    pred = reg.predict(X_test)
    assert pred.shape == (X_test.shape[0], 3)
    train_coverage = (y_train[:, None] <= reg.predict(X_train)).mean(axis=0)
    assert_allclose(train_coverage, quantiles, atol=0.03)
    test_coverage = (y_test[:, None] <= pred).mean(axis=0)
    assert_allclose(test_coverage, quantiles, atol=0.15)


@pytest.mark.parametrize(
    "Estimator, params, data",
    [
        (HistGradientBoostingClassifier, {}, "binary"),
        (HistGradientBoostingRegressor, {"loss": "quantile", "quantile": 0.3}, "reg"),
        (HistGradientBoostingRegressor, {}, "reg"),
    ],
)
def test_vector_leaf_ignored_for_single_output(Estimator, params, data):
    if data == "binary":
        X, y = _multiclass_data(n_samples=500, n_classes=2)
    else:
        X, y = _regression_data(n_samples=500)
    params.update(max_iter=10, early_stopping=False, random_state=0)
    est = Estimator(vector_leaf=False, **params).fit(X, y)
    est_vector = Estimator(vector_leaf=True, **params).fit(X, y)
    assert est_vector.n_trees_per_iteration_ == 1
    assert est_vector._predictors[0][0].n_outputs is None
    assert_array_equal(est.predict(X), est_vector.predict(X))


def _fit_vector_estimators(**params):
    """Yield fitted (estimator, X, y) with vector leaves, for all tasks."""
    X, y = _multiclass_data(n_samples=600)
    yield HistGradientBoostingClassifier(vector_leaf=True, **params).fit(X, y), X, y
    X, Y = _multilabel_data(n_samples=600)
    yield HistGradientBoostingClassifier(vector_leaf=True, **params).fit(X, Y), X, Y
    X, y = _regression_data(n_samples=600)
    reg = HistGradientBoostingRegressor(
        loss="quantile", quantile=[0.2, 0.8], vector_leaf=True, **params
    )
    yield reg.fit(X, y), X, y


def _raw_and_pred(est, X):
    if isinstance(est, HistGradientBoostingClassifier):
        return est.decision_function(X), est.predict_proba(X)
    return est._raw_predict(X), est.predict(X)


@pytest.mark.parametrize("scoring", ["loss", None])
def test_vector_leaf_early_stopping(scoring):
    params = dict(
        max_iter=200,
        early_stopping=True,
        n_iter_no_change=3,
        validation_fraction=0.2,
        scoring=scoring,
        random_state=0,
    )
    estimators = _fit_vector_estimators(**params)
    if scoring is None:
        # The default regression scorer does not support several quantiles.
        estimators = list(itertools.islice(estimators, 2))
        X, y = _regression_data(n_samples=100)
        reg = HistGradientBoostingRegressor(
            loss="quantile", quantile=[0.2, 0.8], **params
        )
        with pytest.raises(ValueError, match="scoring must be 'loss' or a callable"):
            reg.fit(X, y)
    for est, X, y in estimators:
        assert est.n_iter_ < 200
        assert len(est.validation_score_) == est.n_iter_ + 1
        if scoring == "loss":
            # with raw predictions incrementally updated
            assert est.validation_score_[-1] > est.validation_score_[0]


def test_vector_leaf_staged_predictions():
    for est, X, y in _fit_vector_estimators(max_iter=10, early_stopping=False):
        raw, pred = _raw_and_pred(est, X)
        staged_raw = list(est._staged_raw_predict(X))
        assert len(staged_raw) == 10
        assert_allclose(staged_raw[-1], est._raw_predict(X))
        if isinstance(est, HistGradientBoostingClassifier):
            staged_proba = list(est.staged_predict_proba(X))
            assert_allclose(staged_proba[-1], pred)
            assert_array_equal(list(est.staged_predict(X))[-1], est.predict(X))
            assert_allclose(list(est.staged_decision_function(X))[-1], raw)
        else:
            assert_allclose(list(est.staged_predict(X))[-1], pred)


def test_vector_leaf_warm_start():
    params = dict(early_stopping=False, random_state=0)
    for est, X, y in _fit_vector_estimators(max_iter=10, **params):
        est_ws = type(est)(**{**est.get_params(), "max_iter": 5, "warm_start": True})
        est_ws.fit(X, y)
        est_ws.set_params(max_iter=10).fit(X, y)
        assert est_ws.n_iter_ == 10
        assert_allclose(_raw_and_pred(est_ws, X)[0], _raw_and_pred(est, X)[0])

        with pytest.raises(ValueError, match="vector_leaf cannot be changed"):
            est_ws.set_params(max_iter=15, vector_leaf=False).fit(X, y)


def test_vector_leaf_pickle():
    for est, X, y in _fit_vector_estimators(max_iter=5, early_stopping=False):
        est_unpickled = pickle.loads(pickle.dumps(est))
        for a, b in zip(_raw_and_pred(est, X), _raw_and_pred(est_unpickled, X)):
            assert_array_equal(a, b)


def test_vector_leaf_partial_dependence_recursion():
    # With a single feature, the partial dependence computed with the
    # recursion method is exactly the sum of the trees predictions.
    for est, X, y in _fit_vector_estimators(max_iter=5, early_stopping=False):
        X = X[:, :1]
        est.fit(X, y)
        grid = np.sort(X[:20], axis=0)
        pd = est._compute_partial_dependence_recursion(grid, np.array([0]))
        expected = est._raw_predict(grid) - est._baseline_prediction
        assert pd.shape == expected.T.shape
        assert_allclose(pd, expected.T)


@pytest.mark.parametrize("vector_leaf", [False, True])
def test_vector_leaf_missing_values_and_categorical_features(vector_leaf):
    rng = np.random.RandomState(0)
    n_samples = 2000
    X = rng.randn(n_samples, 3)
    categories = rng.randint(0, 6, size=n_samples)
    X[:, 2] = categories
    # class depends on a categorical feature and on missingness
    y = categories % 3
    X[y == 0, 0] = np.nan
    X[rng.rand(n_samples) < 0.1, 1] = np.nan
    clf = HistGradientBoostingClassifier(
        vector_leaf=vector_leaf,
        categorical_features=[2],
        max_iter=20,
        early_stopping=False,
        random_state=0,
    ).fit(X, y)
    assert clf.score(X, y) > 0.99


@pytest.mark.parametrize("vector_leaf", [False, True])
def test_multiple_outputs_monotonic_cst_not_supported(vector_leaf):
    X, Y = _multilabel_data(n_samples=100)
    clf = HistGradientBoostingClassifier(
        vector_leaf=vector_leaf, monotonic_cst=[1] + [0] * 9
    )
    with pytest.raises(ValueError, match="monotonic constraints are not supported"):
        clf.fit(X, Y)
    X, y = _regression_data(n_samples=100)
    reg = HistGradientBoostingRegressor(
        loss="quantile",
        quantile=[0.1, 0.9],
        vector_leaf=vector_leaf,
        monotonic_cst=[1, 0, 0, 0, 0],
    )
    with pytest.raises(ValueError, match="monotonic constraints are not supported"):
        reg.fit(X, y)


def test_multiclass_multioutput_not_supported():
    X = np.random.RandomState(0).randn(20, 2)
    y = np.array([[0, 1, 2, 0]] * 5).T.reshape(-1, 1).repeat(2, axis=1)
    y[::2, 1] = 3
    with pytest.raises(ValueError, match="got a multiclass-multioutput target"):
        HistGradientBoostingClassifier().fit(X, y)


def test_multilabel_validation_data():
    X, Y = _multilabel_data(n_samples=600)
    clf = HistGradientBoostingClassifier(
        vector_leaf=True, max_iter=5, early_stopping=True
    )
    clf.fit(X[:500], Y[:500], X_val=X[500:], y_val=Y[500:])
    assert len(clf.validation_score_) == 6
    with pytest.raises(ValueError, match="binary indicator matrix with 4 columns"):
        clf.fit(X[:500], Y[:500], X_val=X[500:], y_val=Y[500:, :2])


def test_multilabel_sample_and_class_weight():
    X, Y = _multilabel_data(n_samples=600)
    params = dict(max_iter=5, vector_leaf=True, early_stopping=False, random_state=0)
    clf = HistGradientBoostingClassifier(**params).fit(X, Y)
    clf_sw = HistGradientBoostingClassifier(**params).fit(
        X, Y, sample_weight=np.ones(X.shape[0])
    )
    assert_allclose(clf.decision_function(X), clf_sw.decision_function(X))
    HistGradientBoostingClassifier(class_weight="balanced", **params).fit(X, Y)


@pytest.mark.parametrize("vector_leaf", [False, True])
@pytest.mark.parametrize("quantile", [0.3, [0.05, 0.31, 0.5, 0.9]])
def test_quantile_leaves_values_fast_path(monkeypatch, vector_leaf, quantile):
    # The parallel computation of the unweighted leaf quantiles matches the
    # generic line search based on PinballLoss.fit_intercept_only.
    rng = np.random.RandomState(0)
    X = rng.randn(1000, 3)
    y = X[:, 0] + 0.3 * rng.standard_cauchy(1000)
    params = dict(
        loss="quantile",
        quantile=quantile,
        vector_leaf=vector_leaf,
        max_iter=10,
        min_samples_leaf=1,
        random_state=0,
    )
    pred_fast = HistGradientBoostingRegressor(**params).fit(X, y).predict(X)
    monkeypatch.setattr(gb, "_get_pinball_quantiles", lambda loss: None)
    pred_slow = HistGradientBoostingRegressor(**params).fit(X, y).predict(X)
    assert_allclose(pred_fast, pred_slow, rtol=1e-12, atol=1e-12)
