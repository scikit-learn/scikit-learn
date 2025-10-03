import numpy as np
import numpy.testing as npt
import pytest

from sklearn.preprocessing import TargetEncoder


def _fit_pair_np(X_fit, y_fit):
    """Return (fast, slow) TargetEncoder pair fit on the same data.
    fast: small-batch fast path enabled
    slow: vectorized path forced via threshold = -1
    """
    te_fast = TargetEncoder(smooth=5.0).fit(X_fit, y_fit)
    te_fast._small_batch_threshold = 256

    te_slow = TargetEncoder(smooth=5.0).fit(X_fit, y_fit)
    te_slow._small_batch_threshold = -1
    return te_fast, te_slow


def _fit_te_binary_pair_np():
    X_fit = np.array(
        [["u", "x"], ["v", "y"], ["u", "x"], ["w", "x"]],
        dtype=object,
    )
    y = np.array([0, 1, 1, 0])
    return _fit_pair_np(X_fit, y)


def _fit_te_multiclass_pair_np(n_classes=5, n_samples=64, seed=0):
    rng = np.random.RandomState(seed)
    a = rng.choice(list("abcdef"), size=n_samples)
    b = rng.choice(list("wxyz"), size=n_samples)
    X_fit = np.column_stack([a, b]).astype(object)
    y = rng.randint(0, n_classes, size=n_samples)
    return _fit_pair_np(X_fit, y)


def test_norm_keys_cover_nan_nat_and_except_paths_nopandas():
    te, te_ref = _fit_te_binary_pair_np()

    class Weird:
        # Trigger failures in np.array coercion paths if ever attempted
        def __array__(self, *_, **__):
            raise TypeError("cannot array-coerce")

    # Build object array cell-by-cell to avoid any coercion issues
    X = np.empty((5, 2), dtype=object)
    X[0] = ["u", "x"]
    X[1] = [float("nan"), "x"]  # np.nan → _NAN_SENTINEL
    X[2] = [None, "x"]  # None → _NONE_SENTINEL
    X[3] = ["u", np.datetime64("NaT")]  # datetime NaT → _NAT_SENTINEL
    # assign elementwise so NumPy doesn't try to array-coerce Weird()
    X[4, 0] = Weird()
    X[4, 1] = "x"

    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def test_transform_raises_attributeerror_when_not_fitted_nopandas():
    te = TargetEncoder(smooth=5.0)
    X = np.array([["u", "x"]], dtype=object)
    with pytest.raises(AttributeError):
        te.transform(X)


def test_key_collision_disables_feature_fastpath_nopandas(monkeypatch):
    te, te_ref = _fit_te_binary_pair_np()

    def const_key(_):
        return "__COLLIDE__"

    # Force collisions in key normalization; assert behavioral parity only
    monkeypatch.setattr(type(te), "_norm_key", staticmethod(const_key), raising=False)

    X = np.array([["u", "x"], ["v", "x"], ["w", "y"]], dtype=object)
    out_fast = te.transform(X)
    out_slow = te_ref.transform(X)
    npt.assert_allclose(out_fast, out_slow, rtol=0, atol=0)


def test_enc_vec_ndim_reshape_only_nopandas():
    te, te_ref = _fit_te_binary_pair_np()
    # Force 2D encoding vector for feature 0 → fast path should ravel/reshape
    te.encodings_[0] = np.asarray(te.encodings_[0]).reshape(-1, 1)

    X = np.array([["u", "x"], ["zz_unseen", "x"]], dtype=object)
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def test_multiclass_default_scalar_and_vector_nopandas():
    X = np.array(
        [["unseen1", "x"], ["unseen2", "y"], ["a", "z"], ["b", "x"]], dtype=object
    )

    # (a) scalar default on fast; vector filled on slow → parity
    te_fast_a, te_slow_a = _fit_te_multiclass_pair_np(n_classes=5)
    n_classes = len(te_fast_a.classes_)
    scalar = float(np.mean(te_fast_a.target_mean_))
    te_fast_a.target_mean_ = np.array(scalar)  # 0-D ndarray (fast)
    te_slow_a.target_mean_ = np.full((n_classes,), scalar)  # 1-D vector (slow)
    npt.assert_allclose(te_fast_a.transform(X), te_slow_a.transform(X), rtol=0, atol=0)

    # (b) correct-length vector on both → pass-through parity
    te_fast_b, te_slow_b = _fit_te_multiclass_pair_np(n_classes=5)
    vec = np.full((len(te_fast_b.classes_),), float(np.mean(te_fast_b.target_mean_)))
    te_fast_b.target_mean_ = vec
    te_slow_b.target_mean_ = vec
    npt.assert_allclose(te_fast_b.transform(X), te_slow_b.transform(X), rtol=0, atol=0)


def test_threshold_boundary_routes_fast_vs_slow_nopandas():
    # Larger fit set to make threshold check meaningful
    base = np.array([["u", "x"], ["v", "y"], ["u", "x"], ["w", "x"]], dtype=object)
    Xf = np.vstack([base for _ in range(80)])[:300]
    y = np.array(([0, 1, 1, 0] * 80)[:300])

    te_fast, te_slow = _fit_pair_np(Xf, y)
    te_fast._small_batch_threshold = 256

    # boundary: 256 rows (fast path taken)
    X256 = Xf[:256]
    npt.assert_allclose(
        te_fast.transform(X256), te_slow.transform(X256), rtol=0, atol=0
    )

    # just over: 257 rows (slow/vectorized)
    X257 = Xf[:257]
    npt.assert_allclose(
        te_fast.transform(X257), te_slow.transform(X257), rtol=0, atol=0
    )


def test_all_missing_and_all_unseen_column_nopandas():
    te, te_ref = _fit_te_binary_pair_np()

    X_all_missing = np.array(
        [[np.nan, "x"], [None, "x"], [np.datetime64("NaT"), "x"]],
        dtype=object,
    )
    X_all_unseen = np.array(
        [["zzz", "new"], ["yyy", "new"], ["xxx", "new"]],
        dtype=object,
    )

    for X in (X_all_missing, X_all_unseen):
        out_f = te.transform(X)
        out_s = te_ref.transform(X)
        npt.assert_allclose(out_f, out_s, rtol=0, atol=0)
        assert np.isfinite(out_f).all()


def test_datetime_nat_variant_column_nopandas():
    # Fit with datetime dtype in column 1 so transform dtype matches
    X_fit = np.array(
        [
            ["u", np.datetime64("2024-01-01")],
            ["v", np.datetime64("2024-01-03")],
            ["u", np.datetime64("2024-01-05")],
            ["w", np.datetime64("2024-01-07")],
        ],
        dtype=object,
    )
    y = np.array([0, 1, 1, 0])
    te, te_ref = _fit_pair_np(X_fit, y)

    X = np.array(
        [
            ["u", np.datetime64("2024-01-01")],
            ["u", np.datetime64("NaT")],
            ["v", np.datetime64("2024-02-01")],
        ],
        dtype=object,
    )
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def test_timedelta_nat_variant_column_nopandas():
    # Fit with timedelta dtype in column 1 so transform dtype matches
    X_fit = np.array(
        [
            ["u", np.timedelta64(1, "D")],
            ["v", np.timedelta64(2, "D")],
            ["u", np.timedelta64(3, "D")],
            ["w", np.timedelta64(4, "D")],
        ],
        dtype=object,
    )
    y = np.array([0, 1, 1, 0])
    te, te_ref = _fit_pair_np(X_fit, y)

    X = np.array(
        [
            ["u", np.timedelta64(1, "D")],
            ["v", np.timedelta64("NaT")],
            ["u", np.timedelta64(2, "D")],
        ],
        dtype=object,
    )
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def test_numpy_array_input_path_nopandas():
    te, te_ref = _fit_te_binary_pair_np()
    X_np = np.array([["a", "x"], ["nan", "x"], ["new", "y"], ["a", "x"]], dtype=object)
    npt.assert_allclose(te.transform(X_np), te_ref.transform(X_np), rtol=0, atol=0)


def test_output_invariants_binary_and_multiclass_nopandas():
    # binary invariants
    Xb = np.array([["a", "x"], ["b", "y"], ["c", "x"], ["a", "z"]], dtype=object)
    yb = np.array([0, 1, 1, 0])
    te_b, _ = _fit_pair_np(Xb, yb)
    Z = te_b.transform(np.array([["a", "x"], ["new", "x"]], dtype=object))
    assert Z.dtype == float and Z.ndim == 2 and Z.shape[0] == 2

    # multiclass layout: n_features * n_classes columns
    Xm = np.array([["a", "x"], ["b", "y"], ["c", "x"], ["a", "z"]] * 3, dtype=object)
    ym = np.array([0, 1, 2, 3] * 3) % 5
    te_m, _ = _fit_pair_np(Xm, ym)
    Zm = te_m.transform(np.array([["a", "x"]], dtype=object))
    assert Zm.shape[1] == 2 * len(te_m.classes_)  # 2 features * n_classes


def test_lazy_structs_cache_reuse_nopandas():
    X_fit = np.array([["a", "x"], ["b", "y"], ["c", "x"], ["a", "z"]], dtype=object)
    y = np.array([0, 1, 1, 0])
    te, _ = _fit_pair_np(X_fit, y)

    X = np.array([["a", "x"], ["new", "x"], ["a", "x"], ["a", "x"]], dtype=object)
    out1 = te.transform(X)
    out2 = te.transform(X)  # should reuse lazily-built caches
    npt.assert_allclose(out1, out2, rtol=0, atol=0)


def test_zero_row_input_raises_nopandas():
    Xf = np.array([["a", "x"], ["b", "y"], ["c", "x"], ["a", "z"]], dtype=object)
    y = np.array([0, 1, 1, 0])
    te, _ = _fit_pair_np(Xf, y)

    X0 = np.empty((0, 2), dtype=object)
    with pytest.raises(ValueError):
        te.transform(X0)


# ------------------------------------------------------------------------------
# Explicit default-value behavior checks (exercise unseen fallbacks numerically)
# ------------------------------------------------------------------------------


def test_unseen_uses_binary_default_mean_nopandas():
    te, _ = _fit_te_binary_pair_np()
    # Derive expected default for unseen in binary case:
    # encodings_[j] holds per-category target means; unseen should use te.target_mean_
    expected_default = float(te.target_mean_)
    X = np.array([["__totally_new__", "x"]], dtype=object)
    Z = te.transform(X)
    # First feature unseen → first column equals expected_default
    assert Z.shape == (1, 2)
    assert np.isclose(Z[0, 0], expected_default)


def test_unseen_uses_multiclass_default_vector_nopandas():
    te, _ = _fit_te_multiclass_pair_np(n_classes=4, n_samples=80, seed=42)
    # For multiclass, unseen should use target_mean_ vector (length n_classes)
    X = np.array([["__unseen__", "x"]], dtype=object)
    Z = te.transform(X)
    # Output columns: [f0_c0..f0_cK, f1_c0..f1_cK]
    n_classes = len(te.classes_)
    assert Z.shape[1] == 2 * n_classes
    # First feature's block should equal target_mean_ (within machine epsilon)
    npt.assert_allclose(
        Z[0, :n_classes], np.asarray(te.target_mean_, dtype=float), rtol=0, atol=0
    )
