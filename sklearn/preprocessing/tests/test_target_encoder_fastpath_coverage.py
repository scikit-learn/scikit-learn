import numpy as np
import pandas as pd
import numpy.testing as npt
import pytest

from sklearn.preprocessing import TargetEncoder


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def _fit_pair(X_fit, y_fit):
    """Return (fast, slow) TargetEncoder pair fit on the same data.
    fast: small-batch fast path enabled
    slow: vectorized slow path forced via threshold = -1
    """
    te_fast = TargetEncoder(smooth=5.0).fit(X_fit, y_fit)
    te_fast._small_batch_threshold = 256

    te_slow = TargetEncoder(smooth=5.0).fit(X_fit, y_fit)
    te_slow._small_batch_threshold = -1
    return te_fast, te_slow


def _fit_te_binary_pair():
    X_fit = pd.DataFrame({"a": list("uvuw"), "b": list("xyxx")})
    y = np.array([0, 1, 1, 0])
    return _fit_pair(X_fit, y)


def _fit_te_multiclass_pair(n_classes=5, n_samples=64, seed=0):
    rng = np.random.RandomState(seed)
    X_fit = pd.DataFrame({
        "a": rng.choice(list("abcdef"), size=n_samples),
        "b": rng.choice(list("wxyz"), size=n_samples),
    })
    y = rng.randint(0, n_classes, size=n_samples)
    return _fit_pair(X_fit, y)


# ---------------------------------------------------------------------
# Core coverage tests for fast-path branches (parity vs slow path)
# ---------------------------------------------------------------------

def test_norm_keys_cover_nan_pd_na_nat_numpy_nat_and_except_paths():
    te, te_ref = _fit_te_binary_pair()

    class Weird:
        # Make np.isnan/np.isnat raise, exercising the 'except: pass' bodies
        def __array__(self, *args, **kwargs):
            raise TypeError("cannot array-coerce")

    X = pd.DataFrame({
        "a": ["u", float("nan"), pd.NA, "u", Weird()],
        "b": ["x", "x", pd.NaT, np.datetime64("NaT"), "x"],
    }, dtype="object")

    out_fast = te.transform(X)
    out_slow = te_ref.transform(X)
    npt.assert_allclose(out_fast, out_slow, rtol=0, atol=0)


def test_dataframe_shape_mismatch_and_exception_branch():
    te, te_ref = _fit_te_binary_pair()

    class BadShape(pd.DataFrame):
        @property
        def shape(self):  # raises so the fast-path DF branch excepts
            raise RuntimeError("bad shape")

    bad = BadShape(pd.DataFrame({"a": ["u", "v"], "b": ["x", "y"]}))
    X_obj = bad.astype("object")

    out_fast = te.transform(X_obj)     # triggers internal except path
    out_slow = te_ref.transform(X_obj)
    npt.assert_allclose(out_fast, out_slow, rtol=0, atol=0)


def test_transform_raises_attributeerror_when_not_fitted():
    te = TargetEncoder(smooth=5.0)
    with pytest.raises(AttributeError):
        te.transform(pd.DataFrame({"a": ["u"], "b": ["x"]}))


def test_key_collision_disables_feature_fastpath(monkeypatch):
    te, te_ref = _fit_te_binary_pair()

    def const_key(x):
        return "__COLLIDE__"

    monkeypatch.setattr(type(te), "_norm_key", staticmethod(const_key), raising=False)

    X = pd.DataFrame({"a": ["u", "v", "w"], "b": ["x", "x", "y"]})
    out_fast = te.transform(X)
    out_slow = te_ref.transform(X)
    npt.assert_allclose(out_fast, out_slow, rtol=0, atol=0)


def test_enc_vec_ndim_reshape_only():
    te, te_ref = _fit_te_binary_pair()
    te.encodings_[0] = np.asarray(te.encodings_[0]).reshape(-1, 1)  # force 2D on fast
    X = pd.DataFrame({"a": ["u", "zz_unseen"], "b": ["x", "x"]})
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def test_multiclass_block_and_default_vector_two_branches():
    X = pd.DataFrame({"a": ["unseen1", "unseen2", "a", "b"],
                      "b": ["x", "y", "z", "x"]})

    # (a) scalar default on fast; slow gets vector-filled to match expected behavior
    te_fast_a, te_slow_a = _fit_te_multiclass_pair(n_classes=5)
    n_classes = len(te_fast_a.classes_)
    scalar = float(np.mean(te_fast_a.target_mean_))
    te_fast_a.target_mean_ = np.array(scalar)                      # 0-D ndarray (fast)
    te_slow_a.target_mean_ = np.full((n_classes,), scalar)         # 1-D vector (slow)
    npt.assert_allclose(te_fast_a.transform(X), te_slow_a.transform(X), rtol=0, atol=0)

    # (b) correct-length vector on both → pass-through
    te_fast_b, te_slow_b = _fit_te_multiclass_pair(n_classes=5)
    vec = np.full((len(te_fast_b.classes_),), float(np.mean(te_fast_b.target_mean_)))
    te_fast_b.target_mean_ = vec
    te_slow_b.target_mean_ = vec
    npt.assert_allclose(te_fast_b.transform(X), te_slow_b.transform(X), rtol=0, atol=0)



# ---------------------------------------------------------------------
# Extra edge cases to harden behavior and lift patch coverage
# ---------------------------------------------------------------------

def test_threshold_boundary_routes_fast_vs_slow():
    # Build a slightly larger fit set
    Xf = pd.DataFrame({"a": list("uvuw") * 80, "b": list("xyxx") * 80})[:300]
    y = np.array(([0, 1, 1, 0] * 80)[:300])
    te_fast, te_slow = _fit_pair(Xf, y)

    # exactly at threshold → fast path parity with slow
    te_fast._small_batch_threshold = 256
    X256 = Xf.iloc[:256]
    npt.assert_allclose(te_fast.transform(X256), te_slow.transform(X256), rtol=0, atol=0)

    # just over threshold → should match the slow/vectorized output as well
    X257 = Xf.iloc[:257]
    npt.assert_allclose(te_fast.transform(X257), te_slow.transform(X257), rtol=0, atol=0)


def test_pandas_categorical_ordered_and_unordered():
    X_fit = pd.DataFrame({
        "a": pd.Categorical(list("abca"), categories=list("abcd"), ordered=True),
        "b": pd.Categorical(list("xyxz"), categories=list("wxyz"), ordered=False),
    })
    y = np.array([0, 1, 1, 0])

    te, te_ref = _fit_pair(X_fit, y)

    # include unseen categories ('d' and 'w') and NA
    X = pd.DataFrame({
        "a": pd.Categorical(["a", "d", None], categories=list("abcd"), ordered=True),
        "b": pd.Categorical(["w", "y", pd.NA], categories=list("wxyz"), ordered=False),
    })
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def test_all_missing_and_all_unseen_column():
    te, te_ref = _fit_te_binary_pair()

    X_all_missing = pd.DataFrame({"a": [np.nan, pd.NA, pd.NaT], "b": ["x", "x", "x"]}, dtype="object")
    X_all_unseen = pd.DataFrame({"a": ["zzz", "yyy", "xxx"], "b": ["new", "new", "new"]})

    for X in (X_all_missing, X_all_unseen):
        npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def _fit_te_binary_pair_datetime():
    X_fit = pd.DataFrame({
        "a": list("uvuw"),
        "b": np.array(["2024-01-01", "2024-01-03", "2024-01-05", "2024-01-07"],
                      dtype="datetime64[ns]"),
    })
    y = np.array([0, 1, 1, 0])
    return _fit_pair(X_fit, y)

def test_datetime_nat_variant_column():
    te, te_ref = _fit_te_binary_pair_datetime()
    X = pd.DataFrame({
        "a": ["u", "u", "v"],
        "b": np.array(["2024-01-01", "NaT", "2024-02-01"], dtype="datetime64[ns]"),
    })
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def _fit_te_binary_pair_timedelta():
    X_fit = pd.DataFrame({
        "a": list("uvuw"),
        "b": pd.to_timedelta(["1D", "2D", "3D", "4D"]),  # dtype=timedelta64[ns]
    })
    y = np.array([0, 1, 1, 0])
    return _fit_pair(X_fit, y)

def test_timedelta_nat_variant_column():
    te, te_ref = _fit_te_binary_pair_timedelta()
    X = pd.DataFrame({
        "a": ["u", "v", "u"],
        "b": pd.to_timedelta(["1D", "NaT", "2D"]),  # includes NaT safely
    })
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)



def test_timezone_aware_nat():
    X_fit = pd.DataFrame({"a": ["u", "v", "u", "w"], "b": ["x", "y", "x", "y"]})
    y = np.array([0, 1, 1, 0])

    te, te_ref = _fit_pair(X_fit, y)

    X = pd.DataFrame({
        "a": ["u", pd.Timestamp("NaT", tz="UTC")],   # tz-aware NaT
        "b": ["x", "x"],
    }, dtype="object")
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)


def test_numpy_array_input_path():
    X_fit = pd.DataFrame({"a": list("abca"), "b": list("xyxz")})
    y = np.array([0, 1, 1, 0])

    te, te_ref = _fit_pair(X_fit, y)

    X_np = np.array([["a", "x"], ["nan", "x"], ["new", "y"], ["a", "x"]], dtype=object)
    npt.assert_allclose(te.transform(X_np), te_ref.transform(X_np), rtol=0, atol=0)


def test_output_invariants_binary_and_multiclass():
    # binary invariants
    Xb = pd.DataFrame({"a": list("abca"), "b": list("xyxz")})
    yb = np.array([0, 1, 1, 0])
    te_b, _ = _fit_pair(Xb, yb)
    Z = te_b.transform(pd.DataFrame({"a": ["a", "new"], "b": ["x", "x"]}))
    assert Z.dtype == float and Z.ndim == 2 and Z.shape[0] == 2

    # multiclass layout (#cols == n_features * n_classes)
    Xm = pd.DataFrame({"a": list("abca") * 3, "b": list("xyxz") * 3})
    ym = np.array([0, 1, 2, 3] * 3) % 5
    te_m, _ = _fit_pair(Xm, ym)
    Zm = te_m.transform(pd.DataFrame({"a": ["a"], "b": ["x"]}))
    assert Zm.shape[1] == 2 * len(te_m.classes_)  # f0_c0..f0_ck, f1_c0..f1_ck


def test_lazy_structs_cache_reuse():
    X_fit = pd.DataFrame({"a": list("abca"), "b": list("xyxz")})
    y = np.array([0, 1, 1, 0])
    te, _ = _fit_pair(X_fit, y)

    X = pd.DataFrame({"a": ["a", "new", "a", "a"], "b": ["x", "x", "x", "x"]})
    out1 = te.transform(X)
    out2 = te.transform(X)  # should reuse lazily-built caches
    npt.assert_allclose(out1, out2, rtol=0, atol=0)


def test_zero_row_input_raises():
    Xf = pd.DataFrame({"a": list("abca"), "b": list("xyxz")})
    y = np.array([0, 1, 1, 0])
    te, _ = _fit_pair(Xf, y)

    X0 = pd.DataFrame({"a": pd.Series([], dtype=object),
                       "b": pd.Series([], dtype=object)})
    with pytest.raises(ValueError):
        te.transform(X0)


def test_pandas_int64_extension_with_pd_NA():
    from pandas import Int64Dtype

    Xf = pd.DataFrame({
        "a": pd.Series([1, 2, 1, 3], dtype=Int64Dtype()),
        "b": ["x", "y", "x", "y"],
    })
    y = np.array([0, 1, 1, 0])

    te, te_ref = _fit_pair(Xf, y)

    X = pd.DataFrame({
        "a": pd.Series([1, pd.NA, 999], dtype=Int64Dtype()),
        "b": ["x", "x", "x"],
    })
    npt.assert_allclose(te.transform(X), te_ref.transform(X), rtol=0, atol=0)
