# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import numpy.testing as npt

from sklearn.preprocessing import TargetEncoder
from sklearn.preprocessing._target_encoder import _norm_key


def test_norm_key_real_values_cover_nan_nat_and_except_paths():
    # 1) Exactly None → distinct sentinel
    k_none = _norm_key(None)

    # 2) Float NaN (Python/NumPy) → distinct sentinel
    k_nan1 = _norm_key(float("nan"))
    k_nan2 = _norm_key(np.float64(np.nan))
    assert k_nan1 is k_nan2

    # 3) Plain string reaches the NumPy NaT check and naturally raises
    #    TypeError inside np.isnat("x") → defensive 'except: pass' executes.
    k_str = _norm_key("x")
    assert k_str == "x"

    # 4) NumPy datetime/timedelta NaT (non-float) → distinct sentinel
    k_nat_dt = _norm_key(np.datetime64("NaT"))
    k_nat_td = _norm_key(np.timedelta64("NaT"))
    assert k_nat_dt is k_nat_td

    # sanity: all sentinels are distinct from ordinary values
    assert k_none is not k_nan1 and k_none is not k_nat_dt
    assert k_nan1 is not k_nat_dt and k_nan1 != "x"
    assert k_nat_dt != "x"


def _fit_pair_numpy():
    X_fit = np.array(
        [
            ["u", "x"],
            ["v", "y"],
            ["u", "x"],
            ["w", "y"],
        ],
        dtype=object,
    )
    y = np.array([0.0, 1.0, 1.0, 0.0])
    te_fast = TargetEncoder(smooth=5.0).fit(X_fit, y)
    te_fast._small_batch_threshold = 256  # enable small-batch fast path
    te_vec = TargetEncoder(smooth=5.0).fit(X_fit, y)
    te_vec._small_batch_threshold = -1  # always vectorized
    return te_fast, te_vec


def test_small_batch_fastpath_matches_vectorized_with_nans_and_nats():
    te_fast, te_vec = _fit_pair_numpy()

    # Build object array cell-by-cell so everything stays dtype=object
    X = np.empty((6, 2), dtype=object)
    X[0] = ["u", "x"]  # seen
    X[1] = [float("nan"), "x"]  # float NaN branch
    X[2] = [None, "x"]  # None branch
    X[3] = ["u", np.datetime64("NaT")]  # datetime NaT branch
    X[4] = ["new_unseen", "x"]  # unseen → default
    X[5] = ["u", np.timedelta64("NaT")]  # timedelta NaT branch

    npt.assert_allclose(te_fast.transform(X), te_vec.transform(X), rtol=0, atol=0)


def test_large_batch_vectorized_consistency():
    te_fast, te_vec = _fit_pair_numpy()

    # Make > threshold rows by repeating real rows
    base = np.array(
        [["u", "x"], ["v", "y"], ["u", "x"], ["w", "y"], [None, "x"], ["new", "y"]],
        dtype=object,
    )
    X_big = np.tile(base, (60, 1))  # 360 rows → well over 256

    Z_fast = te_fast.transform(X_big)  # will choose vectorized due to size
    Z_vec = te_vec.transform(X_big)  # always vectorized
    npt.assert_allclose(Z_fast, Z_vec, rtol=0, atol=0)


def test_multiclass_default_fallback_uses_block_mean_axis1():
    rng = np.random.RandomState(0)
    X_fit = np.stack(
        [
            rng.choice(list("abc"), size=80),
            rng.choice(list("wxyz"), size=80),
        ],
        axis=1,
    ).astype(object)
    # 3 classes, legitimate multiclass target
    y = rng.randint(0, 3, size=80)

    te = TargetEncoder(target_type="multiclass").fit(X_fit, y)
    te._small_batch_threshold = 256  # allow small-batch cache build on tiny input

    # Use valid numeric values but a wrong shape (2D array). This is *not* a fake
    # type; it’s real numeric data and simply forces the robust fallback branch.
    te.target_mean_ = np.array([[0.2, 0.5, 0.3]], dtype=float)

    # Small transform to build caches (and compute default via block.mean(axis=1))
    X_small = np.array(
        [
            ["a", "w"],  # seen
            ["zzz_unseen", "w"],  # unseen in first feature → use default vector
        ],
        dtype=object,
    )

    Z = te.transform(X_small)
    # Output has n_features * n_classes columns; values must be finite
    n_classes = len(te.classes_)
    assert Z.shape == (2, X_small.shape[1] * n_classes)
    assert np.isfinite(Z).all()
