import numpy as np
import pytest
from sklearn.preprocessing import TargetEncoder
from sklearn.utils._testing import assert_allclose

def _force_slow_transform(enc: TargetEncoder, X):
    # Force reference/slow path by disabling the small-batch route and clearing caches.
    old_thresh = enc._small_batch_threshold
    enc._small_batch_threshold = -1  # never take fast path
    # Clear caches (simulate a fresh instance for fairness)
    enc._te_index_maps_ = None
    enc._te_is_multiclass_ = None
    enc._te_enc_vecs_ = None
    enc._te_enc_blocks_ = None
    enc._te_defaults_ = None
    try:
        out = enc.transform(X)
    finally:
        enc._small_batch_threshold = old_thresh
    return out

@pytest.mark.parametrize("n_small", [1, 2, 8, 32])
@pytest.mark.parametrize("n_cats", [10_000, 50_000])
def test_small_batch_binary_parity_and_missing(n_small, n_cats):
    rng = np.random.default_rng(0)
    X_fit = rng.integers(0, n_cats, size=(150_000, 1)).astype(object)
    # sprinkle missing during fit to learn its encoding too
    X_fit[:100, 0] = None
    y_fit = rng.integers(0, 2, size=(150_000,))
    enc = TargetEncoder(random_state=0).fit(X_fit, y_fit)

    X_small = rng.integers(0, n_cats * 2, size=(n_small, 1)).astype(object)
    # add missing/unseen explicitly
    if n_small >= 2:
        X_small[0, 0] = None
        X_small[-1, 0] = np.nan

    ref = _force_slow_transform(enc, X_small)
    fast = enc.transform(X_small)
    assert_allclose(fast, ref, rtol=0, atol=0)

def test_small_batch_multiclass_parity_and_order():
    rng = np.random.default_rng(1)
    n_cats = 30_000
    n_classes = 5
    X_fit = rng.integers(0, n_cats, size=(200_000, 1)).astype(object)
    y_fit = rng.integers(0, n_classes, size=(200_000,))
    enc = TargetEncoder(random_state=0).fit(X_fit, y_fit)

    X_small = np.array([[0], [n_cats + 123], [42], [None], [np.nan]], dtype=object)
    ref = _force_slow_transform(enc, X_small)
    fast = enc.transform(X_small)

    # 1) numerical parity
    assert_allclose(fast, ref, rtol=0, atol=0)
    # 2) shape and class-order parity
    assert fast.shape[1] % n_classes == 0
    assert ref.shape == fast.shape
