import re
import warnings

import numpy as np
import pytest
import scipy.sparse as sp

from sklearn import clone
from sklearn.preprocessing import KBinsDiscretizer, OneHotEncoder
from sklearn.utils._mask import _get_mask
from sklearn.utils._testing import (
    assert_allclose,
    assert_allclose_dense_sparse,
    assert_array_almost_equal,
    assert_array_equal,
    ignore_warnings,
)

X = [[-2, 1.5, -4, -1], [-1, 2.5, -3, -0.5], [0, 3.5, -2, 0.5], [1, 4.5, -1, 2]]


@pytest.mark.parametrize(
    "strategy, quantile_method, expected, sample_weight",
    [
        (
            "uniform",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [1, 1, 1, 0], [2, 2, 2, 1], [2, 2, 2, 2]],
            None,
        ),
        (
            "kmeans",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2]],
            None,
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2], [2, 2, 2, 2]],
            None,
        ),
        (
            "uniform",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [1, 1, 1, 0], [2, 2, 2, 1], [2, 2, 2, 2]],
            [1, 1, 2, 1],
        ),
        (
            "uniform",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [1, 1, 1, 0], [2, 2, 2, 1], [2, 2, 2, 2]],
            [1, 1, 1, 1],
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2], [2, 2, 2, 2]],
            [1, 1, 2, 1],
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2], [2, 2, 2, 2]],
            [1, 1, 1, 1],
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2]],
            [0, 1, 1, 1],
        ),
        (
            "kmeans",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1], [2, 2, 2, 2]],
            [1, 0, 3, 1],
        ),
        (
            "kmeans",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2]],
            [1, 1, 1, 1],
        ),
    ],
)
def test_fit_transform(strategy, quantile_method, expected, sample_weight):
    est = KBinsDiscretizer(
        n_bins=3, encode="ordinal", strategy=strategy, quantile_method=quantile_method
    )
    with ignore_warnings(category=UserWarning):
        # Ignore the warning on removed small bins.
        est.fit(X, sample_weight=sample_weight)
    assert_array_equal(est.transform(X), expected)


def test_valid_n_bins():
    KBinsDiscretizer(n_bins=2, quantile_method="averaged_inverted_cdf").fit_transform(X)
    KBinsDiscretizer(
        n_bins=np.array([2])[0], quantile_method="averaged_inverted_cdf"
    ).fit_transform(X)
    assert KBinsDiscretizer(n_bins=2, quantile_method="averaged_inverted_cdf").fit(
        X
    ).n_bins_.dtype == np.dtype(int)


def test_invalid_n_bins_array():
    # Bad shape
    n_bins = np.full((2, 4), 2.0)
    est = KBinsDiscretizer(n_bins=n_bins, quantile_method="averaged_inverted_cdf")
    err_msg = r"n_bins must be a scalar or array of shape \(n_features,\)."
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)

    # Incorrect number of features
    n_bins = [1, 2, 2]
    est = KBinsDiscretizer(n_bins=n_bins, quantile_method="averaged_inverted_cdf")
    err_msg = r"n_bins must be a scalar or array of shape \(n_features,\)."
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)

    # Bad bin values
    n_bins = [1, 2, 2, 1]
    est = KBinsDiscretizer(n_bins=n_bins, quantile_method="averaged_inverted_cdf")
    err_msg = (
        "KBinsDiscretizer received an invalid number of bins "
        "at indices 0, 3. Number of bins must be at least 2, "
        "and must be an int."
    )
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)

    # Float bin values
    n_bins = [2.1, 2, 2.1, 2]
    est = KBinsDiscretizer(n_bins=n_bins, quantile_method="averaged_inverted_cdf")
    err_msg = (
        "KBinsDiscretizer received an invalid number of bins "
        "at indices 0, 2. Number of bins must be at least 2, "
        "and must be an int."
    )
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)


@pytest.mark.parametrize(
    "strategy, quantile_method, expected, sample_weight",
    [
        (
            "uniform",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [0, 1, 1, 0], [1, 2, 2, 1], [1, 2, 2, 2]],
            None,
        ),
        (
            "kmeans",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [1, 2, 2, 2]],
            None,
        ),
        (
            "quantile",
            "linear",
            [[0, 0, 0, 0], [0, 1, 1, 1], [1, 2, 2, 2], [1, 2, 2, 2]],
            None,
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [0, 1, 1, 1], [1, 2, 2, 2], [1, 2, 2, 2]],
            None,
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [0, 1, 1, 1], [1, 2, 2, 2], [1, 2, 2, 2]],
            [1, 1, 1, 1],
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [1, 1, 1, 1]],
            [0, 1, 3, 1],
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [[0, 0, 0, 0], [0, 0, 0, 0], [1, 2, 2, 2], [1, 2, 2, 2]],
            [1, 1, 3, 1],
        ),
        (
            "kmeans",
            "averaged_inverted_cdf",  # default
            [[0, 0, 0, 0], [0, 1, 1, 0], [1, 1, 1, 1], [1, 2, 2, 2]],
            [1, 0, 3, 1],
        ),
    ],
)
def test_fit_transform_n_bins_array(strategy, quantile_method, expected, sample_weight):
    est = KBinsDiscretizer(
        n_bins=[2, 3, 3, 3],
        encode="ordinal",
        strategy=strategy,
        quantile_method=quantile_method,
    ).fit(X, sample_weight=sample_weight)
    assert_array_equal(est.transform(X), expected)

    # test the shape of bin_edges_
    n_features = np.array(X).shape[1]
    assert est.bin_edges_.shape == (n_features,)
    for bin_edges, n_bins in zip(est.bin_edges_, est.n_bins_):
        assert bin_edges.shape == (n_bins + 1,)


@pytest.mark.filterwarnings("ignore: Bins whose width are too small")
def test_kbinsdiscretizer_effect_sample_weight():
    """Check the impact of `sample_weight` one computed quantiles."""
    X = np.array([[-2], [-1], [1], [3], [500], [1000]])
    # add a large number of bins such that each sample with a non-null weight
    # will be used as bin edge
    est = KBinsDiscretizer(
        n_bins=10,
        encode="ordinal",
        strategy="quantile",
        quantile_method="averaged_inverted_cdf",
    )
    est.fit(X, sample_weight=[1, 1, 1, 1, 0, 0])
    assert_allclose(est.bin_edges_[0], [-2, -1, 0, 1, 3])
    assert_allclose(est.transform(X), [[0.0], [1.0], [3.0], [3.0], [3.0], [3.0]])


@pytest.mark.parametrize("strategy", ["kmeans", "quantile"])
def test_kbinsdiscretizer_no_mutating_sample_weight(strategy):
    """Make sure that `sample_weight` is not changed in place."""

    if strategy == "quantile":
        est = KBinsDiscretizer(
            n_bins=3,
            encode="ordinal",
            strategy=strategy,
            quantile_method="averaged_inverted_cdf",
        )
    else:
        est = KBinsDiscretizer(n_bins=3, encode="ordinal", strategy=strategy)
    sample_weight = np.array([1, 3, 1, 2], dtype=np.float64)
    sample_weight_copy = np.copy(sample_weight)
    est.fit(X, sample_weight=sample_weight)
    assert_allclose(sample_weight, sample_weight_copy)


@pytest.mark.parametrize("strategy", ["uniform", "kmeans", "quantile"])
def test_same_min_max(strategy):
    warnings.simplefilter("always")
    X = np.array([[1, -2], [1, -1], [1, 0], [1, 1]])
    if strategy == "quantile":
        est = KBinsDiscretizer(
            strategy=strategy,
            n_bins=3,
            encode="ordinal",
            quantile_method="averaged_inverted_cdf",
        )
    else:
        est = KBinsDiscretizer(strategy=strategy, n_bins=3, encode="ordinal")
    warning_message = "Feature 0 is constant and will be replaced with 0."
    with pytest.warns(UserWarning, match=warning_message):
        est.fit(X)
    assert est.n_bins_[0] == 1
    # replace the feature with zeros
    Xt = est.transform(X)
    assert_array_equal(Xt[:, 0], np.zeros(X.shape[0]))


def test_transform_1d_behavior():
    X = np.arange(4)
    est = KBinsDiscretizer(n_bins=2, quantile_method="averaged_inverted_cdf")
    with pytest.raises(ValueError):
        est.fit(X)

    est = KBinsDiscretizer(n_bins=2, quantile_method="averaged_inverted_cdf")
    est.fit(X.reshape(-1, 1))
    with pytest.raises(ValueError):
        est.transform(X)


@pytest.mark.parametrize("i", range(1, 9))
def test_numeric_stability(i):
    X_init = np.array([2.0, 4.0, 6.0, 8.0, 10.0]).reshape(-1, 1)
    Xt_expected = np.array([0, 0, 1, 1, 1]).reshape(-1, 1)

    # Test up to discretizing nano units
    X = X_init / 10**i
    Xt = KBinsDiscretizer(
        n_bins=2, encode="ordinal", quantile_method="averaged_inverted_cdf"
    ).fit_transform(X)
    assert_array_equal(Xt_expected, Xt)


def test_encode_options():
    est = KBinsDiscretizer(
        n_bins=[2, 3, 3, 3], encode="ordinal", quantile_method="averaged_inverted_cdf"
    ).fit(X)
    Xt_1 = est.transform(X)
    est = KBinsDiscretizer(
        n_bins=[2, 3, 3, 3],
        encode="onehot-dense",
        quantile_method="averaged_inverted_cdf",
    ).fit(X)
    Xt_2 = est.transform(X)
    assert not sp.issparse(Xt_2)
    assert_array_equal(
        OneHotEncoder(
            categories=[np.arange(i) for i in [2, 3, 3, 3]], sparse_output=False
        ).fit_transform(Xt_1),
        Xt_2,
    )
    est = KBinsDiscretizer(
        n_bins=[2, 3, 3, 3], encode="onehot", quantile_method="averaged_inverted_cdf"
    ).fit(X)
    Xt_3 = est.transform(X)
    assert sp.issparse(Xt_3)
    assert_array_equal(
        OneHotEncoder(
            categories=[np.arange(i) for i in [2, 3, 3, 3]], sparse_output=True
        )
        .fit_transform(Xt_1)
        .toarray(),
        Xt_3.toarray(),
    )


@pytest.mark.parametrize(
    "strategy, quantile_method, expected_2bins, expected_3bins, expected_5bins",
    [
        (
            "uniform",
            "averaged_inverted_cdf",
            [0, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 2, 2],
            [0, 0, 1, 1, 4, 4],
        ),
        (
            "kmeans",
            "averaged_inverted_cdf",
            [0, 0, 0, 0, 1, 1],
            [0, 0, 1, 1, 2, 2],
            [0, 0, 1, 2, 3, 4],
        ),
        (
            "quantile",
            "averaged_inverted_cdf",
            [0, 0, 0, 1, 1, 1],
            [0, 0, 1, 1, 2, 2],
            [0, 1, 2, 3, 4, 4],
        ),
    ],
)
def test_nonuniform_strategies(
    strategy, quantile_method, expected_2bins, expected_3bins, expected_5bins
):
    X = np.array([0, 0.5, 2, 3, 9, 10]).reshape(-1, 1)

    # with 2 bins
    est = KBinsDiscretizer(
        n_bins=2, strategy=strategy, quantile_method=quantile_method, encode="ordinal"
    )
    Xt = est.fit_transform(X)
    assert_array_equal(expected_2bins, Xt.ravel())

    # with 3 bins
    est = KBinsDiscretizer(
        n_bins=3, strategy=strategy, quantile_method=quantile_method, encode="ordinal"
    )
    Xt = est.fit_transform(X)
    assert_array_equal(expected_3bins, Xt.ravel())

    # with 5 bins
    est = KBinsDiscretizer(
        n_bins=5, strategy=strategy, quantile_method=quantile_method, encode="ordinal"
    )
    Xt = est.fit_transform(X)
    assert_array_equal(expected_5bins, Xt.ravel())


@pytest.mark.parametrize(
    "strategy, expected_inv,quantile_method",
    [
        (
            "uniform",
            [
                [-1.5, 2.0, -3.5, -0.5],
                [-0.5, 3.0, -2.5, -0.5],
                [0.5, 4.0, -1.5, 0.5],
                [0.5, 4.0, -1.5, 1.5],
            ],
            "averaged_inverted_cdf",  # default
        ),
        (
            "kmeans",
            [
                [-1.375, 2.125, -3.375, -0.5625],
                [-1.375, 2.125, -3.375, -0.5625],
                [-0.125, 3.375, -2.125, 0.5625],
                [0.75, 4.25, -1.25, 1.625],
            ],
            "averaged_inverted_cdf",  # default
        ),
        (
            "quantile",
            [
                [-1.5, 2.0, -3.5, -0.75],
                [-0.5, 3.0, -2.5, 0.0],
                [0.5, 4.0, -1.5, 1.25],
                [0.5, 4.0, -1.5, 1.25],
            ],
            "averaged_inverted_cdf",
        ),
    ],
)
@pytest.mark.parametrize("encode", ["ordinal", "onehot", "onehot-dense"])
def test_inverse_transform(strategy, encode, expected_inv, quantile_method):
    kbd = KBinsDiscretizer(
        n_bins=3, strategy=strategy, quantile_method=quantile_method, encode=encode
    )
    Xt = kbd.fit_transform(X)
    Xinv = kbd.inverse_transform(Xt)
    assert_array_almost_equal(expected_inv, Xinv)


@pytest.mark.parametrize("strategy", ["uniform", "kmeans", "quantile"])
def test_transform_outside_fit_range(strategy):
    X = np.array([0, 1, 2, 3])[:, None]

    if strategy == "quantile":
        kbd = KBinsDiscretizer(
            n_bins=4,
            strategy=strategy,
            encode="ordinal",
            quantile_method="averaged_inverted_cdf",
        )
    else:
        kbd = KBinsDiscretizer(n_bins=4, strategy=strategy, encode="ordinal")
    kbd.fit(X)

    X2 = np.array([-2, 5])[:, None]
    X2t = kbd.transform(X2)
    assert_array_equal(X2t.max(axis=0) + 1, kbd.n_bins_)
    assert_array_equal(X2t.min(axis=0), [0])


def test_overwrite():
    X = np.array([0, 1, 2, 3])[:, None]
    X_before = X.copy()

    est = KBinsDiscretizer(
        n_bins=3, quantile_method="averaged_inverted_cdf", encode="ordinal"
    )
    Xt = est.fit_transform(X)
    assert_array_equal(X, X_before)

    Xt_before = Xt.copy()
    Xinv = est.inverse_transform(Xt)
    assert_array_equal(Xt, Xt_before)
    assert_array_equal(Xinv, np.array([[0.5], [1.5], [2.5], [2.5]]))


@pytest.mark.parametrize(
    "strategy, expected_bin_edges, quantile_method",
    [
        ("quantile", [0, 1.5, 3], "averaged_inverted_cdf"),
        ("kmeans", [0, 1.5, 3], "averaged_inverted_cdf"),
    ],
)
def test_redundant_bins(strategy, expected_bin_edges, quantile_method):
    X = [[0], [0], [0], [0], [3], [3]]
    kbd = KBinsDiscretizer(
        n_bins=3, strategy=strategy, quantile_method=quantile_method, subsample=None
    )
    warning_message = "Consider decreasing the number of bins."
    with pytest.warns(UserWarning, match=warning_message):
        kbd.fit(X)

    assert_array_almost_equal(kbd.bin_edges_[0], expected_bin_edges)


def test_percentile_numeric_stability():
    X = np.array([0.05, 0.05, 0.95]).reshape(-1, 1)
    bin_edges = np.array([0.05, 0.23, 0.41, 0.59, 0.77, 0.95])
    Xt = np.array([0, 0, 4]).reshape(-1, 1)
    kbd = KBinsDiscretizer(
        n_bins=10,
        encode="ordinal",
        strategy="quantile",
        quantile_method="linear",
    )
    ## TODO: change to averaged inverted cdf, but that means we only get bin
    ## edges of 0.05 and 0.95 and nothing in between

    warning_message = "Consider decreasing the number of bins."
    with pytest.warns(UserWarning, match=warning_message):
        kbd.fit(X)

    assert_array_almost_equal(kbd.bin_edges_[0], bin_edges)
    assert_array_almost_equal(kbd.transform(X), Xt)


@pytest.mark.parametrize("in_dtype", [np.float16, np.float32, np.float64])
@pytest.mark.parametrize("out_dtype", [None, np.float32, np.float64])
@pytest.mark.parametrize("encode", ["ordinal", "onehot", "onehot-dense"])
def test_consistent_dtype(in_dtype, out_dtype, encode):
    X_input = np.array(X, dtype=in_dtype)
    kbd = KBinsDiscretizer(
        n_bins=3,
        encode=encode,
        quantile_method="averaged_inverted_cdf",
        dtype=out_dtype,
    )
    kbd.fit(X_input)

    # test output dtype
    if out_dtype is not None:
        expected_dtype = out_dtype
    elif out_dtype is None and X_input.dtype == np.float16:
        # wrong numeric input dtype are cast in np.float64
        expected_dtype = np.float64
    else:
        expected_dtype = X_input.dtype
    Xt = kbd.transform(X_input)
    assert Xt.dtype == expected_dtype


@pytest.mark.parametrize("input_dtype", [np.float16, np.float32, np.float64])
@pytest.mark.parametrize("encode", ["ordinal", "onehot", "onehot-dense"])
def test_32_equal_64(input_dtype, encode):
    # TODO this check is redundant with common checks and can be removed
    #  once #16290 is merged
    X_input = np.array(X, dtype=input_dtype)

    # 32 bit output
    kbd_32 = KBinsDiscretizer(
        n_bins=3,
        encode=encode,
        quantile_method="averaged_inverted_cdf",
        dtype=np.float32,
    )
    kbd_32.fit(X_input)
    Xt_32 = kbd_32.transform(X_input)

    # 64 bit output
    kbd_64 = KBinsDiscretizer(
        n_bins=3,
        encode=encode,
        quantile_method="averaged_inverted_cdf",
        dtype=np.float64,
    )
    kbd_64.fit(X_input)
    Xt_64 = kbd_64.transform(X_input)

    assert_allclose_dense_sparse(Xt_32, Xt_64)


def test_kbinsdiscretizer_subsample_default():
    # Since the size of X is small (< 2e5), subsampling will not take place.
    X = np.array([-2, 1.5, -4, -1]).reshape(-1, 1)
    kbd_default = KBinsDiscretizer(
        n_bins=10,
        encode="ordinal",
        strategy="quantile",
        quantile_method="averaged_inverted_cdf",
    )
    kbd_default.fit(X)

    kbd_without_subsampling = clone(kbd_default)
    kbd_without_subsampling.set_params(subsample=None)
    kbd_without_subsampling.fit(X)

    for bin_kbd_default, bin_kbd_with_subsampling in zip(
        kbd_default.bin_edges_[0], kbd_without_subsampling.bin_edges_[0]
    ):
        np.testing.assert_allclose(bin_kbd_default, bin_kbd_with_subsampling)
    assert kbd_default.bin_edges_.shape == kbd_without_subsampling.bin_edges_.shape


@pytest.mark.parametrize(
    "encode, expected_names",
    [
        (
            "onehot",
            [
                f"feat{col_id}_{float(bin_id)}"
                for col_id in range(3)
                for bin_id in range(4)
            ],
        ),
        (
            "onehot-dense",
            [
                f"feat{col_id}_{float(bin_id)}"
                for col_id in range(3)
                for bin_id in range(4)
            ],
        ),
        ("ordinal", [f"feat{col_id}" for col_id in range(3)]),
    ],
)
def test_kbinsdiscrtizer_get_feature_names_out(encode, expected_names):
    """Check get_feature_names_out for different settings.
    Non-regression test for #22731
    """
    X = [[-2, 1, -4], [-1, 2, -3], [0, 3, -2], [1, 4, -1]]

    kbd = KBinsDiscretizer(
        n_bins=4, encode=encode, quantile_method="averaged_inverted_cdf"
    ).fit(X)
    Xt = kbd.transform(X)

    input_features = [f"feat{i}" for i in range(3)]
    output_names = kbd.get_feature_names_out(input_features)
    assert Xt.shape[1] == output_names.shape[0]

    assert_array_equal(output_names, expected_names)


@pytest.mark.parametrize("strategy", ["uniform", "kmeans", "quantile"])
def test_kbinsdiscretizer_subsample(strategy, global_random_seed):
    # Check that the bin edges are almost the same when subsampling is used.
    X = np.random.RandomState(global_random_seed).random_sample((100000, 1)) + 1

    if strategy == "quantile":
        kbd_subsampling = KBinsDiscretizer(
            strategy=strategy,
            subsample=50000,
            random_state=global_random_seed,
            quantile_method="averaged_inverted_cdf",
        )
    else:
        kbd_subsampling = KBinsDiscretizer(
            strategy=strategy, subsample=50000, random_state=global_random_seed
        )
    kbd_subsampling.fit(X)

    kbd_no_subsampling = clone(kbd_subsampling)
    kbd_no_subsampling.set_params(subsample=None)
    kbd_no_subsampling.fit(X)

    # We use a large tolerance because we can't expect the bin edges to be exactly the
    # same when subsampling is used.
    assert_allclose(
        kbd_subsampling.bin_edges_[0], kbd_no_subsampling.bin_edges_[0], rtol=1e-2
    )


def test_invalid_quantile_method_with_sample_weight():
    X = [[-2, 1, -4], [-1, 2, -3], [0, 3, -2], [1, 4, -1]]
    expected_msg = (
        "When fitting with strategy='quantile' and sample weights, "
        "quantile_method should either be set to 'averaged_inverted_cdf' or "
        "'inverted_cdf', got quantile_method='linear' instead."
    )
    with pytest.raises(
        ValueError,
        match=expected_msg,
    ):
        KBinsDiscretizer(strategy="quantile", quantile_method="linear").fit(
            X,
            sample_weight=[1, 1, 2, 2],
        )


@pytest.mark.parametrize("strategy", ["uniform", "quantile", "kmeans"])
@pytest.mark.parametrize("encode", ["onehot", "onehot-dense", "ordinal"])
def test_kbinsdiscretizer_handle_missing_error(strategy, encode):
    """Check that handle_missing='error' raises clear ValueError on NaNs."""
    X_clean = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0], [4.0, 40.0]])
    X_nan = X_clean.copy()
    X_nan[1, 0] = np.nan

    msg = "Input X contains NaN values and `KBinsDiscretizer` is configured to error"
    kbd = KBinsDiscretizer(
        n_bins=2, encode=encode, strategy=strategy, handle_missing="error"
    )
    with pytest.raises(ValueError, match=re.escape(msg)):
        kbd.fit(X_nan)

    kbd.fit(X_clean)
    with pytest.raises(ValueError, match="Input X contains NaN"):
        kbd.transform(X_nan)

    with pytest.raises(ValueError, match="Input contains NaN"):
        kbd.inverse_transform(X_nan)


@pytest.mark.parametrize("strategy", ["uniform", "quantile", "kmeans"])
@pytest.mark.parametrize("encode", ["onehot", "onehot-dense", "ordinal"])
def test_kbinsdiscretizer_handle_missing_zeros(strategy, encode):
    """Check NaNs are ignored in fit and encoded as zeros (onehot) or NaN (ordinal)."""
    rng = np.random.RandomState(42)
    X_clean = rng.uniform(0.0, 10.0, size=(30, 3))
    X_nan = X_clean.copy()
    for row, col in [(0, 0), (5, 0), (12, 1), (20, 1), (29, 2)]:
        X_nan[row, col] = np.nan

    kbd = KBinsDiscretizer(
        n_bins=3,
        encode=encode,
        strategy=strategy,
        handle_missing="zeros",
        subsample=None,
    )
    Xt = kbd.fit_transform(X_nan)
    assert_allclose_dense_sparse(Xt, clone(kbd).fit(X_nan).transform(X_nan))

    # Bin edges must match fitting each feature on its non-NaN values alone.
    for jj in range(X_nan.shape[1]):
        valid = ~np.isnan(X_nan[:, jj])
        kbd_col = (
            clone(kbd).set_params(handle_missing="error").fit(X_nan[valid, jj : jj + 1])
        )
        assert_allclose(kbd.bin_edges_[jj], kbd_col.bin_edges_[0])

    # Transform on new data must use the test-time NaN mask, not the fit-time mask.
    X_test = X_clean[:10].copy()
    X_test[2, 2] = np.nan
    X_test[7, 0] = np.nan
    Xt_test = kbd.transform(X_test)

    nan_mask = _get_mask(X_nan, np.nan)
    if encode == "ordinal":
        assert_array_equal(np.isnan(Xt), nan_mask)
        assert_array_equal(np.isnan(Xt_test), _get_mask(X_test, np.nan))
        for jj in range(X_nan.shape[1]):
            valid_bins = Xt[~nan_mask[:, jj], jj]
            assert np.all((valid_bins >= 0) & (valid_bins < kbd.n_bins_[jj]))
    else:
        Xt_dense = Xt.toarray() if sp.issparse(Xt) else Xt
        offset = 0
        for jj, n_b in enumerate(kbd.n_bins_):
            block = Xt_dense[:, offset : offset + n_b]
            assert_array_equal(block[nan_mask[:, jj]], 0.0)
            assert_allclose(block[~nan_mask[:, jj]].sum(axis=1), 1.0)
            offset += n_b

    # inverse_transform must preserve NaNs and recover bin centers elsewhere.
    X_inv = kbd.inverse_transform(Xt)
    assert_array_equal(np.isnan(X_inv), nan_mask)
    for jj in range(X_nan.shape[1]):
        centers = (kbd.bin_edges_[jj][1:] + kbd.bin_edges_[jj][:-1]) * 0.5
        assert np.all(np.isin(X_inv[~nan_mask[:, jj], jj], centers))


@pytest.mark.parametrize(
    "strategy, quantile_method",
    [
        ("uniform", "averaged_inverted_cdf"),
        ("quantile", "averaged_inverted_cdf"),
        ("quantile", "inverted_cdf"),
        ("kmeans", "averaged_inverted_cdf"),
    ],
)
def test_kbinsdiscretizer_handle_missing_sample_weight_and_subsample(
    strategy, quantile_method, global_random_seed
):
    """Check sample_weight (including zero weights) and subsample with NaNs."""
    rng = np.random.RandomState(global_random_seed)
    X_nan = rng.normal(size=(40, 2))
    X_nan[[0, 3, 15, 27], 0] = np.nan
    X_nan[[1, 8, 19], 1] = np.nan

    # Unit weights must match unweighted fit exactly.
    ones_w = np.ones(X_nan.shape[0])
    kbd_unweighted = KBinsDiscretizer(
        n_bins=4,
        encode="ordinal",
        strategy=strategy,
        quantile_method=quantile_method,
        handle_missing="zeros",
        subsample=None,
    ).fit(X_nan)
    kbd_unit_w = clone(kbd_unweighted).fit(X_nan, sample_weight=ones_w)
    for jj in range(2):
        assert_allclose(kbd_unweighted.bin_edges_[jj], kbd_unit_w.bin_edges_[jj])

    # Non-uniform weights with zeros must match fitting on non-NaN rows with weights.
    weights = rng.choice([0.0, 1.0, 3.0], size=X_nan.shape[0])
    weights[weights == 0.0][:1] = 1.0  # ensure non-zero total weight
    kbd_weighted = clone(kbd_unweighted).fit(X_nan, sample_weight=weights)
    for jj in range(2):
        valid = ~np.isnan(X_nan[:, jj])
        kbd_col = clone(kbd_unweighted).fit(
            X_nan[valid, jj : jj + 1], sample_weight=weights[valid]
        )
        assert_allclose(kbd_weighted.bin_edges_[jj], kbd_col.bin_edges_[0])

    # Subsampling with NaNs and weights must run cleanly and preserve NaN masking.
    kbd_sub = KBinsDiscretizer(
        n_bins=4,
        encode="onehot-dense",
        strategy=strategy,
        quantile_method=quantile_method,
        handle_missing="zeros",
        subsample=20,
        random_state=global_random_seed,
    )
    Xt_sub = kbd_sub.fit_transform(X_nan, sample_weight=weights)
    assert Xt_sub.shape[0] == X_nan.shape[0]
    assert not np.any(np.isnan(Xt_sub))


@pytest.mark.parametrize("strategy", ["uniform", "quantile", "kmeans"])
@pytest.mark.parametrize("encode", ["onehot", "onehot-dense", "ordinal"])
def test_kbinsdiscretizer_handle_missing_all_nan_feature(strategy, encode):
    """Check all-NaN columns emit a warning and fall back to 1-bin encoding."""
    X_all_nan = np.array([[np.nan, 1.0], [np.nan, 2.0], [np.nan, 3.0], [np.nan, 4.0]])
    kbd = KBinsDiscretizer(
        n_bins=2, encode=encode, strategy=strategy, handle_missing="zeros"
    )
    with pytest.warns(UserWarning, match="Feature 0 has only missing values"):
        Xt = kbd.fit_transform(X_all_nan)

    assert kbd.n_bins_[0] == 1
    assert_array_equal(kbd.bin_edges_[0], [-np.inf, np.inf])

    if encode == "ordinal":
        assert np.all(np.isnan(Xt[:, 0]))
        assert not np.any(np.isnan(Xt[:, 1]))
    else:
        Xt_dense = Xt.toarray() if sp.issparse(Xt) else Xt
        assert_array_equal(Xt_dense[:, 0], 0.0)
        assert_allclose(Xt_dense[:, 1:].sum(axis=1), 1.0)

    # A non-NaN test value in the all-NaN feature falls into its single bin (0).
    Xt_new = kbd.transform(np.array([[10.0, 2.5]]))
    if encode == "ordinal":
        assert Xt_new[0, 0] == 0.0
    else:
        Xt_new_dense = Xt_new.toarray() if sp.issparse(Xt_new) else Xt_new
        assert Xt_new_dense[0, 0] == 1.0


@pytest.mark.parametrize("encode", ["onehot-dense", "ordinal"])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_kbinsdiscretizer_handle_missing_pandas_output_and_dtype(encode, dtype):
    """Check pandas set_output, get_feature_names_out, and dtype with NaNs."""
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame(
        {"a": [1.0, np.nan, 3.0, 4.0], "b": [10.0, 20.0, np.nan, 40.0]}, dtype=dtype
    )
    kbd = KBinsDiscretizer(n_bins=2, encode=encode, handle_missing="zeros").set_output(
        transform="pandas"
    )
    Xt = kbd.fit_transform(df)

    assert isinstance(Xt, pd.DataFrame)
    assert list(Xt.columns) == list(kbd.get_feature_names_out())
    assert Xt.dtypes.iloc[0] == dtype
    if encode == "ordinal":
        assert np.isnan(Xt.loc[1, "a"])
        assert np.isnan(Xt.loc[2, "b"])
    else:
        assert_array_equal(Xt.loc[1, ["a_0.0", "a_1.0"]], [0.0, 0.0])
        assert_array_equal(Xt.loc[2, ["b_0.0", "b_1.0"]], [0.0, 0.0])
