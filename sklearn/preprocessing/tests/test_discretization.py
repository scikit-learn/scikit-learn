
import pytest
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
import scipy.sparse as sp
import warnings

from sklearn.preprocessing import KBinsDiscretizer, MDLPDiscretizer
from sklearn.preprocessing import OneHotEncoder


X = [[-2, 1.5, -4, -1],
     [-1, 2.5, -3, -0.5],
     [0, 3.5, -2, 0.5],
     [1, 4.5, -1, 2]]

y = [1, 1, 2, 2]


def check_fit_transform(discretizer, expected_Xt):
    # Check the fit and transform methods
    discretizer = discretizer.fit(X, y)
    Xt = discretizer.transform(X)
    assert_array_equal(expected_Xt, Xt)


@pytest.mark.parametrize(
    "strategy, expected_Xt",
    [("uniform", [[0, 0, 0, 0], [1, 1, 1, 0], [2, 2, 2, 1], [2, 2, 2, 2]]),
     ("kmeans", [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2]]),
     ("quantile", [[0, 0, 0, 0], [1, 1, 1, 1], [2, 2, 2, 2], [2, 2, 2, 2]])])
def test_kbins_fit_transform(strategy, expected_Xt):
    # Test the fit and transform methods for KBinsDiscretizer
    kbd = KBinsDiscretizer(n_bins=3, encode="ordinal", strategy=strategy)
    check_fit_transform(kbd, expected_Xt)


def test_mdlp_fit_transform():
    # Test the fit and transform methods for MDLPDiscretizer
    mdlpd = MDLPDiscretizer(encode="ordinal")
    expected_Xt = [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [1, 1, 1, 1]]
    check_fit_transform(mdlpd, expected_Xt)


def test_valid_n_bins():
    KBinsDiscretizer(n_bins=2).fit_transform(X)
    KBinsDiscretizer(n_bins=np.array([2])[0]).fit_transform(X)
    assert KBinsDiscretizer(n_bins=2).fit(X).n_bins_.dtype == np.dtype(int)


def test_invalid_n_bins():
    est = KBinsDiscretizer(n_bins=1)
    err_msg = ("KBinsDiscretizer received an invalid "
               "number of bins. Received 1, expected at least 2.")
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)

    est = KBinsDiscretizer(n_bins=1.1)
    err_msg = ("KBinsDiscretizer received an invalid "
               "n_bins type. Received float, expected int.")
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)


def test_invalid_n_bins_array():
    # Bad shape
    n_bins = np.full((2, 4), 2.)
    est = KBinsDiscretizer(n_bins=n_bins)
    err_msg = r"n_bins must be a scalar or array of shape \(n_features,\)."
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)

    # Incorrect number of features
    n_bins = [1, 2, 2]
    est = KBinsDiscretizer(n_bins=n_bins)
    err_msg = r"n_bins must be a scalar or array of shape \(n_features,\)."
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)

    # Bad bin values
    n_bins = [1, 2, 2, 1]
    est = KBinsDiscretizer(n_bins=n_bins)
    err_msg = ("KBinsDiscretizer received an invalid number of bins "
               "at indices 0, 3. Number of bins must be at least 2, "
               "and must be an int.")
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)

    # Float bin values
    n_bins = [2.1, 2, 2.1, 2]
    est = KBinsDiscretizer(n_bins=n_bins)
    err_msg = ("KBinsDiscretizer received an invalid number of bins "
               "at indices 0, 2. Number of bins must be at least 2, "
               "and must be an int.")
    with pytest.raises(ValueError, match=err_msg):
        est.fit_transform(X)


@pytest.mark.parametrize(
    'strategy, expected',
    [('uniform', [[0, 0, 0, 0], [0, 1, 1, 0], [1, 2, 2, 1], [1, 2, 2, 2]]),
     ('kmeans', [[0, 0, 0, 0], [0, 0, 0, 0], [1, 1, 1, 1], [1, 2, 2, 2]]),
     ('quantile', [[0, 0, 0, 0], [0, 1, 1, 1], [1, 2, 2, 2], [1, 2, 2, 2]])])
def test_fit_transform_n_bins_array(strategy, expected):
    est = KBinsDiscretizer(n_bins=[2, 3, 3, 3], encode='ordinal',
                           strategy=strategy).fit(X)
    assert_array_equal(expected, est.transform(X))

    # test the shape of bin_edges_
    n_features = np.array(X).shape[1]
    assert est.bin_edges_.shape == (n_features, )
    for bin_edges, n_bins in zip(est.bin_edges_, est.n_bins_):
        assert bin_edges.shape == (n_bins + 1, )


@pytest.mark.parametrize("Discretizer", [MDLPDiscretizer, KBinsDiscretizer])
def test_invalid_n_features(Discretizer):
    # Test that an invalid number of features raises an error
    bad_X = np.arange(25).reshape(5, -1)
    est = Discretizer(encode="ordinal").fit(X, y)

    msg = (f"X has 5 features, but {Discretizer.__name__} "
           "is expecting 4 features as input.")
    with pytest.raises(ValueError, match=msg):
        est.transform(bad_X)


def check_same_min_max(discretizer):
    # Test that constant features raises an error
    warnings.simplefilter("always")

    X = np.array([[1, -2], [1, -1], [1, 0], [1, 1]])

    msg = "Feature 0 is constant and will be replaced with 0."
    with pytest.warns(UserWarning, match=msg):
        discretizer.fit(X, y)

    assert discretizer.n_bins_[0] == 1

    # Assert that the feature is replaced with zeros
    Xt = discretizer.transform(X)
    assert_array_equal(Xt[:, 0], np.zeros(X.shape[0]))


@pytest.mark.parametrize("strategy", ["uniform", "kmeans", "quantile"])
def test_kbins_same_min_max(strategy):
    # Test that KBinsDiscretizer raises an error on constant features
    kbd = KBinsDiscretizer(n_bins=3, encode="ordinal", strategy=strategy)
    check_same_min_max(kbd)


def test_mdlp_same_min_max():
    # Test that MDLPDiscretizer raises an error on constant features
    mdlpd = MDLPDiscretizer(encode="ordinal")
    check_same_min_max(mdlpd)


@pytest.mark.parametrize("Discretizer", [MDLPDiscretizer, KBinsDiscretizer])
def test_transform_1d_behavior(Discretizer):
    # Test that 1-D arrays for the input data raises an error
    X = np.arange(4)
    est = Discretizer(encode="ordinal")

    msg = "Expected 2D array, got 1D array instead"
    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)

    est = Discretizer(encode="ordinal").fit(X.reshape(-1, 1), y)

    msg = "Expected 2D array, got 1D array instead"
    with pytest.raises(ValueError, match=msg):
        est.transform(X)


def check_numeric_stability(discretizer, i):
    # Check up to discretizing nano units
    X_init = np.array([2, 4, 6, 8]).reshape(-1, 1)
    Xt_expected = np.array([0, 0, 1, 1]).reshape(-1, 1)

    X = X_init / 10**i
    Xt = discretizer.fit_transform(X, y)
    assert_array_equal(Xt_expected, Xt)


@pytest.mark.parametrize("i", range(1, 8))
def test_kbins_numeric_stability(i):
    # Test KBinsDiscretizer up to discretizing nano units
    kbd = KBinsDiscretizer(n_bins=2, encode="ordinal")
    check_numeric_stability(kbd, i)


@pytest.mark.parametrize("i", range(1, 8))
def test_mdlp_numeric_stability(i):
    # Test KBinsDiscretizer up to discretizing nano units
    mdlpd = MDLPDiscretizer(encode="ordinal")
    check_numeric_stability(mdlpd, i)


@pytest.mark.parametrize("Discretizer", [MDLPDiscretizer, KBinsDiscretizer])
def test_invalid_encode_option(Discretizer):
    # Test that an invalid encode option raises an error
    est = Discretizer(encode="invalid-encode")
    msg = (r"Valid options for 'encode' are "
           r"\('onehot', 'onehot-dense', 'ordinal'\). "
           r"Got encode='invalid-encode' instead.")

    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)


def check_encode_options(discretizer, n_bins):
    # Check the different encode options
    encoder = OneHotEncoder(categories=[np.arange(i) for i in n_bins])

    Xt_1 = discretizer.set_params(encode="ordinal").fit_transform(X, y)
    Xt_2 = discretizer.set_params(encode="onehot-dense").fit_transform(X, y)
    assert not sp.issparse(Xt_1)
    assert not sp.issparse(Xt_2)

    encoder = encoder.set_params(sparse=False)
    assert_array_equal(encoder.fit_transform(Xt_1), Xt_2)

    Xt_3 = discretizer.set_params(encode="onehot").fit_transform(X, y)
    assert sp.issparse(Xt_3)

    encoder = encoder.set_params(sparse=True)
    assert_array_equal(encoder.fit_transform(Xt_1).toarray(), Xt_3.toarray())


def test_kbins_encode_options():
    # Test the different encode options for KBinsDiscretizer
    n_bins = [2, 3, 3, 3]
    discretizer = KBinsDiscretizer(n_bins=n_bins)
    check_encode_options(discretizer, n_bins)


def test_mdlp_encode_options():
    # Test the different encode options for MDLPDiscretizer
    n_bins = [2, 2, 2, 2]
    discretizer = MDLPDiscretizer()
    check_encode_options(discretizer, n_bins)


def test_invalid_strategy_option():
    est = KBinsDiscretizer(n_bins=[2, 3, 3, 3], strategy='invalid-strategy')
    err_msg = (r"Valid options for 'strategy' are "
               r"\('uniform', 'quantile', 'kmeans'\). "
               r"Got strategy='invalid-strategy' instead.")
    with pytest.raises(ValueError, match=err_msg):
        est.fit(X)


@pytest.mark.parametrize(
    'strategy, expected_2bins, expected_3bins, expected_5bins',
    [('uniform', [0, 0, 0, 0, 1, 1], [0, 0, 0, 0, 2, 2], [0, 0, 1, 1, 4, 4]),
     ('kmeans', [0, 0, 0, 0, 1, 1], [0, 0, 1, 1, 2, 2], [0, 0, 1, 2, 3, 4]),
     ('quantile', [0, 0, 0, 1, 1, 1], [0, 0, 1, 1, 2, 2], [0, 1, 2, 3, 4, 4])])
def test_nonuniform_strategies(
        strategy, expected_2bins, expected_3bins, expected_5bins):
    X = np.array([0, 0.5, 2, 3, 9, 10]).reshape(-1, 1)

    # with 2 bins
    est = KBinsDiscretizer(n_bins=2, strategy=strategy, encode='ordinal')
    Xt = est.fit_transform(X)
    assert_array_equal(expected_2bins, Xt.ravel())

    # with 3 bins
    est = KBinsDiscretizer(n_bins=3, strategy=strategy, encode='ordinal')
    Xt = est.fit_transform(X)
    assert_array_equal(expected_3bins, Xt.ravel())

    # with 5 bins
    est = KBinsDiscretizer(n_bins=5, strategy=strategy, encode='ordinal')
    Xt = est.fit_transform(X)
    assert_array_equal(expected_5bins, Xt.ravel())


def check_inverse_transform(discretizer, expected_inv):
    # Check the inverse_transform method
    Xt = discretizer.fit_transform(X, y)
    Xinv = discretizer.inverse_transform(Xt)
    assert_array_almost_equal(expected_inv, Xinv)


@pytest.mark.parametrize(
    'strategy, expected_inv',
    [('uniform', [[-1.5, 2., -3.5, -0.5], [-0.5, 3., -2.5, -0.5],
                  [0.5, 4., -1.5, 0.5], [0.5, 4., -1.5, 1.5]]),
     ('kmeans', [[-1.375, 2.125, -3.375, -0.5625],
                 [-1.375, 2.125, -3.375, -0.5625],
                 [-0.125, 3.375, -2.125, 0.5625],
                 [0.75, 4.25, -1.25, 1.625]]),
     ('quantile', [[-1.5, 2., -3.5, -0.75], [-0.5, 3., -2.5, 0.],
                   [0.5, 4., -1.5, 1.25], [0.5, 4., -1.5, 1.25]])])
@pytest.mark.parametrize('encode', ['ordinal', 'onehot', 'onehot-dense'])
def test_kbins_inverse_transform(strategy, encode, expected_inv):
    # Test the inverse_transform method for KBinsDiscretizer
    kbd = KBinsDiscretizer(n_bins=3, strategy=strategy, encode=encode)
    check_inverse_transform(kbd, expected_inv)


@pytest.mark.parametrize("encode", ["ordinal", "onehot", "onehot-dense"])
def test_mdlp_inverse_transform(encode):
    # Test the inverse_transform method for MDLPDiscretizer
    mdlpd = MDLPDiscretizer(encode=encode)
    expected_inv = [[-1.25, 2.25, -3.25, -0.5],
                    [-1.25, 2.25, -3.25, -0.5],
                    [0.25, 3.75, -1.75, 1],
                    [0.25, 3.75, -1.75, 1]]

    check_inverse_transform(mdlpd, expected_inv)


def check_transform_outside_fit_range(discretizer):
    # Check the transform method when the data
    # to discretize is outside of fit range
    X = np.array([0, 1, 2, 3])[:, None]
    discretizer.fit(X, y)

    X2 = np.array([-2, 5])[:, None]
    X2t = discretizer.transform(X2)
    assert_array_equal(X2t.max(axis=0) + 1, discretizer.n_bins_)
    assert_array_equal(X2t.min(axis=0), [0])


@pytest.mark.parametrize("strategy", ["uniform", "kmeans", "quantile"])
def test_kbins_transform_outside_fit_range(strategy):
    # Test the transform method for KBinsDiscretizer when
    # the data to discretize is outside of fit range
    kbd = KBinsDiscretizer(n_bins=4, strategy=strategy, encode="ordinal")

    check_transform_outside_fit_range(kbd)


def test_mdlp_transform_outside_fit_range():
    # Test the transform method for MDLPDiscretizer when
    # the data to discretize is outside of fit range
    mdlpd = MDLPDiscretizer(encode="ordinal")

    check_transform_outside_fit_range(mdlpd)


@pytest.mark.parametrize("Discretizer", [MDLPDiscretizer, KBinsDiscretizer])
def test_overwrite(Discretizer):
    # Test that the input data is not overwritten
    discretizer = Discretizer(encode="ordinal")

    X_before = X.copy()
    Xt = discretizer.fit_transform(X, y)
    assert_array_equal(X, X_before)

    Xt_before = Xt.copy()
    discretizer.inverse_transform(Xt)
    assert_array_equal(Xt, Xt_before)


@pytest.mark.parametrize(
    'strategy, expected_bin_edges',
    [('quantile', [0, 1, 3]), ('kmeans', [0, 1.5, 3])])
def test_redundant_bins(strategy, expected_bin_edges):
    X = [[0], [0], [0], [0], [3], [3]]
    kbd = KBinsDiscretizer(n_bins=3, strategy=strategy)
    msg = (r"Bins whose width are too small \(i.e., <= 1e-8\) in feature "
           r"0 are removed. Consider decreasing the number of bins.")
    with pytest.warns(UserWarning, match=msg):
        kbd.fit(X)
    assert_array_almost_equal(kbd.bin_edges_[0], expected_bin_edges)


def test_percentile_numeric_stability():
    X = np.array([0.05, 0.05, 0.95]).reshape(-1, 1)
    bin_edges = np.array([0.05, 0.23, 0.41, 0.59, 0.77, 0.95])
    Xt = np.array([0, 0, 4]).reshape(-1, 1)
    kbd = KBinsDiscretizer(n_bins=10, encode='ordinal',
                           strategy='quantile')
    msg = (r"Bins whose width are too small \(i.e., <= 1e-8\) in feature "
           r"0 are removed. Consider decreasing the number of bins.")
    with pytest.warns(UserWarning, match=msg):
        kbd.fit(X)
    assert_array_almost_equal(kbd.bin_edges_[0], bin_edges)
    assert_array_almost_equal(kbd.transform(X), Xt)


@pytest.mark.parametrize("Discretizer", [MDLPDiscretizer, KBinsDiscretizer])
@pytest.mark.parametrize("in_dtype", [np.float16, np.float32, np.float64])
@pytest.mark.parametrize("out_dtype", [None, np.float16, np.float32,
                                       np.float64])
@pytest.mark.parametrize("encode", ["ordinal", "onehot", "onehot-dense"])
def test_consistent_dtype(Discretizer, in_dtype, out_dtype, encode):
    # Test that the output data type is consistent
    X_input = np.array(X, dtype=in_dtype)
    discretizer = Discretizer(encode=encode, dtype=out_dtype)

    if out_dtype not in [None, np.float32, np.float64]:
        msg = "Valid options for 'dtype' are"
        # An error is raised if a wrong data type is defined
        with pytest.raises(ValueError, match=msg):
            discretizer.fit(X_input, y)
    else:
        discretizer.fit(X_input, y)

        if out_dtype is not None:
            expected_dtype = out_dtype
        elif out_dtype is None and X_input.dtype == np.float16:
            # Wrong numeric data types are cast to np.float64
            expected_dtype = np.float64
        else:
            expected_dtype = X_input.dtype

        # Test the output data type
        Xt = discretizer.transform(X_input)
        assert Xt.dtype == expected_dtype
