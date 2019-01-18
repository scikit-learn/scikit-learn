import numpy as np
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_almost_equal
import pytest

from sklearn._fast_gradient_boosting.types import HISTOGRAM_DTYPE
from sklearn._fast_gradient_boosting.types import Y_DTYPE
from sklearn._fast_gradient_boosting.types import X_BINNED_DTYPE
from sklearn._fast_gradient_boosting.splitting import Splitter


@pytest.mark.parametrize('n_bins', [3, 32, 256])
def test_histogram_split(n_bins):
    rng = np.random.RandomState(42)
    feature_idx = 0
    l2_regularization = 0
    min_hessian_to_split = 1e-3
    min_samples_leaf = 1
    min_gain_to_split = 0.
    X_binned = np.asfortranarray(
        rng.randint(0, n_bins, size=(int(1e4), 2)), dtype=X_BINNED_DTYPE)
    binned_feature = X_binned.T[feature_idx]
    sample_indices = np.arange(binned_feature.shape[0], dtype=np.uint32)
    ordered_hessians = np.ones_like(binned_feature, dtype=Y_DTYPE)
    all_hessians = ordered_hessians
    sum_hessians = all_hessians.sum()

    for true_bin in range(1, n_bins - 1):
        for sign in [-1, 1]:
            ordered_gradients = np.full_like(binned_feature, sign,
                                             dtype=Y_DTYPE)
            ordered_gradients[binned_feature <= true_bin] *= -1
            all_gradients = ordered_gradients
            sum_gradients = all_gradients.sum()

            n_bins_per_feature = np.array([n_bins] * X_binned.shape[1],
                                          dtype=np.uint32)
            splitter = Splitter(X_binned,
                                n_bins,
                                n_bins_per_feature,
                                all_gradients, all_hessians,
                                l2_regularization,
                                min_hessian_to_split,
                                min_samples_leaf, min_gain_to_split)

            histograms = np.zeros(shape=(1, n_bins), dtype=HISTOGRAM_DTYPE)
            split_info = splitter.find_best_split_wrapper(
                feature_idx, sample_indices, histograms, sum_gradients,
                sum_hessians)

            assert split_info.bin_idx == true_bin
            assert split_info.gain >= 0
            assert split_info.feature_idx == feature_idx
            assert (split_info.n_samples_left + split_info.n_samples_right
                    == sample_indices.shape[0])
            # Constant hessian: 1. per sample.
            assert split_info.n_samples_left == split_info.hessian_left


@pytest.mark.parametrize('constant_hessian', [True, False])
def test_split_vs_split_subtraction(constant_hessian):
    # Make sure find_node_split and find_node_split_subtraction return the
    # same results.
    # Should we add a test about computation time to make sure
    # time(subtraction) < time(regular)?
    rng = np.random.RandomState(42)

    n_bins = 10
    n_features = 20
    n_samples = 500
    l2_regularization = 0.
    min_hessian_to_split = 1e-3
    min_samples_leaf = 1
    min_gain_to_split = 0.

    X_binned = rng.randint(0, n_bins, size=(n_samples, n_features),
                           dtype=X_BINNED_DTYPE)
    X_binned = np.asfortranarray(X_binned)
    sample_indices = np.arange(n_samples, dtype=np.uint32)
    all_gradients = rng.randn(n_samples).astype(Y_DTYPE)
    if constant_hessian:
        all_hessians = np.ones(1, dtype=Y_DTYPE)
    else:
        all_hessians = rng.lognormal(size=n_samples).astype(Y_DTYPE)

    n_bins_per_feature = np.array([n_bins] * X_binned.shape[1],
                                  dtype=np.uint32)
    splitter = Splitter(X_binned, n_bins, n_bins_per_feature, all_gradients,
                        all_hessians, l2_regularization, min_hessian_to_split,
                        min_samples_leaf, min_gain_to_split)

    hists_parent = np.zeros(shape=(n_features, n_bins), dtype=HISTOGRAM_DTYPE)
    hists_left = np.zeros(shape=(n_features, n_bins), dtype=HISTOGRAM_DTYPE)
    hists_right = np.zeros(shape=(n_features, n_bins), dtype=HISTOGRAM_DTYPE)
    hists_left_sub = np.zeros(shape=(n_features, n_bins),
                              dtype=HISTOGRAM_DTYPE)
    hists_right_sub = np.zeros(shape=(n_features, n_bins),
                               dtype=HISTOGRAM_DTYPE)

    # first split parent, left and right with classical method
    si_parent = splitter.find_node_split(sample_indices, hists_parent)
    sample_indices_left, sample_indices_right, _ = splitter.split_indices(
        si_parent, sample_indices)
    si_left = splitter.find_node_split(sample_indices_left, hists_left)
    si_right = splitter.find_node_split(sample_indices_right, hists_right)

    # split left with subtraction method
    si_left_sub = splitter.find_node_split_subtraction(
        sample_indices_left, si_parent.gradient_left,
        si_parent.hessian_left, hists_parent, hists_right, hists_left_sub)

    # split right with subtraction method
    si_right_sub = splitter.find_node_split_subtraction(
        sample_indices_right, si_parent.gradient_right,
        si_parent.hessian_right, hists_parent, hists_left, hists_right_sub)

    # make sure histograms from classical and subtraction method are the same
    for hists, hists_sub in ((hists_left, hists_left_sub),
                             (hists_right, hists_right_sub)):
        for hist, hist_sub in zip(hists, hists_sub):
            for key in ('count', 'sum_hessians', 'sum_gradients'):
                assert_array_almost_equal(hist[key], hist_sub[key], decimal=4)

    # make sure split_infos from classical and subtraction method are the same
    for si, si_sub in ((si_left, si_left_sub), (si_right, si_right_sub)):
        assert_almost_equal(si.gain, si_sub.gain, decimal=3)
        assert_almost_equal(si.feature_idx, si_sub.feature_idx, decimal=3)
        assert_almost_equal(si.gradient_left, si_sub.gradient_left, decimal=3)
        assert_almost_equal(si.gradient_right, si_sub.gradient_right,
                            decimal=3)
        assert_almost_equal(si.hessian_right, si_sub.hessian_right, decimal=3)
        assert_almost_equal(si.hessian_left, si_sub.hessian_left, decimal=3)


@pytest.mark.parametrize('constant_hessian', [True, False])
def test_gradient_and_hessian_sanity(constant_hessian):
    # This test checks that the values of gradients and hessians are
    # consistent in different places:
    # - in split_info: si.gradient_left + si.gradient_right must be equal to
    #   the gradient at the node. Same for hessians.
    # - in the histograms: summing 'sum_gradients' over the bins must be
    #   constant across all features, and those sums must be equal to the
    #   node's gradient. Same for hessians.
    #
    # These checks are carried out for split_info and histograms resulting
    # from both find_node_split() and find_node_split_subtraction().
    #
    # The structure of this test is exactly the same as in
    # test_split_vs_split_subtraction() but it's probably best to keep them
    # separate because they're not checking the same things.

    rng = np.random.RandomState(42)

    n_bins = 10
    n_features = 20
    n_samples = 500
    l2_regularization = 0.
    min_hessian_to_split = 1e-3
    min_samples_leaf = 1
    min_gain_to_split = 0.

    X_binned = rng.randint(0, n_bins, size=(n_samples, n_features),
                           dtype=X_BINNED_DTYPE)
    X_binned = np.asfortranarray(X_binned)
    sample_indices = np.arange(n_samples, dtype=np.uint32)
    all_gradients = rng.randn(n_samples).astype(Y_DTYPE)
    if constant_hessian:
        all_hessians = np.ones(1, dtype=Y_DTYPE)
    else:
        all_hessians = rng.lognormal(size=n_samples).astype(Y_DTYPE)

    n_bins_per_feature = np.array([n_bins] * X_binned.shape[1],
                                  dtype=np.uint32)
    splitter = Splitter(X_binned, n_bins,
                        n_bins_per_feature,
                        all_gradients, all_hessians,
                        l2_regularization, min_hessian_to_split,
                        min_samples_leaf, min_gain_to_split)

    hists_parent = np.zeros(shape=(n_features, n_bins), dtype=HISTOGRAM_DTYPE)
    hists_left = np.zeros(shape=(n_features, n_bins), dtype=HISTOGRAM_DTYPE)
    hists_right = np.zeros(shape=(n_features, n_bins), dtype=HISTOGRAM_DTYPE)
    hists_left_sub = np.zeros(shape=(n_features, n_bins),
                              dtype=HISTOGRAM_DTYPE)
    hists_right_sub = np.zeros(shape=(n_features, n_bins),
                               dtype=HISTOGRAM_DTYPE)
    # first split parent, left and right with classical method
    si_parent = splitter.find_node_split(sample_indices, hists_parent)
    sample_indices_left, sample_indices_right, _ = splitter.split_indices(
        si_parent, sample_indices)

    si_left = splitter.find_node_split(sample_indices_left, hists_left)
    si_right = splitter.find_node_split(sample_indices_right, hists_right)

    # split left with subtraction method
    si_left_sub = splitter.find_node_split_subtraction(
        sample_indices_left, si_parent.gradient_left,
        si_parent.hessian_left, hists_parent, hists_right, hists_left_sub)

    # split right with subtraction method
    si_right_sub = splitter.find_node_split_subtraction(
        sample_indices_right, si_parent.gradient_right,
        si_parent.hessian_right, hists_parent, hists_left, hists_right_sub)

    # make sure that si.gradient_left + si.gradient_right have their expected
    # value, same for hessians
    for si, indices in (
            (si_parent, sample_indices),
            (si_left, sample_indices_left),
            (si_left_sub, sample_indices_left),
            (si_right, sample_indices_right),
            (si_right_sub, sample_indices_right)):
        gradient = si.gradient_right + si.gradient_left
        expected_gradient = all_gradients[indices].sum()
        hessian = si.hessian_right + si.hessian_left
        if constant_hessian:
            expected_hessian = indices.shape[0] * all_hessians[0]
        else:
            expected_hessian = all_hessians[indices].sum()

        assert_almost_equal(gradient, expected_gradient, decimal=3)
        assert_almost_equal(hessian, expected_hessian, decimal=3)

    # make sure sum of gradients in histograms are the same for all features,
    # and make sure they're equal to their expected value
    for hists, indices in (
            (hists_parent, sample_indices),
            (hists_left, sample_indices_left),
            (hists_left_sub, sample_indices_left),
            (hists_right, sample_indices_right),
            (hists_right_sub, sample_indices_right)):
        # note: gradients and hessians have shape (n_features,),
        # we're comparing them to *scalars*. This has the benefit of also
        # making sure that all the entries are equal.
        gradients = hists['sum_gradients'].sum(axis=1)  # shape = (n_features,)
        expected_gradient = all_gradients[indices].sum()  # scalar
        hessians = hists['sum_hessians'].sum(axis=1)
        if constant_hessian:
            # 0 is not the actual hessian, but it's not computed in this case
            expected_hessian = 0.
        else:
            expected_hessian = all_hessians[indices].sum()

        assert_almost_equal(gradients, expected_gradient, decimal=4)
        assert_almost_equal(hessians, expected_hessian, decimal=4)


def test_split_indices():
    # Check that split_indices returns the correct splits and that
    # splitter.partition is consistent with what is returned.
    rng = np.random.RandomState(421)

    n_bins = 5
    n_samples = 10
    l2_regularization = 0.
    min_hessian_to_split = 1e-3
    min_samples_leaf = 1
    min_gain_to_split = 0.

    # split will happen on feature 1 and on bin 3
    X_binned = [[0, 0],
                [0, 3],
                [0, 4],
                [0, 0],
                [0, 0],
                [0, 0],
                [0, 0],
                [0, 4],
                [0, 0],
                [0, 4]]
    X_binned = np.asfortranarray(X_binned, dtype=X_BINNED_DTYPE)
    sample_indices = np.arange(n_samples, dtype=np.uint32)
    all_gradients = rng.randn(n_samples).astype(Y_DTYPE)
    all_hessians = np.ones(1, dtype=Y_DTYPE)

    n_bins_per_feature = np.array([n_bins] * X_binned.shape[1],
                                  dtype=np.uint32)
    splitter = Splitter(X_binned, n_bins,
                        n_bins_per_feature,
                        all_gradients, all_hessians,
                        l2_regularization, min_hessian_to_split,
                        min_samples_leaf, min_gain_to_split)

    assert_array_almost_equal(sample_indices, splitter.partition)

    histograms = np.zeros(shape=(2, n_bins), dtype=HISTOGRAM_DTYPE)
    si_root = splitter.find_node_split(sample_indices, histograms)

    # sanity checks for best split
    assert si_root.feature_idx == 1
    assert si_root.bin_idx == 3

    samples_left, samples_right, position_right = splitter.split_indices(
        si_root, splitter.partition)
    assert set(samples_left) == set([0, 1, 3, 4, 5, 6, 8])
    assert set(samples_right) == set([2, 7, 9])

    assert_array_almost_equal(samples_left,
                              splitter.partition[:position_right])
    assert_array_almost_equal(samples_right,
                              splitter.partition[position_right:])

    # Check that the resulting split indices sizes are consistent with the
    # count statistics anticipated when looking for the best split.
    assert samples_left.shape[0] == si_root.n_samples_left
    assert samples_right.shape[0] == si_root.n_samples_right


def test_min_gain_to_split():
    # Try to split a pure node (all gradients are equal, same for hessians)
    # with min_gain_to_split = 0 and make sure that the node is not split (best
    # possible gain = -1). Note: before the strict inequality comparison, this
    # test would fail because the node would be split with a gain of 0.
    rng = np.random.RandomState(42)
    l2_regularization = 0
    min_hessian_to_split = 0
    min_samples_leaf = 1
    min_gain_to_split = 0.
    n_bins = 255
    n_samples = 100
    X_binned = np.asfortranarray(
        rng.randint(0, n_bins, size=(n_samples, 1)), dtype=X_BINNED_DTYPE)
    binned_feature = X_binned[:, 0]
    sample_indices = np.arange(n_samples, dtype=np.uint32)
    all_hessians = np.ones_like(binned_feature, dtype=Y_DTYPE)
    all_gradients = np.ones_like(binned_feature, dtype=Y_DTYPE)

    n_bins_per_feature = np.array([n_bins] * X_binned.shape[1],
                                  dtype=np.uint32)
    splitter = Splitter(X_binned, n_bins, n_bins_per_feature,
                        all_gradients, all_hessians,
                        l2_regularization,
                        min_hessian_to_split,
                        min_samples_leaf, min_gain_to_split)

    histograms = np.zeros(shape=(1, n_bins), dtype=HISTOGRAM_DTYPE)
    split_info = splitter.find_node_split(sample_indices, histograms)
    assert split_info.gain == -1
