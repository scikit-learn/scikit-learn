import numpy as np
import pytest
from pytest import approx

from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower
from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper
from sklearn.ensemble._hist_gradient_boosting.common import X_BINNED_DTYPE
from sklearn.ensemble._hist_gradient_boosting.common import Y_DTYPE
from sklearn.ensemble._hist_gradient_boosting.common import G_H_DTYPE


def _make_training_data(n_bins=256, constant_hessian=True):
    rng = np.random.RandomState(42)
    n_samples = 10000

    # Generate some test data directly binned so as to test the grower code
    # independently of the binning logic.
    X_binned = rng.randint(0, n_bins - 1, size=(n_samples, 2),
                           dtype=X_BINNED_DTYPE)
    X_binned = np.asfortranarray(X_binned)

    def true_decision_function(input_features):
        """Ground truth decision function

        This is a very simple yet asymmetric decision tree. Therefore the
        grower code should have no trouble recovering the decision function
        from 10000 training samples.
        """
        if input_features[0] <= n_bins // 2:
            return -1
        else:
            return -1 if input_features[1] <= n_bins // 3 else 1

    target = np.array([true_decision_function(x) for x in X_binned],
                      dtype=Y_DTYPE)

    # Assume a square loss applied to an initial model that always predicts 0
    # (hardcoded for this test):
    all_gradients = target.astype(G_H_DTYPE)
    shape_hessians = 1 if constant_hessian else all_gradients.shape
    all_hessians = np.ones(shape=shape_hessians, dtype=G_H_DTYPE)

    return X_binned, all_gradients, all_hessians


def _check_children_consistency(parent, left, right):
    # Make sure the samples are correctly dispatched from a parent to its
    # children
    assert parent.left_child is left
    assert parent.right_child is right

    # each sample from the parent is propagated to one of the two children
    assert (len(left.sample_indices) + len(right.sample_indices)
            == len(parent.sample_indices))

    assert (set(left.sample_indices).union(set(right.sample_indices))
            == set(parent.sample_indices))

    # samples are sent either to the left or the right node, never to both
    assert (set(left.sample_indices).intersection(set(right.sample_indices))
            == set())


@pytest.mark.parametrize(
    'n_bins, constant_hessian, stopping_param, shrinkage',
    [
        (11, True, "min_gain_to_split", 0.5),
        (11, False, "min_gain_to_split", 1.),
        (11, True, "max_leaf_nodes", 1.),
        (11, False, "max_leaf_nodes", 0.1),
        (42, True, "max_leaf_nodes", 0.01),
        (42, False, "max_leaf_nodes", 1.),
        (256, True, "min_gain_to_split", 1.),
        (256, True, "max_leaf_nodes", 0.1),
    ]
)
def test_grow_tree(n_bins, constant_hessian, stopping_param, shrinkage):
    X_binned, all_gradients, all_hessians = _make_training_data(
        n_bins=n_bins, constant_hessian=constant_hessian)
    n_samples = X_binned.shape[0]

    if stopping_param == "max_leaf_nodes":
        stopping_param = {"max_leaf_nodes": 3}
    else:
        stopping_param = {"min_gain_to_split": 0.01}

    grower = TreeGrower(X_binned, all_gradients, all_hessians,
                        n_bins=n_bins, shrinkage=shrinkage,
                        min_samples_leaf=1, **stopping_param)

    # The root node is not yet splitted, but the best possible split has
    # already been evaluated:
    assert grower.root.left_child is None
    assert grower.root.right_child is None

    root_split = grower.root.split_info
    assert root_split.feature_idx == 0
    assert root_split.bin_idx == n_bins // 2
    assert len(grower.splittable_nodes) == 1

    # Calling split next applies the next split and computes the best split
    # for each of the two newly introduced children nodes.
    left_node, right_node = grower.split_next()

    # All training samples have ben splitted in the two nodes, approximately
    # 50%/50%
    _check_children_consistency(grower.root, left_node, right_node)
    assert len(left_node.sample_indices) > 0.4 * n_samples
    assert len(left_node.sample_indices) < 0.6 * n_samples

    if grower.min_gain_to_split > 0:
        # The left node is too pure: there is no gain to split it further.
        assert left_node.split_info.gain < grower.min_gain_to_split
        assert left_node in grower.finalized_leaves

    # The right node can still be splitted further, this time on feature #1
    split_info = right_node.split_info
    assert split_info.gain > 1.
    assert split_info.feature_idx == 1
    assert split_info.bin_idx == n_bins // 3
    assert right_node.left_child is None
    assert right_node.right_child is None

    # The right split has not been applied yet. Let's do it now:
    assert len(grower.splittable_nodes) == 1
    right_left_node, right_right_node = grower.split_next()
    _check_children_consistency(right_node, right_left_node, right_right_node)
    assert len(right_left_node.sample_indices) > 0.1 * n_samples
    assert len(right_left_node.sample_indices) < 0.2 * n_samples

    assert len(right_right_node.sample_indices) > 0.2 * n_samples
    assert len(right_right_node.sample_indices) < 0.4 * n_samples

    # All the leafs are pure, it is not possible to split any further:
    assert not grower.splittable_nodes

    grower._apply_shrinkage()

    # Check the values of the leaves:
    assert grower.root.left_child.value == approx(shrinkage)
    assert grower.root.right_child.left_child.value == approx(shrinkage)
    assert grower.root.right_child.right_child.value == approx(-shrinkage,
                                                               rel=1e-3)


def test_predictor_from_grower():
    # Build a tree on the toy 3-leaf dataset to extract the predictor.
    n_bins = 256
    X_binned, all_gradients, all_hessians = _make_training_data(
        n_bins=n_bins)
    grower = TreeGrower(X_binned, all_gradients, all_hessians,
                        n_bins=n_bins, shrinkage=1.,
                        max_leaf_nodes=3, min_samples_leaf=5)
    grower.grow()
    assert grower.n_nodes == 5  # (2 decision nodes + 3 leaves)

    # Check that the node structure can be converted into a predictor
    # object to perform predictions at scale
    predictor = grower.make_predictor()
    assert predictor.nodes.shape[0] == 5
    assert predictor.nodes['is_leaf'].sum() == 3

    # Probe some predictions for each leaf of the tree
    # each group of 3 samples corresponds to a condition in _make_training_data
    input_data = np.array([
        [0, 0],
        [42, 99],
        [128, 254],

        [129, 0],
        [129, 85],
        [254, 85],

        [129, 86],
        [129, 254],
        [242, 100],
    ], dtype=np.uint8)
    missing_values_bin_idx = n_bins - 1
    predictions = predictor.predict_binned(input_data, missing_values_bin_idx)
    expected_targets = [1, 1, 1, 1, 1, 1, -1, -1, -1]
    assert np.allclose(predictions, expected_targets)

    # Check that training set can be recovered exactly:
    predictions = predictor.predict_binned(X_binned, missing_values_bin_idx)
    assert np.allclose(predictions, -all_gradients)


@pytest.mark.parametrize(
    'n_samples, min_samples_leaf, n_bins, constant_hessian, noise',
    [
        (11, 10, 7, True, 0),
        (13, 10, 42, False, 0),
        (56, 10, 255, True, 0.1),
        (101, 3, 7, True, 0),
        (200, 42, 42, False, 0),
        (300, 55, 255, True, 0.1),
        (300, 301, 255, True, 0.1),
    ]
)
def test_min_samples_leaf(n_samples, min_samples_leaf, n_bins,
                          constant_hessian, noise):
    rng = np.random.RandomState(seed=0)
    # data = linear target, 3 features, 1 irrelevant.
    X = rng.normal(size=(n_samples, 3))
    y = X[:, 0] - X[:, 1]
    if noise:
        y_scale = y.std()
        y += rng.normal(scale=noise, size=n_samples) * y_scale
    mapper = _BinMapper(n_bins=n_bins)
    X = mapper.fit_transform(X)

    all_gradients = y.astype(G_H_DTYPE)
    shape_hessian = 1 if constant_hessian else all_gradients.shape
    all_hessians = np.ones(shape=shape_hessian, dtype=G_H_DTYPE)
    grower = TreeGrower(X, all_gradients, all_hessians,
                        n_bins=n_bins, shrinkage=1.,
                        min_samples_leaf=min_samples_leaf,
                        max_leaf_nodes=n_samples)
    grower.grow()
    predictor = grower.make_predictor(
        bin_thresholds=mapper.bin_thresholds_)

    if n_samples >= min_samples_leaf:
        for node in predictor.nodes:
            if node['is_leaf']:
                assert node['count'] >= min_samples_leaf
    else:
        assert predictor.nodes.shape[0] == 1
        assert predictor.nodes[0]['is_leaf']
        assert predictor.nodes[0]['count'] == n_samples


@pytest.mark.parametrize('n_samples, min_samples_leaf', [
                         (99, 50),
                         (100, 50)])
def test_min_samples_leaf_root(n_samples, min_samples_leaf):
    # Make sure root node isn't split if n_samples is not at least twice
    # min_samples_leaf
    rng = np.random.RandomState(seed=0)

    n_bins = 256

    # data = linear target, 3 features, 1 irrelevant.
    X = rng.normal(size=(n_samples, 3))
    y = X[:, 0] - X[:, 1]
    mapper = _BinMapper(n_bins=n_bins)
    X = mapper.fit_transform(X)

    all_gradients = y.astype(G_H_DTYPE)
    all_hessians = np.ones(shape=1, dtype=G_H_DTYPE)
    grower = TreeGrower(X, all_gradients, all_hessians,
                        n_bins=n_bins, shrinkage=1.,
                        min_samples_leaf=min_samples_leaf,
                        max_leaf_nodes=n_samples)
    grower.grow()
    if n_samples >= min_samples_leaf * 2:
        assert len(grower.finalized_leaves) >= 2
    else:
        assert len(grower.finalized_leaves) == 1


def assert_is_stump(grower):
    # To assert that stumps are created when max_depth=1
    for leaf in (grower.root.left_child, grower.root.right_child):
        assert leaf.left_child is None
        assert leaf.right_child is None


@pytest.mark.parametrize('max_depth', [1, 2, 3])
def test_max_depth(max_depth):
    # Make sure max_depth parameter works as expected
    rng = np.random.RandomState(seed=0)

    n_bins = 256
    n_samples = 1000

    # data = linear target, 3 features, 1 irrelevant.
    X = rng.normal(size=(n_samples, 3))
    y = X[:, 0] - X[:, 1]
    mapper = _BinMapper(n_bins=n_bins)
    X = mapper.fit_transform(X)

    all_gradients = y.astype(G_H_DTYPE)
    all_hessians = np.ones(shape=1, dtype=G_H_DTYPE)
    grower = TreeGrower(X, all_gradients, all_hessians, max_depth=max_depth)
    grower.grow()

    depth = max(leaf.depth for leaf in grower.finalized_leaves)
    assert depth == max_depth

    if max_depth == 1:
        assert_is_stump(grower)


def test_input_validation():

    X_binned, all_gradients, all_hessians = _make_training_data()

    X_binned_float = X_binned.astype(np.float32)
    with pytest.raises(NotImplementedError,
                       match="X_binned must be of type uint8"):
        TreeGrower(X_binned_float, all_gradients, all_hessians)

    X_binned_C_array = np.ascontiguousarray(X_binned)
    with pytest.raises(
            ValueError,
            match="X_binned should be passed as Fortran contiguous array"):
        TreeGrower(X_binned_C_array, all_gradients, all_hessians)


def test_init_parameters_validation():
    X_binned, all_gradients, all_hessians = _make_training_data()
    with pytest.raises(ValueError,
                       match="min_gain_to_split=-1 must be positive"):

        TreeGrower(X_binned, all_gradients, all_hessians,
                   min_gain_to_split=-1)

    with pytest.raises(ValueError,
                       match="min_hessian_to_split=-1 must be positive"):
        TreeGrower(X_binned, all_gradients, all_hessians,
                   min_hessian_to_split=-1)


def test_missing_value_predict_only():
    # Make sure that missing values are supported at predict time even if they
    # were not encountered in the training data: the missing values are
    # assigned to whichever child has the most samples.

    rng = np.random.RandomState(0)
    n_samples = 100
    X_binned = rng.randint(0, 256, size=(n_samples, 1), dtype=np.uint8)
    X_binned = np.asfortranarray(X_binned)

    gradients = rng.normal(size=n_samples).astype(G_H_DTYPE)
    hessians = np.ones(shape=1, dtype=G_H_DTYPE)

    grower = TreeGrower(X_binned, gradients, hessians, min_samples_leaf=5,
                        has_missing_values=False)
    grower.grow()

    predictor = grower.make_predictor()

    # go from root to a leaf, always following node with the most samples.
    # That's the path nans are supposed to take
    node = predictor.nodes[0]
    while not node['is_leaf']:
        left = predictor.nodes[node['left']]
        right = predictor.nodes[node['right']]
        node = left if left['count'] > right['count'] else right

    prediction_main_path = node['value']

    # now build X_test with only nans, and make sure all predictions are equal
    # to prediction_main_path
    all_nans = np.full(shape=(n_samples, 1), fill_value=np.nan)
    assert np.all(predictor.predict(all_nans) == prediction_main_path)


def test_split_on_nan_with_infinite_values():
    # Make sure the split on nan situations are respected even when there are
    # samples with +inf values (we set the threshold to +inf when we have a
    # split on nan so this test makes sure this does not introduce edge-case
    # bugs). We need to use the private API so that we can also test
    # predict_binned().

    X = np.array([0, 1, np.inf, np.nan, np.nan]).reshape(-1, 1)
    # the gradient values will force a split on nan situation
    gradients = np.array([0, 0, 0, 100, 100], dtype=G_H_DTYPE)
    hessians = np.ones(shape=1, dtype=G_H_DTYPE)

    bin_mapper = _BinMapper()
    X_binned = bin_mapper.fit_transform(X)

    n_bins_non_missing = 3
    has_missing_values = True
    grower = TreeGrower(X_binned, gradients, hessians,
                        n_bins_non_missing=n_bins_non_missing,
                        has_missing_values=has_missing_values,
                        min_samples_leaf=1)

    grower.grow()

    predictor = grower.make_predictor(
        bin_thresholds=bin_mapper.bin_thresholds_
    )

    # sanity check: this was a split on nan
    assert predictor.nodes[0]['threshold'] == np.inf
    assert predictor.nodes[0]['bin_threshold'] == n_bins_non_missing - 1

    # Make sure in particular that the +inf sample is mapped to the left child
    # Note that lightgbm "fails" here and will assign the inf sample to the
    # right child, even though it's a "split on nan" situation.
    predictions = predictor.predict(X)
    predictions_binned = predictor.predict_binned(
        X_binned, missing_values_bin_idx=bin_mapper.missing_values_bin_idx_)
    np.testing.assert_allclose(predictions, -gradients)
    np.testing.assert_allclose(predictions_binned, -gradients)
