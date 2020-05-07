import numpy as np
from numpy.testing import assert_allclose
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import pytest

from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper
from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower
from sklearn.ensemble._hist_gradient_boosting.predictor import TreePredictor
from sklearn.ensemble._hist_gradient_boosting.common import (
    G_H_DTYPE, PREDICTOR_RECORD_DTYPE, ALMOST_INF, X_BINNED_DTYPE,
    X_BITSET_INNER_DTYPE)


@pytest.mark.parametrize('n_bins', [200, 256])
def test_regression_dataset(n_bins):
    X, y = make_regression(n_samples=500, n_features=10, n_informative=5,
                           random_state=42)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=42)

    mapper = _BinMapper(n_bins=n_bins, random_state=42)
    X_train_binned = mapper.fit_transform(X_train)

    # Init gradients and hessians to that of least squares loss
    gradients = -y_train.astype(G_H_DTYPE)
    hessians = np.ones(1, dtype=G_H_DTYPE)

    min_samples_leaf = 10
    max_leaf_nodes = 30
    grower = TreeGrower(X_train_binned, gradients, hessians,
                        min_samples_leaf=min_samples_leaf,
                        max_leaf_nodes=max_leaf_nodes, n_bins=n_bins,
                        n_bins_non_missing=mapper.n_bins_non_missing_)
    grower.grow()

    predictor = grower.make_predictor(bin_thresholds=mapper.bin_thresholds_)

    assert r2_score(y_train, predictor.predict(X_train)) > 0.82
    assert r2_score(y_test, predictor.predict(X_test)) > 0.67


@pytest.mark.parametrize('threshold, expected_predictions', [
    (-np.inf, [0, 1, 1, 1]),
    (10, [0, 0, 1, 1]),
    (20, [0, 0, 0, 1]),
    (ALMOST_INF, [0, 0, 0, 1]),
    (np.inf, [0, 0, 0, 0]),
])
def test_infinite_values_and_thresholds(threshold, expected_predictions):
    # Make sure infinite values and infinite thresholds are handled properly.
    # In particular, if a value is +inf and the threshold is ALMOST_INF the
    # sample should go to the right child. If the threshold is inf (split on
    # nan), the +inf sample will go to the left child.

    X = np.array([-np.inf, 10, 20,  np.inf]).reshape(-1, 1)
    nodes = np.zeros(3, dtype=PREDICTOR_RECORD_DTYPE)

    # We just construct a simple tree with 1 root and 2 children
    # parent node
    nodes[0]['left'] = 1
    nodes[0]['right'] = 2
    nodes[0]['feature_idx'] = 0
    nodes[0]['threshold'] = threshold

    # left child
    nodes[1]['is_leaf'] = True
    nodes[1]['value'] = 0

    # right child
    nodes[2]['is_leaf'] = True
    nodes[2]['value'] = 1

    predictor = TreePredictor(nodes)
    predictions = predictor.predict(X)

    assert np.all(predictions == expected_predictions)


def _construct_bitset(thresholds):
    output = np.zeros(8, dtype=X_BITSET_INNER_DTYPE)

    for thres in thresholds:
        i1 = thres // 32
        i2 = thres % 32
        output[i1] |= X_BITSET_INNER_DTYPE(1) << X_BITSET_INNER_DTYPE(i2)

    return output


@pytest.mark.parametrize('thresholds, expected_predictions', [
    ([0, 124, 240],  [1, 0, 0, 1, 1, 0]),
    ([0, 4, 60],  [1, 1, 1, 0, 0, 0]),
    ([124, 255],  [0, 0, 0, 1, 0, 1]),
    ([10, 14],  [0, 0, 0, 0, 0, 0]),
])
def test_categorical_predictor(thresholds, expected_predictions):
    # Test predictor outputs are correct with categorical features

    cat_bitset = _construct_bitset(thresholds)
    X_binned = np.array([[0, 4, 60, 124, 240, 255]], dtype=X_BINNED_DTYPE).T
    nodes = np.zeros(3, dtype=PREDICTOR_RECORD_DTYPE)

    # We just construct a simple tree with 1 root and 2 children
    # parent node
    nodes[0]['left'] = 1
    nodes[0]['right'] = 2
    nodes[0]['feature_idx'] = 0
    nodes[0]['is_categorical'] = True
    nodes[0]['cat_bitset'] = cat_bitset

    # left child
    nodes[1]['is_leaf'] = True
    nodes[1]['value'] = 1

    # right child
    nodes[2]['is_leaf'] = True
    nodes[2]['value'] = 0

    predictor = TreePredictor(nodes)
    # missing_values_bin_idx is ignored for categories because it is already
    # encoded in cat_bitset
    predictions = predictor.predict_binned(X_binned,
                                           missing_values_bin_idx=255)

    assert_allclose(predictions, expected_predictions)
