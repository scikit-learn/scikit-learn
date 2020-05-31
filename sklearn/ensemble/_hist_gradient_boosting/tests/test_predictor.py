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
    X_BITSET_INNER_DTYPE, X_DTYPE)
from sklearn.ensemble._hist_gradient_boosting._predictor_bitset import \
    PredictorBitSet
from sklearn.ensemble._hist_gradient_boosting._cat_mapper import \
    CategoryMapper


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

    predictor = grower.make_predictor(
        bin_thresholds=mapper.bin_thresholds_,
        category_mapper=CategoryMapper(n_bins - 1))

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

    predictor = TreePredictor(
        nodes, PredictorBitSet([], np.array([False], dtype=np.uint8)),
        CategoryMapper(3))
    predictions = predictor.predict(X)

    assert np.all(predictions == expected_predictions)


def _construct_bitset(bins_go_left):
    output = np.zeros(8, dtype=X_BITSET_INNER_DTYPE)

    for threshold in bins_go_left:
        i1 = threshold // 32
        i2 = threshold % 32
        output[i1] |= X_BITSET_INNER_DTYPE(1) << X_BITSET_INNER_DTYPE(i2)

    return output


@pytest.mark.parametrize(
    'bins_go_left, expected_predictions', [
        ([0, 3, 4, 6], [1, 0, 0, 1, 1, 0]),
        ([0, 1, 2, 6], [1, 1, 1, 0, 0, 0]),
        ([3, 5, 6], [0, 0, 0, 1, 0, 1])
    ])
def test_categorical_predictor(bins_go_left, expected_predictions):
    # Test predictor outputs are correct with categorical features

    X_binned = np.array([[0, 1, 2, 3, 4, 5]], dtype=X_BINNED_DTYPE).T
    category_bins = np.array([2, 5, 6, 8, 10, 15], dtype=X_DTYPE)
    nodes = np.zeros(3, dtype=PREDICTOR_RECORD_DTYPE)

    # We just construct a simple tree with 1 root and 2 children
    # parent node
    nodes[0]['left'] = 1
    nodes[0]['right'] = 2
    nodes[0]['feature_idx'] = 0
    nodes[0]['is_categorical'] = True
    nodes[0]['missing_go_to_left'] = True

    # left child
    nodes[1]['is_leaf'] = True
    nodes[1]['value'] = 1

    # right child
    nodes[2]['is_leaf'] = True
    nodes[2]['value'] = 0

    cat_bitset = _construct_bitset(bins_go_left)
    predictor_bitset = PredictorBitSet([category_bins],
                                       np.array([True], dtype=np.uint8))
    predictor_bitset.insert_categories_bitset(0, category_bins, cat_bitset)

    category_mapper = CategoryMapper(missing_values_bin_idx=6)
    category_mapper.insert(0, category_bins)
    predictor = TreePredictor(nodes, predictor_bitset, category_mapper)

    # Check binned data gives correct predictions
    prediction_binned = predictor.predict_binned(X_binned,
                                                 missing_values_bin_idx=6)
    assert_allclose(prediction_binned, expected_predictions)

    # Check with un-binned data
    predictions = predictor.predict(category_bins.reshape(-1, 1))
    assert_allclose(predictions, expected_predictions)

    # Check missing goes left because missing_values_bin_idx=6
    X_binned_missing = np.array([[6]], dtype=X_BINNED_DTYPE).T
    predictions = predictor.predict_binned(X_binned_missing,
                                           missing_values_bin_idx=6)
    assert_allclose(predictions, [1])

    # missing and unknown go left
    predictions = predictor.predict(np.array([[np.nan, 17.0]],
                                             dtype=X_DTYPE).T)
    assert_allclose(predictions, [1, 1])
