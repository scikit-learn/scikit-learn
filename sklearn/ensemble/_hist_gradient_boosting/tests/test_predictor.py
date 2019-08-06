import numpy as np
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import pytest

from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper
from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower
from sklearn.ensemble._hist_gradient_boosting.predictor import TreePredictor
from sklearn.ensemble._hist_gradient_boosting.types import (
    G_H_DTYPE, PREDICTOR_RECORD_DTYPE)


@pytest.mark.parametrize('max_bins', [200, 256])
def test_boston_dataset(max_bins):
    boston = load_boston()
    X_train, X_test, y_train, y_test = train_test_split(
        boston.data, boston.target, random_state=42)

    mapper = _BinMapper(max_bins=max_bins, random_state=42)
    X_train_binned = mapper.fit_transform(X_train)

    # Init gradients and hessians to that of least squares loss
    gradients = -y_train.astype(G_H_DTYPE)
    hessians = np.ones(1, dtype=G_H_DTYPE)

    min_samples_leaf = 8
    max_leaf_nodes = 31
    grower = TreeGrower(X_train_binned, gradients, hessians,
                        min_samples_leaf=min_samples_leaf,
                        max_leaf_nodes=max_leaf_nodes, max_bins=max_bins,
                        actual_n_bins=mapper.actual_n_bins_)
    grower.grow()

    predictor = grower.make_predictor(bin_thresholds=mapper.bin_thresholds_)

    assert r2_score(y_train, predictor.predict(X_train)) > 0.85
    assert r2_score(y_test, predictor.predict(X_test)) > 0.70


@pytest.mark.parametrize('threshold, expected_predictions', [
    (-np.inf, [0, 1, 1, 1]),
    (10, [0, 0, 1, 1]),
    (20, [0, 0, 0, 1]),
    (np.inf, [0, 0, 0, 1]),
])
def test_infinite_values_and_thresholds(threshold, expected_predictions):
    # Make sure infinite values and infinite thresholds are handled properly.
    # In paticular, if a value is +inf and the threhsold is +inf, the sample
    # should go to the right child.

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
