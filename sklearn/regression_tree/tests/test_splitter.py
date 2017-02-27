import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true

from sklearn.regression_tree.splitter import Splitter
from sklearn.regression_tree.split_record import SplitRecord
from sklearn.regression_tree.stats_node import StatsNode


def test_split_record():
    # create synthetic data
    # the data are already ordered
    X = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]).T
    y = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    sample_weight = np.ones(y.shape)
    sum_total_weighted_samples = sample_weight.sum()

    root_stats = StatsNode(
        sum_y=0.,
        sum_sq_y=0.,
        n_samples=0,
        sum_weighted_samples=0)
    parent_split_record = SplitRecord()
    parent_split_record.c_stats = root_stats

    feature_idx = 10
    start_idx = 10

    # parameters of the tree
    min_samples_leaf = 2
    min_weight_leaf = 2.

    sp = Splitter(X, y, sample_weight, sum_total_weighted_samples, feature_idx,
                  start_idx, parent_split_record, min_samples_leaf,
                  min_weight_leaf)
    assert_true(sp.X is X)
    assert_true(sp.y is y)
    assert_true(sp.sample_weight is sample_weight)
    assert_equal(sp.sum_total_weighted_samples, sum_total_weighted_samples)
    assert_equal(sp.feature_idx, feature_idx)
    assert_equal(sp.start_idx, start_idx)
    assert_equal(sp.min_samples_leaf, min_samples_leaf)
    assert_equal(sp.min_weight_leaf, min_weight_leaf)
    # have to be copy
    assert_false(sp.split_record is parent_split_record)
    assert_false(sp.best_split_record is parent_split_record)

    # set the real value to check the split later on
    root_stats = StatsNode(
        sum_y=np.sum(np.ravel(y) * sample_weight),
        sum_sq_y=np.sum(np.ravel(y ** 2) * sample_weight),
        n_samples=y.size,
        sum_weighted_samples=sum_total_weighted_samples)
    parent_split_record = SplitRecord()
    parent_split_record.c_stats = root_stats
    feature_idx = 0
    start_idx = 0
    sp.reset(feature_idx=feature_idx, start_idx=start_idx,
             split_record=parent_split_record)

    assert_true(sp.X is X)
    assert_true(sp.y is y)
    assert_true(sp.sample_weight is sample_weight)
    assert_equal(sp.sum_total_weighted_samples, sum_total_weighted_samples)
    assert_equal(sp.feature_idx, feature_idx)
    assert_equal(sp.start_idx, start_idx)
    assert_equal(sp.min_samples_leaf, min_samples_leaf)
    assert_equal(sp.min_weight_leaf, min_weight_leaf)
    # have to be copy
    assert_false(sp.split_record is parent_split_record)
    assert_false(sp.best_split_record is parent_split_record)

    for i in range(y.size):
        sp.node_evaluate_split(i)

    assert_equal(sp.best_split_record.pos, 4)
    assert_equal(sp.best_split_record.feature, 0)
    assert_equal(sp.best_split_record.threshold, 4.5)
    assert_equal(sp.best_split_record.impurity_improvement, 0.25)
    assert_equal(sp.best_split_record.l_stats, StatsNode(0., 0., 5, 5.))
    assert_equal(sp.best_split_record.r_stats, StatsNode(5., 5., 5, 5.))
