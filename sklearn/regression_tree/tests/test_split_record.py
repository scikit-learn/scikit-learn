import numpy as np

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true

from sklearn.regression_tree.stats_node import StatsNode
from sklearn.regression_tree.split_record import SplitRecord


def test_split_record():
    sr = SplitRecord()
    assert_equal(type(sr.feature), int)
    assert_equal(type(sr.pos), int)
    assert_equal(type(sr.threshold), float)
    assert_equal(type(sr.impurity), float)
    assert_equal(type(sr.impurity_improvement), float)
    assert_equal(type(sr.nid), int)
    assert_equal(sr.feature, 0)
    assert_equal(sr.pos, 0)
    assert_true(np.isnan(sr.threshold))
    assert_equal(sr.impurity, np.inf)
    assert_equal(sr.impurity_improvement, -np.inf)
    assert_equal(sr.nid, 0)

    sr.reset(10., 10., 10, 10, 10, 10.,
             StatsNode(5., 5., 5, 5.),
             StatsNode(10., 10., 10, 10.),
             StatsNode(15., 15., 15, 15.))
    assert_equal(type(sr.feature), int)
    assert_equal(type(sr.pos), int)
    assert_equal(type(sr.threshold), float)
    assert_equal(type(sr.impurity), float)
    assert_equal(type(sr.impurity_improvement), float)
    assert_equal(type(sr.nid), int)
    assert_equal(sr.feature, 10)
    assert_equal(sr.pos, 10)
    assert_equal(sr.threshold, 10.)
    assert_equal(sr.impurity, 10.)
    assert_equal(sr.impurity_improvement, 10.)
    assert_equal(sr.nid, 10)

    sr.clear()
    assert_equal(type(sr.feature), int)
    assert_equal(type(sr.pos), int)
    assert_equal(type(sr.threshold), float)
    assert_equal(type(sr.impurity), float)
    assert_equal(type(sr.impurity_improvement), float)
    assert_equal(type(sr.nid), int)
    assert_equal(sr.feature, 0)
    assert_equal(sr.pos, 0)
    assert_true(np.isnan(sr.threshold))
    assert_equal(sr.impurity, np.inf)
    assert_equal(sr.impurity_improvement, -np.inf)
    assert_equal(sr.nid, 0)


def test_split_record_expandable_record():
    sr = SplitRecord()
    sr.reset(10., 10., 10, 10, 10, 10.,
             StatsNode(5., 5., 5, 5.),
             StatsNode(10., 10., 10, 10.),
             StatsNode(15., 15., 15, 15.))
    left_sr, right_sr = sr.expand_record()
    assert_equal(left_sr.c_stats, sr.l_stats)
    assert_equal(right_sr.c_stats, sr.r_stats)
