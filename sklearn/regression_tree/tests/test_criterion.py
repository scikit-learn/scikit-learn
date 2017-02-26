from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_raises

from sklearn.regression_tree.criterion import _impurity_mse
from sklearn.regression_tree.criterion import impurity_improvement
from sklearn.regression_tree.stats_node import StatsNode
from sklearn.regression_tree.split_record import SplitRecord


def test_impurity_mse():
    sn = StatsNode(2., 4., 4, 4.)
    assert_equal(_impurity_mse(sn), 0.75)
    sn.reset(0., 0., 10, 10.)
    assert_equal(_impurity_mse(sn), 0.)


def test_impurity_mse_div_zero():
    sn = StatsNode(0., 0., 0, 0.)
    assert_raises_regex(ZeroDivisionError, "float division by zero",
                        _impurity_mse, sn)


def test_impurity_improvement_mse():
    sr = SplitRecord()
    # [2 2 -2 -2]
    sr.c_stats = StatsNode(0., 16., 4, 4.)
    sr.l_stats = StatsNode(4., 8., 2, 2.)
    sr.r_stats = sr.c_stats - sr.l_stats
    assert_equal(impurity_improvement(sr, 4.), 4.0)


def test_impurity_improvement_error():
    sr = SplitRecord()
    sr.c_stats = StatsNode(2., 4., 4, 4.)
    sr.l_stats = StatsNode(1., 2., 2, 2.)
    sr.r_stats = sr.c_stats - sr.l_stats
    assert_raises(NotImplementedError, impurity_improvement,
                  sr, 10, criterion='wrong')
    assert_raises_regex(ValueError, "The number of samples in the root nodes"
                        " is less than in the current children node.",
                        impurity_improvement, sr, 1.)
