from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_true

from sklearn.regression_tree.stats_node import StatsNode


def test_stats_node():
    sn = StatsNode(10, 10, 10., 10)
    assert_equal(type(sn.sum_y), float)
    assert_equal(type(sn.sum_sq_y), float)
    assert_equal(type(sn.n_samples), int)
    assert_equal(type(sn.sum_weighted_samples), float)
    assert_equal(sn.sum_y, 10)
    assert_equal(sn.sum_sq_y, 10)
    assert_equal(sn.n_samples, 10)
    assert_equal(sn.sum_weighted_samples, 10)

    sn.reset(5, 5, 5., 5)
    assert_equal(type(sn.sum_y), float)
    assert_equal(type(sn.sum_sq_y), float)
    assert_equal(type(sn.n_samples), int)
    assert_equal(type(sn.sum_weighted_samples), float)
    assert_equal(sn.sum_y, 5)
    assert_equal(sn.sum_sq_y, 5)
    assert_equal(sn.n_samples, 5)
    assert_equal(sn.sum_weighted_samples, 5)

    sn.clear()
    assert_equal(type(sn.sum_y), float)
    assert_equal(type(sn.sum_sq_y), float)
    assert_equal(type(sn.n_samples), int)
    assert_equal(type(sn.sum_weighted_samples), float)
    assert_equal(sn.sum_y, 0)
    assert_equal(sn.sum_sq_y, 0)
    assert_equal(sn.n_samples, 0)
    assert_equal(sn.sum_weighted_samples, 0)


def test_stats_node_add():
    sn = StatsNode(1, 2, 3, 4)
    sn2 = StatsNode(4, 4, 4, 4)
    sn_id = id(sn)
    sn += sn2
    assert_false(id(sn) == sn_id)
    assert_equal(type(sn.sum_y), float)
    assert_equal(type(sn.sum_sq_y), float)
    assert_equal(type(sn.n_samples), int)
    assert_equal(type(sn.sum_weighted_samples), float)
    assert_equal(sn.sum_y, 5)
    assert_equal(sn.sum_sq_y, 6)
    assert_equal(sn.n_samples, 7)
    assert_equal(sn.sum_weighted_samples, 8)

    assert_equal(sn.__iadd__((1, 2, 3, 4)), NotImplemented)

    sn_id = id(sn)
    sn = sn + sn2
    assert_false(id(sn) == sn_id)
    assert_equal(type(sn.sum_y), float)
    assert_equal(type(sn.sum_sq_y), float)
    assert_equal(type(sn.n_samples), int)
    assert_equal(type(sn.sum_weighted_samples), float)
    assert_equal(sn.sum_y, 9)
    assert_equal(sn.sum_sq_y, 10)
    assert_equal(sn.n_samples, 11)
    assert_equal(sn.sum_weighted_samples, 12)

    assert_equal(sn.__add__((1, 2, 3, 4)), NotImplemented)


def test_stats_node_sub():
    sn = StatsNode(1, 2, 3, 4)
    sn2 = StatsNode(4, 4, 4, 4)
    sn_id = id(sn)
    sn -= sn2
    assert_false(id(sn) == sn_id)
    assert_equal(type(sn.sum_y), float)
    assert_equal(type(sn.sum_sq_y), float)
    assert_equal(type(sn.n_samples), int)
    assert_equal(type(sn.sum_weighted_samples), float)
    assert_equal(sn.sum_y, -3)
    assert_equal(sn.sum_sq_y, -2)
    assert_equal(sn.n_samples, -1)
    assert_equal(sn.sum_weighted_samples, 0)

    assert_equal(sn.__isub__((1, 2, 3, 4)), NotImplemented)

    sn_id = id(sn)
    sn = sn - sn2
    assert_false(id(sn) == sn_id)
    assert_equal(type(sn.sum_y), float)
    assert_equal(type(sn.sum_sq_y), float)
    assert_equal(type(sn.n_samples), int)
    assert_equal(type(sn.sum_weighted_samples), float)
    assert_equal(sn.sum_y, -7)
    assert_equal(sn.sum_sq_y, -6)
    assert_equal(sn.n_samples, -5)
    assert_equal(sn.sum_weighted_samples, -4)

    assert_equal(sn.__sub__((1, 2, 3, 4)), NotImplemented)


def test_stats_node_repr():
    sn = StatsNode(1, 2, 3, 4)
    assert_equal(sn.__str__(), 'Sum of y: 1.0\nSum of the squared y:'
                 ' 2.0\nNumber of samples: 3\nSum of weights associated'
                 ' to samples: 4.0\n')
