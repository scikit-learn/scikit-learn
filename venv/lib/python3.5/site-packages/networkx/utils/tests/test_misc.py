# -*- encoding: utf-8 -*-
from nose.tools import *
from nose import SkipTest
import networkx as nx
from networkx.utils import *


def test_is_string_like():
    assert_true(is_string_like("aaaa"))
    assert_false(is_string_like(None))
    assert_false(is_string_like(123))


def test_iterable():
    assert_false(iterable(None))
    assert_false(iterable(10))
    assert_true(iterable([1, 2, 3]))
    assert_true(iterable((1, 2, 3)))
    assert_true(iterable({1: "A", 2: "X"}))
    assert_true(iterable("ABC"))


def test_graph_iterable():
    K = nx.complete_graph(10)
    assert_true(iterable(K))
    assert_true(iterable(K.nodes()))
    assert_true(iterable(K.edges()))


def test_is_list_of_ints():
    assert_true(is_list_of_ints([1, 2, 3, 42]))
    assert_false(is_list_of_ints([1, 2, 3, "kermit"]))


def test_random_number_distribution():
    # smoke test only
    z = powerlaw_sequence(20, exponent=2.5)
    z = discrete_sequence(20, distribution=[0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3])


def test_make_str_with_bytes():
    import sys
    PY2 = sys.version_info[0] == 2

    x = "qualité"
    y = make_str(x)
    if PY2:
        assert_true(isinstance(y, unicode))
        # Since file encoding is utf-8, the é will be two bytes.
        assert_true(len(y) == 8)
    else:
        assert_true(isinstance(y, str))
        assert_true(len(y) == 7)


def test_make_str_with_unicode():
    import sys
    PY2 = sys.version_info[0] == 2
    if PY2:
        x = unicode("qualité", encoding='utf-8')
        y = make_str(x)
        assert_true(isinstance(y, unicode))
        assert_true(len(y) == 7)
    else:
        x = "qualité"
        y = make_str(x)
        assert_true(isinstance(y, str))
        assert_true(len(y) == 7)


class TestNumpyArray(object):
    @classmethod
    def setupClass(cls):
        global numpy
        global assert_allclose
        try:
            import numpy
            from numpy.testing import assert_allclose
        except ImportError:
            raise SkipTest('NumPy not available.')

    def test_dict_to_numpy_array1(self):
        d = {'a': 1, 'b': 2}
        a = dict_to_numpy_array1(d, mapping={'a': 0, 'b': 1})
        assert_allclose(a, numpy.array([1, 2]))
        a = dict_to_numpy_array1(d, mapping={'b': 0, 'a': 1})
        assert_allclose(a, numpy.array([2, 1]))

        a = dict_to_numpy_array1(d)
        assert_allclose(a.sum(), 3)

    def test_dict_to_numpy_array2(self):
        d = {'a': {'a': 1, 'b': 2},
             'b': {'a': 10, 'b': 20}}

        mapping = {'a': 1, 'b': 0}
        a = dict_to_numpy_array2(d, mapping=mapping)
        assert_allclose(a, numpy.array([[20, 10], [2, 1]]))

        a = dict_to_numpy_array2(d)
        assert_allclose(a.sum(), 33)

    def test_dict_to_numpy_array_a(self):
        d = {'a': {'a': 1, 'b': 2},
             'b': {'a': 10, 'b': 20}}

        mapping = {'a': 0, 'b': 1}
        a = dict_to_numpy_array(d, mapping=mapping)
        assert_allclose(a, numpy.array([[1, 2], [10, 20]]))

        mapping = {'a': 1, 'b': 0}
        a = dict_to_numpy_array(d, mapping=mapping)
        assert_allclose(a, numpy.array([[20, 10], [2, 1]]))

        a = dict_to_numpy_array2(d)
        assert_allclose(a.sum(), 33)

    def test_dict_to_numpy_array_b(self):
        d = {'a': 1, 'b': 2}

        mapping = {'a': 0, 'b': 1}
        a = dict_to_numpy_array(d, mapping=mapping)
        assert_allclose(a, numpy.array([1, 2]))

        a = dict_to_numpy_array1(d)
        assert_allclose(a.sum(), 3)


def test_pairwise():
    nodes = range(4)
    node_pairs = [(0, 1), (1, 2), (2, 3)]
    node_pairs_cycle = node_pairs + [(3, 0)]
    assert_equal(list(pairwise(nodes)), node_pairs)
    assert_equal(list(pairwise(iter(nodes))), node_pairs)
    assert_equal(list(pairwise(nodes, cyclic=True)), node_pairs_cycle)
    empty_iter = iter(())
    assert_equal(list(pairwise(empty_iter)), [])
    empty_iter = iter(())
    assert_equal(list(pairwise(empty_iter, cyclic=True)), [])


def test_groups():
    many_to_one = dict(zip('abcde', [0, 0, 1, 1, 2]))
    actual = groups(many_to_one)
    expected = {0: {'a', 'b'}, 1: {'c', 'd'}, 2: {'e'}}
    assert_equal(actual, expected)
    assert_equal({}, groups({}))


def test_to_tuple():
    a_list = [1, 2, [1, 3]]
    actual = to_tuple(a_list)
    expected = (1, 2, (1, 3))
    assert_equal(actual, expected)

    a_tuple = (1, 2)
    actual = to_tuple(a_tuple)
    expected = a_tuple
    assert_equal(actual, expected)

    a_mix = (1, 2, [1, 3])
    actual = to_tuple(a_mix)
    expected = (1, 2, (1, 3))
    assert_equal(actual, expected)


def test_create_random_state():
    try:
        import numpy as np
    except ImportError:
        raise SkipTest('numpy not available.')

    rs = np.random.RandomState

    assert_true(isinstance(create_random_state(1), rs))
    assert_true(isinstance(create_random_state(None), rs))
    assert_true(isinstance(create_random_state(np.random), rs))
    assert_true(isinstance(create_random_state(rs(1)), rs))
    assert_raises(ValueError, create_random_state, 'a')

    assert_true(np.all((rs(1).rand(10), create_random_state(1).rand(10))))
