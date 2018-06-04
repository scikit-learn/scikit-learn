from nose.tools import assert_true, assert_false, assert_raises
import networkx as nx


class TestFilterFactory(object):
    def test_no_filter(self):
        nf = nx.filters.no_filter
        assert_true(nf())
        assert_true(nf(1))
        assert_true(nf(2, 1))

    def test_hide_nodes(self):
        f = nx.classes.filters.hide_nodes([1, 2, 3])
        assert_false(f(1))
        assert_false(f(2))
        assert_false(f(3))
        assert_true(f(4))
        assert_true(f(0))
        assert_true(f('a'))
        assert_raises(TypeError, f, 1, 2)
        assert_raises(TypeError, f)

    def test_show_nodes(self):
        f = nx.classes.filters.show_nodes([1, 2, 3])
        assert_true(f(1))
        assert_true(f(2))
        assert_true(f(3))
        assert_false(f(4))
        assert_false(f(0))
        assert_false(f('a'))
        assert_raises(TypeError, f, 1, 2)
        assert_raises(TypeError, f)

    def test_hide_edges(self):
        factory = nx.classes.filters.hide_edges
        f = factory([(1, 2), (3, 4)])
        assert_false(f(1, 2))
        assert_false(f(3, 4))
        assert_false(f(4, 3))
        assert_true(f(2, 3))
        assert_true(f(0, -1))
        assert_true(f('a', 'b'))
        assert_raises(TypeError, f, 1, 2, 3)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2, 3)])

    def test_show_edges(self):
        factory = nx.classes.filters.show_edges
        f = factory([(1, 2), (3, 4)])
        assert_true(f(1, 2))
        assert_true(f(3, 4))
        assert_true(f(4, 3))
        assert_false(f(2, 3))
        assert_false(f(0, -1))
        assert_false(f('a', 'b'))
        assert_raises(TypeError, f, 1, 2, 3)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2, 3)])

    def test_hide_diedges(self):
        factory = nx.classes.filters.hide_diedges
        f = factory([(1, 2), (3, 4)])
        assert_false(f(1, 2))
        assert_false(f(3, 4))
        assert_true(f(4, 3))
        assert_true(f(2, 3))
        assert_true(f(0, -1))
        assert_true(f('a', 'b'))
        assert_raises(TypeError, f, 1, 2, 3)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2, 3)])

    def test_show_diedges(self):
        factory = nx.classes.filters.show_diedges
        f = factory([(1, 2), (3, 4)])
        assert_true(f(1, 2))
        assert_true(f(3, 4))
        assert_false(f(4, 3))
        assert_false(f(2, 3))
        assert_false(f(0, -1))
        assert_false(f('a', 'b'))
        assert_raises(TypeError, f, 1, 2, 3)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2, 3)])

    def test_hide_multiedges(self):
        factory = nx.classes.filters.hide_multiedges
        f = factory([(1, 2, 0), (3, 4, 1), (1, 2, 1)])
        assert_false(f(1, 2, 0))
        assert_false(f(1, 2, 1))
        assert_true(f(1, 2, 2))
        assert_true(f(3, 4, 0))
        assert_false(f(3, 4, 1))
        assert_false(f(4, 3, 1))
        assert_true(f(4, 3, 0))
        assert_true(f(2, 3, 0))
        assert_true(f(0, -1, 0))
        assert_true(f('a', 'b', 0))
        assert_raises(TypeError, f, 1, 2, 3, 4)
        assert_raises(TypeError, f, 1, 2)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2)])
        assert_raises(ValueError, factory, [(1, 2, 3, 4)])

    def test_show_multiedges(self):
        factory = nx.classes.filters.show_multiedges
        f = factory([(1, 2, 0), (3, 4, 1), (1, 2, 1)])
        assert_true(f(1, 2, 0))
        assert_true(f(1, 2, 1))
        assert_false(f(1, 2, 2))
        assert_false(f(3, 4, 0))
        assert_true(f(3, 4, 1))
        assert_true(f(4, 3, 1))
        assert_false(f(4, 3, 0))
        assert_false(f(2, 3, 0))
        assert_false(f(0, -1, 0))
        assert_false(f('a', 'b', 0))
        assert_raises(TypeError, f, 1, 2, 3, 4)
        assert_raises(TypeError, f, 1, 2)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2)])
        assert_raises(ValueError, factory, [(1, 2, 3, 4)])

    def test_hide_multidiedges(self):
        factory = nx.classes.filters.hide_multidiedges
        f = factory([(1, 2, 0), (3, 4, 1), (1, 2, 1)])
        assert_false(f(1, 2, 0))
        assert_false(f(1, 2, 1))
        assert_true(f(1, 2, 2))
        assert_true(f(3, 4, 0))
        assert_false(f(3, 4, 1))
        assert_true(f(4, 3, 1))
        assert_true(f(4, 3, 0))
        assert_true(f(2, 3, 0))
        assert_true(f(0, -1, 0))
        assert_true(f('a', 'b', 0))
        assert_raises(TypeError, f, 1, 2, 3, 4)
        assert_raises(TypeError, f, 1, 2)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2)])
        assert_raises(ValueError, factory, [(1, 2, 3, 4)])

    def test_show_multidiedges(self):
        factory = nx.classes.filters.show_multidiedges
        f = factory([(1, 2, 0), (3, 4, 1), (1, 2, 1)])
        assert_true(f(1, 2, 0))
        assert_true(f(1, 2, 1))
        assert_false(f(1, 2, 2))
        assert_false(f(3, 4, 0))
        assert_true(f(3, 4, 1))
        assert_false(f(4, 3, 1))
        assert_false(f(4, 3, 0))
        assert_false(f(2, 3, 0))
        assert_false(f(0, -1, 0))
        assert_false(f('a', 'b', 0))
        assert_raises(TypeError, f, 1, 2, 3, 4)
        assert_raises(TypeError, f, 1, 2)
        assert_raises(TypeError, f, 1)
        assert_raises(TypeError, f)
        assert_raises(TypeError, factory, [1, 2, 3])
        assert_raises(ValueError, factory, [(1, 2)])
        assert_raises(ValueError, factory, [(1, 2, 3, 4)])
