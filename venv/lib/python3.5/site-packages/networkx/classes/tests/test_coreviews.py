from nose.tools import assert_equal, assert_not_equal, assert_is,\
    assert_is_not, assert_true, assert_false, assert_raises
import tempfile
import pickle

import networkx as nx


class TestAtlasView(object):
    # node->data
    def setup(self):
        self.d = {0: {'color': 'blue', 'weight': 1.2}, 1: {}, 2: {'color': 1}}
        self.av = nx.classes.coreviews.AtlasView(self.d)

    def test_pickle(self):
        view = self.av
        pview = pickle.loads(pickle.dumps(view, -1))
        assert_equal(view, pview)
        assert_equal(view.__slots__, pview.__slots__)

    def test_len(self):
        assert_equal(len(self.av), len(self.d))

    def test_iter(self):
        assert_equal(list(self.av), list(self.d))

    def test_getitem(self):
        assert_is(self.av[1], self.d[1])
        assert_equal(self.av[2]['color'], 1)
        assert_raises(KeyError, self.av.__getitem__, 3)

    def test_copy(self):
        avcopy = self.av.copy()
        assert_equal(avcopy[0], self.av[0])
        assert_equal(avcopy, self.av)
        assert_is_not(avcopy[0], self.av[0])
        assert_is_not(avcopy, self.av)
        avcopy[5] = {}
        assert_not_equal(avcopy, self.av)

        avcopy[0]['ht'] = 4
        assert_not_equal(avcopy[0], self.av[0])
        self.av[0]['ht'] = 4
        assert_equal(avcopy[0], self.av[0])
        del self.av[0]['ht']

        assert_false(hasattr(self.av, '__setitem__'))

    def test_items(self):
        assert_equal(sorted(self.av.items()), sorted(self.d.items()))

    def test_str(self):
        out = str(self.d)
        assert_equal(str(self.av), out)

    def test_repr(self):
        out = "AtlasView(" + str(self.d) + ")"
        assert_equal(repr(self.av), out)


class TestAdjacencyView(object):
    # node->nbr->data
    def setup(self):
        dd = {'color': 'blue', 'weight': 1.2}
        self.nd = {0: dd, 1: {}, 2: {'color': 1}}
        self.adj = {3: self.nd, 0: {3: dd}, 1: {}, 2: {3: {'color': 1}}}
        self.adjview = nx.classes.coreviews.AdjacencyView(self.adj)

    def test_pickle(self):
        view = self.adjview
        pview = pickle.loads(pickle.dumps(view, -1))
        assert_equal(view, pview)
        assert_equal(view.__slots__, pview.__slots__)

    def test_len(self):
        assert_equal(len(self.adjview), len(self.adj))

    def test_iter(self):
        assert_equal(list(self.adjview), list(self.adj))

    def test_getitem(self):
        assert_is_not(self.adjview[1], self.adj[1])
        assert_is(self.adjview[3][0], self.adjview[0][3])
        assert_equal(self.adjview[2][3]['color'], 1)
        assert_raises(KeyError, self.adjview.__getitem__, 4)

    def test_copy(self):
        avcopy = self.adjview.copy()
        assert_equal(avcopy[0], self.adjview[0])
        assert_is_not(avcopy[0], self.adjview[0])

        avcopy[2][3]['ht'] = 4
        assert_not_equal(avcopy[2], self.adjview[2])
        self.adjview[2][3]['ht'] = 4
        assert_equal(avcopy[2], self.adjview[2])
        del self.adjview[2][3]['ht']

        assert_false(hasattr(self.adjview, '__setitem__'))

    def test_items(self):
        view_items = sorted((n, dict(d)) for n, d in self.adjview.items())
        assert_equal(view_items, sorted(self.adj.items()))

    def test_str(self):
        out = str(dict(self.adj))
        assert_equal(str(self.adjview), out)

    def test_repr(self):
        out = self.adjview.__class__.__name__ + "(" + str(self.adj) + ")"
        assert_equal(repr(self.adjview), out)


class TestMultiAdjacencyView(TestAdjacencyView):
    # node->nbr->key->data
    def setup(self):
        dd = {'color': 'blue', 'weight': 1.2}
        self.kd = {0: dd, 1: {}, 2: {'color': 1}}
        self.nd = {3: self.kd, 0: {3: dd}, 1: {0: {}}, 2: {3: {'color': 1}}}
        self.adj = {3: self.nd, 0: {3: {3: dd}}, 1: {}, 2: {3: {8: {}}}}
        self.adjview = nx.classes.coreviews.MultiAdjacencyView(self.adj)

    def test_getitem(self):
        assert_is_not(self.adjview[1], self.adj[1])
        assert_is(self.adjview[3][0][3], self.adjview[0][3][3])
        assert_equal(self.adjview[3][2][3]['color'], 1)
        assert_raises(KeyError, self.adjview.__getitem__, 4)

    def test_copy(self):
        avcopy = self.adjview.copy()
        assert_equal(avcopy[0], self.adjview[0])
        assert_is_not(avcopy[0], self.adjview[0])

        avcopy[2][3][8]['ht'] = 4
        assert_not_equal(avcopy[2], self.adjview[2])
        self.adjview[2][3][8]['ht'] = 4
        assert_equal(avcopy[2], self.adjview[2])
        del self.adjview[2][3][8]['ht']

        assert_false(hasattr(self.adjview, '__setitem__'))


class TestUnionAtlas(object):
    # node->data
    def setup(self):
        self.s = {0: {'color': 'blue', 'weight': 1.2}, 1: {}, 2: {'color': 1}}
        self.p = {3: {'color': 'blue', 'weight': 1.2}, 4: {}, 2: {'watch': 2}}
        self.av = nx.classes.coreviews.UnionAtlas(self.s, self.p)

    def test_pickle(self):
        view = self.av
        pview = pickle.loads(pickle.dumps(view, -1))
        assert_equal(view, pview)
        assert_equal(view.__slots__, pview.__slots__)

    def test_len(self):
        assert_equal(len(self.av), len(self.s) + len(self.p))

    def test_iter(self):
        assert_equal(set(self.av), set(self.s) | set(self.p))

    def test_getitem(self):
        assert_is(self.av[0], self.s[0])
        assert_is(self.av[4], self.p[4])
        assert_equal(self.av[2]['color'], 1)
        assert_raises(KeyError, self.av[2].__getitem__, 'watch')
        assert_raises(KeyError, self.av.__getitem__, 8)

    def test_copy(self):
        avcopy = self.av.copy()
        assert_equal(avcopy[0], self.av[0])
        assert_is_not(avcopy[0], self.av[0])
        assert_is_not(avcopy, self.av)
        avcopy[5] = {}
        assert_not_equal(avcopy, self.av)

        avcopy[0]['ht'] = 4
        assert_not_equal(avcopy[0], self.av[0])
        self.av[0]['ht'] = 4
        assert_equal(avcopy[0], self.av[0])
        del self.av[0]['ht']

        assert_false(hasattr(self.av, '__setitem__'))

    def test_items(self):
        expected = dict(self.p.items())
        expected.update(self.s)
        assert_equal(sorted(self.av.items()), sorted(expected.items()))

    def test_str(self):
        out = str(dict(self.av))
        assert_equal(str(self.av), out)

    def test_repr(self):
        out = "{}({}, {})".format(self.av.__class__.__name__, self.s, self.p)
        assert_equal(repr(self.av), out)


class TestUnionAdjacency(object):
    # node->nbr->data
    def setup(self):
        dd = {'color': 'blue', 'weight': 1.2}
        self.nd = {0: dd, 1: {}, 2: {'color': 1}}
        self.s = {3: self.nd, 0: {}, 1: {}, 2: {3: {'color': 1}}}
        self.p = {3: {}, 0: {3: dd}, 1: {0: {}}, 2: {1: {'color': 1}}}
        self.adjview = nx.classes.coreviews.UnionAdjacency(self.s, self.p)

    def test_pickle(self):
        view = self.adjview
        pview = pickle.loads(pickle.dumps(view, -1))
        assert_equal(view, pview)
        assert_equal(view.__slots__, pview.__slots__)

    def test_len(self):
        assert_equal(len(self.adjview), len(self.s))

    def test_iter(self):
        assert_equal(sorted(self.adjview), sorted(self.s))

    def test_getitem(self):
        assert_is_not(self.adjview[1], self.s[1])
        assert_is(self.adjview[3][0], self.adjview[0][3])
        assert_equal(self.adjview[2][3]['color'], 1)
        assert_raises(KeyError, self.adjview.__getitem__, 4)

    def test_copy(self):
        avcopy = self.adjview.copy()
        assert_equal(avcopy[0], self.adjview[0])
        assert_is_not(avcopy[0], self.adjview[0])

        avcopy[2][3]['ht'] = 4
        assert_not_equal(avcopy[2], self.adjview[2])
        self.adjview[2][3]['ht'] = 4
        assert_equal(avcopy[2], self.adjview[2])
        del self.adjview[2][3]['ht']

        assert_false(hasattr(self.adjview, '__setitem__'))

    def test_str(self):
        out = str(dict(self.adjview))
        assert_equal(str(self.adjview), out)

    def test_repr(self):
        clsname = self.adjview.__class__.__name__
        out = "{}({}, {})".format(clsname, self.s, self.p)
        assert_equal(repr(self.adjview), out)


class TestUnionMultiInner(TestUnionAdjacency):
    # nbr->key->data
    def setup(self):
        dd = {'color': 'blue', 'weight': 1.2}
        self.kd = {7: {}, 'ekey': {}, 9: {'color': 1}}
        self.s = {3: self.kd, 0: {7: dd}, 1: {}, 2: {'key': {'color': 1}}}
        self.p = {3: {}, 0: {3: dd}, 1: {}, 2: {1: {'span': 2}}}
        self.adjview = nx.classes.coreviews.UnionMultiInner(self.s, self.p)

    def test_len(self):
        assert_equal(len(self.adjview), len(self.s) + len(self.p))

    def test_getitem(self):
        assert_is_not(self.adjview[1], self.s[1])
        assert_is(self.adjview[0][7], self.adjview[0][3])
        assert_equal(self.adjview[2]['key']['color'], 1)
        assert_equal(self.adjview[2][1]['span'], 2)
        assert_raises(KeyError, self.adjview.__getitem__, 4)
        assert_raises(KeyError, self.adjview[1].__getitem__, 'key')

    def test_copy(self):
        avcopy = self.adjview.copy()
        assert_equal(avcopy[0], self.adjview[0])
        assert_is_not(avcopy[0], self.adjview[0])

        avcopy[2][1]['width'] = 8
        assert_not_equal(avcopy[2], self.adjview[2])
        self.adjview[2][1]['width'] = 8
        assert_equal(avcopy[2], self.adjview[2])
        del self.adjview[2][1]['width']

        assert_false(hasattr(self.adjview, '__setitem__'))
        assert_true(hasattr(avcopy, '__setitem__'))


class TestUnionMultiAdjacency(TestUnionAdjacency):
    # node->nbr->key->data
    def setup(self):
        dd = {'color': 'blue', 'weight': 1.2}
        self.kd = {7: {}, 8: {}, 9: {'color': 1}}
        self.nd = {3: self.kd, 0: {9: dd}, 1: {8: {}}, 2: {9: {'color': 1}}}
        self.s = {3: self.nd, 0: {3: {7: dd}}, 1: {}, 2: {3: {8: {}}}}
        self.p = {3: {}, 0: {3: {9: dd}}, 1: {}, 2: {1: {8: {}}}}
        self.adjview = nx.classes.coreviews.UnionMultiAdjacency(self.s, self.p)

    def test_getitem(self):
        assert_is_not(self.adjview[1], self.s[1])
        assert_is(self.adjview[3][0][9], self.adjview[0][3][9])
        assert_equal(self.adjview[3][2][9]['color'], 1)
        assert_raises(KeyError, self.adjview.__getitem__, 4)

    def test_copy(self):
        avcopy = self.adjview.copy()
        assert_equal(avcopy[0], self.adjview[0])
        assert_is_not(avcopy[0], self.adjview[0])

        avcopy[2][3][8]['ht'] = 4
        assert_not_equal(avcopy[2], self.adjview[2])
        self.adjview[2][3][8]['ht'] = 4
        assert_equal(avcopy[2], self.adjview[2])
        del self.adjview[2][3][8]['ht']

        assert_false(hasattr(self.adjview, '__setitem__'))
        assert_true(hasattr(avcopy, '__setitem__'))
