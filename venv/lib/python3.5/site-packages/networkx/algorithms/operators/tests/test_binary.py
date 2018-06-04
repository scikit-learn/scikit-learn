from nose.tools import *
import networkx as nx
from networkx import *
from networkx.testing import *


def test_union_attributes():
    g = nx.Graph()
    g.add_node(0, x=4)
    g.add_node(1, x=5)
    g.add_edge(0, 1, size=5)
    g.graph['name'] = 'g'

    h = g.copy()
    h.graph['name'] = 'h'
    h.graph['attr'] = 'attr'
    h.nodes[0]['x'] = 7

    gh = nx.union(g, h, rename=('g', 'h'))
    assert_equal(set(gh.nodes()), set(['h0', 'h1', 'g0', 'g1']))
    for n in gh:
        graph, node = n
        assert_equal(gh.nodes[n], eval(graph).nodes[int(node)])

    assert_equal(gh.graph['attr'], 'attr')
    assert_equal(gh.graph['name'], 'h')  # h graph attributes take precendent


def test_intersection():
    G = nx.Graph()
    H = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4])
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    H.add_nodes_from([1, 2, 3, 4])
    H.add_edge(2, 3)
    H.add_edge(3, 4)
    I = nx.intersection(G, H)
    assert_equal(set(I.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(I.edges()), [(2, 3)])


def test_intersection_attributes():
    g = nx.Graph()
    g.add_node(0, x=4)
    g.add_node(1, x=5)
    g.add_edge(0, 1, size=5)
    g.graph['name'] = 'g'

    h = g.copy()
    h.graph['name'] = 'h'
    h.graph['attr'] = 'attr'
    h.nodes[0]['x'] = 7

    gh = nx.intersection(g, h)
    assert_equal(set(gh.nodes()), set(g.nodes()))
    assert_equal(set(gh.nodes()), set(h.nodes()))
    assert_equal(sorted(gh.edges()), sorted(g.edges()))

    h.remove_node(0)
    assert_raises(nx.NetworkXError, nx.intersection, g, h)


def test_intersection_multigraph_attributes():
    g = nx.MultiGraph()
    g.add_edge(0, 1, key=0)
    g.add_edge(0, 1, key=1)
    g.add_edge(0, 1, key=2)
    h = nx.MultiGraph()
    h.add_edge(0, 1, key=0)
    h.add_edge(0, 1, key=3)
    gh = nx.intersection(g, h)
    assert_equal(set(gh.nodes()), set(g.nodes()))
    assert_equal(set(gh.nodes()), set(h.nodes()))
    assert_equal(sorted(gh.edges()), [(0, 1)])
    assert_equal(sorted(gh.edges(keys=True)), [(0, 1, 0)])


def test_difference():
    G = nx.Graph()
    H = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4])
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    H.add_nodes_from([1, 2, 3, 4])
    H.add_edge(2, 3)
    H.add_edge(3, 4)
    D = nx.difference(G, H)
    assert_equal(set(D.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(D.edges()), [(1, 2)])
    D = nx.difference(H, G)
    assert_equal(set(D.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(D.edges()), [(3, 4)])
    D = nx.symmetric_difference(G, H)
    assert_equal(set(D.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(D.edges()), [(1, 2), (3, 4)])


def test_difference2():
    G = nx.Graph()
    H = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4])
    H.add_nodes_from([1, 2, 3, 4])
    G.add_edge(1, 2)
    H.add_edge(1, 2)
    G.add_edge(2, 3)
    D = nx.difference(G, H)
    assert_equal(set(D.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(D.edges()), [(2, 3)])
    D = nx.difference(H, G)
    assert_equal(set(D.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(D.edges()), [])
    H.add_edge(3, 4)
    D = nx.difference(H, G)
    assert_equal(set(D.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(D.edges()), [(3, 4)])


def test_difference_attributes():
    g = nx.Graph()
    g.add_node(0, x=4)
    g.add_node(1, x=5)
    g.add_edge(0, 1, size=5)
    g.graph['name'] = 'g'

    h = g.copy()
    h.graph['name'] = 'h'
    h.graph['attr'] = 'attr'
    h.nodes[0]['x'] = 7

    gh = nx.difference(g, h)
    assert_equal(set(gh.nodes()), set(g.nodes()))
    assert_equal(set(gh.nodes()), set(h.nodes()))
    assert_equal(sorted(gh.edges()), [])

    h.remove_node(0)
    assert_raises(nx.NetworkXError, nx.intersection, g, h)


def test_difference_multigraph_attributes():
    g = nx.MultiGraph()
    g.add_edge(0, 1, key=0)
    g.add_edge(0, 1, key=1)
    g.add_edge(0, 1, key=2)
    h = nx.MultiGraph()
    h.add_edge(0, 1, key=0)
    h.add_edge(0, 1, key=3)
    gh = nx.difference(g, h)
    assert_equal(set(gh.nodes()), set(g.nodes()))
    assert_equal(set(gh.nodes()), set(h.nodes()))
    assert_equal(sorted(gh.edges()), [(0, 1), (0, 1)])
    assert_equal(sorted(gh.edges(keys=True)), [(0, 1, 1), (0, 1, 2)])


@raises(nx.NetworkXError)
def test_difference_raise():
    G = nx.path_graph(4)
    H = nx.path_graph(3)
    GH = nx.difference(G, H)


def test_symmetric_difference_multigraph():
    g = nx.MultiGraph()
    g.add_edge(0, 1, key=0)
    g.add_edge(0, 1, key=1)
    g.add_edge(0, 1, key=2)
    h = nx.MultiGraph()
    h.add_edge(0, 1, key=0)
    h.add_edge(0, 1, key=3)
    gh = nx.symmetric_difference(g, h)
    assert_equal(set(gh.nodes()), set(g.nodes()))
    assert_equal(set(gh.nodes()), set(h.nodes()))
    assert_equal(sorted(gh.edges()), 3 * [(0, 1)])
    assert_equal(sorted(sorted(e) for e in gh.edges(keys=True)),
                 [[0, 1, 1], [0, 1, 2], [0, 1, 3]])


@raises(nx.NetworkXError)
def test_symmetric_difference_raise():
    G = nx.path_graph(4)
    H = nx.path_graph(3)
    GH = nx.symmetric_difference(G, H)


def test_union_and_compose():
    K3 = complete_graph(3)
    P3 = path_graph(3)

    G1 = nx.DiGraph()
    G1.add_edge('A', 'B')
    G1.add_edge('A', 'C')
    G1.add_edge('A', 'D')
    G2 = nx.DiGraph()
    G2.add_edge('1', '2')
    G2.add_edge('1', '3')
    G2.add_edge('1', '4')

    G = union(G1, G2)
    H = compose(G1, G2)
    assert_edges_equal(G.edges(), H.edges())
    assert_false(G.has_edge('A', 1))
    assert_raises(nx.NetworkXError, nx.union, K3, P3)
    H1 = union(H, G1, rename=('H', 'G1'))
    assert_equal(sorted(H1.nodes()),
                 ['G1A', 'G1B', 'G1C', 'G1D',
                  'H1', 'H2', 'H3', 'H4', 'HA', 'HB', 'HC', 'HD'])

    H2 = union(H, G2, rename=("H", ""))
    assert_equal(sorted(H2.nodes()),
                 ['1', '2', '3', '4',
                  'H1', 'H2', 'H3', 'H4', 'HA', 'HB', 'HC', 'HD'])

    assert_false(H1.has_edge('NB', 'NA'))

    G = compose(G, G)
    assert_edges_equal(G.edges(), H.edges())

    G2 = union(G2, G2, rename=('', 'copy'))
    assert_equal(sorted(G2.nodes()),
                 ['1', '2', '3', '4', 'copy1', 'copy2', 'copy3', 'copy4'])

    assert_equal(sorted(G2.neighbors('copy4')), [])
    assert_equal(sorted(G2.neighbors('copy1')), ['copy2', 'copy3', 'copy4'])
    assert_equal(len(G), 8)
    assert_equal(number_of_edges(G), 6)

    E = disjoint_union(G, G)
    assert_equal(len(E), 16)
    assert_equal(number_of_edges(E), 12)

    E = disjoint_union(G1, G2)
    assert_equal(sorted(E.nodes()), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

    G = nx.Graph()
    H = nx.Graph()
    G.add_nodes_from([(1, {'a1': 1})])
    H.add_nodes_from([(1, {'b1': 1})])
    R = compose(G, H)
    assert_equal(R.nodes, {1: {'a1': 1, 'b1': 1}})


def test_union_multigraph():
    G = nx.MultiGraph()
    G.add_edge(1, 2, key=0)
    G.add_edge(1, 2, key=1)
    H = nx.MultiGraph()
    H.add_edge(3, 4, key=0)
    H.add_edge(3, 4, key=1)
    GH = nx.union(G, H)
    assert_equal(set(GH), set(G) | set(H))
    assert_equal(set(GH.edges(keys=True)),
                 set(G.edges(keys=True)) | set(H.edges(keys=True)))


def test_disjoint_union_multigraph():
    G = nx.MultiGraph()
    G.add_edge(0, 1, key=0)
    G.add_edge(0, 1, key=1)
    H = nx.MultiGraph()
    H.add_edge(2, 3, key=0)
    H.add_edge(2, 3, key=1)
    GH = nx.disjoint_union(G, H)
    assert_equal(set(GH), set(G) | set(H))
    assert_equal(set(GH.edges(keys=True)),
                 set(G.edges(keys=True)) | set(H.edges(keys=True)))


def test_compose_multigraph():
    G = nx.MultiGraph()
    G.add_edge(1, 2, key=0)
    G.add_edge(1, 2, key=1)
    H = nx.MultiGraph()
    H.add_edge(3, 4, key=0)
    H.add_edge(3, 4, key=1)
    GH = nx.compose(G, H)
    assert_equal(set(GH), set(G) | set(H))
    assert_equal(set(GH.edges(keys=True)),
                 set(G.edges(keys=True)) | set(H.edges(keys=True)))
    H.add_edge(1, 2, key=2)
    GH = nx.compose(G, H)
    assert_equal(set(GH), set(G) | set(H))
    assert_equal(set(GH.edges(keys=True)),
                 set(G.edges(keys=True)) | set(H.edges(keys=True)))


@raises(nx.NetworkXError)
def test_mixed_type_union():
    G = nx.Graph()
    H = nx.MultiGraph()
    U = nx.union(G, H)


@raises(nx.NetworkXError)
def test_mixed_type_disjoint_union():
    G = nx.Graph()
    H = nx.MultiGraph()
    U = nx.disjoint_union(G, H)


@raises(nx.NetworkXError)
def test_mixed_type_intersection():
    G = nx.Graph()
    H = nx.MultiGraph()
    U = nx.intersection(G, H)


@raises(nx.NetworkXError)
def test_mixed_type_difference():
    G = nx.Graph()
    H = nx.MultiGraph()
    U = nx.difference(G, H)


@raises(nx.NetworkXError)
def test_mixed_type_symmetric_difference():
    G = nx.Graph()
    H = nx.MultiGraph()
    U = nx.symmetric_difference(G, H)


@raises(nx.NetworkXError)
def test_mixed_type_compose():
    G = nx.Graph()
    H = nx.MultiGraph()
    U = nx.compose(G, H)
