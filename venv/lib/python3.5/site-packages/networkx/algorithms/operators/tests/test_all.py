from nose.tools import *
import networkx as nx
from networkx.testing import *


def test_union_all_attributes():
    g = nx.Graph()
    g.add_node(0, x=4)
    g.add_node(1, x=5)
    g.add_edge(0, 1, size=5)
    g.graph['name'] = 'g'

    h = g.copy()
    h.graph['name'] = 'h'
    h.graph['attr'] = 'attr'
    h.nodes[0]['x'] = 7

    j = g.copy()
    j.graph['name'] = 'j'
    j.graph['attr'] = 'attr'
    j.nodes[0]['x'] = 7

    ghj = nx.union_all([g, h, j], rename=('g', 'h', 'j'))
    assert_equal(set(ghj.nodes()), set(['h0', 'h1', 'g0', 'g1', 'j0', 'j1']))
    for n in ghj:
        graph, node = n
        assert_equal(ghj.nodes[n], eval(graph).nodes[int(node)])

    assert_equal(ghj.graph['attr'], 'attr')
    assert_equal(ghj.graph['name'], 'j')  # j graph attributes take precendent


def test_intersection_all():
    G = nx.Graph()
    H = nx.Graph()
    R = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4])
    G.add_edge(1, 2)
    G.add_edge(2, 3)
    H.add_nodes_from([1, 2, 3, 4])
    H.add_edge(2, 3)
    H.add_edge(3, 4)
    R.add_nodes_from([1, 2, 3, 4])
    R.add_edge(2, 3)
    R.add_edge(4, 1)
    I = nx.intersection_all([G, H, R])
    assert_equal(set(I.nodes()), set([1, 2, 3, 4]))
    assert_equal(sorted(I.edges()), [(2, 3)])


def test_intersection_all_attributes():
    g = nx.Graph()
    g.add_node(0, x=4)
    g.add_node(1, x=5)
    g.add_edge(0, 1, size=5)
    g.graph['name'] = 'g'

    h = g.copy()
    h.graph['name'] = 'h'
    h.graph['attr'] = 'attr'
    h.nodes[0]['x'] = 7

    gh = nx.intersection_all([g, h])
    assert_equal(set(gh.nodes()), set(g.nodes()))
    assert_equal(set(gh.nodes()), set(h.nodes()))
    assert_equal(sorted(gh.edges()), sorted(g.edges()))

    h.remove_node(0)
    assert_raises(nx.NetworkXError, nx.intersection, g, h)


def test_intersection_all_multigraph_attributes():
    g = nx.MultiGraph()
    g.add_edge(0, 1, key=0)
    g.add_edge(0, 1, key=1)
    g.add_edge(0, 1, key=2)
    h = nx.MultiGraph()
    h.add_edge(0, 1, key=0)
    h.add_edge(0, 1, key=3)
    gh = nx.intersection_all([g, h])
    assert_equal(set(gh.nodes()), set(g.nodes()))
    assert_equal(set(gh.nodes()), set(h.nodes()))
    assert_equal(sorted(gh.edges()), [(0, 1)])
    assert_equal(sorted(gh.edges(keys=True)), [(0, 1, 0)])


def test_union_all_and_compose_all():
    K3 = nx.complete_graph(3)
    P3 = nx.path_graph(3)

    G1 = nx.DiGraph()
    G1.add_edge('A', 'B')
    G1.add_edge('A', 'C')
    G1.add_edge('A', 'D')
    G2 = nx.DiGraph()
    G2.add_edge('1', '2')
    G2.add_edge('1', '3')
    G2.add_edge('1', '4')

    G = nx.union_all([G1, G2])
    H = nx.compose_all([G1, G2])
    assert_edges_equal(G.edges(), H.edges())
    assert_false(G.has_edge('A', '1'))
    assert_raises(nx.NetworkXError, nx.union, K3, P3)
    H1 = nx.union_all([H, G1], rename=('H', 'G1'))
    assert_equal(sorted(H1.nodes()),
                 ['G1A', 'G1B', 'G1C', 'G1D',
                  'H1', 'H2', 'H3', 'H4', 'HA', 'HB', 'HC', 'HD'])

    H2 = nx.union_all([H, G2], rename=("H", ""))
    assert_equal(sorted(H2.nodes()),
                 ['1', '2', '3', '4',
                  'H1', 'H2', 'H3', 'H4', 'HA', 'HB', 'HC', 'HD'])

    assert_false(H1.has_edge('NB', 'NA'))

    G = nx.compose_all([G, G])
    assert_edges_equal(G.edges(), H.edges())

    G2 = nx.union_all([G2, G2], rename=('', 'copy'))
    assert_equal(sorted(G2.nodes()),
                 ['1', '2', '3', '4', 'copy1', 'copy2', 'copy3', 'copy4'])

    assert_equal(sorted(G2.neighbors('copy4')), [])
    assert_equal(sorted(G2.neighbors('copy1')), ['copy2', 'copy3', 'copy4'])
    assert_equal(len(G), 8)
    assert_equal(nx.number_of_edges(G), 6)

    E = nx.disjoint_union_all([G, G])
    assert_equal(len(E), 16)
    assert_equal(nx.number_of_edges(E), 12)

    E = nx.disjoint_union_all([G1, G2])
    assert_equal(sorted(E.nodes()), [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

    G1 = nx.DiGraph()
    G1.add_edge('A', 'B')
    G2 = nx.DiGraph()
    G2.add_edge(1, 2)
    G3 = nx.DiGraph()
    G3.add_edge(11, 22)
    G4 = nx.union_all([G1, G2, G3], rename=("G1", "G2", "G3"))
    assert_equal(sorted(G4.nodes()),
                 ['G1A', 'G1B', 'G21', 'G22',
                  'G311', 'G322'])


def test_union_all_multigraph():
    G = nx.MultiGraph()
    G.add_edge(1, 2, key=0)
    G.add_edge(1, 2, key=1)
    H = nx.MultiGraph()
    H.add_edge(3, 4, key=0)
    H.add_edge(3, 4, key=1)
    GH = nx.union_all([G, H])
    assert_equal(set(GH), set(G) | set(H))
    assert_equal(set(GH.edges(keys=True)),
                 set(G.edges(keys=True)) | set(H.edges(keys=True)))


def test_input_output():
    l = [nx.Graph([(1, 2)]), nx.Graph([(3, 4)])]
    U = nx.disjoint_union_all(l)
    assert_equal(len(l), 2)
    C = nx.compose_all(l)
    assert_equal(len(l), 2)
    l = [nx.Graph([(1, 2)]), nx.Graph([(1, 2)])]
    R = nx.intersection_all(l)
    assert_equal(len(l), 2)


@raises(nx.NetworkXError)
def test_mixed_type_union():
    G = nx.Graph()
    H = nx.MultiGraph()
    I = nx.Graph()
    U = nx.union_all([G, H, I])


@raises(nx.NetworkXError)
def test_mixed_type_disjoint_union():
    G = nx.Graph()
    H = nx.MultiGraph()
    I = nx.Graph()
    U = nx.disjoint_union_all([G, H, I])


@raises(nx.NetworkXError)
def test_mixed_type_intersection():
    G = nx.Graph()
    H = nx.MultiGraph()
    I = nx.Graph()
    U = nx.intersection_all([G, H, I])


@raises(nx.NetworkXError)
def test_mixed_type_compose():
    G = nx.Graph()
    H = nx.MultiGraph()
    I = nx.Graph()
    U = nx.compose_all([G, H, I])
