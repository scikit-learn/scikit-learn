
from nose.tools import *
import networkx as nx


class TestTreeRecognition(object):

    graph = nx.Graph
    multigraph = nx.MultiGraph

    def setUp(self):

        self.T1 = self.graph()

        self.T2 = self.graph()
        self.T2.add_node(1)

        self.T3 = self.graph()
        self.T3.add_nodes_from(range(5))
        edges = [(i, i + 1) for i in range(4)]
        self.T3.add_edges_from(edges)

        self.T5 = self.multigraph()
        self.T5.add_nodes_from(range(5))
        edges = [(i, i + 1) for i in range(4)]
        self.T5.add_edges_from(edges)

        self.T6 = self.graph()
        self.T6.add_nodes_from([6, 7])
        self.T6.add_edge(6, 7)

        self.F1 = nx.compose(self.T6, self.T3)

        self.N4 = self.graph()
        self.N4.add_node(1)
        self.N4.add_edge(1, 1)

        self.N5 = self.graph()
        self.N5.add_nodes_from(range(5))

        self.N6 = self.graph()
        self.N6.add_nodes_from(range(3))
        self.N6.add_edges_from([(0, 1), (1, 2), (2, 0)])

        self.NF1 = nx.compose(self.T6, self.N6)

    @raises(nx.NetworkXPointlessConcept)
    def test_null_tree(self):
        nx.is_tree(self.graph())
        nx.is_tree(self.multigraph())

    @raises(nx.NetworkXPointlessConcept)
    def test_null_forest(self):
        nx.is_forest(self.graph())
        nx.is_forest(self.multigraph())

    def test_is_tree(self):
        assert_true(nx.is_tree(self.T2))
        assert_true(nx.is_tree(self.T3))
        assert_true(nx.is_tree(self.T5))

    def test_is_not_tree(self):
        assert_false(nx.is_tree(self.N4))
        assert_false(nx.is_tree(self.N5))
        assert_false(nx.is_tree(self.N6))

    def test_is_forest(self):
        assert_true(nx.is_forest(self.T2))
        assert_true(nx.is_forest(self.T3))
        assert_true(nx.is_forest(self.T5))
        assert_true(nx.is_forest(self.F1))
        assert_true(nx.is_forest(self.N5))

    def test_is_not_forest(self):
        assert_false(nx.is_forest(self.N4))
        assert_false(nx.is_forest(self.N6))
        assert_false(nx.is_forest(self.NF1))


class TestDirectedTreeRecognition(TestTreeRecognition):
    graph = nx.DiGraph
    multigraph = nx.MultiDiGraph


def test_disconnected_graph():
    # https://github.com/networkx/networkx/issues/1144
    G = nx.Graph()
    G.add_edges_from([(0, 1), (1, 2), (2, 0), (3, 4)])
    assert_false(nx.is_tree(G))

    G = nx.DiGraph()
    G.add_edges_from([(0, 1), (1, 2), (2, 0), (3, 4)])
    assert_false(nx.is_tree(G))


def test_dag_nontree():
    G = nx.DiGraph()
    G.add_edges_from([(0, 1), (0, 2), (1, 2)])
    assert_false(nx.is_tree(G))
    assert_true(nx.is_directed_acyclic_graph(G))


def test_multicycle():
    G = nx.MultiDiGraph()
    G.add_edges_from([(0, 1), (0, 1)])
    assert_false(nx.is_tree(G))
    assert_true(nx.is_directed_acyclic_graph(G))


def test_emptybranch():
    G = nx.DiGraph()
    G.add_nodes_from(range(10))
    assert_true(nx.is_branching(G))
    assert_false(nx.is_arborescence(G))


def test_path():
    G = nx.DiGraph()
    nx.add_path(G, range(5))
    assert_true(nx.is_branching(G))
    assert_true(nx.is_arborescence(G))


def test_notbranching1():
    # Acyclic violation.
    G = nx.MultiDiGraph()
    G.add_nodes_from(range(10))
    G.add_edges_from([(0, 1), (1, 0)])
    assert_false(nx.is_branching(G))
    assert_false(nx.is_arborescence(G))


def test_notbranching2():
    # In-degree violation.
    G = nx.MultiDiGraph()
    G.add_nodes_from(range(10))
    G.add_edges_from([(0, 1), (0, 2), (3, 2)])
    assert_false(nx.is_branching(G))
    assert_false(nx.is_arborescence(G))


def test_notarborescence1():
    # Not an arborescence due to not spanning.
    G = nx.MultiDiGraph()
    G.add_nodes_from(range(10))
    G.add_edges_from([(0, 1), (0, 2), (1, 3), (5, 6)])
    assert_true(nx.is_branching(G))
    assert_false(nx.is_arborescence(G))


def test_notarborescence2():
    # Not an arborescence due to in-degree violation.
    G = nx.MultiDiGraph()
    nx.add_path(G, range(5))
    G.add_edge(6, 4)
    assert_false(nx.is_branching(G))
    assert_false(nx.is_arborescence(G))
