from nose.tools import *
import networkx as nx


def test_2d_grid_graph():
    # FC article claims 2d grid graph of size n is (3,3)-connected
    # and (5,9)-connected, but I don't think it is (5,9)-connected
    G = nx.grid_2d_graph(8, 8, periodic=True)
    assert_true(nx.is_kl_connected(G, 3, 3))
    assert_false(nx.is_kl_connected(G, 5, 9))
    (H, graphOK) = nx.kl_connected_subgraph(G, 5, 9, same_as_graph=True)
    assert_false(graphOK)


def test_small_graph():
    G = nx.Graph()
    G.add_edge(1, 2)
    G.add_edge(1, 3)
    G.add_edge(2, 3)
    assert_true(nx.is_kl_connected(G, 2, 2))
    H = nx.kl_connected_subgraph(G, 2, 2)
    (H, graphOK) = nx.kl_connected_subgraph(G, 2, 2,
                                            low_memory=True,
                                            same_as_graph=True)
    assert_true(graphOK)
