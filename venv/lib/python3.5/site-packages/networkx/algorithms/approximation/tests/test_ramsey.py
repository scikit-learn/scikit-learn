from nose.tools import *
import networkx as nx
import networkx.algorithms.approximation as apxa


def test_ramsey():
    # this should only find the complete graph
    graph = nx.complete_graph(10)
    c, i = apxa.ramsey_R2(graph)
    cdens = nx.density(graph.subgraph(c))
    eq_(cdens, 1.0, "clique not found by ramsey!")
    idens = nx.density(graph.subgraph(i))
    eq_(idens, 0.0, "i-set not found by ramsey!")

    # this trival graph has no cliques. should just find i-sets
    graph = nx.trivial_graph(nx.Graph())
    c, i = apxa.ramsey_R2(graph)
    cdens = nx.density(graph.subgraph(c))
    eq_(cdens, 0.0, "clique not found by ramsey!")
    idens = nx.density(graph.subgraph(i))
    eq_(idens, 0.0, "i-set not found by ramsey!")

    graph = nx.barbell_graph(10, 5, nx.Graph())
    c, i = apxa.ramsey_R2(graph)
    cdens = nx.density(graph.subgraph(c))
    eq_(cdens, 1.0, "clique not found by ramsey!")
    idens = nx.density(graph.subgraph(i))
    eq_(idens, 0.0, "i-set not found by ramsey!")
