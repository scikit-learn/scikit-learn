from nose.tools import raises
import networkx as nx

# smoke tests for exceptions


@raises(nx.NetworkXException)
def test_raises_networkxexception():
    raise nx.NetworkXException


@raises(nx.NetworkXError)
def test_raises_networkxerr():
    raise nx.NetworkXError


@raises(nx.NetworkXPointlessConcept)
def test_raises_networkx_pointless_concept():
    raise nx.NetworkXPointlessConcept


@raises(nx.NetworkXAlgorithmError)
def test_raises_networkxalgorithmerr():
    raise nx.NetworkXAlgorithmError


@raises(nx.NetworkXUnfeasible)
def test_raises_networkx_unfeasible():
    raise nx.NetworkXUnfeasible


@raises(nx.NetworkXNoPath)
def test_raises_networkx_no_path():
    raise nx.NetworkXNoPath


@raises(nx.NetworkXUnbounded)
def test_raises_networkx_unbounded():
    raise nx.NetworkXUnbounded
