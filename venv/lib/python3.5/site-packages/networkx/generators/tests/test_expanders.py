# Copyright 2014 "cheebee7i".
# Copyright 2014 "alexbrc".
# Copyright 2014 Jeffrey Finkelstein <jeffrey.finkelstein@gmail.com>.
"""Unit tests for the :mod:`networkx.generators.expanders` module.

"""
try:
    import scipy
    is_scipy_available = True
except:
    is_scipy_available = False

import networkx as nx
from networkx import adjacency_matrix
from networkx import number_of_nodes
from networkx.generators.expanders import chordal_cycle_graph
from networkx.generators.expanders import margulis_gabber_galil_graph

from nose import SkipTest
from nose.tools import assert_equal
from nose.tools import assert_less
from nose.tools import assert_raises
from nose.tools import assert_true


def test_margulis_gabber_galil_graph():
    try:
        # Scipy is required for conversion to an adjacency matrix.
        # We also use scipy for computing the eigenvalues,
        # but this second use could be done using only numpy.
        import numpy as np
        import scipy.linalg
        has_scipy = True
    except ImportError as e:
        has_scipy = False
    for n in 2, 3, 5, 6, 10:
        g = margulis_gabber_galil_graph(n)
        assert_equal(number_of_nodes(g), n * n)
        for node in g:
            assert_equal(g.degree(node), 8)
            assert_equal(len(node), 2)
            for i in node:
                assert_equal(int(i), i)
                assert_true(0 <= i < n)
        if has_scipy:
            # Eigenvalues are already sorted using the scipy eigvalsh,
            # but the implementation in numpy does not guarantee order.
            w = sorted(scipy.linalg.eigvalsh(adjacency_matrix(g).A))
            assert_less(w[-2], 5 * np.sqrt(2))


def test_chordal_cycle_graph():
    """Test for the :func:`networkx.chordal_cycle_graph` function."""
    if not is_scipy_available:
        raise SkipTest('SciPy is not available')
    primes = [3, 5, 7, 11]
    for p in primes:
        G = chordal_cycle_graph(p)
        assert_equal(len(G), p)
        # TODO The second largest eigenvalue should be smaller than a constant,
        # independent of the number of nodes in the graph:
        #
        #     eigs = sorted(scipy.linalg.eigvalsh(adjacency_matrix(G).A))
        #     assert_less(eigs[-2], ...)
        #


def test_margulis_gabber_galil_graph_badinput():
    assert_raises(nx.NetworkXError, margulis_gabber_galil_graph, 3,
                  nx.DiGraph())
    assert_raises(nx.NetworkXError, margulis_gabber_galil_graph, 3,
                  nx.Graph())
