# test_triads.py - unit tests for the triads module
#
# Copyright 2015 NetworkX developers.
# Copyright 2009 Diederik van Liere <diederik.vanliere@rotman.utoronto.ca>.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Unit tests for the :mod:`networkx.algorithms.triads` module."""
from nose.tools import assert_equal

import networkx as nx


def test_triadic_census():
    """Tests the triadic census function."""
    G = nx.DiGraph()
    G.add_edges_from(['01', '02', '03', '04', '05', '12', '16', '51', '56',
                      '65'])
    expected = {'030T': 2, '120C': 1, '210': 0, '120U': 0, '012': 9, '102': 3,
                '021U': 0, '111U': 0, '003': 8, '030C': 0, '021D': 9, '201': 0,
                '111D': 1, '300': 0, '120D': 0, '021C': 2}
    actual = nx.triadic_census(G)
    assert_equal(expected, actual)
