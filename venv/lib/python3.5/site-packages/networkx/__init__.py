"""
NetworkX
========

NetworkX is a Python package for the creation, manipulation,
and study of the structure, dynamics, and functions
of complex networks.

Website (including documentation)::

    http://networkx.github.io

Mailing list::

    https://groups.google.com/forum/#!forum/networkx-discuss

Source::

    https://github.com/networkx/networkx

Bug reports::

    https://github.com/networkx/networkx/issues

Simple example
--------------

Find the shortest path between two nodes in an undirected graph::

    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edge('A', 'B', weight=4)
    >>> G.add_edge('B', 'D', weight=2)
    >>> G.add_edge('A', 'C', weight=3)
    >>> G.add_edge('C', 'D', weight=4)
    >>> nx.shortest_path(G, 'A', 'D', weight='weight')
    ['A', 'B', 'D']

Bugs
----

Please report any bugs that you find `here <https://github.com/networkx/networkx/issues>`_.
Or, even better, fork the repository on GitHub and create a pull request (PR).

License
-------

Released under the 3-Clause BSD license::

   Copyright (C) 2004-2018 NetworkX Developers
   Aric Hagberg <hagberg@lanl.gov>
   Dan Schult <dschult@colgate.edu>
   Pieter Swart <swart@lanl.gov>
"""
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Add platform dependent shared library path to sys.path
#

from __future__ import absolute_import

import sys
if sys.version_info[:2] < (2, 7):
    m = "Python 2.7 or later is required for NetworkX (%d.%d detected)."
    raise ImportError(m % sys.version_info[:2])
del sys

# Release data
from networkx import release

__author__ = '%s <%s>\n%s <%s>\n%s <%s>' % \
    (release.authors['Hagberg'] + release.authors['Schult'] +
        release.authors['Swart'])
__license__ = release.license

__date__ = release.date
__version__ = release.version

__bibtex__ = """@inproceedings{hagberg-2008-exploring,
author = {Aric A. Hagberg and Daniel A. Schult and Pieter J. Swart},
title = {Exploring network structure, dynamics, and function using {NetworkX}},
year = {2008},
month = Aug,
urlpdf = {http://math.lanl.gov/~hagberg/Papers/hagberg-2008-exploring.pdf},
booktitle = {Proceedings of the 7th Python in Science Conference (SciPy2008)},
editors = {G\"{a}el Varoquaux, Travis Vaught, and Jarrod Millman},
address = {Pasadena, CA USA},
pages = {11--15}
}"""

# These are import orderwise
from networkx.exception import *
import networkx.utils

import networkx.classes.filters
import networkx.classes
from networkx.classes import *

import networkx.convert
from networkx.convert import *

import networkx.convert_matrix
from networkx.convert_matrix import *


import networkx.relabel
from networkx.relabel import *

import networkx.generators
from networkx.generators import *

import networkx.readwrite
from networkx.readwrite import *

# Need to test with SciPy, when available
import networkx.algorithms
from networkx.algorithms import *
import networkx.linalg

from networkx.linalg import *
from networkx.tests.test import run as test

import networkx.drawing
from networkx.drawing import *
