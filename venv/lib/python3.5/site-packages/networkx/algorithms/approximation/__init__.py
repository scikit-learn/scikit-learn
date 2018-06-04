# __init__.py - package containing heuristics for optimization problems
#
# Copyright 2016-2018 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Approximations of graph properties and Heuristic functions for optimization
problems.

    .. warning:: The approximation submodule is not imported in the top-level
        ``networkx``.

    These functions can be imported with
    ``from networkx.algorithms import approximation``.

"""

from networkx.algorithms.approximation.clustering_coefficient import *
from networkx.algorithms.approximation.clique import *
from networkx.algorithms.approximation.connectivity import *
from networkx.algorithms.approximation.dominating_set import *
from networkx.algorithms.approximation.kcomponents import *
from networkx.algorithms.approximation.independent_set import *
from networkx.algorithms.approximation.matching import *
from networkx.algorithms.approximation.ramsey import *
from networkx.algorithms.approximation.steinertree import *
from networkx.algorithms.approximation.vertex_cover import *
