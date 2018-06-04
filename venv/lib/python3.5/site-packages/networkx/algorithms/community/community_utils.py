# -*- coding: utf-8 -*-
#
# utils.py - helper functions for community-finding algorithms
#
# Copyright 2011 Ben Edwards <bedwards@cs.unm.edu>.
# Copyright 2011 Aric Hagberg <hagberg@lanl.gov>.
# Copyright 2015 NetworkX developers.
#
# This file is part of NetworkX.
#
# NetworkX is distributed under a BSD license; see LICENSE.txt for more
# information.
"""Helper functions for community-finding algorithms."""

__all__ = ['is_partition']


def is_partition(G, communities):
    """Return True if and only if `communities` is a partition of
    the nodes of `G`.

    A partition of a universe set is a family of pairwise disjoint sets
    whose union is the entire universe set.

    `G` is a NetworkX graph.

    `communities` is an iterable of sets of nodes of `G`. This
    iterable will be consumed multiple times during the execution of
    this function.

    """
    # Alternate implementation:
    #
    #     return (len(G) == sum(len(c) for c in community) and
    #             set(G) == set.union(*community))
    #
    return all(sum(1 if v in c else 0 for c in communities) == 1 for v in G)
