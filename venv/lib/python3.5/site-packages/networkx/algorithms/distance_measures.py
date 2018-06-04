# -*- coding: utf-8 -*-
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: Aric Hagberg (hagberg@lanl.gov)
#          Dan Schult (dschult@colgate.edu)
"""Graph diameter, radius, eccentricity and other properties."""
import networkx

__all__ = ['extrema_bounding', 'eccentricity', 'diameter',
           'radius', 'periphery', 'center']


def extrema_bounding(G, compute="diameter"):
    """Compute requested extreme distance metric of undirected graph G

    Computation is based on smart lower and upper bounds, and in practice
    linear in the number of nodes, rather than quadratic (except for some
    border cases such as complete graphs or circle shaped graphs).

    Parameters
    ----------
    G : NetworkX graph
       An undirected graph

    compute : string denoting the requesting metric
       "diameter" for the maximal eccentricity value,
       "radius" for the minimal eccentricity value,
       "periphery" for the set of nodes with eccentricity equal to the diameter
       "center" for the set of nodes with eccentricity equal to the radius

    Returns
    -------
    value : value of the requested metric
       int for "diameter" and "radius" or
       list of nodes for "center" and "periphery"

    Raises
    ------
    NetworkXError
        If the graph consists of multiple components

    Notes
    -----
    This algorithm was proposed in the following papers:

    F.W. Takes and W.A. Kosters, Determining the Diameter of Small World
    Networks, in Proceedings of the 20th ACM International Conference on
    Information and Knowledge Management (CIKM 2011), pp. 1191-1196, 2011.
    doi: http://dx.doi.org/10.1145/2063576.2063748

    F.W. Takes and W.A. Kosters, Computing the Eccentricity Distribution of
    Large Graphs, Algorithms 6(1): 100-118, 2013.
    doi: http://dx.doi.org/10.3390/a6010100

    M. Borassi, P. Crescenzi, M. Habib, W.A. Kosters, A. Marino and F.W. Takes,
    Fast Graph Diameter and Radius BFS-Based Computation in (Weakly Connected)
    Real-World Graphs, Theoretical Computer Science 586: 59-80, 2015.
    doi: http://dx.doi.org/10.1016/j.tcs.2015.02.033
    """

    # init variables
    degrees = dict(G.degree())  # start with the highest degree node
    minlowernode = max(degrees, key=degrees.get)
    N = len(degrees)  # number of nodes
    # alternate between smallest lower and largest upper bound
    high = False
    # status variables
    ecc_lower = dict.fromkeys(G, 0)
    ecc_upper = dict.fromkeys(G, N)
    candidates = set(G)

    # (re)set bound extremes
    minlower = N
    maxlower = 0
    minupper = N
    maxupper = 0

    # repeat the following until there are no more candidates
    while candidates:
        if high:
            current = maxuppernode  # select node with largest upper bound
        else:
            current = minlowernode  # select node with smallest lower bound
        high = not high

        # get distances from/to current node and derive eccentricity
        dist = dict(networkx.single_source_shortest_path_length(G, current))
        if len(dist) != N:
            msg = ('Cannot compute metric because graph is not connected.')
            raise networkx.NetworkXError(msg)
        current_ecc = max(dist.values())

        # print status update
#        print ("ecc of " + str(current) + " (" + str(ecc_lower[current]) + "/"
#        + str(ecc_upper[current]) + ", deg: " + str(dist[current]) + ") is "
#        + str(current_ecc))
#        print(ecc_upper)

        # (re)set bound extremes
        maxuppernode = None
        minlowernode = None

        # update node bounds
        for i in candidates:
            # update eccentricity bounds
            d = dist[i]
            ecc_lower[i] = low = max(ecc_lower[i], max(d, (current_ecc - d)))
            ecc_upper[i] = upp = min(ecc_upper[i], current_ecc + d)

            # update min/max values of lower and upper bounds
            minlower = min(ecc_lower[i], minlower)
            maxlower = max(ecc_lower[i], maxlower)
            minupper = min(ecc_upper[i], minupper)
            maxupper = max(ecc_upper[i], maxupper)

        # update candidate set
        if compute == 'diameter':
            ruled_out = {i for i in candidates if ecc_upper[i] <= maxlower and
                         2 * ecc_lower[i] >= maxupper}

        elif compute == 'radius':
            ruled_out = {i for i in candidates if ecc_lower[i] >= minupper and
                         ecc_upper[i] + 1 <= 2 * minlower}

        elif compute == 'periphery':
            ruled_out = {i for i in candidates if ecc_upper[i] < maxlower and
                         (maxlower == maxupper or ecc_lower[i] > maxupper)}

        elif compute == 'center':
            ruled_out = {i for i in candidates if ecc_lower[i] > minupper and
                         (minlower == minupper or ecc_upper[i] + 1 < 2 * minlower)}

        elif compute == 'eccentricities':
            ruled_out = {}

        ruled_out.update(i for i in candidates if ecc_lower[i] == ecc_upper[i])
        candidates -= ruled_out

#        for i in ruled_out:
#            print("removing %g: ecc_u: %g maxl: %g ecc_l: %g maxu: %g"%
#                    (i,ecc_upper[i],maxlower,ecc_lower[i],maxupper))
#        print("node %g: ecc_u: %g maxl: %g ecc_l: %g maxu: %g"%
#                    (4,ecc_upper[4],maxlower,ecc_lower[4],maxupper))
#        print("NODE 4: %g"%(ecc_upper[4] <= maxlower))
#        print("NODE 4: %g"%(2 * ecc_lower[4] >= maxupper))
#        print("NODE 4: %g"%(ecc_upper[4] <= maxlower
#                            and 2 * ecc_lower[4] >= maxupper))

        # updating maxuppernode and minlowernode for selection in next round
        for i in candidates:
            if minlowernode is None \
                    or (ecc_lower[i] == ecc_lower[minlowernode]
                        and degrees[i] > degrees[minlowernode]) \
                    or (ecc_lower[i] < ecc_lower[minlowernode]):
                minlowernode = i

            if maxuppernode is None \
                    or (ecc_upper[i] == ecc_upper[maxuppernode]
                        and degrees[i] > degrees[maxuppernode]) \
                    or (ecc_upper[i] > ecc_upper[maxuppernode]):
                maxuppernode = i

        # print status update
#        print (" min=" + str(minlower) + "/" + str(minupper) +
#        " max=" + str(maxlower) + "/" + str(maxupper) +
#        " candidates: " + str(len(candidates)))
#        print("cand:",candidates)
#        print("ecc_l",ecc_lower)
#        print("ecc_u",ecc_upper)
#        wait = input("press Enter to continue")

    # return the correct value of the requested metric
    if compute == 'diameter':
        return maxlower
    elif compute == 'radius':
        return minupper
    elif compute == 'periphery':
        p = [v for v in G if ecc_lower[v] == maxlower]
        return p
    elif compute == 'center':
        c = [v for v in G if ecc_upper[v] == minupper]
        return c
    elif compute == 'eccentricities':
        return ecc_lower
    return None


def eccentricity(G, v=None, sp=None):
    """Return the eccentricity of nodes in G.

    The eccentricity of a node v is the maximum distance from v to
    all other nodes in G.

    Parameters
    ----------
    G : NetworkX graph
       A graph

    v : node, optional
       Return value of specified node

    sp : dict of dicts, optional
       All pairs shortest path lengths as a dictionary of dictionaries

    Returns
    -------
    ecc : dictionary
       A dictionary of eccentricity values keyed by node.
    """
#    if v is None:                # none, use entire graph
#        nodes=G.nodes()
#    elif v in G:               # is v a single node
#        nodes=[v]
#    else:                      # assume v is a container of nodes
#        nodes=v
    order = G.order()

    e = {}
    for n in G.nbunch_iter(v):
        if sp is None:
            length = networkx.single_source_shortest_path_length(G, n)
            L = len(length)
        else:
            try:
                length = sp[n]
                L = len(length)
            except TypeError:
                raise networkx.NetworkXError('Format of "sp" is invalid.')
        if L != order:
            if G.is_directed():
                msg = ('Found infinite path length because the digraph is not'
                       ' strongly connected')
            else:
                msg = ('Found infinite path length because the graph is not'
                       ' connected')
            raise networkx.NetworkXError(msg)

        e[n] = max(length.values())

    if v in G:
        return e[v]  # return single value
    else:
        return e


def diameter(G, e=None, usebounds=False):
    """Return the diameter of the graph G.

    The diameter is the maximum eccentricity.

    Parameters
    ----------
    G : NetworkX graph
       A graph

    e : eccentricity dictionary, optional
      A precomputed dictionary of eccentricities.

    Returns
    -------
    d : integer
       Diameter of graph

    See Also
    --------
    eccentricity
    """
    if usebounds is True and e is None and not G.is_directed():
        return extrema_bounding(G, compute="diameter")
    if e is None:
        e = eccentricity(G)
    return max(e.values())


def periphery(G, e=None, usebounds=False):
    """Return the periphery of the graph G.

    The periphery is the set of nodes with eccentricity equal to the diameter.

    Parameters
    ----------
    G : NetworkX graph
       A graph

    e : eccentricity dictionary, optional
      A precomputed dictionary of eccentricities.

    Returns
    -------
    p : list
       List of nodes in periphery
    """
    if usebounds is True and e is None and not G.is_directed():
        return extrema_bounding(G, compute="periphery")
    if e is None:
        e = eccentricity(G)
    diameter = max(e.values())
    p = [v for v in e if e[v] == diameter]
    return p


def radius(G, e=None, usebounds=False):
    """Return the radius of the graph G.

    The radius is the minimum eccentricity.

    Parameters
    ----------
    G : NetworkX graph
       A graph

    e : eccentricity dictionary, optional
      A precomputed dictionary of eccentricities.

    Returns
    -------
    r : integer
       Radius of graph
    """
    if usebounds is True and e is None and not G.is_directed():
        return extrema_bounding(G, compute="radius")
    if e is None:
        e = eccentricity(G)
    return min(e.values())


def center(G, e=None, usebounds=False):
    """Return the center of the graph G.

    The center is the set of nodes with eccentricity equal to radius.

    Parameters
    ----------
    G : NetworkX graph
       A graph

    e : eccentricity dictionary, optional
      A precomputed dictionary of eccentricities.

    Returns
    -------
    c : list
       List of nodes in center
    """
    if usebounds is True and e is None and not G.is_directed():
        return extrema_bounding(G, compute="center")
    if e is None:
        e = eccentricity(G)
    radius = min(e.values())
    p = [v for v in e if e[v] == radius]
    return p
