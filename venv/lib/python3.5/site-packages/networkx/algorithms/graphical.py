# -*- coding: utf-8 -*-
"""Test sequences for graphiness.
"""
#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
import heapq
import networkx as nx
__author__ = "\n".join(['Aric Hagberg (hagberg@lanl.gov)',
                        'Pieter Swart (swart@lanl.gov)',
                        'Dan Schult (dschult@colgate.edu)'
                        'Joel Miller (joel.c.miller.research@gmail.com)'
                        'Ben Edwards'
                        'Brian Cloteaux <brian.cloteaux@nist.gov>'])

__all__ = ['is_graphical',
           'is_multigraphical',
           'is_pseudographical',
           'is_digraphical',
           'is_valid_degree_sequence_erdos_gallai',
           'is_valid_degree_sequence_havel_hakimi',
           ]


def is_graphical(sequence, method='eg'):
    """Returns True if sequence is a valid degree sequence.

    A degree sequence is valid if some graph can realize it.

    Parameters
    ----------
    sequence : list or iterable container
        A sequence of integer node degrees

    method : "eg" | "hh"
        The method used to validate the degree sequence.
        "eg" corresponds to the Erdős-Gallai algorithm, and
        "hh" to the Havel-Hakimi algorithm.

    Returns
    -------
    valid : bool
        True if the sequence is a valid degree sequence and False if not.

    Examples
    --------
    >>> G = nx.path_graph(4)
    >>> sequence = (d for n, d in G.degree())
    >>> nx.is_graphical(sequence)
    True

    References
    ----------
    Erdős-Gallai
        [EG1960]_, [choudum1986]_

    Havel-Hakimi
        [havel1955]_, [hakimi1962]_, [CL1996]_
    """
    if method == 'eg':
        valid = is_valid_degree_sequence_erdos_gallai(list(sequence))
    elif method == 'hh':
        valid = is_valid_degree_sequence_havel_hakimi(list(sequence))
    else:
        msg = "`method` must be 'eg' or 'hh'"
        raise nx.NetworkXException(msg)
    return valid


def _basic_graphical_tests(deg_sequence):
    # Sort and perform some simple tests on the sequence
    if not nx.utils.is_list_of_ints(deg_sequence):
        raise nx.NetworkXUnfeasible
    p = len(deg_sequence)
    num_degs = [0] * p
    dmax, dmin, dsum, n = 0, p, 0, 0
    for d in deg_sequence:
        # Reject if degree is negative or larger than the sequence length
        if d < 0 or d >= p:
            raise nx.NetworkXUnfeasible
        # Process only the non-zero integers
        elif d > 0:
            dmax, dmin, dsum, n = max(dmax, d), min(dmin, d), dsum + d, n + 1
            num_degs[d] += 1
    # Reject sequence if it has odd sum or is oversaturated
    if dsum % 2 or dsum > n * (n - 1):
        raise nx.NetworkXUnfeasible
    return dmax, dmin, dsum, n, num_degs


def is_valid_degree_sequence_havel_hakimi(deg_sequence):
    r"""Returns True if deg_sequence can be realized by a simple graph.

    The validation proceeds using the Havel-Hakimi theorem.
    Worst-case run time is $O(s)$ where $s$ is the sum of the sequence.

    Parameters
    ----------
    deg_sequence : list
        A list of integers where each element specifies the degree of a node
        in a graph.

    Returns
    -------
    valid : bool
        True if deg_sequence is graphical and False if not.

    Notes
    -----
    The ZZ condition says that for the sequence d if

    .. math::
        |d| >= \frac{(\max(d) + \min(d) + 1)^2}{4*\min(d)}

    then d is graphical.  This was shown in Theorem 6 in [1]_.

    References
    ----------
    .. [1] I.E. Zverovich and V.E. Zverovich. "Contributions to the theory
       of graphic sequences", Discrete Mathematics, 105, pp. 292-303 (1992).

    [havel1955]_, [hakimi1962]_, [CL1996]_

    """
    try:
        dmax, dmin, dsum, n, num_degs = _basic_graphical_tests(deg_sequence)
    except nx.NetworkXUnfeasible:
        return False
    # Accept if sequence has no non-zero degrees or passes the ZZ condition
    if n == 0 or 4 * dmin * n >= (dmax + dmin + 1) * (dmax + dmin + 1):
        return True

    modstubs = [0] * (dmax + 1)
    # Successively reduce degree sequence by removing the maximum degree
    while n > 0:
        # Retrieve the maximum degree in the sequence
        while num_degs[dmax] == 0:
            dmax -= 1
        # If there are not enough stubs to connect to, then the sequence is
        # not graphical
        if dmax > n - 1:
            return False

        # Remove largest stub in list
        num_degs[dmax], n = num_degs[dmax] - 1, n - 1
        # Reduce the next dmax largest stubs
        mslen = 0
        k = dmax
        for i in range(dmax):
            while num_degs[k] == 0:
                k -= 1
            num_degs[k], n = num_degs[k] - 1, n - 1
            if k > 1:
                modstubs[mslen] = k - 1
                mslen += 1
        # Add back to the list any non-zero stubs that were removed
        for i in range(mslen):
            stub = modstubs[i]
            num_degs[stub], n = num_degs[stub] + 1, n + 1
    return True


def is_valid_degree_sequence_erdos_gallai(deg_sequence):
    r"""Returns True if deg_sequence can be realized by a simple graph.

    The validation is done using the Erdős-Gallai theorem [EG1960]_.

    Parameters
    ----------
    deg_sequence : list
        A list of integers

    Returns
    -------
    valid : bool
        True if deg_sequence is graphical and False if not.

    Notes
    -----

    This implementation uses an equivalent form of the Erdős-Gallai criterion.
    Worst-case run time is $O(n)$ where $n$ is the length of the sequence.

    Specifically, a sequence d is graphical if and only if the
    sum of the sequence is even and for all strong indices k in the sequence,

     .. math::

       \sum_{i=1}^{k} d_i \leq k(k-1) + \sum_{j=k+1}^{n} \min(d_i,k)
             = k(n-1) - ( k \sum_{j=0}^{k-1} n_j - \sum_{j=0}^{k-1} j n_j )

    A strong index k is any index where d_k >= k and the value n_j is the
    number of occurrences of j in d.  The maximal strong index is called the
    Durfee index.

    This particular rearrangement comes from the proof of Theorem 3 in [2]_.

    The ZZ condition says that for the sequence d if

    .. math::
        |d| >= \frac{(\max(d) + \min(d) + 1)^2}{4*\min(d)}

    then d is graphical.  This was shown in Theorem 6 in [2]_.

    References
    ----------
    .. [1] A. Tripathi and S. Vijay. "A note on a theorem of Erdős & Gallai",
       Discrete Mathematics, 265, pp. 417-420 (2003).
    .. [2] I.E. Zverovich and V.E. Zverovich. "Contributions to the theory
       of graphic sequences", Discrete Mathematics, 105, pp. 292-303 (1992).

    [EG1960]_, [choudum1986]_
    """
    try:
        dmax, dmin, dsum, n, num_degs = _basic_graphical_tests(deg_sequence)
    except nx.NetworkXUnfeasible:
        return False
    # Accept if sequence has no non-zero degrees or passes the ZZ condition
    if n == 0 or 4 * dmin * n >= (dmax + dmin + 1) * (dmax + dmin + 1):
        return True

    # Perform the EG checks using the reformulation of Zverovich and Zverovich
    k, sum_deg, sum_nj, sum_jnj = 0, 0, 0, 0
    for dk in range(dmax, dmin - 1, -1):
        if dk < k + 1:            # Check if already past Durfee index
            return True
        if num_degs[dk] > 0:
            run_size = num_degs[dk]  # Process a run of identical-valued degrees
            if dk < k + run_size:     # Check if end of run is past Durfee index
                run_size = dk - k     # Adjust back to Durfee index
            sum_deg += run_size * dk
            for v in range(run_size):
                sum_nj += num_degs[k + v]
                sum_jnj += (k + v) * num_degs[k + v]
            k += run_size
            if sum_deg > k * (n - 1) - k * sum_nj + sum_jnj:
                return False
    return True


def is_multigraphical(sequence):
    """Returns True if some multigraph can realize the sequence.

    Parameters
    ----------
    deg_sequence : list
        A list of integers

    Returns
    -------
    valid : bool
        True if deg_sequence is a multigraphic degree sequence and False if not.

    Notes
    -----
    The worst-case run time is $O(n)$ where $n$ is the length of the sequence.

    References
    ----------
    .. [1] S. L. Hakimi. "On the realizability of a set of integers as
       degrees of the vertices of a linear graph", J. SIAM, 10, pp. 496-506
       (1962).
    """
    deg_sequence = list(sequence)
    if not nx.utils.is_list_of_ints(deg_sequence):
        return False
    dsum, dmax = 0, 0
    for d in deg_sequence:
        if d < 0:
            return False
        dsum, dmax = dsum + d, max(dmax, d)
    if dsum % 2 or dsum < 2 * dmax:
        return False
    return True


def is_pseudographical(sequence):
    """Returns True if some pseudograph can realize the sequence.

    Every nonnegative integer sequence with an even sum is pseudographical
    (see [1]_).

    Parameters
    ----------
    sequence : list or iterable container
        A sequence of integer node degrees

    Returns
    -------
    valid : bool
      True if the sequence is a pseudographic degree sequence and False if not.

    Notes
    -----
    The worst-case run time is $O(n)$ where n is the length of the sequence.

    References
    ----------
    .. [1] F. Boesch and F. Harary. "Line removal algorithms for graphs
       and their degree lists", IEEE Trans. Circuits and Systems, CAS-23(12),
       pp. 778-782 (1976).
    """
    s = list(sequence)
    if not nx.utils.is_list_of_ints(s):
        return False
    return sum(s) % 2 == 0 and min(s) >= 0


def is_digraphical(in_sequence, out_sequence):
    r"""Returns True if some directed graph can realize the in- and out-degree
    sequences.

    Parameters
    ----------
    in_sequence : list or iterable container
        A sequence of integer node in-degrees

    out_sequence : list or iterable container
        A sequence of integer node out-degrees

    Returns
    -------
    valid : bool
      True if in and out-sequences are digraphic False if not.

    Notes
    -----
    This algorithm is from Kleitman and Wang [1]_.
    The worst case runtime is $O(s \times \log n)$ where $s$ and $n$ are the
    sum and length of the sequences respectively.

    References
    ----------
    .. [1] D.J. Kleitman and D.L. Wang
       Algorithms for Constructing Graphs and Digraphs with Given Valences
       and Factors, Discrete Mathematics, 6(1), pp. 79-88 (1973)
    """
    in_deg_sequence = list(in_sequence)
    out_deg_sequence = list(out_sequence)
    if not nx.utils.is_list_of_ints(in_deg_sequence):
        return False
    if not nx.utils.is_list_of_ints(out_deg_sequence):
        return False
    # Process the sequences and form two heaps to store degree pairs with
    # either zero or non-zero out degrees
    sumin, sumout, nin, nout = 0, 0, len(in_deg_sequence), len(out_deg_sequence)
    maxn = max(nin, nout)
    maxin = 0
    if maxn == 0:
        return True
    stubheap, zeroheap = [], []
    for n in range(maxn):
        in_deg, out_deg = 0, 0
        if n < nout:
            out_deg = out_deg_sequence[n]
        if n < nin:
            in_deg = in_deg_sequence[n]
        if in_deg < 0 or out_deg < 0:
            return False
        sumin, sumout, maxin = sumin + in_deg, sumout + out_deg, max(maxin, in_deg)
        if in_deg > 0:
            stubheap.append((-1 * out_deg, -1 * in_deg))
        elif out_deg > 0:
            zeroheap.append(-1 * out_deg)
    if sumin != sumout:
        return False
    heapq.heapify(stubheap)
    heapq.heapify(zeroheap)

    modstubs = [(0, 0)] * (maxin + 1)
    # Successively reduce degree sequence by removing the maximum out degree
    while stubheap:
        # Take the first value in the sequence with non-zero in degree
        (freeout, freein) = heapq.heappop(stubheap)
        freein *= -1
        if freein > len(stubheap) + len(zeroheap):
            return False

        # Attach out stubs to the nodes with the most in stubs
        mslen = 0
        for i in range(freein):
            if zeroheap and (not stubheap or stubheap[0][0] > zeroheap[0]):
                stubout = heapq.heappop(zeroheap)
                stubin = 0
            else:
                (stubout, stubin) = heapq.heappop(stubheap)
            if stubout == 0:
                return False
            # Check if target is now totally connected
            if stubout + 1 < 0 or stubin < 0:
                modstubs[mslen] = (stubout + 1, stubin)
                mslen += 1

        # Add back the nodes to the heap that still have available stubs
        for i in range(mslen):
            stub = modstubs[i]
            if stub[1] < 0:
                heapq.heappush(stubheap, stub)
            else:
                heapq.heappush(zeroheap, stub[0])
        if freeout < 0:
            heapq.heappush(zeroheap, freeout)
    return True
