# -*- coding: utf-8 -*-
import random
import networkx as nx
import itertools as it
from networkx.utils import pairwise
from nose.tools import (assert_equal, assert_false, assert_true,
                        assert_greater_equal, assert_less, assert_less_equal,
                        assert_raises)
from networkx.algorithms.connectivity import (
    k_edge_augmentation,
)
from networkx.algorithms.connectivity.edge_augmentation import (
    collapse,
    complement_edges,
    is_locally_k_edge_connected,
    is_k_edge_connected,
    _unpack_available_edges,
)

# This should be set to the largest k for which an efficient algorithm is
# explicitly defined.
MAX_EFFICIENT_K = 2


def tarjan_bridge_graph():
    # graph from tarjan paper
    # RE Tarjan - "A note on finding the bridges of a graph"
    # Information Processing Letters, 1974 - Elsevier
    # doi:10.1016/0020-0190(74)90003-9.
    # define 2-connected components and bridges
    ccs = [(1, 2, 4, 3, 1, 4), (5, 6, 7, 5), (8, 9, 10, 8),
           (17, 18, 16, 15, 17), (11, 12, 14, 13, 11, 14)]
    bridges = [(4, 8), (3, 5), (3, 17)]
    G = nx.Graph(it.chain(*(pairwise(path) for path in ccs + bridges)))
    return G


def test_weight_key():
    G = nx.Graph()
    G.add_nodes_from([
        1, 2, 3, 4, 5, 6, 7, 8, 9])
    G.add_edges_from([(3, 8), (1, 2), (2, 3)])
    impossible = {(3, 6), (3, 9)}
    rng = random.Random(0)
    avail_uv = list(set(complement_edges(G)) - impossible)
    avail = [(u, v, {'cost': rng.random()}) for u, v in avail_uv]

    _augment_and_check(G, k=1)
    _augment_and_check(G, k=1, avail=avail_uv)
    _augment_and_check(G, k=1, avail=avail, weight='cost')

    _check_augmentations(G, avail, weight='cost')


def test_is_locally_k_edge_connected_exceptions():
    assert_raises(nx.NetworkXNotImplemented,
                  is_k_edge_connected,
                  nx.DiGraph(), k=0)
    assert_raises(nx.NetworkXNotImplemented,
                  is_k_edge_connected,
                  nx.MultiGraph(), k=0)
    assert_raises(ValueError, is_k_edge_connected,
                  nx.Graph(), k=0)


def test_is_k_edge_connected():
    G = nx.barbell_graph(10, 0)
    assert_true(is_k_edge_connected(G, k=1))
    assert_false(is_k_edge_connected(G, k=2))

    G = nx.Graph()
    G.add_nodes_from([5, 15])
    assert_false(is_k_edge_connected(G, k=1))
    assert_false(is_k_edge_connected(G, k=2))

    G = nx.complete_graph(5)
    assert_true(is_k_edge_connected(G, k=1))
    assert_true(is_k_edge_connected(G, k=2))
    assert_true(is_k_edge_connected(G, k=3))
    assert_true(is_k_edge_connected(G, k=4))


def test_is_k_edge_connected_exceptions():
    assert_raises(nx.NetworkXNotImplemented,
                  is_locally_k_edge_connected,
                  nx.DiGraph(), 1, 2, k=0)
    assert_raises(nx.NetworkXNotImplemented,
                  is_locally_k_edge_connected,
                  nx.MultiGraph(), 1, 2, k=0)
    assert_raises(ValueError,
                  is_locally_k_edge_connected,
                  nx.Graph(), 1, 2, k=0)


def test_is_locally_k_edge_connected():
    G = nx.barbell_graph(10, 0)
    assert_true(is_locally_k_edge_connected(G, 5, 15, k=1))
    assert_false(is_locally_k_edge_connected(G, 5, 15, k=2))

    G = nx.Graph()
    G.add_nodes_from([5, 15])
    assert_false(is_locally_k_edge_connected(G, 5, 15, k=2))


def test_null_graph():
    G = nx.Graph()
    _check_augmentations(G, max_k=MAX_EFFICIENT_K + 2)


def test_cliques():
    for n in range(1, 10):
        G = nx.complete_graph(n)
        _check_augmentations(G, max_k=MAX_EFFICIENT_K + 2)


def test_clique_and_node():
    for n in range(1, 10):
        G = nx.complete_graph(n)
        G.add_node(n + 1)
        _check_augmentations(G, max_k=MAX_EFFICIENT_K + 2)


def test_point_graph():
    G = nx.Graph()
    G.add_node(1)
    _check_augmentations(G, max_k=MAX_EFFICIENT_K + 2)


def test_edgeless_graph():
    G = nx.Graph()
    G.add_nodes_from([1, 2, 3, 4])
    _check_augmentations(G)


def test_invalid_k():
    G = nx.Graph()
    assert_raises(ValueError, list, k_edge_augmentation(G, k=-1))
    assert_raises(ValueError, list, k_edge_augmentation(G, k=0))


def test_unfeasible():
    G = tarjan_bridge_graph()
    assert_raises(nx.NetworkXUnfeasible, list,
                  k_edge_augmentation(G, k=1, avail=[]))

    assert_raises(nx.NetworkXUnfeasible, list,
                  k_edge_augmentation(G, k=2, avail=[]))

    assert_raises(nx.NetworkXUnfeasible, list,
                  k_edge_augmentation(G, k=2, avail=[(7, 9)]))

    # partial solutions should not error if real solutions are infeasible
    aug_edges = list(k_edge_augmentation(G, k=2, avail=[(7, 9)], partial=True))
    assert_equal(aug_edges, [(7, 9)])

    _check_augmentations(G, avail=[], max_k=MAX_EFFICIENT_K + 2)

    _check_augmentations(G, avail=[(7, 9)], max_k=MAX_EFFICIENT_K + 2)


def test_tarjan():
    G = tarjan_bridge_graph()

    aug_edges = set(_augment_and_check(G, k=2)[0])
    print('aug_edges = {!r}'.format(aug_edges))
    # cant assert edge exactly equality due to non-determenant edge order
    # but we do know the size of the solution must be 3
    assert_equal(len(aug_edges), 3)

    avail = [(9, 7), (8, 5), (2, 10), (6, 13), (11, 18), (1, 17), (2, 3),
             (16, 17), (18, 14), (15, 14)]
    aug_edges = set(_augment_and_check(G, avail=avail, k=2)[0])

    # Can't assert exact length since approximation depnds on the order of a
    # dict traversal.
    assert_less_equal(len(aug_edges), 3 * 2)

    _check_augmentations(G, avail)


def test_configuration():
    seeds = [2718183590, 2470619828, 1694705158, 3001036531, 2401251497]
    for seed in seeds:
        deg_seq = nx.random_powerlaw_tree_sequence(20, seed=seed, tries=5000)
        G = nx.Graph(nx.configuration_model(deg_seq, seed=seed))
        G.remove_edges_from(nx.selfloop_edges(G))
        _check_augmentations(G)


def test_shell():
    # seeds = [2057382236, 3331169846, 1840105863, 476020778, 2247498425]
    seeds = [1840105863]
    for seed in seeds:
        constructor = [(12, 70, 0.8), (15, 40, 0.6)]
        G = nx.random_shell_graph(constructor, seed=seed)
        _check_augmentations(G)


def test_karate():
    G = nx.karate_club_graph()
    _check_augmentations(G)


def test_star():
    G = nx.star_graph(3)
    _check_augmentations(G)

    G = nx.star_graph(5)
    _check_augmentations(G)

    G = nx.star_graph(10)
    _check_augmentations(G)


def test_barbell():
    G = nx.barbell_graph(5, 0)
    _check_augmentations(G)

    G = nx.barbell_graph(5, 2)
    _check_augmentations(G)

    G = nx.barbell_graph(5, 3)
    _check_augmentations(G)

    G = nx.barbell_graph(5, 4)
    _check_augmentations(G)


def test_bridge():
    G = nx.Graph([(2393, 2257), (2393, 2685), (2685, 2257), (1758, 2257)])
    _check_augmentations(G)


def test_gnp_augmentation():
    rng = random.Random(0)
    G = nx.gnp_random_graph(30, 0.005, seed=0)
    # Randomly make edges available
    avail = {(u, v): 1 + rng.random()
             for u, v in complement_edges(G)
             if rng.random() < .25}
    _check_augmentations(G, avail)


def _assert_solution_properties(G, aug_edges, avail_dict=None):
    """ Checks that aug_edges are consistently formatted """
    if avail_dict is not None:
        assert_true(all(e in avail_dict for e in aug_edges),
                    'when avail is specified aug-edges should be in avail')

    unique_aug = set(map(tuple, map(sorted, aug_edges)))
    unique_aug = list(map(tuple, map(sorted, aug_edges)))
    assert_true(len(aug_edges) == len(unique_aug),
                'edges should be unique')

    assert_false(any(u == v for u, v in unique_aug),
                 'should be no self-edges')

    assert_false(any(G.has_edge(u, v) for u, v in unique_aug),
                 'aug edges and G.edges should be disjoint')


def _augment_and_check(G, k, avail=None, weight=None, verbose=False,
                       orig_k=None, max_aug_k=None):
    """
    Does one specific augmentation and checks for properties of the result
    """
    if orig_k is None:
        try:
            orig_k = nx.edge_connectivity(G)
        except nx.NetworkXPointlessConcept:
            orig_k = 0
    info = {}
    try:
        if avail is not None:
            # ensure avail is in dict form
            avail_dict = dict(zip(*_unpack_available_edges(avail,
                                                           weight=weight)))
        else:
            avail_dict = None
        try:
            # Find the augmentation if possible
            generator = nx.k_edge_augmentation(G, k=k, weight=weight,
                                               avail=avail)
            assert_false(isinstance(generator, list),
                         'should always return an iter')
            aug_edges = []
            for edge in generator:
                aug_edges.append(edge)
        except nx.NetworkXUnfeasible:
            infeasible = True
            info['infeasible'] = True
            assert_equal(len(aug_edges), 0,
                         'should not generate anything if unfeasible')

            if avail is None:
                n_nodes = G.number_of_nodes()
                assert_less_equal(n_nodes, k, (
                    'unconstrained cases are only unfeasible if |V| <= k. '
                    'Got |V|={} and k={}'.format(n_nodes, k)
                ))
            else:
                if max_aug_k is None:
                    G_aug_all = G.copy()
                    G_aug_all.add_edges_from(avail_dict.keys())
                    try:
                        max_aug_k = nx.edge_connectivity(G_aug_all)
                    except nx.NetworkXPointlessConcept:
                        max_aug_k = 0

                assert_less(max_aug_k, k, (
                    'avail should only be unfeasible if using all edges '
                    'doesnt acheive k-edge-connectivity'))

            # Test for a partial solution
            partial_edges = list(nx.k_edge_augmentation(
                G, k=k, weight=weight, partial=True, avail=avail))

            info['n_partial_edges'] = len(partial_edges)

            if avail_dict is None:
                assert_equal(set(partial_edges), set(complement_edges(G)), (
                    'unweighted partial solutions should be the complement'))
            elif len(avail_dict) > 0:
                H = G.copy()

                # Find the partial / full augmented connectivity
                H.add_edges_from(partial_edges)
                partial_conn = nx.edge_connectivity(H)

                H.add_edges_from(set(avail_dict.keys()))
                full_conn = nx.edge_connectivity(H)

                # Full connectivity should be no better than our partial
                # solution.
                assert_equal(partial_conn, full_conn,
                             'adding more edges should not increase k-conn')

            # Find the new edge-connectivity after adding the augmenting edges
            aug_edges = partial_edges
        else:
            infeasible = False

        # Find the weight of the augmentation
        num_edges = len(aug_edges)
        if avail is not None:
            total_weight = sum([avail_dict[e] for e in aug_edges])
        else:
            total_weight = num_edges

        info['total_weight'] = total_weight
        info['num_edges'] = num_edges

        # Find the new edge-connectivity after adding the augmenting edges
        G_aug = G.copy()
        G_aug.add_edges_from(aug_edges)
        try:
            aug_k = nx.edge_connectivity(G_aug)
        except nx.NetworkXPointlessConcept:
            aug_k = 0
        info['aug_k'] = aug_k

        # Do checks
        if not infeasible and orig_k < k:
            assert_greater_equal(info['aug_k'], k, (
                'connectivity should increase to k={} or more'.format(k)))

        assert_greater_equal(info['aug_k'], orig_k, (
            'augmenting should never reduce connectivity'))

        _assert_solution_properties(G, aug_edges, avail_dict)

    except Exception:
        info['failed'] = True
        print('edges = {}'.format(list(G.edges())))
        print('nodes = {}'.format(list(G.nodes())))
        print('aug_edges = {}'.format(list(aug_edges)))
        print('info  = {}'.format(info))
        raise
    else:
        if verbose:
            print('info  = {}'.format(info))

    if infeasible:
        aug_edges = None
    return aug_edges, info


def _check_augmentations(G, avail=None, max_k=None, weight=None,
                         verbose=False):
    """ Helper to check weighted/unweighted cases with multiple values of k """
    # Using all available edges, find the maximum edge-connectivity
    try:
        orig_k = nx.edge_connectivity(G)
    except nx.NetworkXPointlessConcept:
        orig_k = 0

    if avail is not None:
        all_aug_edges = _unpack_available_edges(avail, weight=weight)[0]
        G_aug_all = G.copy()
        G_aug_all.add_edges_from(all_aug_edges)
        try:
            max_aug_k = nx.edge_connectivity(G_aug_all)
        except nx.NetworkXPointlessConcept:
            max_aug_k = 0
    else:
        max_aug_k = G.number_of_nodes() - 1

    if max_k is None:
        max_k = min(4, max_aug_k)

    avail_uniform = {e: 1 for e in complement_edges(G)}

    if verbose:
        print('\n=== CHECK_AUGMENTATION ===')
        print('G.number_of_nodes = {!r}'.format(G.number_of_nodes()))
        print('G.number_of_edges = {!r}'.format(G.number_of_edges()))
        print('max_k = {!r}'.format(max_k))
        print('max_aug_k = {!r}'.format(max_aug_k))
        print('orig_k = {!r}'.format(orig_k))

    # check augmentation for multiple values of k
    for k in range(1, max_k + 1):
        if verbose:
            print('---------------')
            print('Checking k = {}'.format(k))

        # Check the unweighted version
        if verbose:
            print('unweighted case')
        aug_edges1, info1 = _augment_and_check(
            G, k=k,  verbose=verbose, orig_k=orig_k)

        # Check that the weighted version with all available edges and uniform
        # weights gives a similar solution to the unweighted case.
        if verbose:
            print('weighted uniform case')
        aug_edges2, info2 = _augment_and_check(
            G, k=k, avail=avail_uniform, verbose=verbose,
            orig_k=orig_k,
            max_aug_k=G.number_of_nodes() - 1)

        # Check the weighted version
        if avail is not None:
            if verbose:
                print('weighted case')
            aug_edges3, info3 = _augment_and_check(
                G, k=k, avail=avail, weight=weight, verbose=verbose,
                max_aug_k=max_aug_k, orig_k=orig_k)

        if aug_edges1 is not None:
            # Check approximation ratios
            if k == 1:
                # when k=1, both solutions should be optimal
                assert_equal(info2['total_weight'], info1['total_weight'])
            if k == 2:
                # when k=2, the weighted version is an approximation
                if orig_k == 0:
                    # the approximation ratio is 3 if G is not connected
                    assert_less_equal(info2['total_weight'],
                                      info1['total_weight'] * 3)
                else:
                    # the approximation ratio is 2 if G is was connected
                    assert_less_equal(info2['total_weight'],
                                      info1['total_weight'] * 2)
                _check_unconstrained_bridge_property(G, info1)


def _check_unconstrained_bridge_property(G, info1):
    # Check Theorem 5 from Eswaran and Tarjan. (1975) Augmentation problems
    import math
    bridge_ccs = list(nx.connectivity.bridge_components(G))
    # condense G into an forest C
    C = collapse(G, bridge_ccs)

    p = len([n for n, d in C.degree() if d == 1])  # leafs
    q = len([n for n, d in C.degree() if d == 0])  # isolated
    if p + q > 1:
        size_target = int(math.ceil(p / 2.0)) + q
        size_aug = info1['num_edges']
        assert_equal(size_aug, size_target, (
            'augmentation size is different from what theory predicts'))
