from itertools import combinations

from nose.tools import assert_equal
from nose.tools import assert_false
from nose.tools import assert_in
from nose.tools import assert_raises
from nose.tools import assert_true
from nose.tools import raises
from nose.tools import ok_

import networkx as nx
from networkx.testing.utils import assert_edges_equal
from networkx.utils import arbitrary_element
from networkx.utils import consume
from networkx.utils import pairwise


class TestDagLongestPath(object):
    """Unit tests computing the longest path in a directed acyclic graph."""

    def test_empty(self):
        G = nx.DiGraph()
        assert_equal(nx.dag_longest_path(G), [])

    def test_unweighted1(self):
        edges = [(1, 2), (2, 3), (2, 4), (3, 5), (5, 6), (3, 7)]
        G = nx.DiGraph(edges)
        assert_equal(nx.dag_longest_path(G), [1, 2, 3, 5, 6])

    def test_unweighted2(self):
        edges = [(1, 2), (2, 3), (3, 4), (4, 5), (1, 3), (1, 5), (3, 5)]
        G = nx.DiGraph(edges)
        assert_equal(nx.dag_longest_path(G), [1, 2, 3, 4, 5])

    def test_weighted(self):
        G = nx.DiGraph()
        edges = [(1, 2, -5), (2, 3, 1), (3, 4, 1), (4, 5, 0), (3, 5, 4),
                 (1, 6, 2)]
        G.add_weighted_edges_from(edges)
        assert_equal(nx.dag_longest_path(G), [2, 3, 5])

    def test_undirected_not_implemented(self):
        G = nx.Graph()
        assert_raises(nx.NetworkXNotImplemented, nx.dag_longest_path, G)

    def test_unorderable_nodes(self):
        """Tests that computing the longest path does not depend on
        nodes being orderable.

        For more information, see issue #1989.

        """
        # TODO In Python 3, instances of the `object` class are
        # unorderable by default, so we wouldn't need to define our own
        # class here, we could just instantiate an instance of the
        # `object` class. However, we still support Python 2; when
        # support for Python 2 is dropped, this test can be simplified
        # by replacing `Unorderable()` by `object()`.
        class Unorderable(object):
            def __lt__(self, other):
                error_msg = "< not supported between instances of " \
                    "{} and {}".format(type(self).__name__, type(other).__name__)
                raise TypeError(error_msg)

        # Create the directed path graph on four nodes in a diamond shape,
        # with nodes represented as (unorderable) Python objects.
        nodes = [Unorderable() for n in range(4)]
        G = nx.DiGraph()
        G.add_edge(nodes[0], nodes[1])
        G.add_edge(nodes[0], nodes[2])
        G.add_edge(nodes[2], nodes[3])
        G.add_edge(nodes[1], nodes[3])

        # this will raise NotImplementedError when nodes need to be ordered
        nx.dag_longest_path(G)


class TestDagLongestPathLength(object):
    """Unit tests for computing the length of a longest path in a
    directed acyclic graph.

    """

    def test_unweighted(self):
        edges = [(1, 2), (2, 3), (2, 4), (3, 5), (5, 6), (5, 7)]
        G = nx.DiGraph(edges)
        assert_equal(nx.dag_longest_path_length(G), 4)

        edges = [(1, 2), (2, 3), (3, 4), (4, 5), (1, 3), (1, 5), (3, 5)]
        G = nx.DiGraph(edges)
        assert_equal(nx.dag_longest_path_length(G), 4)

        # test degenerate graphs
        G = nx.DiGraph()
        G.add_node(1)
        assert_equal(nx.dag_longest_path_length(G), 0)

    def test_undirected_not_implemented(self):
        G = nx.Graph()
        assert_raises(nx.NetworkXNotImplemented, nx.dag_longest_path_length, G)

    def test_weighted(self):
        edges = [(1, 2, -5), (2, 3, 1), (3, 4, 1), (4, 5, 0), (3, 5, 4),
                 (1, 6, 2)]
        G = nx.DiGraph()
        G.add_weighted_edges_from(edges)
        assert_equal(nx.dag_longest_path_length(G), 5)


class TestDAG:

    def setUp(self):
        pass

    def test_topological_sort1(self):
        DG = nx.DiGraph([(1, 2), (1, 3), (2, 3)])

        for algorithm in [nx.topological_sort,
                          nx.lexicographical_topological_sort]:
            assert_equal(tuple(algorithm(DG)), (1, 2, 3))

        DG.add_edge(3, 2)

        for algorithm in [nx.topological_sort,
                          nx.lexicographical_topological_sort]:
            assert_raises(nx.NetworkXUnfeasible, consume, algorithm(DG))

        DG.remove_edge(2, 3)

        for algorithm in [nx.topological_sort,
                          nx.lexicographical_topological_sort]:
            assert_equal(tuple(algorithm(DG)), (1, 3, 2))

        DG.remove_edge(3, 2)

        assert_in(tuple(nx.topological_sort(DG)), {(1, 2, 3), (1, 3, 2)})
        assert_equal(tuple(nx.lexicographical_topological_sort(DG)), (1, 2, 3))

    def test_is_directed_acyclic_graph(self):
        G = nx.generators.complete_graph(2)
        assert_false(nx.is_directed_acyclic_graph(G))
        assert_false(nx.is_directed_acyclic_graph(G.to_directed()))
        assert_false(nx.is_directed_acyclic_graph(nx.Graph([(3, 4), (4, 5)])))
        assert_true(nx.is_directed_acyclic_graph(nx.DiGraph([(3, 4), (4, 5)])))

    def test_topological_sort2(self):
        DG = nx.DiGraph({1: [2], 2: [3], 3: [4],
                         4: [5], 5: [1], 11: [12],
                         12: [13], 13: [14], 14: [15]})
        assert_raises(nx.NetworkXUnfeasible, consume, nx.topological_sort(DG))

        assert_false(nx.is_directed_acyclic_graph(DG))

        DG.remove_edge(1, 2)
        consume(nx.topological_sort(DG))
        assert_true(nx.is_directed_acyclic_graph(DG))

    def test_topological_sort3(self):
        DG = nx.DiGraph()
        DG.add_edges_from([(1, i) for i in range(2, 5)])
        DG.add_edges_from([(2, i) for i in range(5, 9)])
        DG.add_edges_from([(6, i) for i in range(9, 12)])
        DG.add_edges_from([(4, i) for i in range(12, 15)])

        def validate(order):
            ok_(isinstance(order, list))
            assert_equal(set(order), set(DG))
            for u, v in combinations(order, 2):
                assert_false(nx.has_path(DG, v, u))
        validate(list(nx.topological_sort(DG)))

        DG.add_edge(14, 1)
        assert_raises(nx.NetworkXUnfeasible, consume, nx.topological_sort(DG))

    def test_topological_sort4(self):
        G = nx.Graph()
        G.add_edge(1, 2)
        # Only directed graphs can be topologically sorted.
        assert_raises(nx.NetworkXError, consume, nx.topological_sort(G))

    def test_topological_sort5(self):
        G = nx.DiGraph()
        G.add_edge(0, 1)
        assert_equal(list(nx.topological_sort(G)), [0, 1])

    def test_topological_sort6(self):
        for algorithm in [nx.topological_sort,
                          nx.lexicographical_topological_sort]:
            def runtime_error():
                DG = nx.DiGraph([(1, 2), (2, 3), (3, 4)])
                first = True
                for x in algorithm(DG):
                    if first:
                        first = False
                        DG.add_edge(5 - x, 5)

            def unfeasible_error():
                DG = nx.DiGraph([(1, 2), (2, 3), (3, 4)])
                first = True
                for x in algorithm(DG):
                    if first:
                        first = False
                        DG.remove_node(4)

            def runtime_error2():
                DG = nx.DiGraph([(1, 2), (2, 3), (3, 4)])
                first = True
                for x in algorithm(DG):
                    if first:
                        first = False
                        DG.remove_node(2)

            assert_raises(RuntimeError, runtime_error)
            assert_raises(RuntimeError, runtime_error2)
            assert_raises(nx.NetworkXUnfeasible, unfeasible_error)

    def test_ancestors(self):
        G = nx.DiGraph()
        ancestors = nx.algorithms.dag.ancestors
        G.add_edges_from([
            (1, 2), (1, 3), (4, 2), (4, 3), (4, 5), (2, 6), (5, 6)])
        assert_equal(ancestors(G, 6), set([1, 2, 4, 5]))
        assert_equal(ancestors(G, 3), set([1, 4]))
        assert_equal(ancestors(G, 1), set())
        assert_raises(nx.NetworkXError, ancestors, G, 8)

    def test_descendants(self):
        G = nx.DiGraph()
        descendants = nx.algorithms.dag.descendants
        G.add_edges_from([
            (1, 2), (1, 3), (4, 2), (4, 3), (4, 5), (2, 6), (5, 6)])
        assert_equal(descendants(G, 1), set([2, 3, 6]))
        assert_equal(descendants(G, 4), set([2, 3, 5, 6]))
        assert_equal(descendants(G, 3), set())
        assert_raises(nx.NetworkXError, descendants, G, 8)

    def test_transitive_closure(self):
        G = nx.DiGraph([(1, 2), (2, 3), (3, 4)])
        transitive_closure = nx.algorithms.dag.transitive_closure
        solution = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
        assert_edges_equal(transitive_closure(G).edges(), solution)
        G = nx.DiGraph([(1, 2), (2, 3), (2, 4)])
        solution = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)]
        assert_edges_equal(transitive_closure(G).edges(), solution)
        G = nx.Graph([(1, 2), (2, 3), (3, 4)])
        assert_raises(nx.NetworkXNotImplemented, transitive_closure, G)

        # test if edge data is copied
        G = nx.DiGraph([(1, 2, {"a": 3}), (2, 3, {"b": 0}), (3, 4)])
        H = transitive_closure(G)
        for u, v in G.edges():
            assert_equal(G.get_edge_data(u, v), H.get_edge_data(u, v))

        k = 10
        G = nx.DiGraph((i, i + 1, {"foo": "bar", "weight": i}) for i in range(k))
        H = transitive_closure(G)
        for u, v in G.edges():
            assert_equal(G.get_edge_data(u, v), H.get_edge_data(u, v))

    def test_transitive_reduction(self):
        G = nx.DiGraph([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)])
        transitive_reduction = nx.algorithms.dag.transitive_reduction
        solution = [(1, 2), (2, 3), (3, 4)]
        assert_edges_equal(transitive_reduction(G).edges(), solution)
        G = nx.DiGraph([(1, 2), (1, 3), (1, 4), (2, 3), (2, 4)])
        transitive_reduction = nx.algorithms.dag.transitive_reduction
        solution = [(1, 2), (2, 3), (2, 4)]
        assert_edges_equal(transitive_reduction(G).edges(), solution)
        G = nx.Graph([(1, 2), (2, 3), (3, 4)])
        assert_raises(nx.NetworkXNotImplemented, transitive_reduction, G)

    def _check_antichains(self, solution, result):
        sol = [frozenset(a) for a in solution]
        res = [frozenset(a) for a in result]
        assert_true(set(sol) == set(res))

    def test_antichains(self):
        antichains = nx.algorithms.dag.antichains
        G = nx.DiGraph([(1, 2), (2, 3), (3, 4)])
        solution = [[], [4], [3], [2], [1]]
        self._check_antichains(list(antichains(G)), solution)
        G = nx.DiGraph([(1, 2), (2, 3), (2, 4), (3, 5), (5, 6), (5, 7)])
        solution = [[], [4], [7], [7, 4], [6], [6, 4], [6, 7], [6, 7, 4],
                    [5], [5, 4], [3], [3, 4], [2], [1]]
        self._check_antichains(list(antichains(G)), solution)
        G = nx.DiGraph([(1, 2), (1, 3), (3, 4), (3, 5), (5, 6)])
        solution = [[], [6], [5], [4], [4, 6], [4, 5], [3], [2], [2, 6],
                    [2, 5], [2, 4], [2, 4, 6], [2, 4, 5], [2, 3], [1]]
        self._check_antichains(list(antichains(G)), solution)
        G = nx.DiGraph({0: [1, 2], 1: [4], 2: [3], 3: [4]})
        solution = [[], [4], [3], [2], [1], [1, 3], [1, 2], [0]]
        self._check_antichains(list(antichains(G)), solution)
        G = nx.DiGraph()
        self._check_antichains(list(antichains(G)), [[]])
        G = nx.DiGraph()
        G.add_nodes_from([0, 1, 2])
        solution = [[], [0], [1], [1, 0], [2], [2, 0], [2, 1], [2, 1, 0]]
        self._check_antichains(list(antichains(G)), solution)

        def f(x): return list(antichains(x))
        G = nx.Graph([(1, 2), (2, 3), (3, 4)])
        assert_raises(nx.NetworkXNotImplemented, f, G)
        G = nx.DiGraph([(1, 2), (2, 3), (3, 1)])
        assert_raises(nx.NetworkXUnfeasible, f, G)

    def test_lexicographical_topological_sort(self):
        G = nx.DiGraph([(1, 2), (2, 3), (1, 4), (1, 5), (2, 6)])
        assert_equal(list(nx.lexicographical_topological_sort(G)),
                     [1, 2, 3, 4, 5, 6])
        assert_equal(list(nx.lexicographical_topological_sort(
            G, key=lambda x: x)),
            [1, 2, 3, 4, 5, 6])
        assert_equal(list(nx.lexicographical_topological_sort(
            G, key=lambda x: -x)),
            [1, 5, 4, 2, 6, 3])


def test_is_aperiodic_cycle():
    G = nx.DiGraph()
    nx.add_cycle(G, [1, 2, 3, 4])
    assert_false(nx.is_aperiodic(G))


def test_is_aperiodic_cycle2():
    G = nx.DiGraph()
    nx.add_cycle(G, [1, 2, 3, 4])
    nx.add_cycle(G, [3, 4, 5, 6, 7])
    assert_true(nx.is_aperiodic(G))


def test_is_aperiodic_cycle3():
    G = nx.DiGraph()
    nx.add_cycle(G, [1, 2, 3, 4])
    nx.add_cycle(G, [3, 4, 5, 6])
    assert_false(nx.is_aperiodic(G))


def test_is_aperiodic_cycle4():
    G = nx.DiGraph()
    nx.add_cycle(G, [1, 2, 3, 4])
    G.add_edge(1, 3)
    assert_true(nx.is_aperiodic(G))


def test_is_aperiodic_selfloop():
    G = nx.DiGraph()
    nx.add_cycle(G, [1, 2, 3, 4])
    G.add_edge(1, 1)
    assert_true(nx.is_aperiodic(G))


def test_is_aperiodic_raise():
    G = nx.Graph()
    assert_raises(nx.NetworkXError,
                  nx.is_aperiodic,
                  G)


def test_is_aperiodic_bipartite():
    # Bipartite graph
    G = nx.DiGraph(nx.davis_southern_women_graph())
    assert_false(nx.is_aperiodic(G))


def test_is_aperiodic_rary_tree():
    G = nx.full_rary_tree(3, 27, create_using=nx.DiGraph())
    assert_false(nx.is_aperiodic(G))


def test_is_aperiodic_disconnected():
    # disconnected graph
    G = nx.DiGraph()
    nx.add_cycle(G, [1, 2, 3, 4])
    nx.add_cycle(G, [5, 6, 7, 8])
    assert_false(nx.is_aperiodic(G))
    G.add_edge(1, 3)
    G.add_edge(5, 7)
    assert_true(nx.is_aperiodic(G))


def test_is_aperiodic_disconnected2():
    G = nx.DiGraph()
    nx.add_cycle(G, [0, 1, 2])
    G.add_edge(3, 3)
    assert_false(nx.is_aperiodic(G))


class TestDagToBranching(object):
    """Unit tests for the :func:`networkx.dag_to_branching` function."""

    def test_single_root(self):
        """Tests that a directed acyclic graph with a single degree
        zero node produces an arborescence.

        """
        G = nx.DiGraph([(0, 1), (0, 2), (1, 3), (2, 3)])
        B = nx.dag_to_branching(G)
        expected = nx.DiGraph([(0, 1), (1, 3), (0, 2), (2, 4)])
        assert_true(nx.is_arborescence(B))
        assert_true(nx.is_isomorphic(B, expected))

    def test_multiple_roots(self):
        """Tests that a directed acyclic graph with multiple degree zero
        nodes creates an arborescence with multiple (weakly) connected
        components.

        """
        G = nx.DiGraph([(0, 1), (0, 2), (1, 3), (2, 3), (5, 2)])
        B = nx.dag_to_branching(G)
        expected = nx.DiGraph([(0, 1), (1, 3), (0, 2), (2, 4), (5, 6), (6, 7)])
        assert_true(nx.is_branching(B))
        assert_false(nx.is_arborescence(B))
        assert_true(nx.is_isomorphic(B, expected))

    # # Attributes are not copied by this function. If they were, this would
    # # be a good test to uncomment.
    # def test_copy_attributes(self):
    #     """Tests that node attributes are copied in the branching."""
    #     G = nx.DiGraph([(0, 1), (0, 2), (1, 3), (2, 3)])
    #     for v in G:
    #         G.node[v]['label'] = str(v)
    #     B = nx.dag_to_branching(G)
    #     # Determine the root node of the branching.
    #     root = next(v for v, d in B.in_degree() if d == 0)
    #     assert_equal(B.node[root]['label'], '0')
    #     children = B[root]
    #     # Get the left and right children, nodes 1 and 2, respectively.
    #     left, right = sorted(children, key=lambda v: B.node[v]['label'])
    #     assert_equal(B.node[left]['label'], '1')
    #     assert_equal(B.node[right]['label'], '2')
    #     # Get the left grandchild.
    #     children = B[left]
    #     assert_equal(len(children), 1)
    #     left_grandchild = arbitrary_element(children)
    #     assert_equal(B.node[left_grandchild]['label'], '3')
    #     # Get the right grandchild.
    #     children = B[right]
    #     assert_equal(len(children), 1)
    #     right_grandchild = arbitrary_element(children)
    #     assert_equal(B.node[right_grandchild]['label'], '3')

    def test_already_arborescence(self):
        """Tests that a directed acyclic graph that is already an
        arborescence produces an isomorphic arborescence as output.

        """
        A = nx.balanced_tree(2, 2, create_using=nx.DiGraph())
        B = nx.dag_to_branching(A)
        assert_true(nx.is_isomorphic(A, B))

    def test_already_branching(self):
        """Tests that a directed acyclic graph that is already a
        branching produces an isomorphic branching as output.

        """
        T1 = nx.balanced_tree(2, 2, create_using=nx.DiGraph())
        T2 = nx.balanced_tree(2, 2, create_using=nx.DiGraph())
        G = nx.disjoint_union(T1, T2)
        B = nx.dag_to_branching(G)
        assert_true(nx.is_isomorphic(G, B))

    @raises(nx.HasACycle)
    def test_not_acyclic(self):
        """Tests that a non-acyclic graph causes an exception."""
        G = nx.DiGraph(pairwise('abc', cyclic=True))
        nx.dag_to_branching(G)

    @raises(nx.NetworkXNotImplemented)
    def test_undirected(self):
        nx.dag_to_branching(nx.Graph())

    @raises(nx.NetworkXNotImplemented)
    def test_multigraph(self):
        nx.dag_to_branching(nx.MultiGraph())

    @raises(nx.NetworkXNotImplemented)
    def test_multidigraph(self):
        nx.dag_to_branching(nx.MultiDiGraph())
