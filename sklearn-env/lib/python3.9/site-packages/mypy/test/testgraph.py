"""Test cases for graph processing code in build.py."""

import sys
from typing import AbstractSet, Dict, Set, List

from mypy.test.helpers import assert_equal, Suite
from mypy.build import BuildManager, State, BuildSourceSet
from mypy.modulefinder import SearchPaths
from mypy.build import topsort, strongly_connected_components, sorted_components, order_ascc
from mypy.version import __version__
from mypy.options import Options
from mypy.report import Reports
from mypy.plugin import Plugin
from mypy.errors import Errors
from mypy.fscache import FileSystemCache


class GraphSuite(Suite):

    def test_topsort(self) -> None:
        a = frozenset({'A'})
        b = frozenset({'B'})
        c = frozenset({'C'})
        d = frozenset({'D'})
        data: Dict[AbstractSet[str], Set[AbstractSet[str]]] = {a: {b, c}, b: {d}, c: {d}}
        res = list(topsort(data))
        assert_equal(res, [{d}, {b, c}, {a}])

    def test_scc(self) -> None:
        vertices = {"A", "B", "C", "D"}
        edges: Dict[str, List[str]] = {"A": ["B", "C"], "B": ["C"], "C": ["B", "D"], "D": []}
        sccs = set(frozenset(x) for x in strongly_connected_components(vertices, edges))
        assert_equal(sccs,
                     {frozenset({'A'}),
                      frozenset({'B', 'C'}),
                      frozenset({'D'})})

    def _make_manager(self) -> BuildManager:
        errors = Errors()
        options = Options()
        fscache = FileSystemCache()
        search_paths = SearchPaths((), (), (), ())
        manager = BuildManager(
            data_dir='',
            search_paths=search_paths,
            ignore_prefix='',
            source_set=BuildSourceSet([]),
            reports=Reports('', {}),
            options=options,
            version_id=__version__,
            plugin=Plugin(options),
            plugins_snapshot={},
            errors=errors,
            flush_errors=lambda msgs, serious: None,
            fscache=fscache,
            stdout=sys.stdout,
            stderr=sys.stderr,
        )
        return manager

    def test_sorted_components(self) -> None:
        manager = self._make_manager()
        graph = {'a': State('a', None, 'import b, c', manager),
                 'd': State('d', None, 'pass', manager),
                 'b': State('b', None, 'import c', manager),
                 'c': State('c', None, 'import b, d', manager),
                 }
        res = sorted_components(graph)
        assert_equal(res, [frozenset({'d'}), frozenset({'c', 'b'}), frozenset({'a'})])

    def test_order_ascc(self) -> None:
        manager = self._make_manager()
        graph = {'a': State('a', None, 'import b, c', manager),
                 'd': State('d', None, 'def f(): import a', manager),
                 'b': State('b', None, 'import c', manager),
                 'c': State('c', None, 'import b, d', manager),
                 }
        res = sorted_components(graph)
        assert_equal(res, [frozenset({'a', 'd', 'c', 'b'})])
        ascc = res[0]
        scc = order_ascc(graph, ascc)
        assert_equal(scc, ['d', 'c', 'b', 'a'])
