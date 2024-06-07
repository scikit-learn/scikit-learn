"""Helpers for manipulations with graphs."""

from __future__ import annotations

from typing import AbstractSet, Iterable, Iterator, TypeVar

T = TypeVar("T")


def strongly_connected_components(
    vertices: AbstractSet[T], edges: dict[T, list[T]]
) -> Iterator[set[T]]:
    """Compute Strongly Connected Components of a directed graph.

    Args:
      vertices: the labels for the vertices
      edges: for each vertex, gives the target vertices of its outgoing edges

    Returns:
      An iterator yielding strongly connected components, each
      represented as a set of vertices.  Each input vertex will occur
      exactly once; vertices not part of a SCC are returned as
      singleton sets.

    From https://code.activestate.com/recipes/578507/.
    """
    identified: set[T] = set()
    stack: list[T] = []
    index: dict[T, int] = {}
    boundaries: list[int] = []

    def dfs(v: T) -> Iterator[set[T]]:
        index[v] = len(stack)
        stack.append(v)
        boundaries.append(index[v])

        for w in edges[v]:
            if w not in index:
                yield from dfs(w)
            elif w not in identified:
                while index[w] < boundaries[-1]:
                    boundaries.pop()

        if boundaries[-1] == index[v]:
            boundaries.pop()
            scc = set(stack[index[v] :])
            del stack[index[v] :]
            identified.update(scc)
            yield scc

    for v in vertices:
        if v not in index:
            yield from dfs(v)


def prepare_sccs(
    sccs: list[set[T]], edges: dict[T, list[T]]
) -> dict[AbstractSet[T], set[AbstractSet[T]]]:
    """Use original edges to organize SCCs in a graph by dependencies between them."""
    sccsmap = {v: frozenset(scc) for scc in sccs for v in scc}
    data: dict[AbstractSet[T], set[AbstractSet[T]]] = {}
    for scc in sccs:
        deps: set[AbstractSet[T]] = set()
        for v in scc:
            deps.update(sccsmap[x] for x in edges[v])
        data[frozenset(scc)] = deps
    return data


def topsort(data: dict[T, set[T]]) -> Iterable[set[T]]:
    """Topological sort.

    Args:
      data: A map from vertices to all vertices that it has an edge
            connecting it to.  NOTE: This data structure
            is modified in place -- for normalization purposes,
            self-dependencies are removed and entries representing
            orphans are added.

    Returns:
      An iterator yielding sets of vertices that have an equivalent
      ordering.

    Example:
      Suppose the input has the following structure:

        {A: {B, C}, B: {D}, C: {D}}

      This is normalized to:

        {A: {B, C}, B: {D}, C: {D}, D: {}}

      The algorithm will yield the following values:

        {D}
        {B, C}
        {A}

    From https://code.activestate.com/recipes/577413/.
    """
    # TODO: Use a faster algorithm?
    for k, v in data.items():
        v.discard(k)  # Ignore self dependencies.
    for item in set.union(*data.values()) - set(data.keys()):
        data[item] = set()
    while True:
        ready = {item for item, dep in data.items() if not dep}
        if not ready:
            break
        yield ready
        data = {item: (dep - ready) for item, dep in data.items() if item not in ready}
    assert not data, f"A cyclic dependency exists amongst {data!r}"
