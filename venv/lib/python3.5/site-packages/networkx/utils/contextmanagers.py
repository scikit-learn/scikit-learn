from __future__ import absolute_import

from contextlib import contextmanager

__all__ = [
    'reversed',
]


@contextmanager
def reversed(G):
    """A context manager for temporarily reversing a directed graph in place.

    This is a no-op for undirected graphs.

    Parameters
    ----------
    G : graph
        A NetworkX graph.
    """
    directed = G.is_directed()
    if directed:
        G._pred, G._succ = G._succ, G._pred
        G._adj = G._succ

    try:
        yield
    finally:
        if directed:
            # Reverse the reverse.
            G._pred, G._succ = G._succ, G._pred
            G._adj = G._succ
