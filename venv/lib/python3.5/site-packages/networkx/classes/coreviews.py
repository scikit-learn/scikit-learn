#    Copyright (C) 2004-2018 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.
#
# Authors: Aric Hagberg (hagberg@lanl.gov),
#          Pieter Swart (swart@lanl.gov),
#          Dan Schult(dschult@colgate.edu)
"""
"""
from collections import Mapping
import networkx as nx

__all__ = ['AtlasView', 'AdjacencyView', 'MultiAdjacencyView',
           'UnionAtlas', 'UnionAdjacency',
           'UnionMultiInner', 'UnionMultiAdjacency',
           'FilterAtlas', 'FilterAdjacency',
           'FilterMultiInner', 'FilterMultiAdjacency',
           'ReadOnlyGraph',
           ]


class AtlasView(Mapping):
    """An AtlasView is a Read-only Mapping of Mappings.

    It is a View into a dict-of-dict data structure.
    The inner level of dict is read-write. But the
    outer level is read-only.

    See Also
    ========
    AdjacencyView - View into dict-of-dict-of-dict
    MultiAdjacencyView - View into dict-of-dict-of-dict-of-dict
    """
    __slots__ = ('_atlas',)

    def __getstate__(self):
        return {'_atlas': self._atlas}

    def __setstate__(self, state):
        self._atlas = state['_atlas']

    def __init__(self, d):
        self._atlas = d

    def __len__(self):
        return len(self._atlas)

    def __iter__(self):
        return iter(self._atlas)

    def __getitem__(self, key):
        return self._atlas[key]

    def copy(self):
        return {n: self[n].copy() for n in self._atlas}

    def __str__(self):
        return str(self._atlas)  # {nbr: self[nbr] for nbr in self})

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self._atlas)


class AdjacencyView(AtlasView):
    """An AdjacencyView is a Read-only Map of Maps of Maps.

    It is a View into a dict-of-dict-of-dict data structure.
    The inner level of dict is read-write. But the
    outer levels are read-only.

    See Also
    ========
    AtlasView - View into dict-of-dict
    MultiAdjacencyView - View into dict-of-dict-of-dict-of-dict
    """
    __slots__ = ()   # Still uses AtlasView slots names _atlas

    def __getitem__(self, name):
        return AtlasView(self._atlas[name])

    def copy(self):
        return {n: self[n].copy() for n in self._atlas}


class MultiAdjacencyView(AdjacencyView):
    """An MultiAdjacencyView is a Read-only Map of Maps of Maps of Maps.

    It is a View into a dict-of-dict-of-dict-of-dict data structure.
    The inner level of dict is read-write. But the
    outer levels are read-only.

    See Also
    ========
    AtlasView - View into dict-of-dict
    AdjacencyView - View into dict-of-dict-of-dict
    """
    __slots__ = ()   # Still uses AtlasView slots names _atlas

    def __getitem__(self, name):
        return AdjacencyView(self._atlas[name])

    def copy(self):
        return {n: self[n].copy() for n in self._atlas}


class UnionAtlas(Mapping):
    """A read-only union of two atlases (dict-of-dict).

    The two dict-of-dicts represent the inner dict of
    an Adjacency:  `G.succ[node]` and `G.pred[node]`.
    The inner level of dict of both hold attribute key:value
    pairs and is read-write. But the outer level is read-only.

    See Also
    ========
    UnionAdjacency - View into dict-of-dict-of-dict
    UnionMultiAdjacency - View into dict-of-dict-of-dict-of-dict
    """
    __slots__ = ('_succ', '_pred')

    def __getstate__(self):
        return {'_succ': self._succ, '_pred': self._pred}

    def __setstate__(self, state):
        self._succ = state['_succ']
        self._pred = state['_pred']

    def __init__(self, succ, pred):
        self._succ = succ
        self._pred = pred

    def __len__(self):
        return len(self._succ) + len(self._pred)

    def __iter__(self):
        return iter(set(self._succ.keys()) | set(self._pred.keys()))

    def __getitem__(self, key):
        try:
            return self._succ[key]
        except KeyError:
            return self._pred[key]

    def copy(self):
        result = {nbr: dd.copy() for nbr, dd in self._succ.items()}
        for nbr, dd in self._pred.items():
            if nbr in result:
                result[nbr].update(dd)
            else:
                result[nbr] = dd.copy()
        return result

    def __str__(self):
        return str({nbr: self[nbr] for nbr in self})

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self._succ, self._pred)


class UnionAdjacency(Mapping):
    """A read-only union of dict Adjacencies as a Map of Maps of Maps.

    The two input dict-of-dict-of-dicts represent the union of
    `G.succ` and `G.pred`. Return values are UnionAtlas
    The inner level of dict is read-write. But the
    middle and outer levels are read-only.

    succ : a dict-of-dict-of-dict {node: nbrdict}
    pred : a dict-of-dict-of-dict {node: nbrdict}
    The keys for the two dicts should be the same

    See Also
    ========
    UnionAtlas - View into dict-of-dict
    UnionMultiAdjacency - View into dict-of-dict-of-dict-of-dict
    """
    __slots__ = ('_succ', '_pred')

    def __getstate__(self):
        return {'_succ': self._succ, '_pred': self._pred}

    def __setstate__(self, state):
        self._succ = state['_succ']
        self._pred = state['_pred']

    def __init__(self, succ, pred):
        # keys must be the same for two input dicts
        assert(len(set(succ.keys()) ^ set(pred.keys())) == 0)
        self._succ = succ
        self._pred = pred

    def __len__(self):
        return len(self._succ)  # length of each dict should be the same

    def __iter__(self):
        return iter(self._succ)

    def __getitem__(self, nbr):
        return UnionAtlas(self._succ[nbr], self._pred[nbr])

    def copy(self):
        return {n: self[n].copy() for n in self._succ}

    def __str__(self):
        return str({nbr: self[nbr] for nbr in self})

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self._succ, self._pred)


class UnionMultiInner(UnionAtlas):
    """A read-only union of two inner dicts of MultiAdjacencies.

    The two input dict-of-dict-of-dicts represent the union of
    `G.succ[node]` and `G.pred[node]` for MultiDiGraphs.
    Return values are UnionAtlas.
    The inner level of dict is read-write. But the outer levels are read-only.

    See Also
    ========
    UnionAtlas - View into dict-of-dict
    UnionAdjacency - View into dict-of-dict-of-dict
    UnionMultiAdjacency - View into dict-of-dict-of-dict-of-dict
    """
    __slots__ = ()   # Still uses UnionAtlas slots names _succ, _pred

    def __getitem__(self, node):
        in_succ = node in self._succ
        in_pred = node in self._pred
        if in_succ:
            if in_pred:
                return UnionAtlas(self._succ[node], self._pred[node])
            return UnionAtlas(self._succ[node], {})
        return UnionAtlas({}, self._pred[node])

    def copy(self):
        nodes = set(self._succ.keys()) | set(self._pred.keys())
        return {n: self[n].copy() for n in nodes}


class UnionMultiAdjacency(UnionAdjacency):
    """A read-only union of two dict MultiAdjacencies.

    The two input dict-of-dict-of-dict-of-dicts represent the union of
    `G.succ` and `G.pred` for MultiDiGraphs. Return values are UnionAdjacency.
    The inner level of dict is read-write. But the outer levels are read-only.

    See Also
    ========
    UnionAtlas - View into dict-of-dict
    UnionMultiInner - View into dict-of-dict-of-dict
    """
    __slots__ = ()   # Still uses UnionAdjacency slots names _succ, _pred

    def __getitem__(self, node):
        return UnionMultiInner(self._succ[node], self._pred[node])


class ReadOnlyGraph(object):
    """A Mixin Class to mask the write methods of a graph class."""

    def not_allowed(self, *args, **kwds):
        msg = "SubGraph Views are readonly. Mutations not allowed"
        raise nx.NetworkXError(msg)

    add_node = not_allowed
    remove_node = not_allowed
    add_nodes_from = not_allowed
    remove_nodes_from = not_allowed

    add_edge = not_allowed
    remove_edge = not_allowed
    add_edges_from = not_allowed
    add_weighted_edges_from = not_allowed
    remove_edges_from = not_allowed

    clear = not_allowed


class FilterAtlas(Mapping):  # nodedict, nbrdict, keydict
    def __init__(self, d, NODE_OK):
        self._atlas = d
        self.NODE_OK = NODE_OK

    def __len__(self):
        return sum(1 for n in self)

    def __iter__(self):
        if hasattr(self.NODE_OK, 'nodes'):
            return (n for n in self.NODE_OK.nodes if n in self._atlas)
        return (n for n in self._atlas if self.NODE_OK(n))

    def __getitem__(self, key):
        if key in self._atlas and self.NODE_OK(key):
            return self._atlas[key]
        raise KeyError("Key {} not found".format(key))

    def copy(self):
        if hasattr(self.NODE_OK, 'nodes'):
            return {u: self._atlas[u] for u in self.NODE_OK.nodes
                    if u in self._atlas}
        return {u: d for u, d in self._atlas.items()
                if self.NODE_OK(u)}

    def __str__(self):
        return str({nbr: self[nbr] for nbr in self})

    def __repr__(self):
        return '%s(%r, %r)' % (self.__class__.__name__, self._atlas,
                               self.NODE_OK)


class FilterAdjacency(Mapping):   # edgedict
    def __init__(self, d, NODE_OK, EDGE_OK):
        self._atlas = d
        self.NODE_OK = NODE_OK
        self.EDGE_OK = EDGE_OK

    def __len__(self):
        return sum(1 for n in self)

    def __iter__(self):
        if hasattr(self.NODE_OK, 'nodes'):
            return (n for n in self.NODE_OK.nodes if n in self._atlas)
        return (n for n in self._atlas if self.NODE_OK(n))

    def __getitem__(self, node):
        if node in self._atlas and self.NODE_OK(node):
            def new_node_ok(nbr):
                return self.NODE_OK(nbr) and self.EDGE_OK(node, nbr)
            return FilterAtlas(self._atlas[node], new_node_ok)
        raise KeyError("Key {} not found".format(node))

    def copy(self):
        if hasattr(self.NODE_OK, 'nodes'):
            return {u: {v: d for v, d in self._atlas[u].items()
                        if self.NODE_OK(v) if self.EDGE_OK(u, v)}
                    for u in self.NODE_OK.nodes if u in self._atlas}
        return {u: {v: d for v, d in nbrs.items() if self.NODE_OK(v)
                    if self.EDGE_OK(u, v)}
                for u, nbrs in self._atlas.items()
                if self.NODE_OK(u)}

    def __str__(self):
        return str({nbr: self[nbr] for nbr in self})

    def __repr__(self):
        return '%s(%r, %r, %r)' % (self.__class__.__name__, self._atlas,
                                   self.NODE_OK, self.EDGE_OK)


class FilterMultiInner(FilterAdjacency):  # muliedge_seconddict
    def __iter__(self):
        if hasattr(self.NODE_OK, 'nodes'):
            my_nodes = (n for n in self.NODE_OK.nodes if n in self._atlas)
        else:
            my_nodes = (n for n in self._atlas if self.NODE_OK(n))
        for n in my_nodes:
            some_keys_ok = False
            for key in self._atlas[n]:
                if self.EDGE_OK(n, key):
                    some_keys_ok = True
                    break
            if some_keys_ok is True:
                yield n

    def __getitem__(self, nbr):
        if nbr in self._atlas and self.NODE_OK(nbr):
            def new_node_ok(key):
                return self.EDGE_OK(nbr, key)
            return FilterAtlas(self._atlas[nbr], new_node_ok)
        raise KeyError("Key {} not found".format(nbr))

    def copy(self):
        if hasattr(self.NODE_OK, 'nodes'):
            return {v: {k: d for k, d in self._atlas[v].items()
                        if self.EDGE_OK(v, k)}
                    for v in self.NODE_OK.nodes if v in self._atlas}
        return {v: {k: d for k, d in nbrs.items() if self.EDGE_OK(v, k)}
                for v, nbrs in self._atlas.items() if self.NODE_OK(v)}


class FilterMultiAdjacency(FilterAdjacency):  # multiedgedict
    def __getitem__(self, node):
        if node in self._atlas and self.NODE_OK(node):
            def edge_ok(nbr, key):
                return self.NODE_OK(nbr) and self.EDGE_OK(node, nbr, key)
            return FilterMultiInner(self._atlas[node], self.NODE_OK, edge_ok)
        raise KeyError("Key {} not found".format(node))

    def copy(self):
        if hasattr(self.NODE_OK, 'nodes'):
            my_nodes = self.NODE_OK.nodes
            return {u: {v: {k: d for k, d in kd.items()
                            if self.EDGE_OK(u, v, k)}
                        for v, kd in self._atlas[u].items() if v in my_nodes}
                    for u in my_nodes if u in self._atlas}
        return {u: {v: {k: d for k, d in kd.items()
                        if self.EDGE_OK(u, v, k)}
                    for v, kd in nbrs.items() if self.NODE_OK(v)}
                for u, nbrs in self._atlas.items() if self.NODE_OK(u)}
