'''
Minimal directed graph replacement for networkx.DiGraph

This has the sole advantage of being a standalone file that doesn't bring any
dependency with it.
'''

class DiGraph(object):
    def __init__(self):
        # adjacency[i][j] = True means j is a successor of i
        self._adjacency = {}
        self._edges = {}

    def successors(self, node):
        return (n for n in self._adjacency[node])

    def predecessors(self, node):
        return (k for k, v in self._adjacency.items() if node in v)

    def add_node(self, node):
        self._adjacency.setdefault(node, {})

    def add_edge(self, src, dest, **props):
        self.add_node(dest)
        self._adjacency.setdefault(src, {}).setdefault(dest, None)
        self._edges[(src, dest)] = props

    @property
    def edges(self):
        return self._edges

    def has_edge(self, src, dest):
        return dest in self._adjacency[src]

    def remove_edge(self, src, dest):
        self._adjacency[src].pop(dest)
        del self._edges[(src, dest)]

    def __len__(self):
        return len(self._adjacency)

    def __iter__(self):
        return iter(self._adjacency.keys())

    def __contains__(self, value):
        return value in self._adjacency

    def __getitem__(self, node):
        return self._adjacency[node]

class Unfeasible(RuntimeError):
    pass

def has_path(graph, src, dest):
    visited = set()
    worklist = [src]
    while worklist:
        current = worklist.pop()
        if current in visited:
            continue
        visited.add(current)
        if dest in graph.successors(current):
            return True
        worklist.extend(graph.successors(current))
    return False
