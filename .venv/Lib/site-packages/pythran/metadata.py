"""
This module provides a way to pass information between passes as metadata.

* add attaches a metadata to a node
* get retrieves all metadata from a particular class attached to a node
"""

from gast import AST  # so that metadata are walkable as regular ast nodes


class Metadata(AST):

    """ Base class to add information on a node to improve code generation. """

    def __init__(self):
        """ Initialize content of these metadata. """
        self.data = list()
        self._fields = ('data',)
        super(Metadata, self).__init__()

    def __iter__(self):
        """ Enable iteration over every metadata informations. """
        return iter(self.data)

    def append(self, data):
        """ Add a metadata information. """
        self.data.append(data)


class Lazy(AST):

    """ Metadata to mark variable which doesn't need to be evaluated now. """


class Comprehension(AST):
    def __init__(self, *args):  # no positional argument to be deep copyable
        super(Comprehension, self).__init__()
        if args:
            self.target = args[0]


class StaticReturn(AST):

    """ Metadata to mark return with a constant value. """


class Local(AST):
    """ Metadata to mark function as non exported. """


def add(node, data):
    if not hasattr(node, 'metadata'):
        node.metadata = Metadata()
        node._fields += ('metadata',)
    node.metadata.append(data)


def get(node, class_):
    if hasattr(node, 'metadata'):
        return [s for s in node.metadata if isinstance(s, class_)]
    else:
        return []


def clear(node, class_):
    if hasattr(node, 'metadata'):
        node.metadata.data = [s for s in node.metadata
                              if not isinstance(s, class_)]
        if not node.metadata.data:
            del node.metadata
            assert node._fields[-1] == 'metadata'
            node._fields = node._fields[:-1]


def visit(self, node):
    if hasattr(node, 'metadata'):
        self.visit(node.metadata)
