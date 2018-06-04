"""
The ``Parser`` tries to convert the available Python code in an easy to read
format, something like an abstract syntax tree. The classes who represent this
tree, are sitting in the :mod:`parso.tree` module.

The Python module ``tokenize`` is a very important part in the ``Parser``,
because it splits the code into different words (tokens).  Sometimes it looks a
bit messy. Sorry for that! You might ask now: "Why didn't you use the ``ast``
module for this? Well, ``ast`` does a very good job understanding proper Python
code, but fails to work as soon as there's a single line of broken code.

There's one important optimization that needs to be known: Statements are not
being parsed completely. ``Statement`` is just a representation of the tokens
within the statement. This lowers memory usage and cpu time and reduces the
complexity of the ``Parser`` (there's another parser sitting inside
``Statement``, which produces ``Array`` and ``Call``).
"""
from parso import tree
from parso.pgen2.parse import PgenParser


class ParserSyntaxError(Exception):
    """
    Contains error information about the parser tree.

    May be raised as an exception.
    """
    def __init__(self, message, error_leaf):
        self.message = message
        self.error_leaf = error_leaf


class BaseParser(object):
    node_map = {}
    default_node = tree.Node

    leaf_map = {
    }
    default_leaf = tree.Leaf

    def __init__(self, pgen_grammar, start_symbol='file_input', error_recovery=False):
        self._pgen_grammar = pgen_grammar
        self._start_symbol = start_symbol
        self._error_recovery = error_recovery

    def parse(self, tokens):
        start_number = self._pgen_grammar.symbol2number[self._start_symbol]
        self.pgen_parser = PgenParser(
            self._pgen_grammar, self.convert_node, self.convert_leaf,
            self.error_recovery, start_number
        )

        node = self.pgen_parser.parse(tokens)
        # The stack is empty now, we don't need it anymore.
        del self.pgen_parser
        return node

    def error_recovery(self, pgen_grammar, stack, arcs, typ, value, start_pos, prefix,
                       add_token_callback):
        if self._error_recovery:
            raise NotImplementedError("Error Recovery is not implemented")
        else:
            error_leaf = tree.ErrorLeaf('TODO %s' % typ, value, start_pos, prefix)
            raise ParserSyntaxError('SyntaxError: invalid syntax', error_leaf)

    def convert_node(self, pgen_grammar, type_, children):
        # TODO REMOVE symbol, we don't want type here.
        symbol = pgen_grammar.number2symbol[type_]
        try:
            return self.node_map[symbol](children)
        except KeyError:
            return self.default_node(symbol, children)

    def convert_leaf(self, pgen_grammar, type_, value, prefix, start_pos):
        try:
            return self.leaf_map[type_](value, start_pos, prefix)
        except KeyError:
            return self.default_leaf(value, start_pos, prefix)
