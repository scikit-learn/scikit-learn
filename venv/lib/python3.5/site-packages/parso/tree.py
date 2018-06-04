from abc import abstractmethod, abstractproperty
from parso._compatibility import utf8_repr, encoding, py_version


def search_ancestor(node, *node_types):
    """
    Recursively looks at the parents of a node and returns the first found node
    that matches node_types. Returns ``None`` if no matching node is found.

    :param node: The ancestors of this node will be checked.
    :param node_types: type names that are searched for.
    :type node_types: tuple of str
    """
    while True:
        node = node.parent
        if node is None or node.type in node_types:
            return node


class NodeOrLeaf(object):
    """
    The base class for nodes and leaves.
    """
    __slots__ = ()
    type = None
    '''
    The type is a string that typically matches the types of the grammar file.
    '''

    def get_root_node(self):
        """
        Returns the root node of a parser tree. The returned node doesn't have
        a parent node like all the other nodes/leaves.
        """
        scope = self
        while scope.parent is not None:
            scope = scope.parent
        return scope

    def get_next_sibling(self):
        """
        Returns the node immediately following this node in this parent's
        children list. If this node does not have a next sibling, it is None
        """
        # Can't use index(); we need to test by identity
        for i, child in enumerate(self.parent.children):
            if child is self:
                try:
                    return self.parent.children[i + 1]
                except IndexError:
                    return None

    def get_previous_sibling(self):
        """
        Returns the node immediately preceding this node in this parent's
        children list. If this node does not have a previous sibling, it is
        None.
        """
        # Can't use index(); we need to test by identity
        for i, child in enumerate(self.parent.children):
            if child is self:
                if i == 0:
                    return None
                return self.parent.children[i - 1]

    def get_previous_leaf(self):
        """
        Returns the previous leaf in the parser tree.
        Returns `None` if this is the first element in the parser tree.
        """
        node = self
        while True:
            c = node.parent.children
            i = c.index(node)
            if i == 0:
                node = node.parent
                if node.parent is None:
                    return None
            else:
                node = c[i - 1]
                break

        while True:
            try:
                node = node.children[-1]
            except AttributeError:  # A Leaf doesn't have children.
                return node

    def get_next_leaf(self):
        """
        Returns the next leaf in the parser tree.
        Returns None if this is the last element in the parser tree.
        """
        node = self
        while True:
            c = node.parent.children
            i = c.index(node)
            if i == len(c) - 1:
                node = node.parent
                if node.parent is None:
                    return None
            else:
                node = c[i + 1]
                break

        while True:
            try:
                node = node.children[0]
            except AttributeError:  # A Leaf doesn't have children.
                return node

    @abstractproperty
    def start_pos(self):
        """
        Returns the starting position of the prefix as a tuple, e.g. `(3, 4)`.

        :return tuple of int: (line, column)
        """

    @abstractproperty
    def end_pos(self):
        """
        Returns the end position of the prefix as a tuple, e.g. `(3, 4)`.

        :return tuple of int: (line, column)
        """

    @abstractmethod
    def get_start_pos_of_prefix(self):
        """
        Returns the start_pos of the prefix. This means basically it returns
        the end_pos of the last prefix. The `get_start_pos_of_prefix()` of the
        prefix `+` in `2 + 1` would be `(1, 1)`, while the start_pos is
        `(1, 2)`.

        :return tuple of int: (line, column)
        """

    @abstractmethod
    def get_first_leaf(self):
        """
        Returns the first leaf of a node or itself if this is a leaf.
        """

    @abstractmethod
    def get_last_leaf(self):
        """
        Returns the last leaf of a node or itself if this is a leaf.
        """

    @abstractmethod
    def get_code(self, include_prefix=True):
        """
        Returns the code that was input the input for the parser for this node.

        :param include_prefix: Removes the prefix (whitespace and comments) of
            e.g. a statement.
        """


class Leaf(NodeOrLeaf):
    '''
    Leafs are basically tokens with a better API. Leafs exactly know where they
    were defined and what text preceeds them.
    '''
    __slots__ = ('value', 'parent', 'line', 'column', 'prefix')

    def __init__(self, value, start_pos, prefix=''):
        self.value = value
        '''
        :py:func:`str` The value of the current token.
        '''
        self.start_pos = start_pos
        self.prefix = prefix
        '''
        :py:func:`str` Typically a mixture of whitespace and comments. Stuff
        that is syntactically irrelevant for the syntax tree.
        '''
        self.parent = None
        '''
        The parent :class:`BaseNode` of this leaf.
        '''

    @property
    def start_pos(self):
        return self.line, self.column

    @start_pos.setter
    def start_pos(self, value):
        self.line = value[0]
        self.column = value[1]

    def get_start_pos_of_prefix(self):
        previous_leaf = self.get_previous_leaf()
        if previous_leaf is None:
            return self.line - self.prefix.count('\n'), 0  # It's the first leaf.
        return previous_leaf.end_pos

    def get_first_leaf(self):
        return self

    def get_last_leaf(self):
        return self

    def get_code(self, include_prefix=True):
        if include_prefix:
            return self.prefix + self.value
        else:
            return self.value

    @property
    def end_pos(self):
        lines = self.value.split('\n')
        end_pos_line = self.line + len(lines) - 1
        # Check for multiline token
        if self.line == end_pos_line:
            end_pos_column = self.column + len(lines[-1])
        else:
            end_pos_column = len(lines[-1])
        return end_pos_line, end_pos_column

    @utf8_repr
    def __repr__(self):
        value = self.value
        if not value:
            value = self.type
        return "<%s: %s>" % (type(self).__name__, value)


class TypedLeaf(Leaf):
    __slots__ = ('type',)
    def __init__(self, type, value, start_pos, prefix=''):
        super(TypedLeaf, self).__init__(value, start_pos, prefix)
        self.type = type


class BaseNode(NodeOrLeaf):
    """
    The super class for all nodes.
    A node has children, a type and possibly a parent node.
    """
    __slots__ = ('children', 'parent')
    type = None

    def __init__(self, children):
        for c in children:
            c.parent = self
        self.children = children
        """
        A list of :class:`NodeOrLeaf` child nodes.
        """
        self.parent = None
        '''
        The parent :class:`BaseNode` of this leaf.
        None if this is the root node.
        '''

    @property
    def start_pos(self):
        return self.children[0].start_pos

    def get_start_pos_of_prefix(self):
        return self.children[0].get_start_pos_of_prefix()

    @property
    def end_pos(self):
        return self.children[-1].end_pos

    def _get_code_for_children(self, children, include_prefix):
        if include_prefix:
            return "".join(c.get_code() for c in children)
        else:
            first = children[0].get_code(include_prefix=False)
            return first + "".join(c.get_code() for c in children[1:])

    def get_code(self, include_prefix=True):
        return self._get_code_for_children(self.children, include_prefix)

    def get_leaf_for_position(self, position, include_prefixes=False):
        """
        Get the :py:class:`parso.tree.Leaf` at ``position``

        :param tuple position: A position tuple, row, column. Rows start from 1
        :param bool include_prefixes: If ``False``, ``None`` will be returned if ``position`` falls
            on whitespace or comments before a leaf
        :return: :py:class:`parso.tree.Leaf` at ``position``, or ``None``
        """
        def binary_search(lower, upper):
            if lower == upper:
                element = self.children[lower]
                if not include_prefixes and position < element.start_pos:
                    # We're on a prefix.
                    return None
                # In case we have prefixes, a leaf always matches
                try:
                    return element.get_leaf_for_position(position, include_prefixes)
                except AttributeError:
                    return element


            index = int((lower + upper) / 2)
            element = self.children[index]
            if position <= element.end_pos:
                return binary_search(lower, index)
            else:
                return binary_search(index + 1, upper)

        if not ((1, 0) <= position <= self.children[-1].end_pos):
            raise ValueError('Please provide a position that exists within this node.')
        return binary_search(0, len(self.children) - 1)

    def get_first_leaf(self):
        return self.children[0].get_first_leaf()

    def get_last_leaf(self):
        return self.children[-1].get_last_leaf()

    @utf8_repr
    def __repr__(self):
        code = self.get_code().replace('\n', ' ').strip()
        if not py_version >= 30:
            code = code.encode(encoding, 'replace')
        return "<%s: %s@%s,%s>" % \
            (type(self).__name__, code, self.start_pos[0], self.start_pos[1])


class Node(BaseNode):
    """Concrete implementation for interior nodes."""
    __slots__ = ('type',)

    def __init__(self, type, children):
        super(Node, self).__init__(children)
        self.type = type

    def __repr__(self):
        return "%s(%s, %r)" % (self.__class__.__name__, self.type, self.children)


class ErrorNode(BaseNode):
    """
    A node that contains valid nodes/leaves that we're follow by a token that
    was invalid. This basically means that the leaf after this node is where
    Python would mark a syntax error.
    """
    __slots__ = ()
    type = 'error_node'


class ErrorLeaf(Leaf):
    """
    A leaf that is either completely invalid in a language (like `$` in Python)
    or is invalid at that position. Like the star in `1 +* 1`.
    """
    __slots__ = ('original_type',)
    type = 'error_leaf'

    def __init__(self, original_type, value, start_pos, prefix=''):
        super(ErrorLeaf, self).__init__(value, start_pos, prefix)
        self.original_type = original_type

    def __repr__(self):
        return "<%s: %s:%s, %s>" % \
            (type(self).__name__, self.original_type, repr(self.value), self.start_pos)
