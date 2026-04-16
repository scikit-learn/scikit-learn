"""IterTransformation replaces expressions by iterators when possible."""

from pythran.analyses import PotentialIterator, Aliases
from pythran.passmanager import Transformation
from pythran.utils import path_to_attr, path_to_node


EQUIVALENT_ITERATORS = {
    ('builtins', "list"): None,
    ('builtins', "tuple"): None,
    ('numpy', "array"): None,
    ('numpy', "asarray"): None,
    ('numpy', "copy"): None,
}


class IterTransformation(Transformation[PotentialIterator, Aliases]):

    """
    Replaces expressions by iterators when possible.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('''
    ... def foo(l):
    ...     return builtins.sum(l)
    ... def bar(n):
    ...     return foo(builtins.list(n))
    ... ''')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(IterTransformation, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(l):
        return builtins.sum(l)
    def bar(n):
        return foo(n)
    """

    def find_matching_builtin(self, node):
        """
        Return matched keyword.

        If the node alias on a correct keyword (and only it), it matches.
        """
        for path in EQUIVALENT_ITERATORS.keys():
            correct_alias = {path_to_node(path)}
            if self.aliases[node.func] == correct_alias:
                return path

    def visit_Call(self, node):
        """Replace function call by its correct iterator if it is possible."""
        if node in self.potential_iterator:
            matched_path = self.find_matching_builtin(node)
            if matched_path is None:
                return self.generic_visit(node)

            # if any kind of specific (~ with more arg) behavior is required
            if len(node.args) != 1:
                return self.generic_visit(node)

            path = EQUIVALENT_ITERATORS[matched_path]
            if path:
                node.func = path_to_attr(path)
            else:
                node = node.args[0]

            self.update = True
        return self.generic_visit(node)
