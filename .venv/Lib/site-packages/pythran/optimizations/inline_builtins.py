""" Expand some builtins implementation when it is profitable."""

from pythran.analyses import Aliases
from pythran.passmanager import Transformation
from pythran.tables import MODULES
from pythran.intrinsic import FunctionIntr
from pythran.utils import path_to_attr, path_to_node
from pythran.syntax import PythranSyntaxError

from copy import deepcopy
import gast as ast


class InlineBuiltins(Transformation[Aliases]):

    """
    Replace some builtins by their bodies.

    This may trigger some extra optimizations later on!

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> pm = passmanager.PassManager("test")
    >>> node = ast.parse('''
    ... def foo(a):
    ...     return  a + 1
    ... def bar(b):
    ...     return builtins.map(bar, (1, 2))''')
    >>> _, node = pm.apply(InlineBuiltins, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(a):
        return (a + 1)
    def bar(b):
        return [bar(1), bar(2)]
    """

    def inlineBuiltinsXMap(self, node):
        self.update = True

        elts = []
        nelts = min(len(n.elts) for n in node.args[1:])
        for i in range(nelts):
            elts.append([n.elts[i] for n in node.args[1:]])
        return ast.List([ast.Call(node.args[0], elt, []) for elt in elts],
                        ast.Load())

    def inlineBuiltinsMap(self, node):

        if not isinstance(node, ast.Call):
            return node

        func_aliases = self.aliases[node.func]
        if len(func_aliases) != 1:
            return node

        obj = next(iter(func_aliases))
        if obj is not MODULES['builtins']['map']:
            return node

        if not all(isinstance(arg, (ast.List, ast.Tuple))
                   for arg in node.args[1:]):
            return node

        mapped_func_aliases = self.aliases[node.args[0]]
        if len(mapped_func_aliases) != 1:
            return node

        obj = next(iter(mapped_func_aliases))
        if not isinstance(obj, (ast.FunctionDef, FunctionIntr)):
            return node

        # all preconditions are met, do it!
        return self.inlineBuiltinsXMap(node)

    def visit_Call(self, node):
        node = self.generic_visit(node)
        node = self.inlineBuiltinsMap(node)
        return node

    def make_array_index(self, base, size, index):
        if isinstance(base, ast.Constant):
            return ast.Constant(base.value, None)
        if size == 1:
            return deepcopy(base.elts[0])
        return base.elts[index]

    def fixedSizeArray(self, node):
        if isinstance(node, ast.Constant):
            return node, 1

        if isinstance(node, (ast.List, ast.Tuple)):
            return node, len(node.elts)

        if not isinstance(node, ast.Call):
            return None, 0

        func_aliases = self.aliases[node.func]
        if len(func_aliases) != 1:
            return None, 0

        obj = next(iter(func_aliases))
        if obj not in (MODULES['numpy']['array'], MODULES['numpy']['asarray']):
            return None, 0

        if len(node.args) != 1:
            return None, 0

        if isinstance(node.args[0], (ast.List, ast.Tuple)):
            return node.args[0], len(node.args[0].elts)

        return None, 0

    def inlineFixedSizeArrayBinOp(self, node):

        alike = ast.List, ast.Tuple, ast.Constant
        if isinstance(node.left, alike) and isinstance(node.right, alike):
            return node

        lbase, lsize = self.fixedSizeArray(node.left)
        rbase, rsize = self.fixedSizeArray(node.right)
        if not lbase or not rbase:
            return node

        if rsize != 1 and lsize != 1 and rsize != lsize:
            raise PythranSyntaxError("Invalid numpy broadcasting", node)

        self.update = True

        operands = [ast.BinOp(self.make_array_index(lbase, lsize, i),
                              type(node.op)(),
                              self.make_array_index(rbase, rsize, i))
                    for i in range(max(lsize, rsize))]
        res = ast.Call(path_to_attr(('numpy', 'array')),
                       [ast.Tuple(operands, ast.Load())],
                       [])
        self.aliases[res.func] = {path_to_node(('numpy', 'array'))}
        return res

    def visit_BinOp(self, node):
        node = self.generic_visit(node)
        node = self.inlineFixedSizeArrayBinOp(node)
        return node

    def inlineFixedSizeArrayUnaryOp(self, node):

        if isinstance(node.operand, (ast.Constant, ast.List, ast.Tuple)):
            return node

        base, size = self.fixedSizeArray(node.operand)
        if not base:
            return node

        self.update = True

        operands = [ast.UnaryOp(type(node.op)(),
                                self.make_array_index(base, size, i))
                    for i in range(size)]
        res = ast.Call(path_to_attr(('numpy', 'array')),
                       [ast.Tuple(operands, ast.Load())],
                       [])
        self.aliases[res.func] = {path_to_node(('numpy', 'array'))}
        return res

    def visit_UnaryOp(self, node):
        node = self.generic_visit(node)
        node = self.inlineFixedSizeArrayUnaryOp(node)
        return node

    def inlineFixedSizeArrayCompare(self, node):
        if len(node.comparators) != 1:
            return node

        node_right = node.comparators[0]

        alike = ast.Constant, ast.List, ast.Tuple
        if isinstance(node.left, alike) and isinstance(node_right, alike):
            return node

        lbase, lsize = self.fixedSizeArray(node.left)
        rbase, rsize = self.fixedSizeArray(node_right)
        if not lbase or not rbase:
            return node

        if rsize != 1 and lsize != 1 and rsize != lsize:
            raise PythranSyntaxError("Invalid numpy broadcasting", node)

        self.update = True

        operands = [ast.Compare(self.make_array_index(lbase, lsize, i),
                                [type(node.ops[0])()],
                                [self.make_array_index(rbase, rsize, i)])
                    for i in range(max(lsize, rsize))]
        res = ast.Call(path_to_attr(('numpy', 'array')),
                       [ast.Tuple(operands, ast.Load())],
                       [])
        self.aliases[res.func] = {path_to_node(('numpy', 'array'))}
        return res

    def visit_Compare(self, node):
        node = self.generic_visit(node)
        node = self.inlineFixedSizeArrayCompare(node)
        return node
