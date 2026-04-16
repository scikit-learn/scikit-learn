""" Module to manage dependencies between pythran types. """

import gast as ast
import itertools
import os

from pythran.analyses import GlobalDeclarations
from pythran.errors import PythranInternalError
from pythran.passmanager import ModuleAnalysis
from pythran.types.conversion import PYTYPE_TO_CTYPE_TABLE
from pythran.utils import get_variable
from pythran.typing import List, Set, Dict, NDArray, Tuple, Pointer, Fun, Type
from pythran.graph import DiGraph


def pytype_to_deps_hpp(t):
    """python -> pythonic type hpp filename."""
    if t in (list, set, dict, tuple):
        return {f"{t.__name__}.hpp"}
    elif isinstance(t, List):
        return {'list.hpp'}.union(pytype_to_deps_hpp(t.__args__[0]))
    elif isinstance(t, Type):
        return {'type.hpp'}.union(pytype_to_deps_hpp(t.__args__[0]))
    elif isinstance(t, Set):
        return {'set.hpp'}.union(pytype_to_deps_hpp(t.__args__[0]))
    elif isinstance(t, Dict):
        tkey, tvalue = t.__args__
        return {'dict.hpp'}.union(pytype_to_deps_hpp(tkey),
                                  pytype_to_deps_hpp(tvalue))
    elif isinstance(t, Tuple):
        return {'tuple.hpp'}.union(*[pytype_to_deps_hpp(elt)
                                     for elt in t.__args__])
    elif isinstance(t, NDArray):
        out = {'ndarray.hpp'}
        # it's a transpose!
        if t.__args__[1].start == -1:
            out.add('numpy_texpr.hpp')
        return out.union(pytype_to_deps_hpp(t.__args__[0]))
    elif isinstance(t, Pointer):
        return {'pointer.hpp'}.union(pytype_to_deps_hpp(t.__args__[0]))
    elif isinstance(t, Fun):
        return {'cfun.hpp'}.union(*[pytype_to_deps_hpp(a) for a in t.__args__])
    elif t in PYTYPE_TO_CTYPE_TABLE:
        return {'{}.hpp'.format(t.__name__)}
    else:
        raise NotImplementedError("{0}:{1}".format(type(t), t))


def pytype_to_deps(t):
    """ python -> pythonic type header full path. """
    res = set()
    for hpp_dep in pytype_to_deps_hpp(t):
        res.add(os.path.join('pythonic', 'types', hpp_dep))
        res.add(os.path.join('pythonic', 'include', 'types', hpp_dep))
    return res


class TypeDependencies(ModuleAnalysis[GlobalDeclarations]):

    """
    Gathers the callees of each function required for type inference.

    This analyse produces a directed graph with functions as nodes and edges
    between nodes when a function might call another.

    Check usual behavior.
    >>> import gast as ast
    >>> from pythran import passmanager
    >>> pm = passmanager.PassManager("test")
    >>> node = ast.parse('''
    ... def foo(n):
    ...     return 1 if copy(n) else copy(n)
    ... def copy(n):
    ...     return n == 2''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    3

    foo result depend on : NoDeps and copy
    copy result depend on : NoDeps

    Check that content assignment is a dependency.
    >>> node = ast.parse('''
    ... def foo(n):
    ...     n[1] = copy(n)
    ...     return 1 if copy(n) else n
    ... def copy(n):
    ...     return n == 2''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    3

    foo result depend on : NoDeps and copy
    copy result depend on : NoDeps

    Check augassign add a dependencies but don't remove the old one.
    >>> node = ast.parse('''
    ... def bar(n):
    ...     return n
    ... def foo(n):
    ...     n[1] = copy(n)
    ...     n[1] += bar(1)
    ...     return 1 if copy(n) else n
    ... def copy(n):
    ...     return n == 2''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    5

    bar result depend on : NoDeps
    foo result depend on : NoDeps, bar and copy
    copy depend on : NoDeps

    Check a if statement handle both branches
    >>> node = ast.parse('''
    ... def bar(n):
    ...     return n
    ... def foo(n):
    ...     if n:
    ...         n = bar()
    ...     else:
    ...         n = 4
    ...     return 1 or n''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    3

    Check we do not add everything from a conditional statement.
    >>> node = ast.parse('''
    ... def bar(n):
    ...     return n
    ... def foo(n):
    ...     if n:
    ...         n = bar()
    ...         n = 3
    ...     else:
    ...         n = 4
    ...     return 1 or n''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    2

    bar result depend on : NoDeps
    foo result depend on : NoDeps only

    Check dependency on for target variable
    >>> node = ast.parse('''
    ... def bar(n):
    ...     return builtins.range(n)
    ... def foo(n):
    ...     for i in bar(n):
    ...         i = 2
    ...     return i''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    2

    bar result depend on : NoDeps
    foo result depend on : NoDeps

    Check dependency on for target variable with no deps if we don't start
    >>> node = ast.parse('''
    ... def bar(n):
    ...     return builtins.range(n)
    ... def foo(n):
    ...     i = 4
    ...     for i in bar(n):
    ...         pass
    ...     return i''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    3

    bar result depend on : NoDeps
    foo result depend on : NoDeps and bar

    Check dependency on for target variable with deps
    >>> node = ast.parse('''
    ... def bar(n):
    ...     return builtins.range(n)
    ... def foo(n):
    ...     for i in bar(n):
    ...         pass
    ...     return i''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    2

    bar result depend on : NoDeps
    foo result depend on : NoDeps and bar

    Check conditional without else branch.
    >>> node = ast.parse('''
    ... def foo(n):
    ...     res = 3
    ...     if n:
    ...         res = foo(n - 1)
    ...     return res''')
    >>> res = pm.gather(TypeDependencies, node)
    >>> len(res.edges)
    2

    foo result depend on : NoDeps and foo

    FIXME : We should use CFG to perform better function dependencies.

    Check conditional without break
    >> node = ast.parse('''
    .. def bar2(n):
    ..     return builtins.range(n)
    .. def bar(n):
    ..     return builtins.range(n)
    .. def foo(n):
    ..     for i in bar(n):
    ..         if i:
    ..             j = bar(n)
    ..             break
    ..         j = bar2(n)
    ..     return j''')
    >> res = pm.gather(TypeDependencies, node)
    >> len(res.edges)
    4

    bar result depend on : NoDeps
    bar2 result depend on : NoDeps
    foo result depend on : bar ad bar2
    """

    NoDeps = "None"

    ResultType = DiGraph

    def __init__(self):
        """ Create empty result graph and gather global declarations. """
        super().__init__()
        self.current_function = None
        self.naming = dict()  # variable to dependencies for current function.
        # variable to dependencies for current conditional statement
        self.in_cond = dict()

    def prepare(self, node):
        """
        Add nodes for each global declarations in the result graph.

        No edges are added as there are no type builtin type dependencies.
        """
        super(TypeDependencies, self).prepare(node)
        for v in self.global_declarations.values():
            self.result.add_node(v)
        self.result.add_node(TypeDependencies.NoDeps)

    def visit_any_conditionnal(self, node1, node2):
        """
        Set and restore the in_cond variable before visiting subnode.

        Compute correct dependencies on a value as both branch are possible
        path.
        """

        true_naming = false_naming = None

        try:
            tmp = self.naming.copy()
            for expr in node1:
                self.visit(expr)
            true_naming = self.naming
            self.naming = tmp
        except KeyError:
            pass

        try:
            tmp = self.naming.copy()
            for expr in node2:
                self.visit(expr)
            false_naming = self.naming
            self.naming = tmp
        except KeyError:
            pass

        if true_naming and not false_naming:
            self.naming = true_naming

        elif false_naming and not true_naming:
            self.naming = false_naming

        elif true_naming and false_naming:
            self.naming = false_naming
            for k, v in true_naming.items():
                if k not in self.naming:
                    self.naming[k] = v
                else:
                    for dep in v:
                        if dep not in self.naming[k]:
                            self.naming[k].append(dep)

    def visit_FunctionDef(self, node):
        """
        Initialize variable for the current function to add edges from calls.

        We compute variable to call dependencies and add edges when returns
        are reach.
        """
        # Ensure there are no nested functions.
        assert self.current_function is None
        self.current_function = node
        self.naming = dict()
        self.in_cond = False  # True when we are in a if, while or for
        self.generic_visit(node)
        self.current_function = None

    def visit_Return(self, node):
        """
        Add edge from all possible callee to current function.

        Gather all the function call that led to the creation of the
        returned expression and add an edge to each of this function.

        When visiting an expression, one returns a list of frozensets. Each
        element of the list is linked to a possible path, each element of a
        frozenset is linked to a dependency.
        """
        if not node.value:
            # Yielding function can't return values
            return
        for dep_set in self.visit(node.value):
            if dep_set:
                for dep in dep_set:
                    self.result.add_edge(dep, self.current_function)
            else:
                self.result.add_edge(TypeDependencies.NoDeps,
                                     self.current_function)

    visit_Yield = visit_Return

    def visit_Assign(self, node):
        """
        In case of assignment assign value depend on r-value type dependencies.

        It is valid for subscript, `a[i] = foo()` means `a` type depend on
        `foo` return type.
        """
        value_deps = self.visit(node.value)
        for target in node.targets:
            name = get_variable(target)
            if isinstance(name, ast.Name):
                self.naming[name.id] = value_deps

    def visit_AnnAssign(self, node):
        deps = []
        if node.value:
            deps.extend(self.visit(node.value))
        deps.extend(self.visit(node.annotation))
        name = get_variable(node.target)
        if isinstance(name, ast.Name):
            self.naming[name.id] = deps

    def visit_AugAssign(self, node):
        """
        AugAssigned value depend on r-value type dependencies.

        It is valid for subscript, `a[i] += foo()` means `a` type depend on
        `foo` return type and previous a types too.
        """
        args = (self.naming[get_variable(node.target).id],
                self.visit(node.value))
        merge_dep = list({frozenset.union(*x)
                          for x in itertools.product(*args)})
        self.naming[get_variable(node.target).id] = merge_dep

    def visit_For(self, node):
        """
        Handle iterator variable in for loops.

        Iterate variable may be the correct one at the end of the loop.
        """
        body = node.body
        if node.target.id in self.naming:
            body = [ast.Assign(targets=[node.target], value=node.iter,
                               type_comment=None)] + body
            self.visit_any_conditionnal(body, node.orelse)
        else:
            iter_dep = self.visit(node.iter)
            self.naming[node.target.id] = iter_dep
            self.visit_any_conditionnal(body, body + node.orelse)

    def visit_BoolOp(self, node):
        """ Return type may come from any boolop operand. """
        return sum((self.visit(value) for value in node.values), [])

    def visit_BinOp(self, node):
        """ Return type depend from both operand of the binary operation. """
        args = [self.visit(arg) for arg in (node.left, node.right)]
        return list({frozenset.union(*x) for x in itertools.product(*args)})

    def visit_UnaryOp(self, node):
        """ Return type depend on operand only. """
        return self.visit(node.operand)

    @staticmethod
    def visit_Lambda(_):
        """ Lambda have to be remove before. """
        assert False

    def visit_IfExp(self, node):
        """ Return value depend on both if branch. """
        return self.visit(node.body) + self.visit(node.orelse)

    @staticmethod
    def visit_Compare(_):
        """ Comparison return a bool so there are no dependencies. """
        return [frozenset()]

    def visit_Call(self, node):
        """
        Function call depend on all function use in the call.

        >> a = foo(bar(c) or foobar(d))

        Return type depend on [foo, bar] or [foo, foobar]
        """
        args = [self.visit(arg) for arg in node.args]
        func = self.visit(node.func) or []
        params = args + [func]
        return list({frozenset.union(*p) for p in itertools.product(*params)})

    @staticmethod
    def visit_Constant(_):
        """ Return no dependencies on others functions. """
        return [frozenset()]

    @staticmethod
    def visit_Attribute(_):
        """ Return no dependencies on others functions. """
        return [frozenset()]

    def visit_Subscript(self, node):
        """
        Return dependencies of the subscripted value.

        a = foo()[0] means `a` have a dependency on `foo` return type.
        """
        return self.visit(node.value)

    def visit_Name(self, node):
        """
        Return dependencies for given variable.

        It has to be registered first.
        """
        if node.id in self.naming:
            return self.naming[node.id]
        elif node.id in self.global_declarations:
            return [frozenset([self.global_declarations[node.id]])]
        elif isinstance(node.ctx, ast.Param):
            deps = [frozenset()]
            self.naming[node.id] = deps
            return deps
        else:
            raise PythranInternalError("Variable '{}' used before assignment"
                                       "".format(node.id))

    def visit_List(self, node):
        """ List construction depend on each elements type dependency. """
        if node.elts:
            return list(set(sum([self.visit(elt) for elt in node.elts], [])))
        else:
            return [frozenset()]

    visit_Set = visit_List

    def visit_Dict(self, node):
        """ Dict construction depend on each element/value type dependency."""
        if node.keys:
            items = node.keys + node.values
            return list(set(sum([self.visit(item) for item in items], [])))
        else:
            return [frozenset()]

    visit_Tuple = visit_List

    @staticmethod
    def visit_Slice(_):
        """ Slice are not part of return type dependency information. """
        assert False

    @staticmethod
    def visit_Index(_):
        """ Index are not part of return type dependency information. """
        assert False

    def visit_If(self, node):
        """ Both if branches may be evaluate first. """
        return self.visit_any_conditionnal(node.body, node.orelse)

    def visit_While(self, node):
        """ Both while branches may be evaluate first. """
        return self.visit_any_conditionnal(node.body, node.orelse)

    def visit_ExceptHandler(self, node):
        """ Exception may declare a new variable. """
        if node.name:
            self.naming[node.name.id] = [frozenset()]
        for stmt in node.body:
            self.visit(stmt)
