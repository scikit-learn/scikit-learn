""" DeadCodeElimination remove useless code. """

from pythran.analyses import PureExpressions, DefUseChains, Ancestors
from pythran.openmp import OMPDirective
from pythran.passmanager import Transformation
import pythran.metadata as metadata

import gast as ast


class ClumsyOpenMPDependencyHandler(ast.NodeVisitor):

    def __init__(self):
        self.blacklist = set()

    def visit_OMPDirective(self, node):
        for dep in node.deps:
            if isinstance(dep, ast.Name):
                self.blacklist.add(dep.id)
        return node


class DeadCodeElimination(Transformation[PureExpressions, DefUseChains, Ancestors]):
    """
    Remove useless statement like:
        - assignment to unused variables
        - remove alone pure statement
        - remove empty if

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> pm = passmanager.PassManager("test")
    >>> node = ast.parse("def foo(): a = [2, 3]; return 1")
    >>> _, node = pm.apply(DeadCodeElimination, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        pass
        return 1
    >>> node = ast.parse("def foo(): 'a simple string'; return 1")
    >>> _, node = pm.apply(DeadCodeElimination, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        pass
        return 1
    >>> node = ast.parse('''
    ... def bar(a):
    ...     return a
    ... def foo(a):
    ...    bar(a)
    ...    return 1''')
    >>> _, node = pm.apply(DeadCodeElimination, node)
    >>> print(pm.dump(backend.Python, node))
    def bar(a):
        return a
    def foo(a):
        pass
        return 1
    """
    def __init__(self):
        super().__init__()
        self.blacklist = set()

    def used_target(self, node):
        if isinstance(node, ast.Name):
            if node.id in self.blacklist:
                return True
            chain = self.def_use_chains.chains[node]
            return bool(chain.users())
        return True

    def visit_FunctionDef(self, node):
        codh = ClumsyOpenMPDependencyHandler()
        codh.visit(node)
        self.blacklist = codh.blacklist
        return self.generic_visit(node)

    def visit_Pass(self, node):
        ancestor = self.ancestors[node][-1]
        if getattr(ancestor, 'body', ()) == [node]:
            return node
        if getattr(ancestor, 'orelse', ()) == [node]:
            return node
        if metadata.get(node, OMPDirective):
            return node
        return None

    def visit_Assign(self, node):
        targets = [target for target in node.targets
                   if self.used_target(target)]
        if len(targets) == len(node.targets):
            return node
        self.update = True
        if targets:
            node.targets = targets
            return node
        if node.value in self.pure_expressions:
            return ast.Pass()
        else:
            return ast.Expr(value=node.value)

    def visit_AnnAssign(self, node):
        if not node.value:
            return node
        if self.used_target(node.target):
            return node
        self.update = True
        if node.value in self.pure_expressions:
            return ast.Pass()
        else:
            return ast.Expr(value=node.value)

    def visit_Expr(self, node):
        if (node.value in self.pure_expressions and
                not isinstance(node.value, ast.Yield)):
            self.update = True
            return ast.Pass()
        self.generic_visit(node)
        return node

    def visit_If(self, node):
        self.generic_visit(node)

        try:
            if ast.literal_eval(node.test):
                if not metadata.get(node, OMPDirective):
                    self.update = True
                    return node.body
            else:
                if not metadata.get(node, OMPDirective):
                    self.update = True
                    return node.orelse
        except ValueError:
            # not a constant expression
            pass

        have_body = any(not isinstance(x, ast.Pass) for x in node.body)
        have_else = any(not isinstance(x, ast.Pass) for x in node.orelse)
        # If the "body" is empty but "else content" is useful, switch branches
        # and remove else content
        if not have_body and have_else:
            test = ast.UnaryOp(op=ast.Not(), operand=node.test)
            self.update = True
            return ast.If(test=test, body=node.orelse, orelse=list())
        # if neither "if" and "else" are useful, keep test if it is not pure
        elif not have_body:
            self.update = True
            if node.test in self.pure_expressions:
                return ast.Pass()
            else:
                node = ast.Expr(value=node.test)
                self.generic_visit(node)
        return node

    def visit(self, node):
        """ Add OMPDirective from the old node to the new one. """
        old_omp = metadata.get(node, OMPDirective)
        node = super(DeadCodeElimination, self).visit(node)
        if not metadata.get(node, OMPDirective):
            for omp_directive in old_omp:
                metadata.add(node, omp_directive)
        return node
