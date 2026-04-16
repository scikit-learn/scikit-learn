"""
ExpandGlobals replaces globals variables by function call.

It also turn globals assignment in function definition.
"""

from pythran.analyses import LocalNameDeclarations
from pythran.passmanager import Transformation
from pythran.syntax import PythranSyntaxError
from pythran.utils import path_to_attr
from pythran import metadata

import gast as ast

class GlobalTransformer(ast.NodeTransformer):
    '''
    Use assumptions on globals to improve code generation
    '''

    def visit_Call(self, node):
        # because a list can be a call parameter during global init
        return node

    def visit_List(self, node):
        # because global lists in pythran are static lists
        return ast.Call(path_to_attr(('builtins', 'pythran', 'static_list')),
                        [ast.Tuple([self.visit(elt) for elt in node.elts],
                                  ast.Load())],
                        [])


class ExpandGlobals(Transformation):

    """
    Expands all builtins into full paths.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('''
    ... a = 1
    ... def foo():
    ...     return a''')
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ExpandGlobals, node)
    >>> print(pm.dump(backend.Python, node))
    def a():
        return 1
    def foo():
        return a()
    """

    def __init__(self):
        """ Initialize local declaration and constant name to expand. """
        super().__init__()
        self.local_decl = set()
        self.to_expand = set()

    def visit_Module(self, node):
        """Turn globals assignment to functionDef and visit function defs. """
        module_body = list()
        symbols = set()
        # Gather top level assigned variables.
        for stmt in node.body:
            if isinstance(stmt, (ast.Import, ast.ImportFrom)):
                for alias in stmt.names:
                    name = alias.asname or alias.name
                    symbols.add(name)  # no warning here
            elif isinstance(stmt, ast.FunctionDef):
                if stmt.name in symbols:
                    raise PythranSyntaxError(
                        "Multiple top-level definition of %s." % stmt.name,
                        stmt)
                else:
                    symbols.add(stmt.name)

            if not isinstance(stmt, ast.Assign):
                continue

            for target in stmt.targets:
                if not isinstance(target, ast.Name):
                    raise PythranSyntaxError(
                        "Top-level assignment to an expression.",
                        target)
                if target.id in self.to_expand:
                    raise PythranSyntaxError(
                        "Multiple top-level definition of %s." % target.id,
                        target)
                if isinstance(stmt.value, ast.Name):
                    if stmt.value.id in symbols:
                        continue  # create aliasing between top level symbols
                self.to_expand.add(target.id)

        for stmt in node.body:
            if isinstance(stmt, ast.Assign):
                # that's not a global var, but a module/function aliasing
                if all(isinstance(t, ast.Name) and t.id not in self.to_expand
                       for t in stmt.targets):
                    module_body.append(stmt)
                    continue

                self.local_decl = set()
                cst_value = GlobalTransformer().visit(self.visit(stmt.value))
                for target in stmt.targets:
                    assert isinstance(target, ast.Name)
                    module_body.append(
                        ast.FunctionDef(target.id,
                                        ast.arguments([], [], None, [],
                                                      [], None, []),
                                        [ast.Return(value=cst_value)],
                                        [], None, None))
                    metadata.add(module_body[-1].body[0],
                                 metadata.StaticReturn())
            else:
                self.local_decl = self.gather(
                    LocalNameDeclarations, stmt)
                module_body.append(self.visit(stmt))

        self.update |= bool(self.to_expand)

        node.body = module_body
        return node

    def visit_Name(self, node):
        """
        Turn global variable used not shadows to function call.

        We check it is a name from an assignment as import or functions use
        should not be turn into call.
        """
        if (isinstance(node.ctx, ast.Load) and
                node.id not in self.local_decl and
                node.id in self.to_expand):
            self.update = True
            return ast.Call(func=node,
                            args=[], keywords=[])
        return node
