""" ExpandImports replaces imports by their full paths. """

from pythran.passmanager import Transformation
from pythran.utils import path_to_attr, path_to_node
from pythran.conversion import mangle
from pythran.syntax import PythranSyntaxError
from pythran.analyses import Ancestors

import gast as ast


class ExpandImports(Transformation[Ancestors]):

    """
    Expands all imports into full paths.

    Attributes
    ----------
    imports : {str}
        Imported module (python base module name)
    symbols : {str : (str,)}
        Matching between used name and real cxx name.

    Examples
    --------
    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("from math import cos ; cos(2)")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ExpandImports, node)
    >>> print(pm.dump(backend.Python, node))
    import math as __pythran_import_math
    __pythran_import_math.cos(2)
    >>> node = ast.parse("from os.path import join ; join('a', 'b')")
    >>> _, node = pm.apply(ExpandImports, node)
    >>> print(pm.dump(backend.Python, node))
    import os as __pythran_import_os
    __pythran_import_os.path.join('a', 'b')
    """

    def __init__(self):
        super().__init__()
        self.imports = set()
        self.symbols = dict()

    def visit_Module(self, node):
        """
        Visit the whole module and add all import at the top level.

        >> import numpy.linalg

        Becomes

        >> import numpy

        """
        node.body = [k for k in (self.visit(n) for n in node.body) if k]
        imports = [ast.Import([ast.alias(i, mangle(i))]) for i in self.imports]
        node.body = imports + node.body
        ast.fix_missing_locations(node)
        return node

    def visit_Import(self, node):
        """ Register imported modules and usage symbols.  """
        for alias in node.names:
            alias_name = tuple(alias.name.split('.'))
            self.imports.add(alias_name[0])
            if alias.asname:
                self.symbols[alias.asname] = alias_name
            else:
                self.symbols[alias_name[0]] = alias_name[:1]
            self.update = True
        return None

    def visit_ImportFrom(self, node):
        """ Register imported modules and usage symbols.  """
        module_path = tuple(node.module.split('.'))
        self.imports.add(module_path[0])
        for alias in node.names:
            path = module_path + (alias.name,)
            self.symbols[alias.asname or alias.name] = path
        self.update = True
        return None

    def visit_FunctionDef(self, node):
        """
        Update import context using overwriting name information.

        Examples
        --------
        >> import foo
        >> import bar
        >> def foo(bar):
        >>     print(bar)

        In this case, neither bar nor foo can be used in the foo function and
        in future function, foo will not be usable.
        """
        self.symbols.pop(node.name, None)
        gsymbols = self.symbols.copy()
        [self.symbols.pop(arg.id, None) for arg in node.args.args]
        self.generic_visit(node)
        self.symbols = gsymbols
        return node

    def visit_Assign(self, node):
        """
        Update import context using overwriting name information.

        Examples
        --------
        >> import foo
        >> def bar():
        >>     foo = 2
        >>     print(foo)

        In this case, foo can't be used after assign.
        """
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        if isinstance(node.value, ast.Name) and node.value.id in self.symbols:
            symbol = path_to_node(self.symbols[node.value.id])
            if not getattr(symbol, 'isliteral', lambda: False)():
                for target in targets:
                    if not isinstance(target, ast.Name):
                        err = "Unsupported module aliasing"
                        raise PythranSyntaxError(err, target)
                    self.symbols[target.id] = self.symbols[node.value.id]
                return None  # this assignment is no longer needed
        new_node = self.generic_visit(node)
        ntargets = new_node.targets if isinstance(node, ast.Assign) else (new_node.target,)
        # no problem if targets contains a subscript, it is not a new assign.
        [self.symbols.pop(t.id, None)
         for t in ntargets if isinstance(t, ast.Name)]
        return new_node
    visit_AnnAssign = visit_Assign

    def visit_Name(self, node):
        """
        Replace name with full expanded name.

        Examples
        --------
        >> from numpy.linalg import det

        >> det(a)

        Becomes

        >> numpy.linalg.det(a)
        """
        self.generic_visit(node)
        if node.id in self.symbols:
            symbol = path_to_node(self.symbols[node.id])
            if not getattr(symbol, 'isliteral', lambda: False)():
                parent = self.ancestors[node][-1]
                blacklist = (ast.Tuple,
                             ast.List,
                             ast.Set,
                             ast.Return)
                if isinstance(parent, blacklist):
                    raise PythranSyntaxError(
                        "Unsupported module identifier manipulation",
                        node)
            new_node = path_to_attr(self.symbols[node.id])
            new_node.ctx = node.ctx
            ast.copy_location(new_node, node)
            return new_node
        return node
