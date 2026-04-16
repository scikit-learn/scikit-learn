""" NormalizeMethodCalls turns built in method calls into function calls. """

from pythran.analyses import Globals, Ancestors
from pythran.passmanager import Transformation
from pythran.syntax import PythranSyntaxError
from pythran.tables import attributes, functions, methods, MODULES
from pythran.tables import duplicated_methods
from pythran.conversion import mangle, demangle
from pythran.utils import isstr

import gast as ast
from functools import reduce


class NormalizeMethodCalls(Transformation[Globals, Ancestors]):
    '''
    Turns built in method calls into function calls.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("[].append(12)")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(NormalizeMethodCalls, node)
    >>> print(pm.dump(backend.Python, node))
    import __dispatch__ as __pythran_import___dispatch__
    __pythran_import___dispatch__.append([], 12)
    '''

    def __init__(self):
        super().__init__()
        self.imports = {'builtins': 'builtins',
                        mangle('__dispatch__'): '__dispatch__'}
        self.to_import = set()

    def visit_Module(self, node):
        """
            When we normalize call, we need to add correct import for method
            to function transformation.

            a.max()

            for numpy array will become:

            numpy.max(a)

            so we have to import numpy.
        """
        self.skip_functions = True
        self.generic_visit(node)
        self.skip_functions = False
        self.generic_visit(node)
        new_imports = self.to_import - self.globals
        imports = [ast.Import(names=[ast.alias(name=mod[17:], asname=mod)])
                   for mod in new_imports]
        node.body = imports + node.body
        self.update |= bool(imports)
        return node

    def visit_FunctionDef(self, node):
        if self.skip_functions:
            return node
        old_imports = self.imports
        self.imports = old_imports.copy()
        for arg in node.args.args:
            self.imports.pop(arg.id, None)
        self.generic_visit(node)

        self.imports = old_imports

        return node

    def visit_Import(self, node):
        for alias in node.names:
            name = alias.asname or alias.name
            self.imports[name] = name
        return node

    def visit_Assign(self, node):
        # aliasing between modules
        if isinstance(node.value, ast.Name) and node.value.id in self.imports:
            for t in node.targets:
                if isinstance(t, ast.Name):
                    self.imports[t.id] = self.imports[node.value.id]
            return None
        else:
            n = self.generic_visit(node)
            targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
            for t in targets:
                if isinstance(t, ast.Name):
                    self.imports.pop(t.id, None)
            return n
    visit_AnnAssign = visit_Assign

    def visit_For(self, node):
        node.iter = self.visit(node.iter)
        if isinstance(node.target, ast.Name):
            self.imports.pop(node.target.id, None)
        if node.body:
            node.body = [self.visit(n) for n in node.body]
        if node.orelse:
            node.orelse = [self.visit(n) for n in node.orelse]
        return node

    def baseobj(self, obj):
        # Get the most left identifier
        while isinstance(obj, ast.Attribute):
            obj = obj.value

        # Check if it's a module
        if isinstance(obj, ast.Name) and obj.id in self.imports:
            return None
        else:
            return obj

    def keyword_based_disambiguification(self, node):
        assert isinstance(node.func, ast.Attribute)
        if getattr(node.func.value, 'id', None) != mangle('__dispatch__'):
            return
        if not node.keywords:
            return
        if node.func.attr not in duplicated_methods:
            return

        node_keywords = {kw.arg for kw in node.keywords}
        for disamb_path, disamb_node in duplicated_methods[node.func.attr]:
            disamb_args = {arg.id for arg in disamb_node.args.args}
            if all(kw in disamb_args for kw in node_keywords):
                node.func = self.attr_to_func(node.func, disamb_path)
                return


    def attr_to_func(self, node, mod=None):
        if mod is None:
            mod = methods[node.attr][0]
        # Submodules import full module
        self.to_import.add(mangle(mod[0]))
        func = reduce(
            lambda v, o: ast.Attribute(v, o, ast.Load()),
            mod[1:] + (node.attr,),
            ast.Name(mangle(mod[0]), ast.Load(), None, None)
            )
        return func

    def visit_Attribute(self, node):
        node = self.generic_visit(node)
        # method name -> not a getattr
        if node.attr in methods:

            # Make sure parent is'nt a call, it's already handled in visit_Call
            for parent in reversed(self.ancestors.get(node, ())):
                if isinstance(parent, ast.Attribute):
                    continue
                if isinstance(parent, ast.Call):
                    return node
                break

            # we have a bound method which is not a call
            obj = self.baseobj(node)
            if obj is not None:
                self.update = True
                mod = methods[node.attr][0]
                self.to_import.add(mangle(mod[0]))
                self.to_import.add(mangle('functools'))
                func = self.attr_to_func(node)
                z = ast.Call(
                    ast.Attribute(
                        ast.Name(mangle('functools'), ast.Load(), None, None),
                        "partial",
                        ast.Load()
                        ),
                    [func, obj],
                    [])
                return z
            else:
                return node
        # imported module -> not a getattr
        elif (isinstance(node.value, ast.Name) and
              node.value.id in self.imports):
            module_id = self.imports[node.value.id]
            if node.attr not in MODULES[self.renamer(module_id, MODULES)[1]]:
                msg = ("`" + node.attr + "' is not a member of " +
                       demangle(module_id) + " or Pythran does not support it")
                raise PythranSyntaxError(msg, node)
            node.value.id = module_id  # patch module aliasing
            self.update = True
            return node
        # not listed as attributed -> not a getattr
        elif node.attr not in attributes:
            return node
        # A getattr !
        else:
            self.update = True
            call = ast.Call(
                ast.Attribute(
                    ast.Name('builtins', ast.Load(), None, None),
                    'getattr',
                    ast.Load()),
                [node.value, ast.Constant(node.attr, None)],
                [])
            if isinstance(node.ctx, ast.Store):
                # the only situation where this arises is for real/imag of
                # a ndarray. As a call is not valid for a store, add a slice
                # to ends up with a valid lhs
                if node.attr not in ('real', 'imag'):
                    raise PythranSyntaxError(
                            "Unsupported store to attribute {}".format(node.attr),
                            node)

                return ast.Subscript(call,
                                     ast.Slice(None, None, None),
                                     node.ctx)
            else:
                return call

    @staticmethod
    def renamer(v, cur_module):
        """
        Rename function path to fit Pythonic naming.
        """
        mname = demangle(v)
        return v, mname

    def visit_Call(self, node):
        """
        Transform call site to have normal function call.

        Examples
        --------
        For methods:
        >> a = [1, 2, 3]

        >> a.append(1)

        Becomes

        >> __list__.append(a, 1)


        For functions:
        >> builtins.dict.fromkeys([1, 2, 3])

        Becomes

        >> builtins.__dict__.fromkeys([1, 2, 3])
        """
        node = self.generic_visit(node)

        # Only attributes function can be Pythonic and should be normalized
        if isinstance(node.func, ast.Attribute):
            if node.func.attr in methods:
                # Check object targeted by methods
                if self.baseobj(node.func) is not None:
                    self.update = True
                    # As it was a methods call, push targeted object as first
                    # arguments and add correct module prefix
                    node.args.insert(0, node.func.value)
                    mod = methods[node.func.attr][0]
                    # Submodules import full module
                    self.to_import.add(mangle(mod[0]))
                    node.func = self.attr_to_func(node.func)
                # else methods have been called using function syntax

            ismethod = node.func.attr in methods
            isfunction = node.func.attr in functions
            if ismethod or isfunction:  # i.e. if it's callable
                if not ismethod and self.baseobj(node.func.value):
                    raise PythranSyntaxError(
                            "'{}' called as a method while it's a function"
                            .format(node.func.attr),
                            node)

                # Now, methods and function have both function syntax
                def rec(path, cur_module):
                    """
                    Recursively rename path content looking in matching module.

                    Prefers __module__ to module if it exists.
                    This recursion is done as modules are visited top->bottom
                    while attributes have to be visited bottom->top.
                    """
                    err = "Function path is chained attributes and name"
                    assert isinstance(path, (ast.Name, ast.Attribute)), err
                    if isinstance(path, ast.Attribute):
                        new_node, cur_module = rec(path.value, cur_module)
                        new_id, mname = self.renamer(path.attr, cur_module)
                        return (ast.Attribute(new_node, new_id, ast.Load()),
                                cur_module[mname])
                    else:
                        new_id, mname = self.renamer(path.id, cur_module)
                        return (ast.Name(new_id, ast.Load(), None, None),
                                cur_module[mname])

                # Rename module path to avoid naming issue.
                node.func.value, _ = rec(node.func.value, MODULES)
                self.update = True
            self.keyword_based_disambiguification(node)

        return node

    def visit_BinOp(self, node):
        # replace "str" % (...) by builtins.str.__mod__(...)
        # the reason why we do this is that % formatting is handled by
        # a third party library that's relatively costly to load, so using a
        # function name instead of an operator overload makes it possible to
        # load it only when needed. The drawback is that % formatting is no
        # longer supported when lhs is not a literal
        self.generic_visit(node)
        if isinstance(node.op, ast.Mod) and isstr(node.left):
            self.update = True
            return ast.Call(
                ast.Attribute(
                    ast.Attribute(
                        ast.Name('builtins', ast.Load(), None, None),
                        'str',
                        ast.Load()),
                    '__mod__',
                    ast.Load()),
                [node.left, node.right],
                [])
        return node
