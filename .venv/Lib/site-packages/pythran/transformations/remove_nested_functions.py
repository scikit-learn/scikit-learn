""" RemoveNestedFunctions turns nested function into top-level functions. """

from pythran.analyses import GlobalDeclarations, NonlocalDeclarations, ImportedIds
from pythran.passmanager import Transformation
from pythran.tables import MODULES
from pythran.conversion import mangle

import pythran.metadata as metadata

import gast as ast


class _NestedFunctionRemover(ast.NodeTransformer):
    def __init__(self, parent):
        ast.NodeTransformer.__init__(self)
        self.parent = parent
        self.identifiers = set(self.global_declarations.keys())
        self.boxes = {}
        self.nonlocal_boxes = {}

    def __getattr__(self, attr):
        return getattr(self.parent, attr)

    def visit_Nonlocal(self, node):
        return ast.Pass()

    def visit_FunctionDef(self, node):
        self.update = True
        if MODULES['functools'] not in self.global_declarations.values():
            import_ = ast.Import([ast.alias('functools', mangle('functools'))])
            self.ctx.module.body.insert(0, import_)
            functools_module = MODULES['functools']
            self.global_declarations[mangle('functools')] = functools_module

        self.ctx.module.body.append(node)

        former_name = node.name
        seed = 0
        new_name = "pythran_{}{}"

        while new_name.format(former_name, seed) in self.identifiers:
            seed += 1

        new_name = new_name.format(former_name, seed)
        self.identifiers.add(new_name)

        ii = self.gather(ImportedIds, node)
        sii = sorted(ii)
        binded_args = [ast.Name('__pythran_boxed_' + iin, ast.Load(), None,
                                              None)
                       for iin in sii]
        node.args.args = ([ast.Name(iin, ast.Param(), None, None)
                           for iin in sii] +
                          node.args.args)

        unboxing = []
        nonlocal_boxes = {}
        for iin in sii:
            if iin in self.nonlocal_declarations[node]:
                nonlocal_boxes.setdefault(iin, None)
                self.nonlocal_boxes.setdefault(iin, None)
            else:
                unboxing.append(ast.Assign([ast.Name(iin, ast.Store(), None, None)],
                               ast.Subscript(ast.Name('__pythran_boxed_args_' + iin, ast.Load(), None,
                                                      None),
                                             ast.Constant(None, None),
                                             ast.Load())))
                self.boxes.setdefault(iin, None)
        BoxArgsInserter(nonlocal_boxes).visit(node)

        for arg in node.args.args:
            if arg.id in ii:
                arg.id = '__pythran_boxed_args_' + arg.id

        node.body = unboxing + node.body

        metadata.add(node, metadata.Local())

        class Renamer(ast.NodeTransformer):
            def visit_Call(self, node):
                self.generic_visit(node)
                if (isinstance(node.func, ast.Name) and
                        node.func.id == former_name):
                    node.func.id = new_name
                    node.args = (
                        [ast.Name('__pythran_boxed_args_' + iin, ast.Load(), None, None)
                         for iin in sii] +
                        node.args
                        )
                return node

        Renamer().visit(node)

        node.name = new_name
        self.global_declarations[node.name] = node
        proxy_call = ast.Name(new_name, ast.Load(), None, None)

        new_node = ast.Assign(
            [ast.Name(former_name, ast.Store(), None, None)],
            ast.Call(
                ast.Attribute(
                    ast.Name(mangle('functools'), ast.Load(), None, None),
                    "partial",
                    ast.Load()
                    ),
                [proxy_call] + binded_args,
                [],
                ),
            None)

        nfr = _NestedFunctionRemover(self)
        nfr.remove_nested(node)

        return new_node

    def remove_nested(self, node):
        node.body = [self.visit(stmt) for stmt in node.body]
        if self.update:
            boxes = []
            arg_ids = {arg.id for arg in node.args.args}
            all_boxes = list(self.boxes)
            all_boxes.extend(b for b in self.nonlocal_boxes if b not in
                             self.boxes)
            for i in all_boxes:
                if i in arg_ids:
                    box_value = ast.Dict([ast.Constant(None, None)], [ast.Name(i,
                                                                               ast.Load(),
                                                                               None,
                                                                               None)])
                else:
                    box_value = ast.Dict([], [])
                box = ast.Assign([ast.Name('__pythran_boxed_' + i, ast.Store(),
                                          None, None)], box_value)
                boxes.append(box)

            pre_boxes = self.boxes.copy()
            for k in self.nonlocal_boxes:
                pre_boxes.pop(k, None)

            BoxPreInserter(pre_boxes).visit(node)
            BoxInserter(self.nonlocal_boxes).visit(node)
            node.body = boxes + node.body
        return self.update

class BoxPreInserter(ast.NodeTransformer):

    def __init__(self, insertion_points):
        self.insertion_points = insertion_points

    def insert_target(self, target):
        if getattr(target, 'id', None) not in self.insertion_points:
            return None
        return ast.Assign(
                [ast.Subscript(
                    ast.Name('__pythran_boxed_' + target.id, ast.Load(), None, None),
                             ast.Constant(None, None),
                             ast.Store())],
                ast.Name(target.id, ast.Load(), None, None))


    def visit_Assign(self, node):
        extras = []
        for t in node.targets:
            extra = self.insert_target(t)
            if extra:
                extras.append(extra)
        if extras:
            return [node] + extras
        else:
            return node

    def visit_AugAssign(self, node):
        extra = self.insert_target(node.target)
        if extra:
            return [node, extra]
        else:
            return node


class BoxInserter(ast.NodeTransformer):

    def __init__(self, insertion_points):
        self.insertion_points = insertion_points
        self.prefix = '__pythran_boxed_'

    def visit_Name(self, node):
        if node.id not in self.insertion_points:
            return node
        if isinstance(node.ctx, ast.Param):
            return node
        return ast.Subscript(
                ast.Name(self.prefix + node.id,
                        ast.Load(), None, None),
                ast.Constant(None, None),
                node.ctx)

class BoxArgsInserter(BoxInserter):

    def __init__(self, insertion_points):
        super().__init__(insertion_points)
        self.prefix += 'args_'


class RemoveNestedFunctions(Transformation[GlobalDeclarations,
                                           NonlocalDeclarations]):

    """
    Replace nested function by top-level functions.

    Also add a call to a bind intrinsic that
    generates a local function with some arguments binded.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(x):\\n def bar(y): return x+y\\n bar(12)")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RemoveNestedFunctions, node)
    >>> print(pm.dump(backend.Python, node))
    import functools as __pythran_import_functools
    def foo(x):
        __pythran_boxed_x = {None: x}
        bar = __pythran_import_functools.partial(pythran_bar0, __pythran_boxed_x)
        bar(12)
    def pythran_bar0(__pythran_boxed_args_x, y):
        x = __pythran_boxed_args_x[None]
        return (x + y)
    """

    def visit_Module(self, node):
        # keep original node as it's updated by _NestedFunctionRemover
        for stmt in node.body:
            self.visit(stmt)
        return node

    def visit_FunctionDef(self, node):
        nfr = _NestedFunctionRemover(self)
        self.update |= nfr.remove_nested(node)
        return node
