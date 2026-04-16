
""" Comprehension patterns transforms list comprehension into intrinsics.  """

from pythran.analyses import OptimizableComprehension
from pythran.passmanager import Transformation
from pythran.transformations.normalize_tuples import ConvertToTuple
from pythran.conversion import mangle
from pythran.utils import attr_to_path, path_to_attr

import gast as ast


class ComprehensionPatterns(Transformation[OptimizableComprehension]):
    '''
    Transforms list comprehension into intrinsics.
    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(y) : return (x for x in y)")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ComprehensionPatterns, node)
    >>> 'map' in pm.dump(backend.Python, node)
    True

    >>> node = ast.parse("def foo(y) : return [0 for _ in builtins.range(y)]")
    >>> _, node = pm.apply(ComprehensionPatterns, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(y):
        return ([0] * builtins.len(builtins.range(y)))
    '''

    def visit_Module(self, node):
        self.use_itertools = False
        self.generic_visit(node)
        if self.use_itertools:
            import_alias = ast.alias(name='itertools',
                                     asname=mangle('itertools'))
            importIt = ast.Import(names=[import_alias])
            node.body.insert(0, importIt)
        return node

    def make_Iterator(self, gen):
        if gen.ifs:
            ldFilter = ast.Lambda(
                ast.arguments([ast.Name(gen.target.id, ast.Param(),
                                        None, None)],
                              [], None, [], [], None, []),
                ast.BoolOp(ast.And(), gen.ifs)
                if len(gen.ifs) > 1 else gen.ifs[0])
            ifilterName = ast.Attribute(
                value=ast.Name(id='builtins',
                               ctx=ast.Load(),
                               annotation=None, type_comment=None),
                attr='filter', ctx=ast.Load())
            return ast.Call(ifilterName, [ldFilter, gen.iter], [])
        else:
            return gen.iter

    def visitComp(self, node, make_attr):

        if node in self.optimizable_comprehension:
            self.update = True
            self.generic_visit(node)

            iters = [self.make_Iterator(gen) for gen in node.generators]
            variables = [ast.Name(gen.target.id, ast.Param(), None, None)
                         for gen in node.generators]

            # If dim = 1, product is useless
            if len(iters) == 1:
                iterAST = iters[0]
                varAST = ast.arguments([variables[0]], [],
                                       None, [], [], None, [])
            else:
                self.use_itertools = True
                prodName = ast.Attribute(
                    value=ast.Name(id=mangle('itertools'),
                                   ctx=ast.Load(),
                                   annotation=None, type_comment=None),
                    attr='product', ctx=ast.Load())

                varid = variables[0].id  # retarget this id, it's free
                renamings = {v.id: (i,) for i, v in enumerate(variables)}
                node.elt = ConvertToTuple(varid, renamings).visit(node.elt)
                iterAST = ast.Call(prodName, iters, [])
                varAST = ast.arguments([ast.Name(varid, ast.Param(),
                                                 None, None)],
                                       [], None, [], [], None, [])

            ldBodymap = node.elt
            ldmap = ast.Lambda(varAST, ldBodymap)

            return make_attr(ldmap, iterAST)

        else:
            return self.generic_visit(node)

    def visit_ListComp(self, node):
        def makeattr(*args):
            r = ast.Attribute(
                value=ast.Name(id='builtins',
                               ctx=ast.Load(),
                               annotation=None,
                               type_comment=None),
                attr='map', ctx=ast.Load())
            r = ast.Call(r, list(args), [])
            r = ast.Call(ast.Attribute(ast.Name('builtins', ast.Load(),
                                                None, None),
                                       'list', ast.Load()),
                         [r], [])
            return r

        if isinstance(node.elt, ast.Constant) and len(node.generators) == 1:
            gen = node.generators[0]
            if not gen.ifs and isinstance(gen.iter, ast.Call):
                try:
                    path = attr_to_path(gen.iter.func)[1]
                    range_path = 'pythonic', 'builtins', 'functor', 'range'
                    if path == range_path and len(gen.iter.args) == 1:
                        self.update = True
                        return ast.BinOp(
                            ast.List([node.elt], ast.Load()),
                            ast.Mult(),
                            ast.Call(path_to_attr(('builtins', 'len')),
                                     [gen.iter],
                                     []))
                except TypeError:
                    pass

        return self.visitComp(node, makeattr)

    def visit_GeneratorExp(self, node):
        def makeattr(*args):
            return ast.Call(ast.Attribute(
                value=ast.Name(id='builtins',
                               ctx=ast.Load(),
                               annotation=None, type_comment=None),
                attr='map', ctx=ast.Load()), list(args), [])
        return self.visitComp(node, makeattr)
