""" RemoveLambdas turns lambda into regular functions.  """

from pythran.analyses import GlobalDeclarations, ImportedIds
from pythran.analyses import Check
from pythran.analyses import ExtendedDefUseChains
from pythran.passmanager import Transformation
from pythran.tables import MODULES
from pythran.conversion import mangle

import pythran.metadata as metadata

from copy import copy, deepcopy
import gast as ast

unaryops = {
    ast.Invert: 'invert',
    ast.Not: 'not_',
    ast.UAdd: 'pos',
    ast.USub: 'neg',
}

binaryops = {
    ast.Add: "add",
    ast.Sub: "sub",
    ast.Mult: "mul",
    ast.Div: "truediv",
    ast.Mod: "mod",
    ast.Pow: "pow",
    ast.LShift: "lshift",
    ast.RShift: "rshift",
    ast.BitOr: "or_",
    ast.BitXor: "xor",
    ast.BitAnd: "and_",
    ast.MatMult: "matmul",
    ast.FloorDiv: "floordiv",
}


def issimpleoperator(node):
    if node.args.defaults:
        return None
    body = node.body
    args = node.args.args

    if isinstance(body, ast.UnaryOp) and len(args) == 1:
        if not isinstance(body.operand, ast.Name):
            return None
        return unaryops[type(body.op)]
    if isinstance(body, ast.BinOp) and len(args) == 2:
        if not all(isinstance(op, ast.Name) for op in (body.left, body.right)):
            return None
        if body.left.id != args[0].id or body.right.id != args[1].id:
            return None
        return binaryops[type(body.op)]


def issamelambda(pattern, f1):
    f0, duc = pattern
    if len(f0.args.args) != len(f1.args.args):
        return False
    for arg0, arg1 in zip(f0.args.args, f1.args.args):
        arg0.id = arg1.id
        for u in duc.chains[arg0].users():
            u.node.id = arg1.id
    return Check(f0, {}).visit(f1)


class _LambdaRemover(ast.NodeTransformer):

    def __init__(self, parent, prefix):
        super(_LambdaRemover, self).__init__()
        self.prefix = prefix
        self.parent = parent
        self.patterns = parent.patterns

    def __getattr__(self, attr):
        return getattr(self.parent, attr)

    def visit_Lambda(self, node):
        op = issimpleoperator(node)
        if op is not None:
            if mangle('operator') not in self.global_declarations:
                import_ = ast.Import([ast.alias('operator', mangle('operator'))])
                self.imports.append(import_)
                operator_module = MODULES['operator']
                self.global_declarations[mangle('operator')] = operator_module
            return ast.Attribute(ast.Name(mangle('operator'), ast.Load(),
                                          None, None),
                                 op, ast.Load())

        self.generic_visit(node)
        forged_name = "{0}_lambda{1}".format(
            self.prefix,
            len(self.lambda_functions))


        ii = self.gather(ImportedIds, node)
        ii.difference_update(self.lambda_functions)  # remove current lambdas

        binded_args = [ast.Name(iin, ast.Load(), None, None)
                       for iin in sorted(ii)]
        node.args.args = ([ast.Name(iin, ast.Param(), None, None)
                           for iin in sorted(ii)] +
                          node.args.args)
        for patternname, pattern in self.patterns.items():
            if issamelambda(pattern, node):
                proxy_call = ast.Name(patternname, ast.Load(), None, None)
                break
        else:
            duc = ExtendedDefUseChains()
            nodepattern = deepcopy(node)
            duc.visit(ast.Module([ast.Expr(nodepattern)], []))
            self.patterns[forged_name] = nodepattern, duc

            forged_fdef = ast.FunctionDef(
                forged_name,
                copy(node.args),
                [ast.Return(node.body)],
                [], None, None)
            metadata.add(forged_fdef, metadata.Local())
            self.lambda_functions.append(forged_fdef)
            self.global_declarations[forged_name] = forged_fdef
            proxy_call = ast.Name(forged_name, ast.Load(), None, None)

        if binded_args:
            if MODULES['functools'] not in self.global_declarations.values():
                import_ = ast.Import([ast.alias('functools', mangle('functools'))])
                self.imports.append(import_)
                functools_module = MODULES['functools']
                self.global_declarations[mangle('functools')] = functools_module

            return ast.Call(
                ast.Attribute(
                    ast.Name(mangle('functools'), ast.Load(), None, None),
                    "partial",
                    ast.Load()
                    ),
                [proxy_call] + binded_args,
                [])
        else:
            return proxy_call


class RemoveLambdas(Transformation[GlobalDeclarations]):

    """
    Turns lambda into top-level functions.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(y): lambda x:y+x")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(RemoveLambdas, node)
    >>> print(pm.dump(backend.Python, node))
    import functools as __pythran_import_functools
    def foo(y):
        __pythran_import_functools.partial(foo_lambda0, y)
    def foo_lambda0(y, x):
        return (y + x)
    """

    def visit_Module(self, node):
        self.lambda_functions = list()
        self.patterns = {}
        self.imports = list()
        self.generic_visit(node)
        node.body = self.imports + node.body + self.lambda_functions
        self.update |= bool(self.imports) or bool(self.lambda_functions)
        return node

    def visit_FunctionDef(self, node):
        lr = _LambdaRemover(self, node.name)
        node.body = [lr.visit(n) for n in node.body]
        return node
