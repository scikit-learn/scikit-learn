""" Inlining inline functions body. """

from pythran.analyses import Inlinable, Aliases
from pythran.passmanager import Transformation

import gast as ast
import copy


class Inlining(Transformation[Inlinable, Aliases]):

    """
    Inline one line functions.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> pm = passmanager.PassManager("test")
    >>> node = ast.parse('''
    ... def foo(a, b):
    ...     return b + b * a
    ... def bar(b):
    ...     return foo(2 * b, b) * foo(b, b)''')
    >>> _, node = pm.apply(Inlining, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(a, b):
        return (b + (b * a))
    def bar(b):
        __pythran_inlinefooa0 = (2 * b)
        __pythran_inlinefoob0 = b
        __pythran_inlinefooa1 = b
        __pythran_inlinefoob1 = b
        return ((__pythran_inlinefoob0 + (__pythran_inlinefoob0 * \
__pythran_inlinefooa0)) * (__pythran_inlinefoob1 + \
(__pythran_inlinefoob1 * __pythran_inlinefooa1)))
    """

    def __init__(self):
        """ fun : Function {name :body} for inlinable functions. """
        super().__init__()
        self.update = False
        self.defs = list()
        self.call_count = 0

    def visit_Stmt(self, node):
        """ Add new variable definition before the Statement. """
        save_defs, self.defs = self.defs or list(), list()
        self.generic_visit(node)
        new_defs, self.defs = self.defs, save_defs
        return new_defs + [node]

    visit_Return = visit_Stmt
    visit_Assign = visit_Stmt
    visit_AnnAssign = visit_Stmt
    visit_AugAssign = visit_Stmt
    visit_Print = visit_Stmt
    visit_For = visit_Stmt
    visit_If = visit_Stmt
    visit_With = visit_Stmt
    visit_Assert = visit_Stmt
    visit_Expr = visit_Stmt

    def visit_While(self, node):
        # FIXME: we're only preventing inlining within test because it's
        # difficult to compute predecessors of the test while also handling
        # exceptions. We could use the cfg analysis for this but I'm a bit lazy
        # and it's not a critical optimization.
        test, node.test = node.test, None
        self.generic_visit(node)
        node.test = test
        return node

    def visit_Call(self, node):
        """
        Replace function call by inlined function's body.

        We can inline if it aliases on only one function.
        """
        func_aliases = self.aliases[node.func]
        if len(func_aliases) == 1:
            function_def = next(iter(func_aliases))
            if (isinstance(function_def, ast.FunctionDef) and
                    function_def.name in self.inlinable):
                self.update = True
                to_inline = copy.deepcopy(self.inlinable[function_def.name])
                arg_to_value = dict()
                values = node.args
                values += to_inline.args.defaults[len(node.args) -
                                                  len(to_inline.args.args):]
                for arg_fun, arg_call in zip(to_inline.args.args, values):
                    v_name = "__pythran_inline{}{}{}".format(function_def.name,
                                                             arg_fun.id,
                                                             self.call_count)
                    new_var = ast.Name(id=v_name,
                                       ctx=ast.Store(),
                                       annotation=None, type_comment=None)
                    self.defs.append(ast.Assign(targets=[new_var],
                                                value=arg_call,
                                                type_comment=None))
                    arg_to_value[arg_fun.id] = ast.Name(id=v_name,
                                                        ctx=ast.Load(),
                                                        annotation=None,
                                                        type_comment=None)

                self.call_count += 1
                return Inliner(arg_to_value).visit(to_inline.body[0])
        return node


class Inliner(ast.NodeTransformer):

    """ Helper transform that performed inlined body transformation. """

    def __init__(self, match):
        """ match : {original_variable_name : Arguments use on call}. """
        self.match = match
        super(Inliner, self).__init__()

    def visit_Name(self, node):
        """ Transform name from match values if available. """
        return self.match.get(node.id, node)

    def visit_Return(self, node):
        """ Remove return keyword after inline. """
        return self.visit(node.value)
