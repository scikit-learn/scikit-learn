""" Aliases gather aliasing informations. """

from pythran.analyses.global_declarations import GlobalDeclarations
from pythran.intrinsic import Intrinsic, Class, UnboundValue
from pythran.passmanager import ModuleAnalysis
from pythran.tables import functions, methods, MODULES
from pythran.unparse import Unparser
from pythran.conversion import demangle
import pythran.metadata as md
from pythran.utils import isnum
from pythran.syntax import PythranSyntaxError

import gast as ast
from copy import deepcopy
from itertools import product
import io
import numpy as np


IntrinsicAliases = dict()

class UnboundIdentifierError(RuntimeError):
    pass


class DictLayer(object):
    def __init__(self, base, layer=None):
        self.base = base
        self.layer = dict() if layer is None else layer

    def __getitem__(self, key):
        return (self.layer if key in self.layer else self.base)[key]

    def __contains__(self, key):
        return (key in self.layer) or (key in self.base)

    def __setitem__(self, key, value):
        self.layer[key] = value

    def get(self, key, default=None):
        return self.layer.get(key, self.base.get(key, default))


class ContainerOf(object):
    '''
    Represents a container of something

    We just know that if indexed by the integer value `index',
    we get `containees'
    '''
    UnknownIndex = float('nan')

    __slots__ = 'index', 'containees'

    def __init__(self, containees, index=UnknownIndex):
        self.index = index
        self.containees = containees


def save_intrinsic_alias(module):
    """ Recursively save default aliases for pythonic functions. """
    for v in module.values():
        if isinstance(v, dict):  # Submodules case
            save_intrinsic_alias(v)
        else:
            IntrinsicAliases[v] = frozenset((v,))
            if isinstance(v, Class):
                save_intrinsic_alias(v.fields)


for module in MODULES.values():
    save_intrinsic_alias(module)


class Aliases(ModuleAnalysis[GlobalDeclarations]):
    '''
    Gather aliasing informations across nodes

    As a result, each node from the module is associated to a set of node or
    Intrinsic to which it *may* alias to.
    '''

    RetId = '@'

    ResultType = dict

    def __init__(self):
        super().__init__()
        self.aliases = None

    @staticmethod
    def dump(result, filter=None):
        def pp(n):
            output = io.StringIO()
            if isinstance(n, ContainerOf):
                if n.index == n.index:  # always valid except for UnknownIndex
                    output.write('[{}]='.format(n.index))
                containees = sorted(map(pp, n.containees))
                output.write(', '.join(map("|{}|".format, containees)))
            else:
                Unparser(n, output)
            return output.getvalue().strip()

        if isinstance(result, dict):
            for k, v in result.items():
                if (filter is None) or isinstance(k, filter):
                    print('{} => {}'.format(pp(k), sorted(map(pp, v))))
        elif isinstance(result, (frozenset, set)):
            print(sorted(map(pp, result)))

    def get_unbound_value_set(self):
        return {UnboundValue}

    @staticmethod
    def access_path(node):
        if isinstance(node, ast.Name):
            return MODULES.get(demangle(node.id), node.id)
        elif isinstance(node, ast.Attribute):
            attr_key = demangle(node.attr)
            value_dict = Aliases.access_path(node.value)
            if attr_key not in value_dict:
                raise PythranSyntaxError(
                    "Unsupported attribute '{}' for this object"
                    .format(attr_key),
                    node.value)
            return value_dict[attr_key]
        elif isinstance(node, ast.FunctionDef):
            return node.name
        else:
            return node

    # aliasing created by expressions
    def add(self, node, values=None):
        if values is None:  # no given target for the alias
            if isinstance(node, Intrinsic):
                values = {node}  # an Intrinsic always aliases to itself
            else:
                values = self.get_unbound_value_set()
        self.result[node] = values
        return values

    def visit_BoolOp(self, node):
        '''
        Resulting node may alias to either operands:

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b): return a or b')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.BoolOp)
        (a or b) => ['a', 'b']

        Note that a literal does not create any alias

        >>> module = ast.parse('def foo(a, b): return a or 0')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.BoolOp)
        (a or 0) => ['<unbound-value>', 'a']
        '''
        return self.add(node, set.union(*[self.visit(n) for n in node.values]))

    def visit_UnaryOp(self, node):
        '''
        Resulting node does not alias to anything

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a): return -a')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.UnaryOp)
        (- a) => ['<unbound-value>']
        '''
        self.generic_visit(node)
        return self.add(node)

    visit_BinOp = visit_UnaryOp
    visit_Compare = visit_UnaryOp

    def visit_IfExp(self, node):
        '''
        Resulting node alias to either branch

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b, c): return a if c else b')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.IfExp)
        (a if c else b) => ['a', 'b']
        '''
        self.visit(node.test)
        rec = [self.visit(n) for n in (node.body, node.orelse)]
        return self.add(node, set.union(*rec))

    def visit_Dict(self, node):
        '''
        A dict is abstracted as an unordered container of its values

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b): return {0: a, 1: b}')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Dict)
        {0: a, 1: b} => ['|a|, |b|']

        where the |id| notation means something that may contain ``id``.
        '''
        if node.keys:
            elts_aliases = set()
            for key, val in zip(node.keys, node.values):
                self.visit(key)  # res ignored, just to fill self.aliases
                elt_aliases = self.visit(val)
                elts_aliases.update(elt_aliases)
            aliases = {ContainerOf(elts_aliases)}
        else:
            aliases = None
        return self.add(node, aliases)

    def visit_Set(self, node):
        '''
        A set is abstracted as an unordered container of its elements

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b): return {a, b}')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Set)
        {a, b} => ['|a|, |b|']

        where the |id| notation means something that may contain ``id``.
        '''
        if node.elts:
            elts_aliases = {ContainerOf({alias
                            for elt in node.elts
                            for alias in self.visit(elt)})}
        else:
            elts_aliases = None
        return self.add(node, elts_aliases)

    def visit_Return(self, node):
        '''
        A side effect of computing aliases on a Return is that it updates the
        ``return_alias`` field of current function

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b): return a')
        >>> result = pm.gather(Aliases, module)
        >>> module.body[0].return_alias # doctest: +ELLIPSIS
        <function ...merge_return_aliases at...>

        This field is a function that takes as many nodes as the function
        argument count as input and returns an expression based on
        these arguments if the function happens to create aliasing
        between its input and output. In our case:

        >>> f = module.body[0].return_alias
        >>> Aliases.dump(f([ast.Name('A', ast.Load(), None, None),
        ...                 ast.Constant(1, None)]))
        ['A']

        This also works if the relationship between input and output
        is more complex:

        >>> module = ast.parse('def foo(a, b): return a or b[0]')
        >>> result = pm.gather(Aliases, module)
        >>> f = module.body[0].return_alias
        >>> List = ast.List([ast.Name('L0', ast.Load(), None, None)],
        ...                 ast.Load())
        >>> Aliases.dump(f([ast.Name('B', ast.Load(), None, None), List]))
        ['B', '[L0][0]']

        Which actually means that when called with two arguments ``B`` and
        the single-element list ``[L[0]]``, ``foo`` may returns either the
        first argument, or the first element of the second argument.
        '''
        if not node.value:
            return
        ret_aliases = self.visit(node.value)
        if Aliases.RetId in self.aliases:
            ret_aliases = ret_aliases.union(self.aliases[Aliases.RetId])
        self.aliases[Aliases.RetId] = ret_aliases

    def call_return_alias(self, node):

        def interprocedural_aliases(func, args):
            arg_aliases = [self.result[arg] or {arg} for arg in args]
            return_aliases = set()
            for args_combination in product(*arg_aliases):
                for ra in func.return_alias(args_combination):
                    if isinstance(ra, ast.Subscript):
                        if isinstance(ra.value, ContainerOf):
                            if np.isnan(ra.value.index) or (isnum(ra.slice) and
                                                            ra.slice.value ==
                                                            ra.value.index):
                                return_aliases.update(ra.value.containees)
                            continue
                    return_aliases.add(ra)
            return return_aliases

        def full_args(func, call):
            args = call.args
            if isinstance(func, ast.FunctionDef):
                extra = len(func.args.args) - len(args)
                if extra:
                    tail = [deepcopy(n) for n in func.args.defaults[extra:]]
                    for arg in tail:
                        self.visit(arg)
                    args = args + tail
            return args

        func = node.func
        aliases = set()

        if node.keywords:
            # too soon, we don't support keywords in interprocedural_aliases
            pass
        elif isinstance(func, ast.Attribute):
            _, signature = methods.get(func.attr,
                                       functions.get(func.attr,
                                                     [(None, None)])[0])
            if signature:
                args = full_args(signature, node)
                aliases = interprocedural_aliases(signature, args)

        elif isinstance(func, ast.Name):
            func_aliases = self.result[func]
            for func_alias in func_aliases:
                if hasattr(func_alias, 'return_alias'):
                    args = full_args(func_alias, node)
                    aliases.update(interprocedural_aliases(func_alias, args))
                else:
                    pass  # better thing to do ?
        for a in aliases:
            if a not in self.result:
                self.add(a)
        return aliases or self.get_unbound_value_set()

    def visit_Call(self, node):
        '''
        Resulting node alias to the return_alias of called function,
        if the function is already known by Pythran (i.e. it's an Intrinsic)
        or if Pythran already computed it's ``return_alias`` behavior.

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> fun = """
        ... def f(a): return a
        ... def foo(b): c = f(b)"""
        >>> module = ast.parse(fun)

        The ``f`` function create aliasing between
        the returned value and its first argument.

        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Call)
        f(b) => ['b']

        This also works with intrinsics, e.g ``dict.setdefault`` which
        may create alias between its third argument and the return value.

        >>> fun = 'def foo(a, d): builtins.dict.setdefault(d, 0, a)'
        >>> module = ast.parse(fun)
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Call)
        builtins.dict.setdefault(d, 0, a) => ['<unbound-value>', 'a']

        Note that complex cases can arise, when one of the formal parameter
        is already known to alias to various values:

        >>> fun = """
        ... def f(a, b): return a and b
        ... def foo(A, B, C, D): return f(A or B, C or D)"""
        >>> module = ast.parse(fun)
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Call)
        f((A or B), (C or D)) => ['A', 'B', 'C', 'D']
        '''
        self.generic_visit(node)
        f = node.func
        # special handler for bind functions
        if isinstance(f, ast.Attribute) and f.attr == "partial":
            return self.add(node, {node})
        else:
            return_alias = self.call_return_alias(node)
            # expand collected aliases
            all_aliases = set()
            for value in return_alias:
                # no translation
                if isinstance(value, (ContainerOf, ast.FunctionDef,
                                      Intrinsic)):
                    all_aliases.add(value)
                elif value in self.result:
                    all_aliases.update(self.result[value])
                else:
                    try:
                        ap = Aliases.access_path(value)
                        all_aliases.update(self.aliases.get(ap, ()))
                    except NotImplementedError:
                        # should we do something better here?
                        all_aliases.add(value)
            return self.add(node, all_aliases)

    visit_Constant = visit_UnaryOp

    def visit_Attribute(self, node):
        return self.add(node, {Aliases.access_path(node)})

    def visit_Subscript(self, node):
        '''
        Resulting node alias stores the subscript relationship if we don't know
        anything about the subscripted node.

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a): return a[0]')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Subscript)
        a[0] => ['a[0]']

        If we know something about the container, e.g. in case of a list, we
        can use this information to get more accurate informations:

        >>> module = ast.parse('def foo(a, b, c): return [a, b][c]')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Subscript)
        [a, b][c] => ['a', 'b']

        Also work in case of a dict:

        >>> module = ast.parse('def foo(a, b, c): return {a:b}[c]')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Subscript)
        {a: b}[c] => ['b']

        Even when built in several statements:
        >>> module = ast.parse('def foo(a, b, c): d = {} ; d[a] = b; return d[c]')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Subscript)
        d[a] => ['b']
        d[c] => ['b']


        Moreover, in case of a tuple indexed by a constant value, we can
        further refine the aliasing information:

        >>> fun = """
        ... def f(a, b): return a, b
        ... def foo(a, b): return f(a, b)[0]"""
        >>> module = ast.parse(fun)
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Subscript)
        f(a, b)[0] => ['a']

        Nothing is done for slices, even if the indices are known :-/

        >>> module = ast.parse('def foo(a, b, c): return [a, b, c][1:]')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Subscript)
        [a, b, c][1:] => ['<unbound-value>']
        '''
        if isinstance(node.slice, ast.Tuple):
            # could be enhanced through better handling of containers
            self.visit(node.value)
            for elt in node.slice.elts:
                self.visit(elt)
            aliases = None
        else:
            aliases = set()
            self.visit(node.slice)
            value_aliases = self.visit(node.value)
            for alias in value_aliases:
                if isinstance(alias, ContainerOf):
                    if isinstance(node.slice, ast.Slice):
                        continue
                    if isnum(node.slice):
                        if node.slice.value != alias.index:
                            continue
                    # FIXME: what if the index is a slice variable...
                    aliases.update(alias.containees)
                elif isinstance(getattr(alias, 'ctx', None), (ast.Param,
                                                              ast.Store)):
                    aliases.add(ast.Subscript(alias, node.slice, node.ctx))
            if not aliases:
                aliases = None
        return self.add(node, aliases)

    def visit_OMPDirective(self, node):
        '''
        omp directive may introduce new variables, just register them
        '''
        for dep in node.deps:
            self.add(dep)

    def visit_Name(self, node):
        if node.id not in self.aliases:
            raise UnboundIdentifierError(node.id)
        return self.add(node, self.aliases[node.id])

    def visit_Tuple(self, node):
        '''
        A tuple is abstracted as an ordered container of its values

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b): return a, b')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Tuple)
        (a, b) => ['[0]=|a|', '[1]=|b|']

        where the |[i]=id| notation means something that
        may contain ``id`` at index ``i``.
        '''
        if node.elts:
            aliases = set()
            for i, elt in enumerate(node.elts):
                elt_aliases = self.visit(elt)
                aliases.add(ContainerOf(elt_aliases, i))
        else:
            aliases = None
        return self.add(node, aliases)

    visit_List = visit_Set

    def visit_comprehension(self, node):
        self.aliases[node.target.id] = {node.target}
        self.generic_visit(node)

    def visit_ListComp(self, node):
        '''
        A comprehension is not abstracted in any way

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b): return [a for i in b]')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.ListComp)
        [a for i in b] => ['<unbound-value>']
        '''
        for generator in node.generators:
            self.visit_comprehension(generator)
        self.visit(node.elt)
        return self.add(node)

    visit_SetComp = visit_ListComp

    visit_GeneratorExp = visit_ListComp

    def visit_DictComp(self, node):
        '''
        A comprehension is not abstracted in any way

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a, b): return {i: i for i in b}')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.DictComp)
        {i: i for i in b} => ['<unbound-value>']
        '''
        for generator in node.generators:
            self.visit_comprehension(generator)
        self.visit(node.key)
        self.visit(node.value)
        return self.add(node)

    # aliasing created by statements

    def init_function_alias(self, node):
        "each argument is bound to a different identifier"
        self.aliases.update((arg.id, {arg})
                            for arg in node.args.args)

    def visit_FunctionDef(self, node):
        '''
        Initialise aliasing default value before visiting.

        Add aliasing values for :
            - Pythonic
            - globals declarations
            - current function arguments
        '''
        self.aliases = IntrinsicAliases.copy()

        self.aliases.update((k, {v})
                            for k, v in self.global_declarations.items())

        self.init_function_alias(node)

        self.generic_visit(node)
        if Aliases.RetId in self.aliases:
            # parametrize the expression
            def parametrize(exp):
                # constant or global -> no change
                if isinstance(exp, (ast.Constant, Intrinsic, ast.FunctionDef)):
                    return lambda _: {exp}
                elif isinstance(exp, ContainerOf):
                    index = exp.index
                    return lambda args: {
                        ContainerOf({pc for containee in exp.containees for pc
                                     in parametrize(containee)(args)}, index)
                    }
                elif isinstance(exp, ast.Name):
                    try:
                        w = node.args.args.index(exp)

                        def return_alias(args):
                            if w < len(args):
                                return {args[w]}
                            else:
                                return {node.args.defaults[w - len(args)]}
                        return return_alias
                    except ValueError:
                        return lambda _: self.get_unbound_value_set()
                elif isinstance(exp, ast.Subscript):
                    values = parametrize(exp.value)
                    slices = parametrize(exp.slice)
                    return lambda args: {
                        ast.Subscript(value, slice, ast.Load())
                        for value in values(args)
                        for slice in slices(args)}
                else:
                    return lambda _: self.get_unbound_value_set()

            # this is a little tricky: for each returned alias,
            # parametrize builds a function that, given a list of args,
            # returns the alias
            # then as we may have multiple returned alias, we compute the union
            # of these returned aliases
            return_aliases = [parametrize(ret_alias)
                              for ret_alias
                              in self.aliases[Aliases.RetId]]

            def merge_return_aliases(args):
                return {ra
                        for return_alias in return_aliases
                        for ra in return_alias(args)}

            node.return_alias = merge_return_aliases

    def visit_Assert(self, node):
        self.generic_visit(node)

        if not isinstance(node.test, ast.Compare):
            return
        if len(node.test.ops) != 1:
            return
        op = node.test.ops[0]
        comparator = node.test.comparators[0]

        if not isinstance(node.test.left, ast.Name):
            return
        if not isinstance(comparator, ast.Name):
            return

        if isinstance(op, ast.IsNot):
            left_aliases = self.aliases[node.test.left.id]
            right_aliases = self.aliases[comparator.id]
            self.aliases[node.test.left.id] = self.aliases[node.test.left.id].difference(right_aliases)
            right_aliases.difference_update(left_aliases)

    def visit_Assign(self, node):
        r'''
        Assignment creates aliasing between lhs and rhs

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse('def foo(a): c = a ; d = e = c ; {c, d, e}')
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Set)
        {c, d, e} => ['|a|']

        Everyone points to the formal parameter 'a' \o/
        '''
        md.visit(self, node)
        value_aliases = self.visit(node.value) if node.value else {}
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for t in targets:
            if isinstance(t, ast.Name):
                self.aliases[t.id] = set(value_aliases) or {t}
                for alias in list(value_aliases):
                    if isinstance(alias, ast.Name):
                        a_id = alias.id
                        self.aliases[a_id] = self.aliases[a_id].union((t,))
                self.add(t, self.aliases[t.id])
            elif isinstance(t, ast.Subscript):
                def wrap(t, aliases):
                    if isinstance(t, ast.Subscript):
                        alias, wrapped = wrap(t.value, aliases)
                        if isnum(t.slice):
                            return alias, {ContainerOf(wrapped, t.slice.value)}
                        elif isinstance(t.slice, ast.Slice):
                            return alias, wrapped
                        else:
                            return alias, {ContainerOf(wrapped)}
                    elif isinstance(t, ast.Name):
                        return t, aliases
                    else:
                        raise NotImplementedError
                try:
                    alias, wrapped = wrap(t, value_aliases)
                    self.aliases[alias.id] = self.aliases[alias.id].union(wrapped)
                except NotImplementedError:
                    ...
                self.visit(t)
                self.add(t, value_aliases)
            else:
                self.visit(t)

    def visit_AnnAssign(self, node):
        self.visit_Assign(node)
        self.visit(node.annotation)

    def visit_For(self, node):
        '''
        For loop creates aliasing between the target
        and the content of the iterator

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> module = ast.parse("""
        ... def foo(a):
        ...     for i in a:
        ...         {i}""")
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Set)
        {i} => ['|i|']

        Not very useful, unless we know something about the iterated container

        >>> module = ast.parse("""
        ... def foo(a, b):
        ...     for i in [a, b]:
        ...         {i}""")
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Set)
        {i} => ['|a|, |b|']
        '''

        iter_aliases = self.visit(node.iter)
        if all(isinstance(x, ContainerOf) for x in iter_aliases):
            target_aliases = set().union(*[iter_alias.containees
                                           for iter_alias in
                              iter_aliases])
        else:
            target_aliases = {node.target}

        self.add(node.target, target_aliases)
        self.aliases[node.target.id] = self.result[node.target]

        self.generic_visit(node)
        self.generic_visit(node)

    def visit_While(self, node):
        '''

        While statement evaluation is somehow equivalent to the evaluation of a
        sequence, except the fact that in some subtle cases, the first rounds
        of analyse fails because we do not follow the regular execution order

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> fun = """
        ... def foo(a):
        ...     while(a):
        ...         if a == 1: builtins.print(b)
        ...         else: b = a"""
        >>> module = ast.parse(fun)
        >>> result = pm.gather(Aliases, module)
        '''
        self.generic_visit(node)
        self.generic_visit(node)

    def visit_If(self, node):
        '''
        After an if statement, the values from both branches are merged,
        potentially creating more aliasing:

        >>> from pythran import passmanager
        >>> pm = passmanager.PassManager('demo')
        >>> fun = """
        ... def foo(a, b):
        ...     if a: c=a
        ...     else: c=b
        ...     return {c}"""
        >>> module = ast.parse(fun)
        >>> result = pm.gather(Aliases, module)
        >>> Aliases.dump(result, filter=ast.Set)
        {c} => ['|a|, |b|']
        '''

        md.visit(self, node)
        self.visit(node.test)
        true_aliases = false_aliases = None

        # first try the true branch
        try:
            self.aliases = DictLayer(self.aliases)
            for stmt in node.body:
                self.visit(stmt)
            true_aliases = self.aliases
        except UnboundIdentifierError:
            pass
        finally:
            self.aliases = self.aliases.base

        # then try the false branch
        try:
            self.aliases = DictLayer(self.aliases)
            for stmt in node.orelse:
                self.visit(stmt)
            false_aliases = self.aliases
        except UnboundIdentifierError:
            pass
        finally:
            self.aliases = self.aliases.base

        if true_aliases and not false_aliases:
            self.aliases = DictLayer(true_aliases)
            for stmt in node.orelse:
                self.visit(stmt)
            false_aliases = self.aliases
            self.aliases = self.aliases.base

        if false_aliases and not true_aliases:
            self.aliases = DictLayer(false_aliases)
            for stmt in node.body:
                self.visit(stmt)
            true_aliases = self.aliases
            self.aliases = self.aliases.base

        # merge the results from true and false branches
        for k, v in true_aliases.layer.items():
            if k in false_aliases.layer:
                self.aliases[k] = v.union(false_aliases[k])
            elif k in self.aliases:
                self.aliases[k] = v.union(self.aliases[k])
            else:
                self.aliases[k] = v
        for k, v in false_aliases.layer.items():
            if k in true_aliases.layer:
                pass  # already done
            elif k in self.aliases:
                self.aliases[k] = v.union(self.aliases[k])
            else:
                self.aliases[k] = v

    def visit_ExceptHandler(self, node):
        if node.name:
            self.aliases[node.name.id] = {node.name}
        self.generic_visit(node)


class StrictAliases(Aliases):
    """
    Gather aliasing informations across nodes,
    without adding unsure aliases.
    """

    def get_unbound_value_set(self):
        return set()


class InterproceduralAliases(Aliases):
    """
    Gather aliases while assuming two different parameters can point to the same
    value
    """

    def init_function_alias(self, node):
        self.aliases.update((arg.id, set(node.args.args))
                            for arg in node.args.args)
