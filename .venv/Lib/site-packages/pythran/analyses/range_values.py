""" Module Analysing code to extract positive subscripts from code.  """
# TODO check bound of while and if for more accurate values.

import gast as ast
from collections import defaultdict
from functools import reduce

from pythran.analyses import Aliases, CFG
from pythran.analyses.use_omp import UseOMP
from pythran.intrinsic import Intrinsic
from pythran.passmanager import ModuleAnalysis
from pythran.interval import Interval, IntervalTuple, UNKNOWN_RANGE
from pythran.tables import MODULES, attributes


class UnsupportedExpression(NotImplementedError):
    pass


class RangeValueTooCostly(RuntimeError):
    pass


def combine(op, node0, node1):
    key = '__{}__'.format(op.__class__.__name__.lower())
    try:
        return getattr(type(node0), key)(node0, node1)
    except AttributeError:
        return UNKNOWN_RANGE


def negate(node):
    if isinstance(node, ast.Name):
        # Not type info, could be anything :(
        raise UnsupportedExpression()

    if isinstance(node, ast.UnaryOp):
        # !~x <> ~x == 0 <> x == ~0 <> x == -1
        if isinstance(node.op, ast.Invert):
            return ast.Compare(node.operand,
                               [ast.Eq()],
                               [ast.Constant(-1, None)])
        # !!x <> x
        if isinstance(node.op, ast.Not):
            return node.operand
        # !+x <> +x == 0 <> x == 0 <> !x
        if isinstance(node.op, ast.UAdd):
            return node.operand
        # !-x <> -x == 0 <> x == 0 <> !x
        if isinstance(node.op, ast.USub):
            return node.operand

    if isinstance(node, ast.BoolOp):
        new_values = [ast.UnaryOp(ast.Not(), v) for v in node.values]
        # !(x or y) <> !x and !y
        if isinstance(node.op, ast.Or):
            return ast.BoolOp(ast.And(), new_values)
        # !(x and y) <> !x or !y
        if isinstance(node.op, ast.And):
            return ast.BoolOp(ast.Or(), new_values)

    if isinstance(node, ast.Compare):
        cmps = [ast.Compare(x, [negate(o)], [y])
                for x, o, y
                in zip([node.left] + node.comparators[:-1], node.ops,
                       node.comparators)]
        if len(cmps) == 1:
            return cmps[0]
        return ast.BoolOp(ast.Or(), cmps)

    if isinstance(node, ast.Eq):
        return ast.NotEq()
    if isinstance(node, ast.NotEq):
        return ast.Eq()
    if isinstance(node, ast.Gt):
        return ast.LtE()
    if isinstance(node, ast.GtE):
        return ast.Lt()
    if isinstance(node, ast.Lt):
        return ast.GtE()
    if isinstance(node, ast.LtE):
        return ast.Gt()
    if isinstance(node, ast.In):
        return ast.NotIn()
    if isinstance(node, ast.NotIn):
        return ast.In()

    if isinstance(node, ast.Attribute):
        if node.attr == 'False':
            return ast.Constant(True, None)
        if node.attr == 'True':
            return ast.Constant(False, None)

    raise UnsupportedExpression()


def bound_range(mapping, aliases, node, modified=None):
    """
    Bound the identifier in `mapping' with the expression in `node'.
    `aliases' is the result of aliasing analysis and `modified' is
    updated with the set of identifiers possibly `bounded' as the result
    of the call.

    Returns `modified' or a fresh set of modified identifiers.

    """

    if modified is None:
        modified = set()

    if isinstance(node, ast.Name):
        # could be anything not just an integral
        pass

    elif isinstance(node, ast.UnaryOp):
        try:
            negated = negate(node.operand)
            bound_range(mapping, aliases, negated, modified)
        except UnsupportedExpression:
            pass

    elif isinstance(node, ast.BoolOp):
        if isinstance(node.op, ast.And):
            for value in node.values:
                bound_range(mapping, aliases, value, modified)
        elif isinstance(node.op, ast.Or):
            mappings = [mapping.copy() for _ in node.values]
            for value, mapping_cpy in zip(node.values, mappings):
                bound_range(mapping_cpy, aliases, value, modified)
            for k in modified:
                mapping[k] = reduce(lambda x, y: x.union(y[k]),
                                    mappings[1:],
                                    mappings[0][k])

    elif isinstance(node, ast.Compare):
        left = node.left
        if isinstance(node.left, ast.Name):
            modified.add(node.left.id)

        for op, right in zip(node.ops, node.comparators):
            if isinstance(right, ast.Name):
                modified.add(right.id)

            if isinstance(left, ast.Name):
                left_interval = mapping[left.id]
            else:
                left_interval = mapping[left]

            if isinstance(right, ast.Name):
                right_interval = mapping[right.id]
            else:
                right_interval = mapping[right]

            l_l, l_h = left_interval.low, left_interval.high
            r_l, r_h = right_interval.low, right_interval.high

            r_i = l_i = None

            if isinstance(op, ast.Eq):
                low, high = max(l_l, r_l), min(l_h, r_h)
                if low <= high:
                    l_i = r_i = Interval(max(l_l, r_l), min(l_h, r_h))
            elif isinstance(op, ast.Lt):
                # l < r => l.low < r.high & l.high < r.high
                l_i = Interval(min(l_l, r_h - 1), min(l_h, r_h - 1))
                # l < r => r.low < l.low & r.high < l.low
                r_i = Interval(max(r_l, l_l + 1), max(r_h, l_l + 1))
            elif isinstance(op, ast.LtE):
                # l <= r => l.low <= r.high & l.high <= r.high
                l_i = Interval(min(l_l, r_h), min(l_h, r_h))
                # l <= r => r.low <= l.low & r.high <= l.low
                r_i = Interval(max(r_l, l_l), max(r_h, l_l))
            elif isinstance(op, ast.Gt):
                # l > r => l.low > r.low & l.high > r.low
                l_i = Interval(max(l_l, r_l + 1), max(l_h, r_l + 1))
                # l > r => r.low > l.high & r.high > l.high
                r_i = Interval(min(r_l, l_h - 1), min(r_h, l_h - 1))
            elif isinstance(op, ast.GtE):
                # l >= r => l.high >= r.low & l.low >= r.low
                l_i = Interval(max(l_l, r_l), max(l_h, r_l))
                # l >= r => r.low > l.high & r.high >= l.high
                r_i = Interval(min(r_l, l_h), min(r_h, l_h))
            elif isinstance(op, ast.In):
                if isinstance(right, (ast.List, ast.Tuple, ast.Set)):
                    if right.elts:
                        low = min(mapping[elt].low for elt in right.elts)
                        high = max(mapping[elt].high for elt in right.elts)
                        l_i = Interval(low, high)
                elif isinstance(right, ast.Call):
                    for alias in aliases[right.func]:
                        if not hasattr(alias, 'return_range_content'):
                            l_i = None
                            break
                        rrc = alias.return_range_content([mapping[arg] for arg
                                                          in right.args])
                        if l_i is None:
                            l_i = rrc
                        else:
                            l_i = l_i.union(alias.return_range(right))

            if l_i is not None and isinstance(left, ast.Name):
                mapping[left.id] = l_i
            if r_i is not None and isinstance(right, ast.Name):
                mapping[right.id] = r_i

            left = right


class RangeValuesBase(ModuleAnalysis[Aliases, CFG, UseOMP]):

    ResultHolder = object()
    ResultType = lambda: defaultdict(lambda: UNKNOWN_RANGE)

    def __init__(self):
        """Initialize instance variable and gather globals name information."""
        super().__init__()
        self.parent = self

    def add(self, variable, range_):
        """
        Add a new low and high bound for a variable.

        As it is flow insensitive, it compares it with old values and update it
        if needed.
        """
        if variable not in self.result:
            self.result[variable] = range_
        else:
            self.result[variable] = self.result[variable].union(range_)
        return self.result[variable]

    def unionify(self, other):
        for k, v in other.items():
            if k in self.result:
                self.result[k] = self.result[k].union(v)
            else:
                self.result[k] = v

    def widen(self, curr, other):
        self.result = curr
        for k, v in other.items():
            w = self.result.get(k, None)
            if w is None:
                self.result[k] = v
            elif v is not w:
                self.result[k] = w.widen(v)

    def visit_BoolOp(self, node):
        """ Merge right and left operands ranges.

        TODO : We could exclude some operand with this range information...

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = 2
        ...     c = 3
        ...     d = a or c''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['d']
        Interval(low=2, high=3)
        """
        res = list(zip(*[self.visit(elt).bounds() for elt in node.values]))
        return self.add(node, Interval(min(res[0]), max(res[1])))

    def visit_BinOp(self, node):
        """ Combine operands ranges for given operator.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = 2
        ...     c = 3
        ...     d = a - c''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['d']
        Interval(low=-1, high=-1)
        """
        res = combine(node.op, self.visit(node.left), self.visit(node.right))
        return self.add(node, res)

    def visit_UnaryOp(self, node):
        """ Update range with given unary operation.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = 2
        ...     c = -a
        ...     d = ~a
        ...     f = +a
        ...     e = not a''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['f']
        Interval(low=2, high=2)
        >>> res['c']
        Interval(low=-2, high=-2)
        >>> res['d']
        Interval(low=-3, high=-3)
        >>> res['e']
        Interval(low=0, high=1)
        """
        res = self.visit(node.operand)
        if isinstance(node.op, ast.Not):
            res = Interval(0, 1)
        elif(isinstance(node.op, ast.Invert) and
             isinstance(res.high, int) and
             isinstance(res.low, int)):
            res = Interval(~res.high, ~res.low)
        elif isinstance(node.op, ast.UAdd):
            pass
        elif isinstance(node.op, ast.USub):
            res = Interval(-res.high, -res.low)
        else:
            res = UNKNOWN_RANGE
        return self.add(node, res)

    def visit_IfExp(self, node):
        """ Use worst case for both possible values.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = 2 or 3
        ...     b = 4 or 5
        ...     c = a if a else b''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['c']
        Interval(low=2, high=5)
        """
        self.visit(node.test)
        body_res = self.visit(node.body)
        orelse_res = self.visit(node.orelse)
        return self.add(node, orelse_res.union(body_res))

    def visit_Compare(self, node):
        """ Boolean are possible index.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = 2 or 3
        ...     b = 4 or 5
        ...     c = a < b
        ...     d = b < 3
        ...     e = b == 4''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['c']
        Interval(low=1, high=1)
        >>> res['d']
        Interval(low=0, high=0)
        >>> res['e']
        Interval(low=0, high=1)
        """
        if any(isinstance(op, (ast.In, ast.NotIn, ast.Is, ast.IsNot))
               for op in node.ops):
            self.generic_visit(node)
            return self.add(node, Interval(0, 1))

        curr = self.visit(node.left)
        res = []
        for op, comparator in zip(node.ops, node.comparators):
            comparator = self.visit(comparator)
            fake = ast.Compare(ast.Name('x', ast.Load(), None, None),
                               [op],
                               [ast.Name('y', ast.Load(), None, None)])
            fake = ast.Expression(fake)
            ast.fix_missing_locations(fake)
            expr = compile(ast.gast_to_ast(fake), '<range_values>', 'eval')
            res.append(eval(expr, {'x': curr, 'y': comparator}))
        if all(res):
            return self.add(node, Interval(1, 1))
        elif any(r.low == r.high == 0 for r in res):
            return self.add(node, Interval(0, 0))
        else:
            return self.add(node, Interval(0, 1))

    def visit_Call(self, node):
        """ Function calls are not handled for now.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = builtins.range(10)''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=-inf, high=inf)
        """
        for alias in self.aliases[node.func]:
            if alias is MODULES['builtins']['getattr']:
                attr_name = node.args[-1].value
                attribute = attributes[attr_name][-1]
                self.add(node, attribute.return_range(None))
            elif isinstance(alias, Intrinsic):
                alias_range = alias.return_range(
                    [self.visit(n) for n in node.args])
                self.add(node, alias_range)
            elif isinstance(alias, ast.FunctionDef):
                if alias not in self.result:
                    state = self.save_state()
                    self.parent.visit(alias)
                    return_range = self.result[alias]
                    self.restore_state(state)
                else:
                    return_range = self.result[alias]
                self.add(node, return_range)
            else:
                self.result.pop(node, None)
                return self.generic_visit(node)
        return self.result[node]

    def visit_Constant(self, node):
        """ Handle literals integers values. """
        if isinstance(node.value, (bool, int)):
            return self.add(node, Interval(node.value, node.value))
        return UNKNOWN_RANGE

    def visit_Name(self, node):
        """ Get range for parameters for examples or false branching. """
        return self.add(node, self.result[node.id])

    def visit_Tuple(self, node):
        return self.add(node,
                        IntervalTuple(self.visit(elt) for elt in node.elts))

    def visit_Index(self, node):
        return self.add(node, self.visit(node.value))

    def visit_Subscript(self, node):
        if isinstance(node.value, ast.Call):
            for alias in self.aliases[node.value.func]:
                if alias is MODULES['builtins']['getattr']:
                    attr_name = node.value.args[-1].value
                    attribute = attributes[attr_name][-1]
                    self.add(node, attribute.return_range_content(None))
                elif isinstance(alias, Intrinsic):
                    self.add(node,
                             alias.return_range_content(
                                 [self.visit(n) for n in node.value.args]))
                else:
                    return self.generic_visit(node)
            if not self.aliases[node.value.func]:
                return self.generic_visit(node)
            self.visit(node.slice)
            return self.result[node]
        else:
            value = self.visit(node.value)
            slice = self.visit(node.slice)
            return self.add(node, value[slice])

    def visit_FunctionDef(self, node):
        """ Set default range value for globals and attributes.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(a, b): pass")
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=-inf, high=inf)
        """
        if node in self.result:
            return

        if self.use_omp:
            return

        self.result[node] = UNKNOWN_RANGE

        # Set this prematurely to avoid infinite callgraph loop
        prev_result = self.result.get(RangeValuesBase.ResultHolder, None)

        self.function_visitor(node)

        del self.result[node]
        self.add(node, self.result[RangeValuesBase.ResultHolder])

        if prev_result is not None:
            self.result[RangeValuesBase.ResultHolder] = prev_result
        else:
            del self.result[RangeValuesBase.ResultHolder]


class RangeValuesSimple(RangeValuesBase):

    """
    This analyse extract positive subscripts from code.

    It is flow sensitive and aliasing is not taken into account as integer
    doesn't create aliasing in Python.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('''
    ... def foo(a):
    ...     for i in builtins.range(1, 10):
    ...         c = i // 2''')
    >>> pm = passmanager.PassManager("test")
    >>> res = pm.gather(RangeValuesSimple, node)
    >>> res['c'], res['i']
    (Interval(low=0, high=5), Interval(low=1, high=10))
    """

    def __init__(self, parent=None):
        if parent is not None:
            self.parent = parent
            self.ctx = parent.ctx
            self.deps = parent.deps
            self.result = parent.result
            self.aliases = parent.aliases
            self.passmanager = parent.passmanager
        else:
            super(RangeValuesSimple, self).__init__()

    def generic_visit(self, node):
        """ Other nodes are not known and range value neither. """
        super(RangeValuesSimple, self).generic_visit(node)
        return self.add(node, UNKNOWN_RANGE)

    def save_state(self):
        return self.aliases,

    def restore_state(self, state):
        self.aliases, = state

    def function_visitor(self, node):
        for stmt in node.body:
            self.visit(stmt)

    def visit_Return(self, node):
        if node.value:
            return_range = self.visit(node.value)
            return self.add(RangeValues.ResultHolder, return_range)
        else:
            return self.generic_visit(node)

    def visit_Assert(self, node):
        """
        Constraint the range of variables

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(a): assert a >= 1; b = a + 1")
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['a']
        Interval(low=1, high=inf)
        >>> res['b']
        Interval(low=2, high=inf)
        """
        self.generic_visit(node)
        bound_range(self.result, self.aliases, node.test)

    def visit_Assign(self, node):
        """
        Set range value for assigned variable.

        We do not handle container values.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(): a = b = 2")
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['a']
        Interval(low=2, high=2)
        >>> res['b']
        Interval(low=2, high=2)
        """
        if not node.value:
            return
        assigned_range = self.visit(node.value)
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for target in targets:
            if isinstance(target, ast.Name):
                # Make sure all Interval doesn't alias for multiple variables.
                self.add(target.id, assigned_range)
            else:
                self.visit(target)

    visit_AnnAssign = visit_Assign

    def visit_AugAssign(self, node):
        """ Update range value for augassigned variables.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(): a = 2; a -= 1")
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['a']
        Interval(low=1, high=1)
        """
        self.generic_visit(node)
        if isinstance(node.target, ast.Name):
            name = node.target.id
            res = combine(node.op,
                          self.result[name],
                          self.result[node.value])
            self.result[name] = res

    def visit_For(self, node):
        """ Handle iterate variable in for loops.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = b = c = 2
        ...     for i in builtins.range(1):
        ...         a -= 1
        ...         b += 1''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['a']
        Interval(low=-inf, high=2)
        >>> res['b']
        Interval(low=2, high=inf)
        >>> res['c']
        Interval(low=2, high=2)

        >>> node = ast.parse('''
        ... def foo():
        ...     for i in (1, 2, 4):
        ...         a = i''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['a']
        Interval(low=1, high=4)
        """
        assert isinstance(node.target, ast.Name), "For apply on variables."
        self.visit(node.iter)
        if isinstance(node.iter, ast.Call):
            for alias in self.aliases[node.iter.func]:
                if isinstance(alias, Intrinsic):
                    self.add(node.target.id,
                             alias.return_range_content(
                                 [self.visit(n) for n in node.iter.args]))

        self.visit_loop(node,
                        ast.Compare(node.target, [ast.In()], [node.iter]))

    def visit_loop(self, node, cond=None):
        """ Handle incremented variables in loop body.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = b = c = 2
        ...     while a > 0:
        ...         a -= 1
        ...         b += 1''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['a']
        Interval(low=0, high=2)
        >>> res['b']
        Interval(low=2, high=inf)
        >>> res['c']
        Interval(low=2, high=2)
        """

        if cond is not None:
            init_range = self.result
            self.result = self.result.copy()
            bound_range(self.result, self.aliases, cond)

        # visit once to gather newly declared vars
        for stmt in node.body:
            self.visit(stmt)

        # freeze current state
        old_range = self.result.copy()

        # extra round
        for stmt in node.body:
            self.visit(stmt)

        # widen any change
        for expr, range_ in old_range.items():
            self.result[expr] = self.result[expr].widen(range_)

        # propagate the new informations again
        if cond is not None:
            bound_range(self.result, self.aliases, cond)
            for stmt in node.body:
                self.visit(stmt)
            self.unionify(init_range)
            self.visit(cond)

        for stmt in node.orelse:
            self.visit(stmt)

    def visit_While(self, node):
        self.visit(node.test)
        return self.visit_loop(node, node.test)

    def visit_If(self, node):
        """ Handle iterate variable across branches

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> pm = passmanager.PassManager("test")

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if a > 1: b = 1
        ...     else: b = 3''')

        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['b']
        Interval(low=1, high=3)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if a > 1: b = a
        ...     else: b = 3''')
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['b']
        Interval(low=2, high=inf)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if 0 < a < 4: b = a
        ...     else: b = 3''')
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['b']
        Interval(low=1, high=3)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if (0 < a) and (a < 4): b = a
        ...     else: b = 3''')
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['b']
        Interval(low=1, high=3)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if (a == 1) or (a == 2): b = a
        ...     else: b = 3''')
        >>> res = pm.gather(RangeValuesSimple, node)
        >>> res['b']
        Interval(low=1, high=3)
        """
        self.visit(node.test)
        old_range = self.result

        self.result = old_range.copy()
        bound_range(self.result, self.aliases, node.test)

        for stmt in node.body:
            self.visit(stmt)
        body_range = self.result

        self.result = old_range.copy()
        for stmt in node.orelse:
            self.visit(stmt)
        orelse_range = self.result

        self.result = body_range
        self.unionify(orelse_range)

    def visit_Try(self, node):
        init_range = self.result

        self.result = init_range.copy()
        for stmt in node.body:
            self.visit(stmt)
        self.unionify(init_range)

        init_range = self.result.copy()

        for handler in node.handlers:
            self.result, prev_state = init_range.copy(), self.result
            for stmt in handler.body:
                self.visit(stmt)
            self.unionify(prev_state)

        self.result, prev_state = init_range, self.result
        for stmt in node.orelse:
            self.visit(stmt)
        self.unionify(prev_state)

        for stmt in node.finalbody:
            self.visit(stmt)


class RangeValues(RangeValuesBase):

    """
    This analyse extract positive subscripts from code.

    It is flow sensitive and aliasing is not taken into account as integer
    doesn't create aliasing in Python.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse('''
    ... def foo(a):
    ...     for i in builtins.range(1, 10):
    ...         c = i // 2
    ...     return''')
    >>> pm = passmanager.PassManager("test")
    >>> res = pm.gather(RangeValues, node)
    >>> res['c'], res['i']
    (Interval(low=0, high=5), Interval(low=1, high=10))
    """

    def __init__(self):
        super(RangeValues, self).__init__()

    def generic_visit(self, node):
        """ Other nodes are not known and range value neither. """
        super(RangeValues, self).generic_visit(node)

        if isinstance(node, ast.stmt):
            if node in self.cfg:
                return self.cfg.successors(node)
        else:
            return self.add(node, UNKNOWN_RANGE)

    def cfg_visit(self, node, skip=None):
        successors = [node]
        visited = set() if skip is None else skip.copy()
        while successors:
            successor = successors.pop()
            if successor in visited:
                continue
            visited.add(successor)
            nexts = self.visit(successor)
            if nexts:
                successors.extend((n for n in nexts if n is not CFG.NIL))

    def save_state(self):
        return (self.cfg, self.aliases, self.use_omp, self.no_backward,
                self.no_if_split, self.result.copy())

    def restore_state(self, state):
        (self.cfg, self.aliases, self.use_omp, self.no_backward,
         self.no_if_split, self.result) = state

    def function_visitor(self, node):
        parent_result = self.result
        self.result = defaultdict(
            lambda: UNKNOWN_RANGE,
            [(k, v) for k, v in parent_result.items()
             if isinstance(k, ast.FunctionDef)])

        # try to visit the cfg, it's greedy but more accurate
        try:
            self.no_backward = 0
            self.no_if_split = 0
            self.cfg_visit(next(self.cfg.successors(node)))
            for k, v in self.result.items():
                parent_result[k] = v
            self.result = parent_result

        # too greedy? Never mind, we know how to be fast and simple :-)
        except RangeValueTooCostly:
            self.result = parent_result
            rvs = RangeValuesSimple(self)
            rvs.visit(node)

    def visit_Return(self, node):
        if node.value:
            return_range = self.visit(node.value)
            self.add(RangeValues.ResultHolder, return_range)
        return self.cfg.successors(node)

    def visit_Assert(self, node):
        """
        Constraint the range of variables

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(a): assert a >= 1; b = a + 1")
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=1, high=inf)
        >>> res['b']
        Interval(low=2, high=inf)
        """
        self.visit(node.test)
        bound_range(self.result, self.aliases, node.test)
        return self.cfg.successors(node)

    def visit_Assign(self, node):
        """
        Set range value for assigned variable.

        We do not handle container values.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(): a = b = 2")
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=2, high=2)
        >>> res['b']
        Interval(low=2, high=2)
        """
        assigned_range = self.visit(node.value)
        for target in node.targets:
            if isinstance(target, ast.Name):
                # Make sure all Interval doesn't alias for multiple variables.
                self.result[target.id] = assigned_range
            else:
                self.visit(target)
        return self.cfg.successors(node)

    def visit_AugAssign(self, node):
        """ Update range value for augassigned variables.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse("def foo(): a = 2; a -= 1")
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=1, high=1)
        """
        self.generic_visit(node)
        if isinstance(node.target, ast.Name):
            name = node.target.id
            res = combine(node.op,
                          self.result[name],
                          self.result[node.value])
            self.result[name] = res
        return self.cfg.successors(node)

    def visit_loop_successor(self, node):
        for successor in self.cfg.successors(node):
            if successor is not node.body[0]:
                if isinstance(node, ast.While):
                    bound_range(self.result, self.aliases,
                                ast.UnaryOp(ast.Not(), node.test))
                return [successor]

    def visit_For(self, node):
        """ Handle iterate variable in for loops.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = b = c = 2
        ...     for i in builtins.range(1):
        ...         a -= 1
        ...         b += 1
        ...     return''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=-inf, high=2)
        >>> res['b']
        Interval(low=2, high=inf)
        >>> res['c']
        Interval(low=2, high=2)

        >>> node = ast.parse('''
        ... def foo():
        ...     for i in (1, 2, 4):
        ...         a = i
        ...     return''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=1, high=4)
        """
        assert isinstance(node.target, ast.Name), "For apply on variables."
        self.visit(node.iter)
        init_state = self.result.copy()

        bound_range(self.result, self.aliases, ast.Compare(node.target,
                                                           [ast.In()],
                                                           [node.iter]))

        # visit body
        skip = {x for x in self.cfg.successors(node) if x is not node.body[0]}
        skip.add(node)
        next_ = self.cfg_visit(node.body[0], skip=skip)

        if self.no_backward:
            return self.visit_loop_successor(node)
        else:
            pass #self.no_backward += 1

        prev_state = self.result
        self.result = prev_state.copy()

        self.cfg_visit(node.body[0], skip=skip)

        self.widen(self.result, prev_state)
        self.cfg_visit(node.body[0], skip=skip)
        self.unionify(init_state)
        pass #self.no_backward -= 1

        return self.visit_loop_successor(node)

    def visit_While(self, node):
        """ Handle incremented variables in loop body.

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> node = ast.parse('''
        ... def foo():
        ...     a = b = c = 10
        ...     while a > 0:
        ...         a -= 1
        ...         b += 1
        ...     return''')
        >>> pm = passmanager.PassManager("test")
        >>> res = pm.gather(RangeValues, node)
        >>> res['a']
        Interval(low=-inf, high=0)
        >>> res['b']
        Interval(low=11, high=inf)
        >>> res['c']
        Interval(low=10, high=10)
        """
        test_range = self.visit(node.test)
        init_state = self.result.copy()
        skip = {x for x in self.cfg.successors(node) if x is not node.body[0]}
        skip.add(node)

        # if the test may be false, visit the tail
        if 0 in test_range:
            for successor in list(self.cfg.successors(node)):
                if successor is not node.body[0]:
                    self.cfg_visit(successor, skip=skip)

        bound_range(self.result, self.aliases, node.test)

        # visit body
        self.cfg_visit(node.body[0], skip=skip)

        if self.no_backward:
            if 0 in test_range:
                self.unionify(init_state)
            return self.visit_loop_successor(node)
        else:
            pass #self.no_backward += 1

        prev_state = self.result
        self.result = prev_state.copy()

        self.cfg_visit(node.body[0], skip=skip)

        self.widen(self.result, prev_state)

        # propagate the result of the widening
        self.cfg_visit(node.body[0], skip=skip)

        if 0 in test_range:
            self.unionify(init_state)
        else:
            self.unionify(prev_state)

        self.visit(node.test)

        pass #self.no_backward -= 1

        # exit from the while test
        return self.visit_loop_successor(node)

    def visit_If(self, node):
        """ Handle iterate variable across branches

        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> pm = passmanager.PassManager("test")

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if a > 1: b = 1
        ...     else: b = 3
        ...     pass''')

        >>> res = pm.gather(RangeValues, node)
        >>> res['b']
        Interval(low=1, high=3)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if a > 1: b = a
        ...     else: b = 3
        ...     pass''')
        >>> res = pm.gather(RangeValues, node)
        >>> res['b']
        Interval(low=2, high=inf)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if 0 < a < 4: b = a
        ...     else: b = 3
        ...     pass''')
        >>> res = pm.gather(RangeValues, node)
        >>> res['b']
        Interval(low=1, high=3)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if (0 < a) and (a < 4): b = a
        ...     else: b = 3
        ...     pass''')
        >>> res = pm.gather(RangeValues, node)
        >>> res['b']
        Interval(low=1, high=3)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if (a == 1) or (a == 2): b = a
        ...     else: b = 3
        ...     pass''')
        >>> res = pm.gather(RangeValues, node)
        >>> res['b']
        Interval(low=1, high=3)

        >>> node = ast.parse('''
        ... def foo(a):
        ...     b = 5
        ...     if a > 0: b = a
        ...     pass''')
        >>> res = pm.gather(RangeValues, node)
        >>> res['a'], res['b']
        (Interval(low=-inf, high=inf), Interval(low=1, high=inf))

        >>> node = ast.parse('''
        ... def foo(a):
        ...     if a > 3: b = 1
        ...     else: b = 2
        ...     if a > 1: b = 2
        ...     pass''')
        >>> res = pm.gather(RangeValues, node)
        >>> res['b']
        Interval(low=2, high=2)
        """
        # handling each branch becomes too costly, opt for a simpler,
        # less accurate algorithm.
        if self.no_if_split == 4:
            raise RangeValueTooCostly()

        self.no_if_split += 1

        test_range = self.visit(node.test)
        init_state = self.result.copy()

        if 1 in test_range:
            bound_range(self.result, self.aliases, node.test)
            self.cfg_visit(node.body[0])

        visited_successors = {node.body[0]}

        if node.orelse:
            if 0 in test_range:
                prev_state = self.result
                self.result = init_state.copy()
                bound_range(self.result, self.aliases,
                            ast.UnaryOp(ast.Not(), node.test))
                self.cfg_visit(node.orelse[0])
                self.unionify(prev_state)
            visited_successors.add(node.orelse[0])

        elif 0 in test_range:
            successors = self.cfg.successors(node)
            for successor in list(successors):
                # no else branch
                if successor not in visited_successors:
                    self.result, prev_state = init_state.copy(), self.result
                    bound_range(self.result, self.aliases,
                                ast.UnaryOp(ast.Not(), node.test))
                    self.cfg_visit(successor)
                    self.unionify(prev_state)

        self.no_if_split -= 1

    def visit_Try(self, node):
        init_range = self.result
        self.result = init_range.copy()
        self.cfg_visit(node.body[0])
        self.unionify(init_range)

        init_range = self.result.copy()

        for handler in node.handlers:
            self.result, prev_state = init_range.copy(), self.result
            self.cfg_visit(handler.body[0])
            self.unionify(prev_state)


# Un comment the line below to test RangeValuesSimple
# RangeValues = RangeValuesSimple
# RangeValues.__name__ = 'RangeValues'
