""" ConstantFolding performs some kind of partial evaluation.  """

from pythran.analyses import ConstantExpressions, ASTMatcher
from pythran.passmanager import Transformation
from pythran.tables import MODULES
from pythran.conversion import to_ast, ConversionError, ToNotEval, mangle
from pythran.analyses.ast_matcher import DamnTooLongPattern
from pythran.syntax import PythranSyntaxError
from pythran.utils import isintegral, isnum
from pythran.config import cfg

import builtins
import gast as ast
from copy import deepcopy
import logging
import numpy
import sys

logger = logging.getLogger('pythran')

class FunctionDefWrapper(object):

    def __init__(self, evaluator, func):
        self.evaluator = evaluator
        self.func = func

    def __call__(self, *args):
        missing_args = len(args) - len(self.func.args.args)
        if missing_args:
            defaults = tuple(self.evaluator.visit(default)
                        for default in self.func.args.defaults[missing_args:])
        else:
            defaults = ()
        locals = {arg.id: argv
                  for arg, argv in zip(self.func.args.args, args + defaults)}
        curr_locals, self.evaluator.locals = self.evaluator.locals, locals
        try:
            for stmt in self.func.body:
                self.evaluator.visit(stmt)
            res = locals.get('@', None)
            return res
        finally:
            self.evaluator.locals = curr_locals

class DispatchProxy(object):

    def __getattribute__(self, attr):
        return lambda obj, *args: getattr(obj, attr)(*args)

class PythranBuiltins(object):

    @staticmethod
    def static_list(*args):
        return list(*args)

    @staticmethod
    def static_if(cond, true_br, false_br):
        return true_br if cond else false_br

    @staticmethod
    def StaticIfReturn(val):
        return (None, val, None)

    @staticmethod
    def StaticIfNoReturn(val):
        return (None, val, None)

    @staticmethod
    def StaticIfCont(val):
        return (None, val, None)

    @staticmethod
    def StaticIfBreak(val):
        return (None, val, None)

    @staticmethod
    def is_none(val):
        return val is None

    @staticmethod
    def make_shape(*args):
        return args

    @staticmethod
    def restrict_assign(expr, value):
        numpy.copyto(expr, value)

class BreakLoop(Exception):
    pass

class ContinueLoop(Exception):
    pass

class ConstEval(ast.NodeVisitor):

    def __init__(self, globals):
        self.locals = {}
        self.globals = globals
        self.counter = 0
        self.counter_max = cfg.getint('pythran', 'fold_max_steps')

    def visit(self, node):
        self.counter += 1
        if self.counter == self.counter_max:
            raise ToNotEval()
        return getattr(self, 'visit_' + type(node).__name__)(node)

    # stmt
    def visit_Return(self, node):
        if '@' not in self.locals:
            self.locals['@'] = node.value and self.visit(node.value)

    def visit_Delete(self, node):
        if isinstance(node, ast.Name):
            self.locals.pop(node.id)

    def visit_Assign(self, node):
        if not node.value:
            return
        value = self.visit(node.value)
        targets = node.targets if isinstance(node, ast.Assign) else (node.target,)
        for target in targets:
            if isinstance(target, ast.Name):
                self.locals[target.id] = value
            elif isinstance(target, ast.Subscript):
                self.visit(target.value)[self.visit(target.slice)] = value
            else:
                raise NotImplementedError("assign")

    visit_AnnAssign = visit_Assign

    def visit_AugAssign(self, node):
        value = self.visit(node.value)
        ty = type(node.op)
        if isinstance(node.target, ast.Name):
            if ty is ast.Add:
                self.locals[node.target.id] += value
            elif ty is ast.Sub:
                self.locals[node.target.id] -= value
            elif ty is ast.Mult:
                self.locals[node.target.id] *= value
            elif ty is ast.MatMult:
                self.locals[node.target.id] @= value
            elif ty is ast.Div:
                self.locals[node.target.id] /= value
            elif ty is ast.Mod:
                self.locals[node.target.id] %= value
            elif ty is ast.Pow:
                self.locals[node.target.id] **= value
            elif ty is ast.LShift:
                self.locals[node.target.id] <<= value
            elif ty is ast.RShift:
                self.locals[node.target.id] >>= value
            elif ty is ast.BitOr:
                self.locals[node.target.id] |= value
            elif ty is ast.BitXor:
                self.locals[node.target.id] ^= value
            elif ty is ast.BitAnd:
                self.locals[node.target.id] &= value
            elif ty is ast.FloorDiv:
                self.locals[node.target.id] //= value
            else:
                raise ValueError("invalid binary op")
        elif isinstance(node.target, ast.Subscript):
            subscript = self.visit(node.target.slice)
            if ty is ast.Add:
                self.visit(node.target.value)[subscript] += value
            elif ty is ast.Sub:
                self.visit(node.target.value)[subscript] -= value
            elif ty is ast.Mult:
                self.visit(node.target.value)[subscript] *= value
            elif ty is ast.MatMult:
                self.visit(node.target.value)[subscript] @= value
            elif ty is ast.Div:
                self.visit(node.target.value)[subscript] /= value
            elif ty is ast.Mod:
                self.visit(node.target.value)[subscript] %= value
            elif ty is ast.Pow:
                self.visit(node.target.value)[subscript] **= value
            elif ty is ast.LShift:
                self.visit(node.target.value)[subscript] <<= value
            elif ty is ast.RShift:
                self.visit(node.target.value)[subscript] >>= value
            elif ty is ast.BitOr:
                self.visit(node.target.value)[subscript] |= value
            elif ty is ast.BitXor:
                self.visit(node.target.value)[subscript] ^= value
            elif ty is ast.BitAnd:
                self.visit(node.target.value)[subscript] &= value
            elif ty is ast.FloorDiv:
                self.visit(node.target.value)[subscript] //= value
            else:
                raise ValueError("invalid binary op")
        else:
            raise NotImplementedError("assign")

    def visit_For(self, node):
        iter = self.visit(node.iter)
        for elt in iter:
            if isinstance(node.target, ast.Name):
                self.locals[node.target.id] = elt
            else:
                raise ValueError("invalid loop target")
            try:
                for stmt in node.body:
                    self.visit(stmt)
            except ContinueLoop:
                continue
            except BreakLoop:
                break
        else:
            for stmt in node.orelse:
                self.visit(stmt)

    def visit_While(self, node):
        raise ToNotEval
        while self.visit(node.test):
            try:
                for stmt in node.body:
                    self.visit(stmt)
            except ContinueLoop:
                continue
            except BreakLoop:
                break
        else:
            for stmt in node.orelse:
                self.visit(stmt)

    def visit_If(self, node):
        if self.visit(node.test):
            for stmt in node.body:
                self.visit(stmt)
        else:
            for stmt in node.orelse:
                self.visit(stmt)

    def visit_Expr(self, node):
        self.visit(node.value)

    def visit_Break(self, node):
        raise BreakLoop

    def visit_Continue(self, node):
        raise ContinueLoop

    def visit_Pass(self, node):
        pass

    def visit_Yield(self, node):
        raise ToNotEval

    # expr

    def visit_BoolOp(self, node):
        values = (self.visit(value) for value in node.values)
        if type(node.op) is ast.And:
            for value in values:
                if not value:
                    return value
            return value
        elif type(node.op) is ast.Or:
            for value in values:
                if value:
                    return value
            return value
        else:
            raise ValueError("invalid bool op")

    def visit_BinOp(self, node):
        left = self.visit(node.left)
        right = self.visit(node.right)
        if type (node.op) is ast.Add:
            return left + right
        elif type (node.op) is ast.Sub:
            return left - right
        elif type (node.op) is ast.Mult:
            return left * right
        elif type (node.op) is ast.MatMult:
            return left @ right
        elif type (node.op) is ast.Div:
            return left / right
        elif type (node.op) is ast.Mod:
            return left % right
        elif type (node.op) is ast.Pow:
            return left ** right
        elif type (node.op) is ast.LShift:
            return left << right
        elif type (node.op) is ast.RShift:
            return left >> right
        elif type (node.op) is ast.BitOr:
            return left | right
        elif type (node.op) is ast.BitXor:
            return left ^ right
        elif type (node.op) is ast.BitAnd:
            return left & right
        elif type (node.op) is ast.FloorDiv:
            return left // right
        else:
            raise ValueError("invalid binary op")

    def visit_UnaryOp(self, node):
        value = self.visit(node.operand)
        if type(node.op) is ast.Invert:
            return ~value
        elif type(node.op) is ast.Not:
            return not value
        elif type(node.op) is ast.UAdd:
            return +value
        elif type(node.op) is ast.USub:
            return -value
        else:
            raise ValueError("invalid unary op")

    def visit_IfExp(self, node):
        test = self.visit(node.test)
        if test:
            return self.visit(node.body)
        else:
            return self.visit(node.orelse)

    def visit_Dict(self, node):
        return {self.visit(key):self.visit(value)
                for key, value in
                zip(node.keys, node.values)}

    def visit_Set(self, node):
        return {self.visit(elt) for elt in node.elts}

    def visit_Compare(self, node):
        curr = left = self.visit(node.left)
        for op, comparator in zip(node.ops, node.comparators):
            right = self.visit(comparator)
            if type(op) is ast.Eq:
                cond = curr == right
            elif type(op) is ast.NotEq:
                cond = curr != right
            elif type(op) is ast.Lt:
                cond = curr < right
            elif type(op) is ast.LtE:
                cond = curr <= right
            elif type(op) is ast.Gt:
                cond = curr > right
            elif type(op) is ast.GtE:
                cond = curr >= right
            elif type(op) is ast.Is:
                cond = curr is right
            elif type(op) is ast.IsNot:
                cond = curr is not right
            elif type(op) is ast.In:
                cond = curr in right
            elif type(op) is ast.NotIn:
                cond = curr not in right
            else:
                raise ValueError("invalid compare op")
            if not cond:
                return False
            curr = right
        return True

    def visit_Call(self, node):
        func = self.visit(node.func)
        args = [self.visit(arg) for arg in node.args]
        return func(*args)

    def visit_Constant(self, node):
        return node.value

    def visit_Attribute(self, node):
        value = self.visit(node.value)
        return getattr(value, node.attr)

    def visit_Subscript(self, node):
        value = self.visit(node.value)
        slice = self.visit(node.slice)
        return value[slice]

    def visit_Name(self, node):
        if node.id in self.locals:
            return self.locals[node.id]
        return self.globals[node.id]

    def visit_List(self, node):
        return [self.visit(elt) for elt in node.elts]

    def visit_Tuple(self, node):
        return tuple(self.visit(elt) for elt in node.elts)

    # slice

    def visit_Slice(self, node):
        lower = node.lower and self.visit(node.lower)
        upper = node.upper and self.visit(node.upper)
        step = node.step and self.visit(node.step)
        return slice(lower, upper, step)

    def visit_ExtSlice(self, node):
        return tuple(self.visit(dim) for dim in node.dims)

    def visit_Index(self, node):
        return self.visit(node.value)


    def __call__(self, node):
        self.counter = 0
        try:
            return self.visit(node)
        except:
            self.locals.clear()
            raise
        finally:
            assert not self.locals


class ConstantFolding(Transformation[ConstantExpressions]):

    """
    Replace constant expression by their evaluation.

    >>> import gast as ast
    >>> from pythran import passmanager, backend
    >>> node = ast.parse("def foo(): return 1+3")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(ConstantFolding, node)
    >>> print(pm.dump(backend.Python, node))
    def foo():
        return 4
    """

    def prepare(self, node):
        assert isinstance(node, ast.Module)
        self.env = {'builtins': builtins}
        self.consteval = ConstEval(self.env)

        if sys.implementation.name == 'pypy':
            self.env['__builtins__'] = self.env['builtins']

        for module_name in MODULES:
            # __dispatch__ is the only fake top-level module
            alias_module_name = mangle(module_name)
            try:
                if module_name == '__dispatch__':
                    self.env[alias_module_name] = DispatchProxy()
                else:
                    self.env[alias_module_name] = __import__(module_name)
            except ImportError:
                pass

        # we need to parse the whole code to be able to apply user-defined pure
        # function
        for stmt in node.body:
            if isinstance(stmt, ast.FunctionDef):
                self.env[stmt.name] = FunctionDefWrapper(self.consteval, stmt)

        super(ConstantFolding, self).prepare(node)

    def run(self, node):
        builtins.pythran = PythranBuiltins()
        try:
            return super(ConstantFolding, self).run(node)
        finally:
            # do not pollute builtin namespace outside of this pass
            del builtins.pythran

    def skip(self, node):
        return node

    visit_Constant = visit_Name = skip
    visit_Attribute = skip
    def visit_Attribute(self, node):
        # special case for numpy.Inf => numpy.inf because some filesystem don't
        # make the difference between inf.hpp and Inf.hpp :-/
        if node.attr == 'Inf':
            self.update = True
            node.attr = 'inf'
        return node

    visit_List = visit_Set = Transformation.generic_visit
    visit_Dict = visit_Tuple = Transformation.generic_visit

    def visit_Call(self, node):
        if isinstance(node.func, ast.Attribute):
            if isinstance(node.func.value, ast.Attribute):
                if node.func.value.attr == 'pythran':
                    return Transformation.generic_visit(self, node)
        return self.generic_visit(node)

    def generic_visit(self, node):
        if isinstance(node, ast.expr) and node in self.constant_expressions:
            try:
                value = self.consteval(node)
                new_node = to_ast(value)
                if not ASTMatcher(node).match(new_node):
                    self.update = True
                    return new_node
            except DamnTooLongPattern as e:
                logger.info(str(e) + ", skipping constant folding.")
            except ConversionError as e:
                print('error in constant folding: ', e)
                raise
            except ToNotEval:
                pass
            except Exception as e:
                if not cfg.getboolean('pythran', 'ignore_fold_error'):
                    msg = 'when folding expression, pythran met the following '\
                          'runtime exception:\n>>> {}'
                    raise PythranSyntaxError(msg.format(e), node)

        return Transformation.generic_visit(self, node)


class PartialConstantFolding(Transformation[ConstantExpressions]):

    """
    Replace partially constant expression by their evaluation.

    >>> import gast as ast
    >>> from pythran import passmanager, backend

    >>> node = ast.parse("def foo(n): return [n] * 2")
    >>> pm = passmanager.PassManager("test")
    >>> _, node = pm.apply(PartialConstantFolding, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(n):
        return [n, n]

    >>> node = ast.parse("def foo(n): return 2 * (n,)")
    >>> _, node = pm.apply(PartialConstantFolding, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(n):
        return (n, n)

    >>> node = ast.parse("def foo(n): return [n] + [n]")
    >>> _, node = pm.apply(PartialConstantFolding, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(n):
        return [n, n]

    >>> node = ast.parse("def foo(n, m): return (n,) + (m, n)")
    >>> _, node = pm.apply(PartialConstantFolding, node)
    >>> print(pm.dump(backend.Python, node))
    def foo(n, m):
        return (n, m, n)
    """

    def fold_mult_left(self, node):
        if not isinstance(node.left, (ast.List, ast.Tuple)):
            return False
        if not isnum(node.right):
            return False
        # FIXME: remove that check once we have a proper type inference engine
        if not isintegral(node.right):
            raise PythranSyntaxError("Multiplying a sequence by a float", node)
        return isinstance(node.op, ast.Mult)

    def fold_mult_right(self, node):
        if not isinstance(node.right, (ast.List, ast.Tuple)):
            return False
        if not isnum(node.left):
            return False
        # FIXME: remove that check once we have a proper type inference engine
        if not isintegral(node.left):
            raise PythranSyntaxError("Multiplying a sequence by a float", node)
        return isinstance(node.op, ast.Mult)

    def fold_add(self, node, ty):
        if not isinstance(node.left, ty):
            return False
        if not isinstance(node.right, ty):
            return False
        return isinstance(node.op, ast.Add)

    def visit_BinOp(self, node):
        if node in self.constant_expressions:
            return node

        node = self.generic_visit(node)
        if self.fold_mult_left(node):
            self.update = True
            node.left.elts = [deepcopy(elt)
                              for _ in range(node.right.value)
                              for elt in node.left.elts]
            return node.left

        if self.fold_mult_right(node):
            self.update = True
            node.left, node.right = node.right, node.left
            return self.visit(node)

        for ty in (ast.List, ast.Tuple):
            if self.fold_add(node, ty):
                self.update = True
                node.left.elts += node.right.elts
                return node.left

        return node

    def visit_Subscript(self, node):
        """
        >>> import gast as ast
        >>> from pythran import passmanager, backend
        >>> pm = passmanager.PassManager("test")

        >>> node = ast.parse("def foo(a): a[1:][3]")
        >>> _, node = pm.apply(PartialConstantFolding, node)
        >>> _, node = pm.apply(ConstantFolding, node)
        >>> print(pm.dump(backend.Python, node))
        def foo(a):
            a[4]

        >>> node = ast.parse("def foo(a): a[::2][3]")
        >>> _, node = pm.apply(PartialConstantFolding, node)
        >>> _, node = pm.apply(ConstantFolding, node)
        >>> print(pm.dump(backend.Python, node))
        def foo(a):
            a[6]

        >>> node = ast.parse("def foo(a): a[-4:][5]")
        >>> _, node = pm.apply(PartialConstantFolding, node)
        >>> _, node = pm.apply(ConstantFolding, node)
        >>> print(pm.dump(backend.Python, node))
        def foo(a):
            a[1]
        """
        self.generic_visit(node)
        if not isinstance(node.value, ast.Subscript):
            return node
        if not isinstance(node.value.slice, ast.Slice):
            return node

        if not isintegral(node.slice):
            return node

        slice_ = node.value.slice
        index = node.slice
        node = node.value

        lower = slice_.lower or ast.Constant(0, None)
        step = slice_.step or ast.Constant(1, None)
        node.slice = ast.BinOp(lower,
                               ast.Add(),
                               ast.BinOp(index,
                                         ast.Mult(),
                                         step))
        self.update = True
        return node
