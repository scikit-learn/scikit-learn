""" Optimization for Python costly pattern. """

from pythran.conversion import mangle
from pythran.analyses import Check, Placeholder, AST_or
from pythran.passmanager import Transformation

from copy import deepcopy
import gast as ast
from inspect import isclass


class Pattern(object):

    def match(self, node):
        self.check = Check(node, dict())
        return self.check.visit(self.pattern)

    def replace(self):
        return PlaceholderReplace(self.check.placeholders).visit(self.sub())

    def imports(self):
        return deepcopy(getattr(self, 'extra_imports', []))


class LenSetPattern(Pattern):
    # builtins.len(builtins.set(X)) => builtins.pythran.len_set(X)
    pattern = ast.Call(func=ast.Attribute(value=ast.Name('builtins',
                                                         ast.Load(),
                                                         None, None),
                                          attr="len", ctx=ast.Load()),
                       args=[ast.Call(
                           func=ast.Attribute(
                               value=ast.Name('builtins', ast.Load(),
                                              None, None),
                               attr="set", ctx=ast.Load()),
                           args=[Placeholder(0)],
                           keywords=[])],
                       keywords=[])

    @staticmethod
    def sub():
        return ast.Call(
            func=ast.Attribute(
                value=ast.Attribute(value=ast.Name('builtins', ast.Load(),
                                                   None, None),
                                    attr="pythran", ctx=ast.Load()),
                attr="len_set", ctx=ast.Load()),
            args=[Placeholder(0)], keywords=[])


class LenRangePattern(Pattern):
    # builtins.len(builtins.range(X)) => max(0, X)
    pattern = ast.Call(func=ast.Attribute(value=ast.Name('builtins',
                                                         ast.Load(),
                                                         None, None),
                                          attr="len", ctx=ast.Load()),
                       args=[ast.Call(
                           func=ast.Attribute(
                               value=ast.Name('builtins', ast.Load(),
                                              None, None),
                               attr="range", ctx=ast.Load()),
                           args=[Placeholder(0)],
                           keywords=[])],
                       keywords=[])

    @staticmethod
    def sub():
        return ast.Call(
            func=ast.Attribute(value=ast.Name('builtins', ast.Load(),
                                              None, None),
                               attr="max", ctx=ast.Load()),
            args=[ast.Constant(0, None), Placeholder(0)], keywords=[])


class TupleListPattern(Pattern):
    # builtins.tuple(builtins.list(X)) => builtins.tuple(X)

    pattern = ast.Call(func=ast.Attribute(value=ast.Name('builtins',
                                                         ast.Load(),
                                                         None, None),
                                          attr="tuple", ctx=ast.Load()),
                       args=[ast.Call(
                           func=ast.Attribute(
                               value=ast.Name('builtins',
                                              ast.Load(), None, None),
                               attr="list", ctx=ast.Load()),
                           args=[Placeholder(0)],
                           keywords=[])],
                       keywords=[])

    @staticmethod
    def sub():
        return ast.Call(
            func=ast.Attribute(value=ast.Name(id='builtins',
                                              ctx=ast.Load(),
                                              annotation=None,
                                              type_comment=None),
                               attr="tuple", ctx=ast.Load()),
            args=[Placeholder(0)], keywords=[])


class AbsSqrPattern(Pattern):
    # builtins.abs(X) ** 2 => builtins.pythran.abssqr(X)

    pattern = ast.Call(func=ast.Attribute(value=ast.Name(id=mangle('numpy'),
                                                         ctx=ast.Load(),
                                                         annotation=None,
                                                         type_comment=None),
                                          attr="square", ctx=ast.Load()),
                       args=[ast.Call(func=ast.Attribute(
                           value=ast.Name(id='builtins',
                                          ctx=ast.Load(),
                                          annotation=None,
                                          type_comment=None),
                           attr="abs",
                           ctx=ast.Load()),
                           args=[Placeholder(0)],
                           keywords=[])],
                       keywords=[])

    @staticmethod
    def sub():
        return ast.Call(
            func=ast.Attribute(
                value=ast.Attribute(value=ast.Name(id='builtins',
                                                   ctx=ast.Load(),
                                                   annotation=None,
                                                   type_comment=None),
                                    attr="pythran", ctx=ast.Load()),
                attr="abssqr", ctx=ast.Load()),
            args=[Placeholder(0)], keywords=[])


class AbsSqrPatternNumpy(AbsSqrPattern):
    # numpy.abs(X) ** 2 => builtins.pythran.abssqr(X)

    pattern = ast.Call(func=ast.Attribute(value=ast.Name(id=mangle('numpy'),
                                                         ctx=ast.Load(),
                                                         annotation=None,
                                                         type_comment=None),
                                          attr="square", ctx=ast.Load()),
                       args=[ast.Call(func=ast.Attribute(
                           value=ast.Name(id=mangle('numpy'),
                                          ctx=ast.Load(),
                                          annotation=None,
                                          type_comment=None),
                           attr="abs",
                           ctx=ast.Load()),
                           args=[Placeholder(0)],
                           keywords=[])],
                       keywords=[])


class PowFuncPattern(Pattern):
    # builtins.pow(X, Y) => X ** Y

    pattern = ast.Call(func=ast.Attribute(
        value=ast.Name(id=mangle('builtins'), ctx=ast.Load(),
                       annotation=None,
                       type_comment=None),
        attr='pow', ctx=ast.Load()),
        args=[Placeholder(0), Placeholder(1)],
        keywords=[])

    @staticmethod
    def sub():
        return ast.BinOp(Placeholder(0), ast.Pow(), Placeholder(1))


class SqrtPattern(Pattern):
    # X ** .5 => numpy.sqrt(X)

    pattern = ast.BinOp(Placeholder(0), ast.Pow(), ast.Constant(0.5, None))

    @staticmethod
    def sub():
        return ast.Call(
            func=ast.Attribute(value=ast.Name(id=mangle('numpy'),
                                              ctx=ast.Load(),
                                              annotation=None,
                                              type_comment=None),
                               attr="sqrt", ctx=ast.Load()),
            args=[Placeholder(0)], keywords=[])

    extra_imports = [ast.Import([ast.alias('numpy', mangle('numpy'))])]


class CbrtPattern(Pattern):
    # X ** .33333 => numpy.cbrt(X)
    pattern = ast.BinOp(Placeholder(0), ast.Pow(), ast.Constant(1./3., None))

    @staticmethod
    def sub():
        return ast.Call(
            func=ast.Attribute(value=ast.Name(id=mangle('numpy'),
                                              ctx=ast.Load(),
                                              annotation=None,
                                              type_comment=None),
                               attr="cbrt", ctx=ast.Load()),
            args=[Placeholder(0)], keywords=[])

    extra_imports = [ast.Import([ast.alias('numpy', mangle('numpy'))])]


class TuplePattern(Pattern):
    # builtins.tuple([X, ..., Z]) => (X, ..., Z)
    pattern = ast.Call(func=ast.Attribute(value=ast.Name(id='builtins',
                                                         ctx=ast.Load(),
                                                         annotation=None,
                                                         type_comment=None),
                                          attr="tuple", ctx=ast.Load()),
                       args=[ast.List(Placeholder(0), ast.Load())],
                       keywords=[])

    @staticmethod
    def sub():
        return ast.Tuple(Placeholder(0), ast.Load())


class ReversedRangePattern(Pattern):
    # builtins.reversed(builtins.range(X)) =>
    # builtins.range(X-1, -1, -1)
    # FIXME : We should do it even when begin/end/step are given
    pattern = ast.Call(func=ast.Attribute(value=ast.Name(id='builtins',
                                                         ctx=ast.Load(),
                                                         annotation=None,
                                                         type_comment=None),
                                          attr="reversed", ctx=ast.Load()),
                       args=[ast.Call(
                           func=ast.Attribute(
                               value=ast.Name(id='builtins',
                                              ctx=ast.Load(), annotation=None,
                                              type_comment=None),
                               attr='range', ctx=ast.Load()),
                           args=[Placeholder(0)],
                           keywords=[])],
                       keywords=[])

    @staticmethod
    def sub():
        return ast.Call(
            func=ast.Attribute(value=ast.Name(id='builtins',
                                              ctx=ast.Load(), annotation=None,
                                              type_comment=None),
                               attr='range', ctx=ast.Load()),
            args=[ast.BinOp(left=Placeholder(0), op=ast.Sub(),
                            right=ast.Constant(1, None)),
                  ast.Constant(-1, None),
                  ast.Constant(-1, None)],
            keywords=[])


class SqrPattern(Pattern):
    # X * X => X ** 2
    pattern = ast.BinOp(left=Placeholder(0),
                        op=ast.Mult(),
                        right=Placeholder(0))

    @staticmethod
    def sub():
        return ast.BinOp(left=Placeholder(0), op=ast.Pow(),
                         right=ast.Constant(2, None))


class StrJoinPattern(Pattern):
    # a + "..." + b => "...".join((a, b))
    pattern = ast.BinOp(left=ast.BinOp(left=Placeholder(0),
                                       op=ast.Add(),
                                       right=ast.Constant(Placeholder(1, str),
                                                          None)),
                        op=ast.Add(),
                        right=Placeholder(2))

    @staticmethod
    def sub():
        return ast.Call(func=ast.Attribute(
            ast.Attribute(
                ast.Name('builtins', ast.Load(), None, None),
                'str',
                ast.Load()),
            'join', ast.Load()),
            args=[ast.Constant(Placeholder(1), None),
                  ast.Tuple([Placeholder(0), Placeholder(2)], ast.Load())],
            keywords=[])


class ArgminmaxPattern(Pattern):
    # numpy.argminmax(X * positive_cst) => numpy.argminmax(X)

    pattern = ast.Call(
            func=ast.Attribute(
                value=ast.Name(id=mangle('numpy'),
                               ctx=ast.Load(),
                               annotation=None,
                               type_comment=None),
                attr=Placeholder(0,
                                 constraint=lambda s: s in ("argmax", "argmin")),
                ctx=ast.Load()),
            args=[AST_or(
                ast.BinOp(
                    ast.Constant(
                        Placeholder(2, constraint=lambda n: n > 0),
                        None),
                    ast.Mult(),
                    Placeholder(1)),
                ast.BinOp(
                    Placeholder(1),
                    ast.Mult(),
                    ast.Constant(
                        Placeholder(2, constraint=lambda n: n > 0),
                        None))
                    )
                  ],
             keywords=[])

    @staticmethod
    def sub():
        return ast.Call(
                func=ast.Attribute(value=ast.Name(id=mangle('numpy'),
                                                   ctx=ast.Load(),
                                                   annotation=None,
                                                   type_comment=None),
                                    attr=Placeholder(0), ctx=ast.Load()),
            args=[Placeholder(1)], keywords=[])


know_patterns = [x for x in globals().values()
                 if isclass(x) and issubclass(x, Pattern) and x is not Pattern]


class PlaceholderReplace(Transformation):

    """ Helper class to replace the placeholder once value is collected. """

    def __init__(self, placeholders):
        """ Store placeholders value collected. """
        super().__init__()
        self.placeholders = placeholders

    def visit(self, node):
        """ Replace the placeholder if it is one or continue. """
        if isinstance(node, Placeholder):
            return self.placeholders[node.id]
        else:
            return super(PlaceholderReplace, self).visit(node)

class PatternTransform(Transformation):

    """
    Replace all known pattern by pythran function call.

    Based on BaseMatcher to search correct pattern.
    """

    def visit_Module(self, node):
        self.extra_imports = []
        self.generic_visit(node)
        node.body = self.extra_imports + node.body
        return node

    def apply_patterns(self, node, patterns):
        for pattern in patterns:
            matcher = pattern()
            if matcher.match(node):
                self.extra_imports.extend(matcher.imports())
                node = matcher.replace()
                self.update = True
        return self.generic_visit(node)

    CallPatterns = ()
    def visit_Call(self, node):
        return self.apply_patterns(node, PatternTransform.CallPatterns)

    BinOpPatterns = ()
    def visit_BinOp(self, node):
        return self.apply_patterns(node, PatternTransform.BinOpPatterns)


for known_pattern in know_patterns:
    pattern_selector = type(known_pattern.pattern).__name__
    attr_name = "{}Patterns".format(pattern_selector)
    setattr(PatternTransform, attr_name,
            getattr(PatternTransform, attr_name) + (known_pattern,))

