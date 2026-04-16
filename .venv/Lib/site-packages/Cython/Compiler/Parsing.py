# cython: auto_cpdef=True, infer_types=True, py2_import=True
#
#   Parser
#


# This should be done automatically
import cython
cython.declare(Nodes=object, ExprNodes=object, EncodedString=object,
               bytes_literal=object, StringEncoding=object,
               FileSourceDescriptor=object, lookup_unicodechar=object,
               Future=object, Options=object, error=object, warning=object,
               Builtin=object, ModuleNode=object, Utils=object, _unicode=object, _bytes=object,
               re=object, _parse_escape_sequences=object, _parse_escape_sequences_raw=object,
               partial=object, reduce=object,
               _CDEF_MODIFIERS=tuple, COMMON_BINOP_MISTAKES=dict)

from io import StringIO
import re
from unicodedata import lookup as lookup_unicodechar
from functools import partial, reduce

from .Scanning import PyrexScanner, FileSourceDescriptor, tentatively_scan
from . import Nodes
from . import ExprNodes
from . import MatchCaseNodes
from . import Builtin
from . import StringEncoding
from .StringEncoding import EncodedString, bytes_literal
from .ModuleNode import ModuleNode
from .Errors import error, warning, CompileError
from .. import Utils
from . import Future
from . import Options


_CDEF_MODIFIERS = ('inline', 'nogil', 'api')
statement_terminators = cython.declare(frozenset, frozenset((
    ';', 'NEWLINE', 'EOF')))

class Ctx:
    #  Parsing context
    level = 'other'
    visibility = 'private'
    cdef_flag = False
    typedef_flag = False
    api = False
    overridable = False
    nogil = False
    namespace = None
    templates = None
    allow_struct_enum_decorator = False

    def __init__(self, **kwds):
        self.__dict__.update(kwds)

    def __call__(self, **kwds):
        ctx = Ctx()
        d = ctx.__dict__
        d.update(self.__dict__)
        d.update(kwds)
        return ctx


@cython.cfunc
def p_ident(s: PyrexScanner, message="Expected an identifier"):
    if s.sy == 'IDENT':
        name = s.context.intern_ustring(s.systring)
        s.next()
        return name
    else:
        s.error(message)


@cython.cfunc
def p_ident_list(s: PyrexScanner):
    names = []
    while s.sy == 'IDENT':
        names.append(s.context.intern_ustring(s.systring))
        s.next()
        if s.sy != ',':
            break
        s.next()
    return names

#------------------------------------------
#
#   Expressions
#
#------------------------------------------

@cython.cfunc
def p_binop_operator(s: PyrexScanner) -> tuple:
    pos = s.position()
    op = s.sy
    s.next()
    return op, pos


# signature is currently overridden in pxd file
def p_binop_expr(s: PyrexScanner, ops, p_sub_expr):
    n1 = p_sub_expr(s)
    while s.sy in ops:
        op, pos = p_binop_operator(s)
        n2 = p_sub_expr(s)
        n1 = ExprNodes.binop_node(pos, op, n1, n2)
        if op == '/':
            if Future.division in s.context.future_directives:
                n1.truedivision = True
            else:
                n1.truedivision = None  # unknown
    return n1


#lambdef: 'lambda' [varargslist] ':' test

@cython.cfunc
def p_lambdef(s: PyrexScanner):
    # s.sy == 'lambda'
    pos = s.position()
    s.next()
    if s.sy == ':':
        args = []
        star_arg = starstar_arg = None
    else:
        args, star_arg, starstar_arg = p_varargslist(
            s, terminator=':', annotated=False)
    s.expect(':')
    expr = p_test(s)
    return ExprNodes.LambdaNode(
        pos, args = args,
        star_arg = star_arg, starstar_arg = starstar_arg,
        result_expr = expr)


#test: or_test ['if' or_test 'else' test] | lambdef

@cython.cfunc
def p_test(s: PyrexScanner):
    # The check for a following ':=' is only for error reporting purposes.
    # It simply changes a
    #   expected ')', found ':='
    # message into something a bit more descriptive.
    # It is close to what the PEG parser does in CPython, where an expression has
    # a lookahead assertion that it isn't followed by ':='
    expr = p_test_allow_walrus_after(s)
    if s.sy == ':=':
        s.error("invalid syntax: assignment expression not allowed in this context")
    return expr


@cython.cfunc
def p_test_allow_walrus_after(s: PyrexScanner):
    if s.sy == 'lambda':
        return p_lambdef(s)
    pos = s.position()
    expr = p_or_test(s)
    if s.sy == 'if':
        s.next()
        test = p_or_test(s)
        s.expect('else')
        other = p_test(s)
        return ExprNodes.CondExprNode(pos, test=test, true_val=expr, false_val=other)
    else:
        return expr


@cython.cfunc
def p_namedexpr_test(s: PyrexScanner):
    # defined in the LL parser as
    #  namedexpr_test: test [':=' test]
    # The requirement that the LHS is a name is not enforced in the grammar.
    # For comparison the PEG parser does:
    #  1. look for "name :=", if found it's definitely a named expression
    #     so look for expression
    #  2. Otherwise, look for expression
    lhs = p_test_allow_walrus_after(s)
    if s.sy == ':=':
        position = s.position()
        if not lhs.is_name:
            s.error("Left-hand side of assignment expression must be an identifier", fatal=False)
        s.next()
        rhs = p_test(s)
        return ExprNodes.AssignmentExpressionNode(position, lhs=lhs, rhs=rhs)
    return lhs


#or_test: and_test ('or' and_test)*

COMMON_BINOP_MISTAKES = {'||': 'or', '&&': 'and'}

@cython.cfunc
def p_or_test(s: PyrexScanner):
    return p_rassoc_binop_expr(s, 'or', p_and_test)


# signature is currently overridden in pxd file
def p_rassoc_binop_expr(s: PyrexScanner, op, p_subexpr):
    n1 = p_subexpr(s)
    if s.sy == op:
        pos = s.position()
        op = s.sy
        s.next()
        n2 = p_rassoc_binop_expr(s, op, p_subexpr)
        n1 = ExprNodes.binop_node(pos, op, n1, n2)
    elif s.sy in COMMON_BINOP_MISTAKES and COMMON_BINOP_MISTAKES[s.sy] == op:
        # Only report this for the current operator since we pass through here twice for 'and' and 'or'.
        warning(s.position(),
                "Found the C operator '%s', did you mean the Python operator '%s'?" % (s.sy, op),
                level=1)
    return n1


#and_test: not_test ('and' not_test)*

@cython.cfunc
def p_and_test(s: PyrexScanner):
    #return p_binop_expr(s, ('and',), p_not_test)
    return p_rassoc_binop_expr(s, 'and', p_not_test)


#not_test: 'not' not_test | comparison

@cython.cfunc
def p_not_test(s: PyrexScanner):
    if s.sy == 'not':
        pos = s.position()
        s.next()
        return ExprNodes.NotNode(pos, operand = p_not_test(s))
    else:
        return p_comparison(s)


#comparison: expr (comp_op expr)*
#comp_op: '<'|'>'|'=='|'>='|'<='|'<>'|'!='|'in'|'not' 'in'|'is'|'is' 'not'

@cython.cfunc
def p_comparison(s: PyrexScanner):
    n1 = p_starred_expr(s)
    if s.sy in comparison_ops:
        pos = s.position()
        op = p_cmp_op(s)
        n2 = p_starred_expr(s)
        n1 = ExprNodes.PrimaryCmpNode(pos,
            operator = op, operand1 = n1, operand2 = n2)
        if s.sy in comparison_ops:
            n1.cascade = p_cascaded_cmp(s)
    return n1


@cython.cfunc
def p_test_or_starred_expr(s: PyrexScanner):
    if s.sy == '*':
        return p_starred_expr(s)
    else:
        return p_test(s)


@cython.cfunc
def p_namedexpr_test_or_starred_expr(s: PyrexScanner):
    if s.sy == '*':
        return p_starred_expr(s)
    else:
        return p_namedexpr_test(s)


@cython.cfunc
def p_starred_expr(s: PyrexScanner):
    pos = s.position()
    if s.sy == '*':
        starred = True
        s.next()
    else:
        starred = False
    expr = p_bit_expr(s)
    if starred:
        expr = ExprNodes.StarredUnpackingNode(pos, expr)
    return expr


@cython.cfunc
def p_cascaded_cmp(s: PyrexScanner):
    pos = s.position()
    op = p_cmp_op(s)
    n2 = p_starred_expr(s)
    result = ExprNodes.CascadedCmpNode(pos,
        operator = op, operand2 = n2)
    if s.sy in comparison_ops:
        result.cascade = p_cascaded_cmp(s)
    return result


@cython.cfunc
def p_cmp_op(s: PyrexScanner):
    if s.sy == 'not':
        s.next()
        s.expect('in')
        op = 'not_in'
    elif s.sy == 'is':
        s.next()
        if s.sy == 'not':
            s.next()
            op = 'is_not'
        else:
            op = 'is'
    else:
        op = s.sy
        s.next()
    if op == '<>':
        op = '!='
    return op


comparison_ops = cython.declare(frozenset, frozenset((
    '<', '>', '==', '>=', '<=', '<>', '!=',
    'in', 'is', 'not'
)))


#expr: xor_expr ('|' xor_expr)*

@cython.cfunc
def p_bit_expr(s: PyrexScanner):
    return p_binop_expr(s, ('|',), p_xor_expr)


#xor_expr: and_expr ('^' and_expr)*

@cython.cfunc
def p_xor_expr(s: PyrexScanner):
    return p_binop_expr(s, ('^',), p_and_expr)


#and_expr: shift_expr ('&' shift_expr)*

@cython.cfunc
def p_and_expr(s: PyrexScanner):
    return p_binop_expr(s, ('&',), p_shift_expr)


#shift_expr: arith_expr (('<<'|'>>') arith_expr)*

@cython.cfunc
def p_shift_expr(s: PyrexScanner):
    return p_binop_expr(s, ('<<', '>>'), p_arith_expr)


#arith_expr: term (('+'|'-') term)*

@cython.cfunc
def p_arith_expr(s: PyrexScanner):
    return p_binop_expr(s, ('+', '-'), p_term)


#term: factor (('*'|'@'|'/'|'%'|'//') factor)*

@cython.cfunc
def p_term(s: PyrexScanner):
    return p_binop_expr(s, ('*', '@', '/', '%', '//'), p_factor)


#factor: ('+'|'-'|'~'|'&'|typecast|sizeof) factor | power

@cython.cfunc
def p_factor(s: PyrexScanner):
    # little indirection for C-ification purposes
    return _p_factor(s)


@cython.cfunc
def _p_factor(s: PyrexScanner):
    sy = s.sy
    if sy in ('+', '-', '~'):
        op = s.sy
        pos = s.position()
        s.next()
        return ExprNodes.unop_node(pos, op, p_factor(s))
    elif not s.in_python_file:
        if sy == '&':
            pos = s.position()
            s.next()
            arg = p_factor(s)
            return ExprNodes.AmpersandNode(pos, operand = arg)
        elif sy == "<":
            return p_typecast(s)
        elif sy == 'IDENT' and s.systring == "sizeof":
            return p_sizeof(s)
    return p_power(s)


@cython.cfunc
def p_typecast(s: PyrexScanner):
    # s.sy == "<"
    pos = s.position()
    s.next()
    base_type = p_c_base_type(s)
    is_memslice = isinstance(base_type, Nodes.MemoryViewSliceTypeNode)
    is_other_unnamed_type = isinstance(base_type, (
        Nodes.TemplatedTypeNode,
        Nodes.CConstOrVolatileTypeNode,
        Nodes.CTupleBaseTypeNode,
    ))
    if not (is_memslice or is_other_unnamed_type) and base_type.name is None:
        s.error("Unknown type")
    declarator = p_c_declarator(s, empty=True)
    if s.sy == '?':
        s.next()
        typecheck = True
    else:
        typecheck = False
    s.expect(">")
    operand = p_factor(s)
    if is_memslice:
        return ExprNodes.CythonArrayNode(pos, base_type_node=base_type, operand=operand)

    return ExprNodes.TypecastNode(pos,
        base_type = base_type,
        declarator = declarator,
        operand = operand,
        typecheck = typecheck)


@cython.cfunc
def p_sizeof(s: PyrexScanner):
    # s.sy == ident "sizeof"
    pos = s.position()
    s.next()
    s.expect('(')
    # Here we decide if we are looking at an expression or type
    # If it is actually a type, but parsable as an expression,
    # we treat it as an expression here.
    if looking_at_expr(s):
        operand = p_test(s)
        node = ExprNodes.SizeofVarNode(pos, operand = operand)
    else:
        base_type = p_c_base_type(s)
        declarator = p_c_declarator(s, empty=True)
        node = ExprNodes.SizeofTypeNode(pos,
            base_type = base_type, declarator = declarator)
    s.expect(')')
    return node


@cython.cfunc
def p_yield_expression(s: PyrexScanner, statement_terminators: frozenset = statement_terminators):
    # s.sy == "yield"
    pos = s.position()
    s.next()
    is_yield_from = False
    if s.sy == 'from':
        is_yield_from = True
        s.next()
    if s.sy != ')' and s.sy not in statement_terminators:
        # "yield from" does not support implicit tuples, but "yield" does ("yield 1,2")
        arg = p_test(s) if is_yield_from else p_testlist(s)
    else:
        if is_yield_from:
            s.error("'yield from' requires a source argument",
                    pos=pos, fatal=False)
        arg = None
    if is_yield_from:
        return ExprNodes.YieldFromExprNode(pos, arg=arg)
    else:
        return ExprNodes.YieldExprNode(pos, arg=arg)


@cython.cfunc
def p_yield_statement(s: PyrexScanner):
    # s.sy == "yield"
    yield_expr = p_yield_expression(s)
    return Nodes.ExprStatNode(yield_expr.pos, expr=yield_expr)


@cython.cfunc
def p_async_statement(s: PyrexScanner, ctx, decorators):
    # s.sy >> 'async' ...
    if s.sy == 'def':
        # 'async def' statements aren't allowed in pxd files
        if 'pxd' in ctx.level:
            s.error('def statement not allowed here')
        s.level = ctx.level
        return p_def_statement(s, decorators, is_async_def=True)
    elif decorators:
        s.error("Decorators can only be followed by functions or classes")
    elif s.sy == 'for':
        return p_for_statement(s, is_async=True)
    elif s.sy == 'with':
        s.next()
        return p_with_items(s, is_async=True)
    else:
        s.error("expected one of 'def', 'for', 'with' after 'async'")


#power: atom_expr ('**' factor)*
#atom_expr: ['await'] atom trailer*

@cython.cfunc
def p_power(s: PyrexScanner):
    if s.systring == 'new' and s.peek()[0] == 'IDENT':
        return p_new_expr(s)
    await_pos = None
    if s.sy == 'await':
        await_pos = s.position()
        s.next()
    n1 = p_atom(s)
    while s.sy in ('(', '[', '.'):
        n1 = p_trailer(s, n1)
    if await_pos:
        n1 = ExprNodes.AwaitExprNode(await_pos, arg=n1)
    if s.sy == '**':
        pos = s.position()
        s.next()
        n2 = p_factor(s)
        n1 = ExprNodes.binop_node(pos, '**', n1, n2)
    return n1


@cython.cfunc
def p_new_expr(s: PyrexScanner):
    # s.systring == 'new'.
    pos = s.position()
    s.next()
    cppclass = p_c_base_type(s)
    return p_call(s, ExprNodes.NewExprNode(pos, cppclass = cppclass))


#trailer: '(' [arglist] ')' | '[' subscriptlist ']' | '.' NAME

@cython.cfunc
def p_trailer(s: PyrexScanner, node1):
    pos = s.position()
    if s.sy == '(':
        return p_call(s, node1)
    elif s.sy == '[':
        return p_index(s, node1)
    else:  # s.sy == '.'
        s.next()
        name = p_ident(s)
        return ExprNodes.AttributeNode(pos,
            obj=node1, attribute=name)


# arglist:  argument (',' argument)* [',']
# argument: [test '='] test       # Really [keyword '='] test

# since PEP 448:
# argument: ( test [comp_for] |
#             test '=' test |
#             '**' expr |
#             star_expr )

@cython.cfunc
def p_call_parse_args(s: PyrexScanner, allow_genexp: cython.bint = True) -> tuple:
    # s.sy == '('
    s.next()
    positional_args = []
    keyword_args = []
    starstar_seen = False
    last_was_tuple_unpack = False
    while s.sy != ')':
        if s.sy == '*':
            if starstar_seen:
                s.error("Non-keyword arg following keyword arg", pos=s.position())
            s.next()
            positional_args.append(p_test(s))
            last_was_tuple_unpack = True
        elif s.sy == '**':
            s.next()
            keyword_args.append(p_test(s))
            starstar_seen = True
        else:
            arg = p_namedexpr_test(s)
            if s.sy == '=':
                s.next()
                if not arg.is_name:
                    s.error("Expected an identifier before '='",
                            pos=arg.pos)
                encoded_name = s.context.intern_ustring(arg.name)
                keyword = ExprNodes.IdentifierStringNode(
                    arg.pos, value=encoded_name)
                arg = p_test(s)
                keyword_args.append((keyword, arg))
            else:
                if keyword_args:
                    s.error("Non-keyword arg following keyword arg", pos=arg.pos)
                if positional_args and not last_was_tuple_unpack:
                    positional_args[-1].append(arg)
                else:
                    positional_args.append([arg])
                last_was_tuple_unpack = False
        if s.sy != ',':
            break
        s.next()

    if s.sy in ('for', 'async') and allow_genexp:
        if not keyword_args and not last_was_tuple_unpack:
            if len(positional_args) == 1 and len(positional_args[0]) == 1:
                positional_args = [[p_genexp(s, positional_args[0][0])]]
    s.expect(')')
    return positional_args or [[]], keyword_args


@cython.cfunc
def p_call_build_packed_args(pos, positional_args, keyword_args) -> tuple:
    keyword_dict = None

    subtuples = [
        ExprNodes.TupleNode(pos, args=arg) if isinstance(arg, list) else ExprNodes.AsTupleNode(pos, arg=arg)
        for arg in positional_args
    ]
    # TODO: implement a faster way to join tuples than creating each one and adding them
    arg_tuple = reduce(partial(ExprNodes.binop_node, pos, '+'), subtuples)

    if keyword_args:
        kwargs = []
        dict_items = []
        for item in keyword_args:
            if isinstance(item, tuple):
                key, value = item
                dict_items.append(ExprNodes.DictItemNode(pos=key.pos, key=key, value=value))
            elif item.is_dict_literal:
                # unpack "**{a:b}" directly
                dict_items.extend(item.key_value_pairs)
            else:
                if dict_items:
                    kwargs.append(ExprNodes.DictNode(
                        dict_items[0].pos, key_value_pairs=dict_items, reject_duplicates=True))
                    dict_items = []
                kwargs.append(item)

        if dict_items:
            kwargs.append(ExprNodes.DictNode(
                dict_items[0].pos, key_value_pairs=dict_items, reject_duplicates=True))

        if kwargs:
            if len(kwargs) == 1 and kwargs[0].is_dict_literal:
                # only simple keyword arguments found -> one dict
                keyword_dict = kwargs[0]
            else:
                # at least one **kwargs
                keyword_dict = ExprNodes.MergedDictNode(pos, keyword_args=kwargs)

    return arg_tuple, keyword_dict


@cython.cfunc
def p_call(s: PyrexScanner, function):
    # s.sy == '('
    pos = s.position()
    positional_args, keyword_args = p_call_parse_args(s)

    if not keyword_args and len(positional_args) == 1 and isinstance(positional_args[0], list):
        return ExprNodes.SimpleCallNode(pos, function=function, args=positional_args[0])
    else:
        arg_tuple, keyword_dict = p_call_build_packed_args(pos, positional_args, keyword_args)
        return ExprNodes.GeneralCallNode(
            pos, function=function, positional_args=arg_tuple, keyword_args=keyword_dict)


#lambdef: 'lambda' [varargslist] ':' test

#subscriptlist: subscript (',' subscript)* [',']

@cython.cfunc
def p_index(s: PyrexScanner, base):
    # s.sy == '['
    pos = s.position()
    s.next()
    subscripts, is_single_value = p_subscript_list(s)
    if is_single_value and len(subscripts[0]) == 2:
        start, stop = subscripts[0]
        result = ExprNodes.SliceIndexNode(pos,
            base = base, start = start, stop = stop)
    else:
        indexes = make_slice_nodes(pos, subscripts)
        if is_single_value:
            index = indexes[0]
        else:
            index = ExprNodes.TupleNode(pos, args = indexes)
        result = ExprNodes.IndexNode(pos,
            base = base, index = index)
    s.expect(']')
    return result


@cython.cfunc
def p_subscript_list(s: PyrexScanner) -> tuple:
    is_single_value = True
    items = [p_subscript(s)]
    while s.sy == ',':
        is_single_value = False
        s.next()
        if s.sy == ']':
            break
        items.append(p_subscript(s))
    return items, is_single_value


#subscript: '.' '.' '.' | test | [test] ':' [test] [':' [test]]

@cython.cfunc
def p_subscript(s: PyrexScanner) -> list:
    # Parse a subscript and return a list of
    # 1, 2 or 3 ExprNodes, depending on how
    # many slice elements were encountered.
    start = p_slice_element(s, (':',))
    if s.sy != ':':
        return [start]
    s.next()
    stop = p_slice_element(s, (':', ',', ']'))
    if s.sy != ':':
        return [start, stop]
    s.next()
    step = p_slice_element(s, (':', ',', ']'))
    return [start, stop, step]


@cython.cfunc
def p_slice_element(s: PyrexScanner, follow_set):
    # Simple expression which may be missing iff
    # it is followed by something in follow_set.
    if s.sy not in follow_set:
        return p_test(s)
    else:
        return None


@cython.cfunc
def expect_ellipsis(s: PyrexScanner):
    s.expect('...')


@cython.cfunc
def make_slice_nodes(pos, subscripts) -> list:
    # Convert a list of subscripts as returned
    # by p_subscript_list into a list of ExprNodes,
    # creating SliceNodes for elements with 2 or
    # more components.
    result = []
    for subscript in subscripts:
        if len(subscript) == 1:
            result.append(subscript[0])
        else:
            result.append(make_slice_node(pos, *subscript))
    return result


@cython.ccall
def make_slice_node(pos, start, stop = None, step = None):
    if not start:
        start = ExprNodes.NoneNode(pos)
    if not stop:
        stop = ExprNodes.NoneNode(pos)
    if not step:
        step = ExprNodes.NoneNode(pos)
    return ExprNodes.SliceNode(pos,
        start = start, stop = stop, step = step)


#atom: '(' [yield_expr|testlist_comp] ')' | '[' [listmaker] ']' | '{' [dict_or_set_maker] '}' | '`' testlist '`' | NAME | NUMBER | STRING+

@cython.cfunc
def p_atom(s: PyrexScanner):
    pos = s.position()
    sy = s.sy
    if sy == '(':
        s.next()
        if s.sy == ')':
            result = ExprNodes.TupleNode(pos, args = [])
        elif s.sy == 'yield':
            result = p_yield_expression(s)
        else:
            result = p_testlist_comp(s)
        s.expect(')')
        return result
    elif sy == '[':
        return p_list_maker(s)
    elif sy == '{':
        return p_dict_or_set_maker(s)
    elif sy == '`':
        return p_backquote_expr(s)
    elif sy == '...':
        expect_ellipsis(s)
        return ExprNodes.EllipsisNode(pos)
    elif sy == 'INT':
        return p_int_literal(s)
    elif sy == 'FLOAT':
        value = s.systring
        s.next()
        return ExprNodes.FloatNode(pos, value = value)
    elif sy == 'IMAG':
        value = s.systring[:-1]
        s.next()
        return ExprNodes.ImagNode(pos, value = value)
    elif sy == 'BEGIN_STRING' or sy == 'BEGIN_FT_STRING':
        return p_atom_string(s)
    elif sy == 'IDENT':
        result = p_atom_ident_constants(s)
        if result is None:
            result = p_name(s, s.systring)
            s.next()
        return result
    else:
        s.error("Expected an identifier or literal")


@cython.cfunc
def p_atom_string(s: PyrexScanner):
    # s.sy == 'BEGIN_STRING' or s.sy == 'BEGIN_FT_STRING'
    pos = s.position()
    kind, bytes_value, unicode_value = p_cat_string_literal(s)
    if not kind:
        return ExprNodes.UnicodeNode(pos, value=unicode_value, bytes_value=bytes_value)
    kind_char: cython.Py_UCS4 = kind
    if kind_char == 'c':
        return ExprNodes.CharNode(pos, value=bytes_value)
    elif kind_char == 'u':
        return ExprNodes.UnicodeNode(pos, value=unicode_value, bytes_value=bytes_value)
    elif kind_char == 'b':
        return ExprNodes.BytesNode(pos, value=bytes_value)
    elif kind_char == 'f':
        return ExprNodes.JoinedStrNode(pos, values=unicode_value)
    elif kind_char == 't':
        # TODO
        return ExprNodes.TemplateStringNode(pos, values=unicode_value)
    else:
        # This is actually prevented by the scanner (Lexicon.py).
        s.error(f"invalid string kind '{kind}'")


@cython.cfunc
def p_atom_ident_constants(s: PyrexScanner):
    """
    Returns None if it isn't a special-cased named constant.
    Only calls s.next() if it successfully matches a named constant.
    """
    # s.sy == 'IDENT'
    pos = s.position()
    name = s.systring
    if name == "None":
        result = ExprNodes.NoneNode(pos)
    elif name == "True":
        result = ExprNodes.BoolNode(pos, value=True)
    elif name == "False":
        result = ExprNodes.BoolNode(pos, value=False)
    elif name == "NULL" and not s.in_python_file:
        result = ExprNodes.NullNode(pos)
    else:
        return None
    s.next()
    return result


@cython.cfunc
def p_int_literal(s: PyrexScanner):
    pos = s.position()
    value: str = cython.cast(str, s.systring)
    s.next()
    unsigned = ""
    longness = ""
    while value[-1] in "UuLl":
        if value[-1] in "Ll":
            longness += "L"
        else:
            unsigned += "U"
        value = value[:-1]
    # '3L' is ambiguous in Py2 but not in Py3.  '3U' and '3LL' are
    # illegal in Py2 Python files.  All suffixes are illegal in Py3
    # Python files.
    is_c_literal = None
    if unsigned:
        is_c_literal = True
    elif longness:
        if longness == 'LL' or s.context.language_level >= 3:
            is_c_literal = True
    if s.in_python_file:
        if is_c_literal:
            error(pos, "illegal integer literal syntax in Python source file")
        is_c_literal = False
    return ExprNodes.IntNode(pos,
                             is_c_literal = is_c_literal,
                             value = value,
                             unsigned = unsigned,
                             longness = longness)


@cython.cfunc
def p_name(s: PyrexScanner, name):
    pos = s.position()
    if not s.compile_time_expr and name in s.compile_time_env:
        value = s.compile_time_env.lookup_here(name)
        node = wrap_compile_time_constant(pos, value)
        if node is not None:
            return node
    return ExprNodes.NameNode(pos, name=name)


@cython.cfunc
def wrap_compile_time_constant(pos, value):
    if value is None:
        return ExprNodes.NoneNode(pos)
    elif value is Ellipsis:
        return ExprNodes.EllipsisNode(pos)
    elif isinstance(value, bool):
        return ExprNodes.BoolNode(pos, value=value)
    elif isinstance(value, int):
        return ExprNodes.IntNode(pos, value=repr(value), constant_result=value)
    elif isinstance(value, float):
        return ExprNodes.FloatNode(pos, value=repr(value), constant_result=value)
    elif isinstance(value, complex):
        node = ExprNodes.ImagNode(pos, value=repr(value.imag), constant_result=complex(0.0, value.imag))
        if value.real:
            # FIXME: should we care about -0.0 ?
            # probably not worth using the '-' operator for negative imag values
            node = ExprNodes.binop_node(
                pos, '+', ExprNodes.FloatNode(pos, value=repr(value.real), constant_result=value.real), node,
                constant_result=value)
        return node
    elif isinstance(value, str):
        return ExprNodes.UnicodeNode(pos, value=EncodedString(value))
    elif isinstance(value, bytes):
        bvalue = bytes_literal(value, 'ascii')  # actually: unknown encoding, but BytesLiteral requires one
        return ExprNodes.BytesNode(pos, value=bvalue, constant_result=value)
    elif isinstance(value, tuple):
        args = [wrap_compile_time_constant(pos, arg) for arg in value]
        if None in args:
            # error already reported
            return None
        return ExprNodes.TupleNode(pos, args=args)

    error(pos, "Invalid type for compile-time constant: %r (type %s)"
               % (value, value.__class__.__name__))
    return None


@cython.cfunc
def p_cat_string_literal(s: PyrexScanner) -> tuple:
    # A sequence of one or more adjacent string literals.
    # Returns (kind, bytes_value, unicode_value)
    # where kind in ('b', 'c', 'u', 'f', 't', '')
    pos = s.position()
    kind, bytes_value, unicode_value = p_string_literal(s)
    if kind == 'c' or (s.sy != 'BEGIN_STRING' and s.sy != 'BEGIN_FT_STRING'):
        return kind, bytes_value, unicode_value
    bstrings, ustrings, positions = [bytes_value], [unicode_value], [pos]
    bytes_value = unicode_value = None
    while s.sy == 'BEGIN_STRING' or s.sy == 'BEGIN_FT_STRING':
        pos = s.position()
        next_kind, next_bytes_value, next_unicode_value = p_string_literal(s)
        if next_kind == 'c':
            error(pos, "Cannot concatenate char literal with another string or char literal")
            continue
        elif next_kind != kind:
            # concatenating f strings and normal strings is allowed and leads to an f string
            if {kind, next_kind} in ({'f', 'u'}, {'f', ''}):
                kind = 'f'
            elif kind == 't' or next_kind == 't':
                error(pos, "cannot mix t-string literals with string or bytes literals")
                continue
            else:
                error(pos, "Cannot mix string literals of different types, expected %s'', got %s''" % (
                    kind, next_kind))
                continue
        bstrings.append(next_bytes_value)
        ustrings.append(next_unicode_value)
        positions.append(pos)
    # join and rewrap the partial literals
    if kind in ('b', 'c', '') or kind == 'u' and None not in bstrings:
        # Py3 enforced unicode literals are parsed as bytes/unicode combination
        bytes_value = bytes_literal(b''.join(bstrings), s.source_encoding)
    if kind in ('u', ''):
        unicode_value = EncodedString(''.join([u for u in ustrings if u is not None]))
    if kind == 'f':
        unicode_value = []
        for u, pos in zip(ustrings, positions):
            if isinstance(u, list):
                unicode_value += u
            else:
                # non-f-string concatenated into the f-string
                unicode_value.append(ExprNodes.UnicodeNode(pos, value=EncodedString(u)))
    if kind == 't':
        unicode_value = []
        for u in ustrings:
            unicode_value.extend(u)
    return kind, bytes_value, unicode_value


@cython.cfunc
def p_opt_string_literal(s: PyrexScanner, required_type: str = 'u'):
    if s.sy != 'BEGIN_STRING':
        return None
    pos = s.position()
    kind, bytes_value, unicode_value = p_string_literal(s, required_type)
    if required_type == 'u':
        if kind == 'f':
            s.error("f-string not allowed here", pos)
        return unicode_value
    elif required_type == 'b':
        return bytes_value
    else:
        s.error("internal parser configuration error")


@cython.cfunc
def check_for_non_ascii_characters(string) -> cython.bint:
    s = cython.cast(str, string)  # EncodedString
    for c in s:
        if c >= '\x80':
            return True
    return False


@cython.cfunc
def p_string_literal_shared_read(
        s: PyrexScanner, pos, chars, kind,
        is_raw: cython.bint):
    """
    Returns a string of non-escaped characters (if handled) or none.
    If passed an escape sequence returns an empty string.
    """
    sy = s.sy
    systr = s.systring
    result = systr
    is_python3_source: cython.bint = s.context.language_level >= 3
    # print "p_string_literal: sy =", sy, repr(s.systring) ###
    if sy == 'CHARS':
        chars.append(systr)
    elif sy == 'ESCAPE':
        # in Py2, 'ur' raw unicode strings resolve unicode escapes but nothing else
        if is_raw and (is_python3_source or kind != 'u' or len(systr) < 2 or systr[1] not in 'Uu'):
            chars.append(systr)
        else:
            result = ""
            _append_escape_sequence(kind, chars, systr, s)
    elif sy == 'NEWLINE':
        chars.append('\n')
    elif sy == 'EOF':
        s.error("Unclosed string literal", pos=pos)
    else:
        return None
    return result

@cython.cfunc
def _validate_kind_string(pos, systring: str) -> str:
    kind_string = systring.rstrip('"\'').lower()
    if len(kind_string) <= 1 or (len(kind_string) == 2 and kind_string in "rbrurfrtr"):
        return kind_string
    # Otherwise an error of some sort
    unique_string_prefixes = set(kind_string)
    if len(unique_string_prefixes) != len(kind_string):
        error(pos, 'Duplicate string prefix character')
    unique_string_prefixes.discard('r')
    unique_string_prefixes = sorted(unique_string_prefixes)
    if len(unique_string_prefixes) >= 2:
        error(pos, f'String prefixes {unique_string_prefixes[0]} and {unique_string_prefixes[1]} cannot be combined')
    else:
        error(pos, f'Invalid string prefix {kind_string}')
    return ''

@cython.cfunc
def p_string_literal(s: PyrexScanner, kind_override=None) -> tuple:
    # A single string or char literal.  Returns (kind, bvalue, uvalue)
    # where kind in ('b', 'c', 'u', 'f', '').  The 'bvalue' is the source
    # code byte sequence of the string literal, 'uvalue' is the
    # decoded Unicode string.  Either of the two may be None depending
    # on the 'kind' of string, only unprefixed strings have both
    # representations. In f-strings, the uvalue is a list of the Unicode
    # strings and f-string expressions that make up the f-string.
    # s.sy == 'BEGIN_STRING' or s.sy == 'BEGIN_FT_STRING'
    if s.sy == 'BEGIN_FT_STRING':
        assert kind_override is None
        return p_ft_string_literal(s)
    pos = s.position()
    is_python3_source: cython.bint = s.context.language_level >= 3
    has_non_ascii_literal_characters = False
    kind_string = _validate_kind_string(pos, s.systring)

    is_raw: cython.bint = 'r' in kind_string

    if 'c' in kind_string:
        # this should never happen, since the lexer does not allow combining c
        # with other prefix characters
        if len(kind_string) != 1:
            error(pos, 'Invalid string prefix for character literal')
        kind = 'c'
    elif 'b' in kind_string:
        kind = 'b'
    elif 'u' in kind_string:
        kind = 'u'
    else:
        kind = ''

    if kind == '' and kind_override is None and Future.unicode_literals in s.context.future_directives:
        chars = StringEncoding.StrLiteralBuilder(s.source_encoding)
        kind = 'u'
    else:
        if kind_override is not None and kind_override in 'ub':
            kind = kind_override
        if kind in ('u', 'f'):  # f-strings are scanned exactly like Unicode literals, but are parsed further later
            chars = StringEncoding.UnicodeLiteralBuilder()
        elif kind == '':
            chars = StringEncoding.StrLiteralBuilder(s.source_encoding)
        else:
            chars = StringEncoding.BytesLiteralBuilder(s.source_encoding)

    while 1:
        s.next()
        handled_chars = p_string_literal_shared_read(
            s, pos, chars, kind,
            is_raw=is_raw)
        if handled_chars is not None:
            if (not has_non_ascii_literal_characters and
                    is_python3_source and Future.unicode_literals in s.context.future_directives):
                has_non_ascii_literal_characters = check_for_non_ascii_characters(handled_chars)
            continue
        if s.sy == 'END_STRING':
            break
        else:
            s.error("Unexpected token %r:%r in string literal" % (
                s.sy, s.systring))

    if kind == 'c':
        unicode_value = None
        bytes_value = chars.getchar()
        if len(bytes_value) != 1:
            error(pos, "invalid character literal: %r" % bytes_value)
    else:
        bytes_value, unicode_value = chars.getstrings()
        if (has_non_ascii_literal_characters
                and is_python3_source and Future.unicode_literals in s.context.future_directives):
            # Python 3 forbids literal non-ASCII characters in byte strings
            if kind == 'b':
                s.error("bytes can only contain ASCII literal characters.", pos=pos)
            bytes_value = None
    s.next()
    return (kind, bytes_value, unicode_value)


@cython.cfunc
def p_read_ft_string_expression(s: PyrexScanner) -> str:
    strings = []
    while True:
        s.next()
        sy = s.sy
        if sy in ["END_FT_STRING_EXPR",
                    # probably an error, but handle it elsewhere
                   "EOF", None]:
            if sy == "END_FT_STRING_EXPR":
                s.next()
            return ''.join(strings)
        strings.append(s.systring)


@cython.cfunc
def p_ft_string_replacement_field(s: PyrexScanner,
                                is_raw: cython.bint, is_single_quoted: cython.bint,
                                tf_string_kind: cython.Py_UCS4) -> list:
    result = []
    conversion_char = format_spec = expr = None
    t_string_expression = None
    self_documenting = False

    bracket_pos = s.position()
    expr_pos = (bracket_pos[0], bracket_pos[1], bracket_pos[2]+1)
    expr_string = p_read_ft_string_expression(s)
    if not expr_string.strip():
        error(bracket_pos,
              f"empty expression not allowed in {tf_string_kind}-string")
        result = []
    else:
        original_scanner = s
        s = PyrexScanner(
            StringIO(expr_string),
            bracket_pos[0],
            parent_scanner=s,
            source_encoding=s.source_encoding,
            initial_pos=expr_pos
        )
        s.bracket_nesting_level += 1
        if s.sy == "INDENT":
            s.next()
        if s.sy == 'yield':
            expr = p_yield_expression(
                s,
                statement_terminators=statement_terminators | {':', '}', '!'})
        else:
            expr = p_testlist_star_expr(s)

        if s.sy == "=":
            self_documenting = True
            s.next()

        if s.sy == "!":
            # format conversion
            previous_pos = s.position()
            s.next()
            conversion_char = s.systring
            # validate the conversion char
            if conversion_char in ['}', ':', '']:
                error(s.position(), "missing conversion character")
            elif not ExprNodes.FormattedValueNode.find_conversion_func(conversion_char):
                error(s.position(), "invalid conversion character '%s'" % conversion_char)
                s.next()
            elif s.position()[2] != (previous_pos[2] + 1):
                error(s.position(), "f-string: conversion type must come right after the exclamation mark")
                s.next()
            else:
                s.next()

        if self_documenting or tf_string_kind == 't':
            if conversion_char is not None:
                expr_string, _ = expr_string.rsplit('!', 1)
            if tf_string_kind == 't':
                t_string_expression = ExprNodes.UnicodeNode(
                    pos=expr_pos,
                    value=StringEncoding.EncodedString(expr_string.rstrip().rstrip('=').rstrip())
                )
            if self_documenting:
                result.append(
                    ExprNodes.UnicodeNode(
                        pos=expr_pos,
                        value=StringEncoding.EncodedString(expr_string)
                    )
                )

        # Validate that the expression string has actually ended
        while s.sy == "NEWLINE" or s.sy == "DEDENT":
            s.next()
        if s.sy != "EOF":
            error(
                s.position(),
                f"Unexpected characters after {tf_string_kind}-string expression: {s.systring}")

        s = original_scanner

    if s.sy == ":":
        # full format spec
        pos = s.position()
        # Contents of format spec are handled closer to an f-string than a t-string
        # (even for t-strings).
        format_spec_contents = p_ft_string_middles(s, is_raw, is_single_quoted, is_format_string=True, tf_string_kind='f')
        format_spec = ExprNodes.JoinedStrNode(
            pos,
            values=format_spec_contents
        )
    if self_documenting and conversion_char is None and format_spec is None:
        conversion_char = 'r'

    if conversion_char is not None:
        conversion_char = StringEncoding.EncodedString(conversion_char)
    if tf_string_kind == 't':
        result.append(ExprNodes.TStringInterpolationNode(
            bracket_pos, value=expr, conversion_char=conversion_char,
            format_spec=format_spec, expression_str=t_string_expression
        ))
    else:
        result.append(ExprNodes.FormattedValueNode(
            bracket_pos, value=expr, conversion_char=conversion_char,
            format_spec=format_spec
        ))
    return result

@cython.cfunc
def p_ft_string_middles(s: PyrexScanner,
                        is_raw: cython.bint, is_single_quoted: cython.bint,
                        is_format_string: cython.bint,
                        tf_string_kind: cython.Py_UCS4) -> list:
    middles: list = []
    builder = StringEncoding.UnicodeLiteralBuilder()
    pos = s.position()
    while True:
        s.next()
        sy = s.sy

        handled_chars = p_string_literal_shared_read(
            s, pos, builder, "u",
            is_raw=is_raw)
        if handled_chars is not None:
            continue

        if builder.chars:
            middles.append(ExprNodes.UnicodeNode(pos, value=builder.getstring()))
            builder = StringEncoding.UnicodeLiteralBuilder()
        if sy == "{":
            fields = p_ft_string_replacement_field(
                s, is_raw, is_single_quoted, tf_string_kind=tf_string_kind)
            middles.extend(fields)
            if not s.sy == '}':
                s.expected('}')
            continue
        elif sy == "END_FT_STRING":
            break
        elif s.sy == '}':
            if is_format_string:
                break
            # otherwise it's an error, but the scanner has reported it
        else:
            error(
                s.position(),
                "Unexpected token %r:%r in %s-string literal" % (
                s.sy, s.systring, tf_string_kind))
    return middles

@cython.cfunc
def p_ft_string_literal(s: PyrexScanner) -> tuple:
    # s.sy == BEGIN_FT_STRING
    kind_string = _validate_kind_string(s.position(), s.systring)
    tf_string_kind: cython.Py_UCS4 = 't' if 't' in kind_string else 'f'
    is_raw: cython.bint = 'r' in kind_string
    quotes = s.systring.lstrip("rRbBuUfFtT")
    is_single_quoted: cython.bint = len(quotes) != 3
    middles = p_ft_string_middles(s, is_raw, is_single_quoted, is_format_string=False, tf_string_kind=tf_string_kind)
    if s.sy != "END_FT_STRING":
        s.expected(quotes)
    s.next()
    return tf_string_kind, None, middles


@cython.cfunc
def _append_escape_sequence(kind, builder, escape_sequence: str, s: PyrexScanner):
    if len(escape_sequence) < 2:
        builder.append("\\")  # invalid escape sequence, warned earlier
        return
    c = escape_sequence[1]
    if c in "01234567":
        builder.append_charval(int(escape_sequence[1:], 8))
    elif c in "'\"\\":
        builder.append(c)
    elif c in "abfnrtv":
        builder.append(StringEncoding.char_from_escape_sequence(escape_sequence))
    elif c == '\n':
        pass  # line continuation
    elif c == 'x':  # \xXX
        if len(escape_sequence) == 4:
            builder.append_charval(int(escape_sequence[2:], 16))
        else:
            s.error("Invalid hex escape '%s'" % escape_sequence, fatal=False)
    elif c in 'NUu' and kind in ('u', 'f', ''):  # \uxxxx, \Uxxxxxxxx, \N{...}
        chrval = -1
        if c == 'N':
            uchar = None
            try:
                uchar = lookup_unicodechar(escape_sequence[3:-1])
                chrval = ord(uchar)
            except KeyError:
                s.error("Unknown Unicode character name %s" %
                        repr(escape_sequence[3:-1]).lstrip('u'), fatal=False)
        elif len(escape_sequence) in (6, 10):
            chrval = int(escape_sequence[2:], 16)
            if chrval > 1114111:  # sys.maxunicode:
                s.error("Invalid unicode escape '%s'" % escape_sequence)
                chrval = -1
        else:
            s.error("Invalid unicode escape '%s'" % escape_sequence, fatal=False)
        if chrval >= 0:
            builder.append_uescape(chrval, escape_sequence)
    else:
        builder.append(escape_sequence)


# since PEP 448:
# list_display  ::=     "[" [listmaker] "]"
# listmaker     ::=     (named_test|star_expr) ( comp_for | (',' (named_test|star_expr))* [','] )
# comp_iter     ::=     comp_for | comp_if
# comp_for      ::=     ["async"] "for" expression_list "in" testlist [comp_iter]
# comp_if       ::=     "if" test [comp_iter]

@cython.cfunc
def p_list_maker(s: PyrexScanner):
    # s.sy == '['
    pos = s.position()
    s.next()
    if s.sy == ']':
        s.expect(']')
        return ExprNodes.ListNode(pos, args=[])

    expr = p_namedexpr_test_or_starred_expr(s)
    if s.sy in ('for', 'async'):
        if expr.is_starred:
            s.error("iterable unpacking cannot be used in comprehension")
        append = ExprNodes.ComprehensionAppendNode(pos, expr=expr)
        loop = p_comp_for(s, append)
        s.expect(']')
        return ExprNodes.ComprehensionNode(
            pos, loop=loop, append=append, type=Builtin.list_type,
            # list comprehensions leak their loop variable in Py2
            has_local_scope=s.context.language_level >= 3)

    # (merged) list literal
    if s.sy == ',':
        s.next()
        exprs = p_namedexpr_test_or_starred_expr_list(s, expr)
    else:
        exprs = [expr]
    s.expect(']')
    return ExprNodes.ListNode(pos, args=exprs)


@cython.cfunc
def p_comp_iter(s: PyrexScanner, body):
    if s.sy in ('for', 'async'):
        return p_comp_for(s, body)
    elif s.sy == 'if':
        return p_comp_if(s, body)
    else:
        # insert the 'append' operation into the loop
        return body


@cython.cfunc
def p_comp_for(s: PyrexScanner, body):
    pos = s.position()
    # [async] for ...
    is_async = False
    if s.sy == 'async':
        is_async = True
        s.next()

    # s.sy == 'for'
    s.expect('for')
    kw = p_for_bounds(s, allow_testlist=False, is_async=is_async)
    kw.update(else_clause=None, body=p_comp_iter(s, body), is_async=is_async)
    return Nodes.ForStatNode(pos, **kw)


@cython.cfunc
def p_comp_if(s: PyrexScanner, body):
    # s.sy == 'if'
    pos = s.position()
    s.next()
    # Note that Python 3.9+ is actually more restrictive here and Cython now follows
    # the Python 3.9+ behaviour: https://github.com/python/cpython/issues/86014
    # On Python <3.9 `[i for i in range(10) if lambda: i if True else 1]` was disallowed
    # but `[i for i in range(10) if lambda: i]` was allowed.
    # On Python >=3.9 they're both disallowed.
    test = p_or_test(s)
    return Nodes.IfStatNode(pos,
        if_clauses = [Nodes.IfClauseNode(pos, condition = test,
                                         body = p_comp_iter(s, body))],
        else_clause = None )


# since PEP 448:
#dictorsetmaker: ( ((test ':' test | '**' expr)
#                   (comp_for | (',' (test ':' test | '**' expr))* [','])) |
#                  ((test | star_expr)
#                   (comp_for | (',' (test | star_expr))* [','])) )

@cython.cfunc
def p_dict_or_set_maker(s: PyrexScanner):
    # s.sy == '{'
    pos = s.position()
    s.next()
    if s.sy == '}':
        s.next()
        return ExprNodes.DictNode(pos, key_value_pairs=[])

    parts = []
    target_type: cython.int = 0
    last_was_simple_item = False
    while True:
        if s.sy in ('*', '**'):
            # merged set/dict literal
            if target_type == 0:
                target_type = 1 if s.sy == '*' else 2  # 'stars'
            elif target_type != len(s.sy):
                s.error("unexpected %sitem found in %s literal" % (
                    s.sy, 'set' if target_type == 1 else 'dict'))
            s.next()
            if s.sy == '*':
                s.error("expected expression, found '*'")
            item = p_starred_expr(s)
            parts.append(item)
            last_was_simple_item = False
        else:
            item = p_test(s)
            if target_type == 0:
                target_type = 2 if s.sy == ':' else 1  # dict vs. set
            if target_type == 2:
                # dict literal
                s.expect(':')
                key = item
                value = p_test(s)
                item = ExprNodes.DictItemNode(key.pos, key=key, value=value)
            if last_was_simple_item:
                parts[-1].append(item)
            else:
                parts.append([item])
                last_was_simple_item = True

        if s.sy == ',':
            s.next()
            if s.sy == '}':
                break
        else:
            break

    if s.sy in ('for', 'async'):
        # dict/set comprehension
        if len(parts) == 1 and isinstance(parts[0], list) and len(parts[0]) == 1:
            item = parts[0][0]
            if target_type == 2:
                assert isinstance(item, ExprNodes.DictItemNode), type(item)
                comprehension_type = Builtin.dict_type
                append = ExprNodes.DictComprehensionAppendNode(
                    item.pos, key_expr=item.key, value_expr=item.value)
            else:
                comprehension_type = Builtin.set_type
                append = ExprNodes.ComprehensionAppendNode(item.pos, expr=item)
            loop = p_comp_for(s, append)
            s.expect('}')
            return ExprNodes.ComprehensionNode(pos, loop=loop, append=append, type=comprehension_type)
        else:
            # syntax error, try to find a good error message
            if len(parts) == 1 and not isinstance(parts[0], list):
                s.error("iterable unpacking cannot be used in comprehension")
            else:
                # e.g. "{1,2,3 for ..."
                s.expect('}')
            return ExprNodes.DictNode(pos, key_value_pairs=[])

    s.expect('}')
    if target_type == 1:
        # (merged) set literal
        items = []
        set_items = []
        for part in parts:
            if isinstance(part, list):
                set_items.extend(part)
            else:
                if set_items:
                    items.append(ExprNodes.SetNode(set_items[0].pos, args=set_items))
                    set_items = []
                items.append(part)
        if set_items:
            items.append(ExprNodes.SetNode(set_items[0].pos, args=set_items))
        if len(items) == 1 and items[0].is_set_literal:
            return items[0]
        return ExprNodes.MergedSequenceNode(pos, args=items, type=Builtin.set_type)
    else:
        # (merged) dict literal
        items = []
        dict_items = []
        for part in parts:
            if isinstance(part, list):
                dict_items.extend(part)
            else:
                if dict_items:
                    items.append(ExprNodes.DictNode(dict_items[0].pos, key_value_pairs=dict_items))
                    dict_items = []
                items.append(part)
        if dict_items:
            items.append(ExprNodes.DictNode(dict_items[0].pos, key_value_pairs=dict_items))
        if len(items) == 1 and items[0].is_dict_literal:
            return items[0]
        return ExprNodes.MergedDictNode(pos, keyword_args=items, reject_duplicates=False)


# NOTE: no longer in Py3 :)
@cython.cfunc
def p_backquote_expr(s: PyrexScanner):
    # s.sy == '`'
    pos = s.position()
    s.next()
    args = [p_test(s)]
    while s.sy == ',':
        s.next()
        args.append(p_test(s))
    s.expect('`')
    if len(args) == 1:
        arg = args[0]
    else:
        arg = ExprNodes.TupleNode(pos, args = args)
    return ExprNodes.BackquoteNode(pos, arg = arg)


@cython.cfunc
def p_simple_expr_list(s: PyrexScanner, expr=None) -> list:
    exprs: list = [expr] if expr is not None else []
    while s.sy not in expr_terminators:
        exprs.append( p_test(s) )
        if s.sy != ',':
            break
        s.next()
    return exprs


@cython.cfunc
def p_test_or_starred_expr_list(s: PyrexScanner, expr=None) -> list:
    exprs: list = [expr] if expr is not None else []
    while s.sy not in expr_terminators:
        exprs.append(p_test_or_starred_expr(s))
        if s.sy != ',':
            break
        s.next()
    return exprs


@cython.cfunc
def p_namedexpr_test_or_starred_expr_list(s: PyrexScanner, expr=None) -> list:
    exprs: list = [expr] if expr is not None else []
    while s.sy not in expr_terminators:
        exprs.append(p_namedexpr_test_or_starred_expr(s))
        if s.sy != ',':
            break
        s.next()
    return exprs


#testlist: test (',' test)* [',']

@cython.cfunc
def p_testlist(s: PyrexScanner):
    pos = s.position()
    expr = p_test(s)
    if s.sy == ',':
        s.next()
        exprs = p_simple_expr_list(s, expr)
        return ExprNodes.TupleNode(pos, args = exprs)
    else:
        return expr


# testlist_star_expr: (test|star_expr) ( comp_for | (',' (test|star_expr))* [','] )

@cython.cfunc
def p_testlist_star_expr(s: PyrexScanner):
    pos = s.position()
    expr = p_test_or_starred_expr(s)
    if s.sy == ',':
        s.next()
        exprs = p_test_or_starred_expr_list(s, expr)
        return ExprNodes.TupleNode(pos, args = exprs)
    else:
        return expr


# testlist_comp: (test|star_expr) ( comp_for | (',' (test|star_expr))* [','] )

@cython.cfunc
def p_testlist_comp(s: PyrexScanner):
    pos = s.position()
    expr = p_namedexpr_test_or_starred_expr(s)
    if s.sy == ',':
        s.next()
        exprs = p_namedexpr_test_or_starred_expr_list(s, expr)
        return ExprNodes.TupleNode(pos, args = exprs)
    elif s.sy in ('for', 'async'):
        return p_genexp(s, expr)
    else:
        return expr


@cython.cfunc
def p_genexp(s: PyrexScanner, expr):
    # s.sy == 'async' | 'for'
    loop = p_comp_for(s, Nodes.ExprStatNode(
        expr.pos, expr = ExprNodes.YieldExprNode(expr.pos, arg=expr)))
    return ExprNodes.GeneratorExpressionNode(expr.pos, loop=loop)


expr_terminators = cython.declare(frozenset, frozenset((
    ')', ']', '}', ':', '=', 'NEWLINE', 'EOF')))


#-------------------------------------------------------
#
#   Statements
#
#-------------------------------------------------------

@cython.cfunc
def p_global_statement(s: PyrexScanner):
    # assume s.sy == 'global'
    pos = s.position()
    s.next()
    names = p_ident_list(s)
    return Nodes.GlobalNode(pos, names = names)


@cython.cfunc
def p_nonlocal_statement(s: PyrexScanner):
    pos = s.position()
    s.next()
    names = p_ident_list(s)
    return Nodes.NonlocalNode(pos, names = names)


@cython.cfunc
def p_expression_or_assignment(s: PyrexScanner):
    expr = p_testlist_star_expr(s)
    has_annotation = False
    if s.sy == ':' and (expr.is_name or expr.is_subscript or expr.is_attribute):
        has_annotation = True
        s.next()
        expr.annotation = p_annotation(s)

    if s.sy == '=' and expr.is_starred:
        # This is a common enough error to make when learning Cython to let
        # it fail as early as possible and give a very clear error message.
        s.error("a starred assignment target must be in a list or tuple"
                " - maybe you meant to use an index assignment: var[0] = ...",
                pos=expr.pos)

    expr_list = [expr]
    while s.sy == '=':
        s.next()
        if s.sy == 'yield':
            expr = p_yield_expression(s)
        else:
            expr = p_testlist_star_expr(s)
        expr_list.append(expr)
    if len(expr_list) == 1:
        if re.match(r"([-+*/%^&|]|<<|>>|\*\*|//|@)=", s.sy):
            lhs = expr_list[0]
            if isinstance(lhs, ExprNodes.SliceIndexNode):
                # implementation requires IndexNode
                lhs = ExprNodes.IndexNode(
                    lhs.pos,
                    base=lhs.base,
                    index=make_slice_node(lhs.pos, lhs.start, lhs.stop))
            elif not isinstance(lhs, (ExprNodes.AttributeNode, ExprNodes.IndexNode, ExprNodes.NameNode)):
                error(lhs.pos, "Illegal operand for inplace operation.")
            operator = s.sy[:-1]
            s.next()
            if s.sy == 'yield':
                rhs = p_yield_expression(s)
            else:
                rhs = p_testlist(s)
            return Nodes.InPlaceAssignmentNode(lhs.pos, operator=operator, lhs=lhs, rhs=rhs)
        expr = expr_list[0]
        return Nodes.ExprStatNode(expr.pos, expr=expr)

    rhs = expr_list[-1]
    if len(expr_list) == 2:
        return Nodes.SingleAssignmentNode(rhs.pos, lhs=expr_list[0], rhs=rhs, first=has_annotation)
    else:
        return Nodes.CascadedAssignmentNode(rhs.pos, lhs_list=expr_list[:-1], rhs=rhs)


@cython.cfunc
def p_print_statement(s: PyrexScanner):
    # s.sy == 'print'
    pos = s.position()
    ends_with_comma: cython.bint = False
    s.next()
    if s.sy == '>>':
        s.next()
        stream = p_test(s)
        if s.sy == ',':
            s.next()
            ends_with_comma = s.sy in ('NEWLINE', 'EOF')
    else:
        stream = None
    args = []
    if s.sy not in ('NEWLINE', 'EOF'):
        args.append(p_test(s))
        while s.sy == ',':
            s.next()
            if s.sy in ('NEWLINE', 'EOF'):
                ends_with_comma = True
                break
            args.append(p_test(s))
    arg_tuple = ExprNodes.TupleNode(pos, args=args)
    return Nodes.PrintStatNode(pos,
        arg_tuple=arg_tuple, stream=stream,
        append_newline=not ends_with_comma)


@cython.cfunc
def p_exec_statement(s: PyrexScanner):
    # s.sy == 'exec'
    pos = s.position()
    s.next()
    code = p_bit_expr(s)
    if isinstance(code, ExprNodes.TupleNode):
        # Py3 compatibility syntax
        tuple_variant = True
        args = code.args
        if len(args) not in (2, 3):
            s.error("expected tuple of length 2 or 3, got length %d" % len(args),
                    pos=pos, fatal=False)
            args = [code]
    else:
        tuple_variant = False
        args = [code]
    if s.sy == 'in':
        if tuple_variant:
            s.error("tuple variant of exec does not support additional 'in' arguments",
                    fatal=False)
        s.next()
        args.append(p_test(s))
        if s.sy == ',':
            s.next()
            args.append(p_test(s))
    return Nodes.ExecStatNode(pos, args=args)


@cython.cfunc
def p_del_statement(s: PyrexScanner):
    # s.sy == 'del'
    pos = s.position()
    s.next()
    # FIXME: 'exprlist' in Python
    args = p_simple_expr_list(s)
    return Nodes.DelStatNode(pos, args = args)


@cython.cfunc
def p_pass_statement(s: PyrexScanner, with_newline: cython.bint = False):
    pos = s.position()
    s.expect('pass')
    if with_newline:
        s.expect_newline("Expected a newline", ignore_semicolon=True)
    return Nodes.PassStatNode(pos)


@cython.cfunc
def p_break_statement(s: PyrexScanner):
    # s.sy == 'break'
    pos = s.position()
    s.next()
    return Nodes.BreakStatNode(pos)


@cython.cfunc
def p_continue_statement(s: PyrexScanner):
    # s.sy == 'continue'
    pos = s.position()
    s.next()
    return Nodes.ContinueStatNode(pos)


@cython.cfunc
def p_return_statement(s: PyrexScanner):
    # s.sy == 'return'
    pos = s.position()
    s.next()
    if s.sy not in statement_terminators:
        value = p_testlist(s)
    else:
        value = None
    return Nodes.ReturnStatNode(pos, value = value)


@cython.cfunc
def p_raise_statement(s: PyrexScanner):
    # s.sy == 'raise'
    pos = s.position()
    s.next()
    exc_type = None
    exc_value = None
    exc_tb = None
    cause = None
    if s.sy not in statement_terminators:
        exc_type = p_test(s)
        if s.sy == ',':
            s.next()
            exc_value = p_test(s)
            if s.sy == ',':
                s.next()
                exc_tb = p_test(s)
        elif s.sy == 'from':
            s.next()
            cause = p_test(s)
    if exc_type or exc_value or exc_tb:
        return Nodes.RaiseStatNode(pos,
            exc_type = exc_type,
            exc_value = exc_value,
            exc_tb = exc_tb,
            cause = cause)
    else:
        return Nodes.ReraiseStatNode(pos)


@cython.cfunc
def p_import_statement(s: PyrexScanner):
    # s.sy in ('import', 'cimport')
    pos = s.position()
    kind = s.sy
    s.next()
    items = [p_dotted_name(s, as_allowed=True)]
    while s.sy == ',':
        s.next()
        items.append(p_dotted_name(s, as_allowed=True))
    stats = []
    is_absolute = Future.absolute_import in s.context.future_directives
    for pos, target_name, dotted_name, as_name in items:
        if kind == 'cimport':
            stat = Nodes.CImportStatNode(
                pos,
                module_name=dotted_name,
                as_name=as_name,
                is_absolute=is_absolute)
        else:
            stat = Nodes.SingleAssignmentNode(
                pos,
                lhs=ExprNodes.NameNode(pos, name=as_name or target_name),
                rhs=ExprNodes.ImportNode(
                    pos,
                    module_name=ExprNodes.IdentifierStringNode(pos, value=dotted_name),
                    is_import_as_name=bool(as_name),
                    level=0 if is_absolute else None,
                    imported_names=None))
        stats.append(stat)
    return Nodes.StatListNode(pos, stats=stats)


@cython.cfunc
def p_from_import_statement(s: PyrexScanner, first_statement: cython.bint = 0):
    # s.sy == 'from'
    pos = s.position()
    s.next()
    if s.sy in ('.', '...'):
        # count relative import level
        level = 0
        while s.sy in ('.', '...'):
            level += len(s.sy)
            s.next()
    else:
        level = None
    if level is not None and s.sy in ('import', 'cimport'):
        # we are dealing with "from .. import foo, bar"
        dotted_name_pos, dotted_name = s.position(), s.context.intern_ustring('')
    else:
        if level is None and Future.absolute_import in s.context.future_directives:
            level = 0
        (dotted_name_pos, _, dotted_name, _) = p_dotted_name(s, as_allowed=False)
    if s.sy not in ('import', 'cimport'):
        s.error("Expected 'import' or 'cimport'")
    kind = s.sy
    s.next()

    is_cimport = kind == 'cimport'
    is_parenthesized = False
    if s.sy == '*':
        imported_names = [(s.position(), s.context.intern_ustring("*"), None)]
        s.next()
    else:
        if s.sy == '(':
            is_parenthesized = True
            s.next()
        imported_names = [p_imported_name(s)]
    while s.sy == ',':
        s.next()
        if is_parenthesized and s.sy == ')':
            break
        imported_names.append(p_imported_name(s))
    if is_parenthesized:
        s.expect(')')
    if dotted_name == '__future__':
        if not first_statement:
            s.error("from __future__ imports must occur at the beginning of the file")
        elif level:
            s.error("invalid syntax")
        else:
            for (name_pos, name, as_name) in imported_names:
                if name == "braces":
                    s.error("not a chance", name_pos)
                    break
                try:
                    directive = getattr(Future, name)
                except AttributeError:
                    s.error("future feature %s is not defined" % name, name_pos)
                    break
                s.context.future_directives.add(directive)
        return Nodes.PassStatNode(pos)
    elif is_cimport:
        return Nodes.FromCImportStatNode(
            pos, module_name=dotted_name,
            relative_level=level,
            imported_names=imported_names)
    else:
        imported_name_strings = []
        items = []
        for (name_pos, name, as_name) in imported_names:
            imported_name_strings.append(
                ExprNodes.IdentifierStringNode(name_pos, value=name))
            items.append(
                (name, ExprNodes.NameNode(name_pos, name=as_name or name)))
        return Nodes.FromImportStatNode(pos,
            module = ExprNodes.ImportNode(dotted_name_pos,
                module_name = ExprNodes.IdentifierStringNode(pos, value = dotted_name),
                is_import_as_name = False,
                level = level,
                imported_names = imported_name_strings),
            items = items)


@cython.cfunc
def p_imported_name(s: PyrexScanner) -> tuple:
    pos = s.position()
    name = p_ident(s)
    as_name = p_as_name(s)
    return (pos, name, as_name)


@cython.cfunc
def p_dotted_name(s: PyrexScanner, as_allowed: cython.bint) -> tuple:
    pos = s.position()
    target_name = p_ident(s)
    as_name = None
    names = [target_name]
    while s.sy == '.':
        s.next()
        names.append(p_ident(s))
    if as_allowed:
        as_name = p_as_name(s)
    return (pos, target_name, s.context.intern_ustring('.'.join(names)), as_name)


@cython.cfunc
def p_as_name(s: PyrexScanner):
    if s.sy == 'IDENT' and s.systring == 'as':
        s.next()
        return p_ident(s)
    else:
        return None


@cython.cfunc
def p_assert_statement(s: PyrexScanner):
    # s.sy == 'assert'
    pos = s.position()
    s.next()
    cond = p_test(s)
    if s.sy == ',':
        s.next()
        value = p_test(s)
    else:
        value = None
    return Nodes.AssertStatNode(pos, condition=cond, value=value)


@cython.cfunc
def p_if_statement(s: PyrexScanner):
    # s.sy == 'if'
    pos = s.position()
    s.next()
    if_clauses = [p_if_clause(s)]
    while s.sy == 'elif':
        s.next()
        if_clauses.append(p_if_clause(s))
    else_clause = p_else_clause(s)
    return Nodes.IfStatNode(pos,
        if_clauses = if_clauses, else_clause = else_clause)


@cython.cfunc
def p_if_clause(s: PyrexScanner):
    pos = s.position()
    test = p_namedexpr_test(s)
    body = p_suite(s)
    return Nodes.IfClauseNode(pos,
        condition = test, body = body)


@cython.cfunc
def p_else_clause(s: PyrexScanner):
    if s.sy == 'else':
        s.next()
        return p_suite(s)
    else:
        return None


@cython.cfunc
def p_while_statement(s: PyrexScanner):
    # s.sy == 'while'
    pos = s.position()
    s.next()
    test = p_namedexpr_test(s)
    body = p_suite(s)
    else_clause = p_else_clause(s)
    return Nodes.WhileStatNode(pos,
        condition = test, body = body,
        else_clause = else_clause)


@cython.cfunc
def p_for_statement(s: PyrexScanner, is_async: cython.bint = False):
    # s.sy == 'for'
    pos = s.position()
    s.next()
    kw = p_for_bounds(s, allow_testlist=True, is_async=is_async)
    body = p_suite(s)
    else_clause = p_else_clause(s)
    kw.update(body=body, else_clause=else_clause, is_async=is_async)
    return Nodes.ForStatNode(pos, **kw)


@cython.cfunc
def p_for_bounds(s: PyrexScanner, allow_testlist: cython.bint = True, is_async: cython.bint = False) -> dict:
    target = p_for_target(s)
    if s.sy == 'in':
        s.next()
        iterator = p_for_iterator(s, allow_testlist, is_async=is_async)
        return dict(target=target, iterator=iterator)
    elif not s.in_python_file and not is_async:
        if s.sy == 'from':
            s.next()
            bound1 = p_bit_expr(s)
        else:
            # Support shorter "for a <= x < b" syntax
            bound1, target = target, None
        rel1 = p_for_from_relation(s)
        name2_pos = s.position()
        name2 = p_ident(s)
        rel2_pos = s.position()
        rel2 = p_for_from_relation(s)
        bound2 = p_bit_expr(s)
        step = p_for_from_step(s)
        if target is None:
            target = ExprNodes.NameNode(name2_pos, name = name2)
        else:
            if not target.is_name:
                error(target.pos,
                    "Target of for-from statement must be a variable name")
            elif name2 != target.name:
                error(name2_pos,
                    "Variable name in for-from range does not match target")
        if rel1[0] != rel2[0]:
            error(rel2_pos,
                "Relation directions in for-from do not match")
        return dict(target = target,
                    bound1 = bound1,
                    relation1 = rel1,
                    relation2 = rel2,
                    bound2 = bound2,
                    step = step,
                    )
    else:
        s.expect('in')
        return {}


@cython.cfunc
def p_for_from_relation(s: PyrexScanner):
    if s.sy in inequality_relations:
        op = s.sy
        s.next()
        return op
    else:
        s.error("Expected one of '<', '<=', '>' '>='")


@cython.cfunc
def p_for_from_step(s: PyrexScanner):
    if s.sy == 'IDENT' and s.systring == 'by':
        s.next()
        step = p_bit_expr(s)
        return step
    else:
        return None


inequality_relations = cython.declare(frozenset, frozenset((
    '<', '<=', '>', '>=')))


@cython.cfunc
def p_target(s: PyrexScanner, terminator: str):
    pos = s.position()
    expr = p_starred_expr(s)
    if s.sy == ',':
        s.next()
        exprs = [expr]
        while s.sy != terminator:
            exprs.append(p_starred_expr(s))
            if s.sy != ',':
                break
            s.next()
        return ExprNodes.TupleNode(pos, args = exprs)
    else:
        return expr


@cython.cfunc
def p_for_target(s: PyrexScanner):
    return p_target(s, 'in')


@cython.cfunc
def p_for_iterator(s: PyrexScanner, allow_testlist: cython.bint = True, is_async: cython.bint = False):
    pos = s.position()
    if allow_testlist:
        expr = p_testlist(s)
    else:
        expr = p_or_test(s)
    return (ExprNodes.AsyncIteratorNode if is_async else ExprNodes.IteratorNode)(pos, sequence=expr)


@cython.cfunc
def p_try_statement(s: PyrexScanner):
    # s.sy == 'try'
    pos = s.position()
    s.next()
    body = p_suite(s)
    except_clauses = []
    else_clause = None
    if s.sy in ('except', 'else'):
        while s.sy == 'except':
            except_clauses.append(p_except_clause(s))
        if s.sy == 'else':
            s.next()
            else_clause = p_suite(s)
        body = Nodes.TryExceptStatNode(pos,
            body = body, except_clauses = except_clauses,
            else_clause = else_clause)
        if s.sy != 'finally':
            return body
        # try-except-finally is equivalent to nested try-except/try-finally
    if s.sy == 'finally':
        s.next()
        finally_clause = p_suite(s)
        return Nodes.TryFinallyStatNode(pos,
            body = body, finally_clause = finally_clause)
    else:
        s.error("Expected 'except' or 'finally'")


@cython.cfunc
def p_except_clause(s: PyrexScanner):
    # s.sy == 'except'
    pos = s.position()
    s.next()
    exc_type = None
    exc_value = None
    is_except_as = False
    if s.sy != ':':
        exc_type = p_test(s)
        # normalise into list of single exception tests
        if isinstance(exc_type, ExprNodes.TupleNode):
            exc_type = exc_type.args
        else:
            exc_type = [exc_type]
        if s.sy == ',' or (s.sy == 'IDENT' and s.systring == 'as'
                           and s.context.language_level == 2):
            s.next()
            exc_value = p_test(s)
        elif s.sy == 'IDENT' and s.systring == 'as':
            # Py3 syntax requires a name here
            s.next()
            pos2 = s.position()
            name = p_ident(s)
            exc_value = ExprNodes.NameNode(pos2, name = name)
            is_except_as = True
    body = p_suite(s)
    return Nodes.ExceptClauseNode(pos,
        pattern = exc_type, target = exc_value,
        body = body, is_except_as=is_except_as)


@cython.cfunc
def p_include_statement(s: PyrexScanner, ctx):
    pos = s.position()
    s.next()  # 'include'
    unicode_include_file_name = p_string_literal(s, 'u')[2]
    s.expect_newline("Syntax error in include statement")
    if s.compile_time_eval:
        include_file_name = unicode_include_file_name
        include_file_path = s.context.find_include_file(include_file_name, pos)
        if include_file_path:
            s.included_files.append(include_file_name)
            source_desc = FileSourceDescriptor(include_file_path)
            with source_desc.get_file_object() as f:
                s2 = PyrexScanner(f, source_desc, s, source_encoding=f.encoding, parse_comments=s.parse_comments)
                tree = p_statement_list(s2, ctx)
            return tree
        else:
            return None
    else:
        return Nodes.PassStatNode(pos)


@cython.cfunc
def p_with_statement(s: PyrexScanner):
    s.next()  # 'with'
    if s.systring == 'template' and not s.in_python_file:
        node = p_with_template(s)
    else:
        node = p_with_items(s)
    return node


@cython.cfunc
def p_with_items(s: PyrexScanner, is_async: cython.bint = False):
    """
    Copied from CPython:
    | 'with' '(' a[asdl_withitem_seq*]=','.with_item+ ','? ')' ':' b=block {
        _PyAST_With(a, b, NULL, EXTRA) }
    | 'with' a[asdl_withitem_seq*]=','.with_item+ ':' tc=[TYPE_COMMENT] b=block {
        _PyAST_With(a, b, NEW_TYPE_COMMENT(p, tc), EXTRA) }
    Therefore the first thing to try is the bracket-enclosed
    version and if that fails try the regular version
    """
    brackets_succeeded = False
    items = ()  # unused, but static analysis fails to track that below
    if s.sy == '(':
        with tentatively_scan(s) as errors:
            s.next()
            items = p_with_items_list(s, is_async)
            s.expect(")")
            if s.sy != ":":
                # Fail - the message doesn't matter because we'll try the
                # non-bracket version so it'll never be shown
                s.error("")
        brackets_succeeded = not errors
    if not brackets_succeeded:
        # try the non-bracket version
        items = p_with_items_list(s, is_async)
    body = p_suite(s)
    for cls, pos, kwds in reversed(items):
        # construct the actual nodes now that we know what the body is
        body = cls(pos, body=body, **kwds)
    return body


@cython.cfunc
def p_with_items_list(s: PyrexScanner, is_async: cython.bint) -> list:
    items = []
    while True:
        items.append(p_with_item(s, is_async))
        if s.sy != ",":
            break
        s.next()
        if s.sy == ")":
            # trailing commas allowed
            break
    return items


@cython.cfunc
def p_with_item(s: PyrexScanner, is_async: cython.bint) -> tuple:
    # In contrast to most parsing functions, this returns a tuple of
    #  class, pos, kwd_dict
    # This is because GILStatNode does a reasonable amount of initialization in its
    # constructor, and requires "body" to be set, which we don't currently have
    pos = s.position()
    if not s.in_python_file and s.sy == 'IDENT' and s.systring in ('nogil', 'gil'):
        if is_async:
            s.error("with gil/nogil cannot be async")
        state = s.systring
        s.next()

        # support conditional gil/nogil
        condition = None
        if s.sy == '(':
            s.next()
            condition = p_test(s)
            s.expect(')')

        return Nodes.GILStatNode, pos, {"state": state, "condition": condition}
    else:
        manager = p_test(s)
        target = None
        if s.sy == 'IDENT' and s.systring == 'as':
            s.next()
            target = p_starred_expr(s)
        return Nodes.WithStatNode, pos, {"manager": manager, "target": target, "is_async": is_async}


@cython.cfunc
def p_with_template(s: PyrexScanner):
    pos = s.position()
    templates = []
    s.next()
    s.expect('[')
    templates.append(s.systring)
    s.next()
    while s.systring == ',':
        s.next()
        templates.append(s.systring)
        s.next()
    s.expect(']')
    if s.sy == ':':
        s.next()
        s.expect_newline("Syntax error in template function declaration")
        s.expect_indent()
        body_ctx = Ctx()
        body_ctx.templates = templates
        func_or_var = p_c_func_or_var_declaration(s, pos, body_ctx)
        s.expect_dedent()
        return func_or_var
    else:
        error(pos, "Syntax error in template function declaration")


@cython.cfunc
def p_simple_statement(s: PyrexScanner, first_statement: cython.bint = 0):
    #print "p_simple_statement:", s.sy, s.systring ###
    if s.sy == 'global':
        node = p_global_statement(s)
    elif s.sy == 'nonlocal':
        node = p_nonlocal_statement(s)
    elif s.sy == 'print':
        node = p_print_statement(s)
    elif s.sy == 'exec':
        node = p_exec_statement(s)
    elif s.sy == 'del':
        node = p_del_statement(s)
    elif s.sy == 'break':
        node = p_break_statement(s)
    elif s.sy == 'continue':
        node = p_continue_statement(s)
    elif s.sy == 'return':
        node = p_return_statement(s)
    elif s.sy == 'raise':
        node = p_raise_statement(s)
    elif s.sy in ('import', 'cimport'):
        node = p_import_statement(s)
    elif s.sy == 'from':
        node = p_from_import_statement(s, first_statement = first_statement)
    elif s.sy == 'yield':
        node = p_yield_statement(s)
    elif s.sy == 'assert':
        node = p_assert_statement(s)
    elif s.sy == 'pass':
        node = p_pass_statement(s)
    else:
        node = p_expression_or_assignment(s)
    return node


@cython.cfunc
def p_simple_statement_list(s: PyrexScanner, ctx, first_statement: cython.bint = 0):
    # Parse a series of simple statements on one line
    # separated by semicolons.
    stat = p_simple_statement(s, first_statement = first_statement)
    pos = stat.pos
    stats = []
    if not isinstance(stat, Nodes.PassStatNode):
        stats.append(stat)
    while s.sy == ';':
        #print "p_simple_statement_list: maybe more to follow" ###
        s.next()
        if s.sy in ('NEWLINE', 'EOF'):
            break
        stat = p_simple_statement(s, first_statement = first_statement)
        if isinstance(stat, Nodes.PassStatNode):
            continue
        stats.append(stat)
        first_statement = False

    if not stats:
        stat = Nodes.PassStatNode(pos)
    elif len(stats) == 1:
        stat = stats[0]
    else:
        stat = Nodes.StatListNode(pos, stats = stats)

    if s.sy not in ('NEWLINE', 'EOF'):
        # provide a better error message for users who accidentally write Cython code in .py files
        if isinstance(stat, Nodes.ExprStatNode):
            if stat.expr.is_name and stat.expr.name == 'cdef':
                s.error("The 'cdef' keyword is only allowed in Cython files (pyx/pxi/pxd)", pos)
    s.expect_newline("Syntax error in simple statement list")

    return stat


@cython.cfunc
def p_compile_time_expr(s: PyrexScanner):
    old = s.compile_time_expr
    s.compile_time_expr = 1
    expr = p_testlist(s)
    s.compile_time_expr = old
    return expr


@cython.cfunc
def p_DEF_statement(s: PyrexScanner):
    pos = s.position()
    denv = s.compile_time_env
    s.next()  # 'DEF'
    name = p_ident(s)
    s.expect('=')
    expr = p_compile_time_expr(s)
    if s.compile_time_eval:
        value = expr.compile_time_value(denv)
        #print "p_DEF_statement: %s = %r" % (name, value) ###
        denv.declare(name, value)
    s.expect_newline("Expected a newline", ignore_semicolon=True)
    return Nodes.PassStatNode(pos)


@cython.cfunc
def p_IF_statement(s: PyrexScanner, ctx):
    pos = s.position()
    saved_eval = s.compile_time_eval
    current_eval = saved_eval
    denv = s.compile_time_env
    result = None
    while 1:
        s.next()  # 'IF' or 'ELIF'
        expr = p_compile_time_expr(s)
        s.compile_time_eval = current_eval and bool(expr.compile_time_value(denv))
        body = p_suite(s, ctx)
        if s.compile_time_eval:
            result = body
            current_eval = 0
        if s.sy != 'ELIF':
            break
    if s.sy == 'ELSE':
        s.next()
        s.compile_time_eval = current_eval
        body = p_suite(s, ctx)
        if current_eval:
            result = body
    if not result:
        result = Nodes.PassStatNode(pos)
    s.compile_time_eval = saved_eval
    return result


@cython.cfunc
def p_statement(s: PyrexScanner, ctx, first_statement: cython.bint = False):
    cdef_flag: cython.bint = ctx.cdef_flag
    pos = s.position()
    decorators = None
    if s.sy == 'ctypedef':
        if ctx.level not in ('module', 'module_pxd'):
            s.error("ctypedef statement not allowed here")
        #if ctx.api:
        #    error(pos, "'api' not allowed with 'ctypedef'")
        return p_ctypedef_statement(s, ctx)
    elif s.sy == 'DEF':
        # We used to dep-warn about this but removed the warning again since
        # we don't have a good answer yet for all use cases.
        if s.context.compiler_directives.get("warn.deprecated.DEF", False):
            warning(pos,
                    "The 'DEF' statement  will be removed in a future Cython version. "
                    "Consider using global variables, constants, and in-place literals instead. "
                    "See https://github.com/cython/cython/issues/4310", level=1)
        return p_DEF_statement(s)
    elif s.sy == 'IF':
        if s.context.compiler_directives.get("warn.deprecated.IF", True):
            warning(pos,
                    "The 'IF' statement is deprecated and will be removed in a future Cython version. "
                    "Consider using runtime conditions or C macros instead. "
                    "See https://github.com/cython/cython/issues/4310", level=1)
        return p_IF_statement(s, ctx)
    elif s.sy == '@':
        if ctx.level not in ('module', 'class', 'c_class', 'function', 'property', 'module_pxd', 'c_class_pxd', 'other'):
            s.error('decorator not allowed here')
        s.level = ctx.level
        decorators = p_decorators(s)
        if not ctx.allow_struct_enum_decorator and s.sy not in ('def', 'cdef', 'cpdef', 'class', 'async'):
            if s.sy == 'IDENT' and s.systring == 'async':
                pass  # handled below
            else:
                s.error("Decorators can only be followed by functions or classes")
    elif s.sy == 'pass' and cdef_flag:
        # empty cdef block
        return p_pass_statement(s, with_newline=True)

    overridable = False
    if s.sy == 'cdef':
        cdef_flag = True
        s.next()
    elif s.sy == 'cpdef':
        cdef_flag = True
        overridable = True
        s.next()
    if cdef_flag:
        if ctx.level not in ('module', 'module_pxd', 'function', 'c_class', 'c_class_pxd'):
            s.error('cdef statement not allowed here')
        s.level = ctx.level
        node = p_cdef_statement(s, pos, ctx(overridable=overridable))
        if decorators is not None:
            tup = (Nodes.CFuncDefNode, Nodes.CVarDefNode, Nodes.CClassDefNode)
            if ctx.allow_struct_enum_decorator:
                tup += (Nodes.CStructOrUnionDefNode, Nodes.CEnumDefNode)
            if not isinstance(node, tup):
                s.error("Decorators can only be followed by functions or classes")
            node.decorators = decorators
        return node
    else:
        if ctx.api:
            s.error("'api' not allowed with this statement", fatal=False)
        elif s.sy == 'def':
            # def statements aren't allowed in pxd files, except
            # as part of a cdef class
            if ('pxd' in ctx.level) and (ctx.level != 'c_class_pxd'):
                s.error('def statement not allowed here')
            s.level = ctx.level
            return p_def_statement(s, decorators)
        elif s.sy == 'class':
            if ctx.level not in ('module', 'function', 'class', 'other'):
                s.error("class definition not allowed here")
            return p_class_statement(s, decorators)
        elif s.sy == 'include':
            if ctx.level not in ('module', 'module_pxd'):
                s.error("include statement not allowed here")
            return p_include_statement(s, ctx)
        elif ctx.level == 'c_class' and s.sy == 'IDENT' and s.systring == 'property':
            return p_property_decl(s)
        elif s.sy == 'pass' and ctx.level != 'property':
            return p_pass_statement(s, with_newline=True)
        else:
            if ctx.level in ('c_class_pxd', 'property'):
                node = p_ignorable_statement(s)
                if node is not None:
                    return node
                s.error("Executable statement not allowed here")
            if s.sy == 'if':
                return p_if_statement(s)
            elif s.sy == 'while':
                return p_while_statement(s)
            elif s.sy == 'for':
                return p_for_statement(s)
            elif s.sy == 'try':
                return p_try_statement(s)
            elif s.sy == 'with':
                return p_with_statement(s)
            elif s.sy == 'async':
                s.next()
                return p_async_statement(s, ctx, decorators)
            else:
                if s.sy == 'IDENT' and s.systring == 'async':
                    ident_name = s.systring
                    ident_pos = s.position()
                    # PEP 492 enables the async/await keywords when it spots "async def ..."
                    s.next()
                    if s.sy == 'def':
                        return p_async_statement(s, ctx, decorators)
                    elif decorators:
                        s.error("Decorators can only be followed by functions or classes")
                    s.put_back('IDENT', ident_name, ident_pos)  # re-insert original token
                if s.sy == 'IDENT' and s.systring == 'match':
                    # p_match_statement returns None on a "soft" initial failure
                    match_statement = p_match_statement(s, ctx)
                    if match_statement is not None:
                        return match_statement
                return p_simple_statement_list(s, ctx, first_statement=first_statement)


@cython.cfunc
def p_statement_list(s: PyrexScanner, ctx, first_statement: cython.bint = 0):
    # Parse a series of statements separated by newlines.
    pos = s.position()
    stats = []
    while s.sy not in ('DEDENT', 'EOF'):
        stat = p_statement(s, ctx, first_statement = first_statement)
        if isinstance(stat, Nodes.PassStatNode):
            continue
        stats.append(stat)
        first_statement = False
    if not stats:
        return Nodes.PassStatNode(pos)
    elif len(stats) == 1:
        return stats[0]
    else:
        return Nodes.StatListNode(pos, stats = stats)


@cython.cfunc
def p_suite(s: PyrexScanner, ctx=Ctx()):
    return p_suite_with_docstring(s, ctx, with_doc_only=False)[1]


@cython.cfunc
def p_suite_with_docstring(s: PyrexScanner, ctx, with_doc_only: cython.bint = False) -> tuple:
    s.expect(':')
    doc = None
    if s.sy == 'NEWLINE':
        s.next()
        s.expect_indent()
        if with_doc_only:
            doc = p_doc_string(s)
        body = p_statement_list(s, ctx)
        s.expect_dedent()
    else:
        if ctx.api:
            s.error("'api' not allowed with this statement", fatal=False)
        if ctx.level in ('module', 'class', 'function', 'other'):
            body = p_simple_statement_list(s, ctx)
        else:
            body = p_pass_statement(s)
            s.expect_newline("Syntax error in declarations", ignore_semicolon=True)
    if not with_doc_only:
        doc, body = _extract_docstring(body)
    return doc, body


@cython.cfunc
def p_positional_and_keyword_args(s: PyrexScanner, end_sy_set, templates = None) -> tuple:
    """
    Parses positional and keyword arguments. end_sy_set
    should contain any s.sy that terminate the argument list.
    Argument expansion (* and **) are not allowed.

    Returns: (positional_args, keyword_args)
    """
    positional_args = []
    keyword_args = []
    pos_idx = 0

    while s.sy not in end_sy_set:
        if s.sy == '*' or s.sy == '**':
            s.error('Argument expansion not allowed here.', fatal=False)

        parsed_type = False
        if s.sy == 'IDENT' and s.peek()[0] == '=':
            ident = s.systring
            s.next()  # s.sy is '='
            s.next()
            if looking_at_expr(s):
                arg = p_test(s)
            else:
                base_type = p_c_base_type(s, templates = templates)
                declarator = p_c_declarator(s, empty=True)
                arg = Nodes.CComplexBaseTypeNode(base_type.pos,
                    base_type = base_type, declarator = declarator)
                parsed_type = True
            keyword_node = ExprNodes.IdentifierStringNode(arg.pos, value=ident)
            keyword_args.append((keyword_node, arg))
            was_keyword = True

        else:
            if looking_at_expr(s):
                arg = p_test(s)
            else:
                base_type = p_c_base_type(s, templates = templates)
                declarator = p_c_declarator(s, empty=True)
                arg = Nodes.CComplexBaseTypeNode(base_type.pos,
                    base_type = base_type, declarator = declarator)
                parsed_type = True
            positional_args.append(arg)
            pos_idx += 1
            if len(keyword_args) > 0:
                s.error("Non-keyword arg following keyword arg",
                        pos=arg.pos)

        if s.sy != ',':
            if s.sy not in end_sy_set:
                if parsed_type:
                    s.error("Unmatched %s" % " or ".join(end_sy_set))
            break
        s.next()
    return positional_args, keyword_args


@cython.ccall
def p_c_base_type(s: PyrexScanner, nonempty: cython.bint = False, templates=None):
    if s.sy == '(':
        return p_c_complex_base_type(s, templates = templates)
    else:
        return p_c_simple_base_type(s, nonempty=nonempty, templates=templates)


@cython.cfunc
def p_calling_convention(s: PyrexScanner):
    if s.sy == 'IDENT' and s.systring in calling_convention_words:
        result = s.systring
        s.next()
        return result
    else:
        return EncodedString("")


calling_convention_words = cython.declare(frozenset, frozenset((
    "__stdcall", "__cdecl", "__fastcall")))


@cython.cfunc
def p_c_complex_base_type(s: PyrexScanner, templates = None):
    # s.sy == '('
    pos = s.position()
    s.next()
    base_type = p_c_base_type(s, templates=templates)
    declarator = p_c_declarator(s, empty=True)
    type_node = Nodes.CComplexBaseTypeNode(
        pos, base_type=base_type, declarator=declarator)
    if s.sy == ',':
        components = [type_node]
        while s.sy == ',':
            s.next()
            if s.sy == ')':
                break
            base_type = p_c_base_type(s, templates=templates)
            declarator = p_c_declarator(s, empty=True)
            components.append(Nodes.CComplexBaseTypeNode(
                pos, base_type=base_type, declarator=declarator))
        type_node = Nodes.CTupleBaseTypeNode(pos, components = components)

    s.expect(')')
    if s.sy == '[':
        if is_memoryviewslice_access(s):
            type_node = p_memoryviewslice_access(s, type_node)
        else:
            type_node = p_buffer_or_template(s, type_node, templates)
    return type_node


@cython.cfunc
def p_c_simple_base_type(s: PyrexScanner, nonempty: cython.bint, templates=None):
    is_basic = False
    signed = 1
    longness = 0
    complex = False
    module_path = []
    pos = s.position()

    # Handle const/volatile
    is_const = is_volatile = False
    while s.sy == 'IDENT':
        if s.systring == 'const':
            if is_const: error(pos, "Duplicate 'const'")
            is_const = True
        elif s.systring == 'volatile':
            if is_volatile: error(pos, "Duplicate 'volatile'")
            is_volatile = True
        else:
            break
        s.next()
    if is_const or is_volatile:
        base_type = p_c_base_type(s, nonempty=nonempty, templates=templates)
        if isinstance(base_type, Nodes.MemoryViewSliceTypeNode):
            # reverse order to avoid having to write "(const int)[:]"
            base_type.base_type_node = Nodes.CConstOrVolatileTypeNode(pos,
                base_type=base_type.base_type_node, is_const=is_const, is_volatile=is_volatile)
            return base_type
        return Nodes.CConstOrVolatileTypeNode(pos,
            base_type=base_type, is_const=is_const, is_volatile=is_volatile)

    if s.sy != 'IDENT':
        error(pos, "Expected an identifier, found '%s'" % s.sy)
    if looking_at_base_type(s):
        #print "p_c_simple_base_type: looking_at_base_type at", s.position()
        is_basic = True
        if s.sy == 'IDENT' and s.systring in special_basic_c_types:
            signed, longness = special_basic_c_types[s.systring]
            name = s.systring
            s.next()
        else:
            signed, longness = p_sign_and_longness(s)
            if s.sy == 'IDENT' and s.systring in basic_c_type_names:
                name = s.systring
                s.next()
            else:
                name = 'int'  # long [int], short [int], long [int] complex, etc.
        if s.sy == 'IDENT' and s.systring == 'complex':
            complex = True
            s.next()
    elif looking_at_dotted_name(s):
        #print "p_c_simple_base_type: looking_at_type_name at", s.position()
        name = s.systring
        s.next()
        while s.sy == '.':
            module_path.append(name)
            s.next()
            name = p_ident(s)
    else:
        name = s.systring
        name_pos = s.position()
        s.next()
        if nonempty and s.sy != 'IDENT':
            # Make sure this is not a declaration of a variable or function.
            if s.sy == '(':
                old_pos = s.position()
                s.next()
                if (s.sy == '*' or s.sy == '**' or s.sy == '&'
                        or (s.sy == 'IDENT' and s.systring in calling_convention_words)):
                    s.put_back('(', '(', old_pos)
                else:
                    s.put_back('(', '(', old_pos)
                    s.put_back('IDENT', name, name_pos)
                    name = None
            elif s.sy not in ('*', '**', '[', '&'):
                s.put_back('IDENT', name, name_pos)
                name = None

    type_node = Nodes.CSimpleBaseTypeNode(pos,
        name = name, module_path = module_path,
        is_basic_c_type = is_basic, signed = signed,
        complex = complex, longness = longness,
        templates = templates)

    #    declarations here.
    if s.sy == '[':
        if is_memoryviewslice_access(s):
            type_node = p_memoryviewslice_access(s, type_node)
        else:
            type_node = p_buffer_or_template(s, type_node, templates)

    if s.sy == '.':
        s.next()
        name = p_ident(s)
        type_node = Nodes.CNestedBaseTypeNode(pos, base_type = type_node, name = name)

    return type_node


@cython.cfunc
def p_buffer_or_template(s: PyrexScanner, base_type_node, templates):
    # s.sy == '['
    pos = s.position()
    s.next()
    # Note that buffer_positional_options_count=1, so the only positional argument is dtype.
    # For templated types, all parameters are types.
    positional_args, keyword_args = (
        p_positional_and_keyword_args(s, (']',), templates)
    )
    s.expect(']')

    if s.sy == '[':
        base_type_node = p_buffer_or_template(s, base_type_node, templates)

    keyword_dict = ExprNodes.DictNode(pos,
        key_value_pairs = [
            ExprNodes.DictItemNode(pos=key.pos, key=key, value=value)
            for key, value in keyword_args
        ])
    result = Nodes.TemplatedTypeNode(pos,
        positional_args = positional_args,
        keyword_args = keyword_dict,
        base_type_node = base_type_node)
    return result


@cython.cfunc
def is_memoryviewslice_access(s: PyrexScanner) -> cython.bint:
    # s.sy == '['
    # a memoryview slice declaration is distinguishable from a buffer access
    # declaration by the first entry in the bracketed list.  The buffer will
    # not have an unnested colon in the first entry; the memoryview slice will.
    saved = [(s.sy, s.systring, s.position())]
    s.next()
    retval = False
    if s.systring == ':':
        retval = True
    elif s.sy == 'INT':
        saved.append((s.sy, s.systring, s.position()))
        s.next()
        if s.sy == ':':
            retval = True

    for sv in saved[::-1]:
        s.put_back(*sv)

    return retval


@cython.cfunc
def p_memoryviewslice_access(s: PyrexScanner, base_type_node):
    # s.sy == '['
    pos = s.position()
    s.next()
    subscripts, _ = p_subscript_list(s)
    # make sure each entry in subscripts is a slice
    for subscript in subscripts:
        if len(subscript) < 2:
            s.error("An axis specification in memoryview declaration does not have a ':'.")
    s.expect(']')
    indexes = make_slice_nodes(pos, subscripts)
    result = Nodes.MemoryViewSliceTypeNode(pos,
            base_type_node = base_type_node,
            axes = indexes)
    return result


@cython.cfunc
def looking_at_name(s: PyrexScanner) -> cython.bint:
    return s.sy == 'IDENT' and s.systring not in calling_convention_words


@cython.cfunc
def looking_at_expr(s: PyrexScanner) -> cython.bint:
    if s.systring in base_type_start_words:
        return False
    elif s.sy == 'IDENT':
        is_type = False
        name = s.systring
        name_pos = s.position()
        dotted_path = []
        s.next()

        while s.sy == '.':
            s.next()
            dotted_path.append((s.systring, s.position()))
            s.expect('IDENT')

        saved = s.sy, s.systring, s.position()
        if s.sy == 'IDENT':
            is_type = True
        elif s.sy == '*' or s.sy == '**':
            s.next()
            is_type = s.sy in (')', ']')
            s.put_back(*saved)
        elif s.sy == '(':
            s.next()
            is_type = s.sy == '*'
            s.put_back(*saved)
        elif s.sy == '[':
            s.next()
            is_type = s.sy == ']' or not looking_at_expr(s)  # could be a nested template type
            s.put_back(*saved)

        dotted_path.reverse()
        for p in dotted_path:
            s.put_back('IDENT', *p)
            s.put_back('.', '.', p[1])  # gets the position slightly wrong

        s.put_back('IDENT', name, name_pos)
        return not is_type and saved[0]
    else:
        return True


@cython.cfunc
def looking_at_base_type(s: PyrexScanner) -> cython.bint:
    #print "looking_at_base_type?", s.sy, s.systring, s.position()
    return s.sy == 'IDENT' and s.systring in base_type_start_words


@cython.cfunc
def looking_at_dotted_name(s: PyrexScanner) -> cython.bint:
    if s.sy == 'IDENT':
        name = s.systring
        name_pos = s.position()
        s.next()
        result: cython.bint = s.sy == '.'
        s.put_back('IDENT', name, name_pos)
        return result
    else:
        return False


basic_c_type_names = cython.declare(frozenset, frozenset((
    "void", "char", "int", "float", "double", "bint")))

special_basic_c_types = cython.declare(dict, {
    # name : (signed, longness)
    "Py_UNICODE" : (0, 0),
    "Py_UCS4"    : (0, 0),
    "Py_hash_t"  : (2, 0),
    "Py_ssize_t" : (2, 0),
    "ssize_t"    : (2, 0),
    "size_t"     : (0, 0),
    "ptrdiff_t"  : (2, 0),
    "Py_tss_t"   : (1, 0),
})

sign_and_longness_words = cython.declare(frozenset, frozenset((
    "short", "long", "signed", "unsigned")))

base_type_start_words = cython.declare(
    frozenset,
    basic_c_type_names
    | sign_and_longness_words
    | frozenset(special_basic_c_types))

struct_enum_union = cython.declare(frozenset, frozenset((
    "struct", "union", "enum", "packed")))


@cython.cfunc
def p_sign_and_longness(s: PyrexScanner) -> tuple:
    signed = 1
    longness = 0
    while s.sy == 'IDENT' and s.systring in sign_and_longness_words:
        if s.systring == 'unsigned':
            signed = 0
        elif s.systring == 'signed':
            signed = 2
        elif s.systring == 'short':
            longness = -1
        elif s.systring == 'long':
            longness += 1
        s.next()
    return signed, longness


@cython.cfunc
def p_opt_cname(s: PyrexScanner):
    literal = p_opt_string_literal(s, 'u')
    if literal is not None:
        cname = EncodedString(literal)
        cname.encoding = s.source_encoding
    else:
        cname = None
    return cname


@cython.ccall
def p_c_declarator(s: PyrexScanner, ctx = Ctx(),
                   empty: cython.bint = False, is_type: cython.bint = False, cmethod_flag: cython.bint = False,
                   assignable: cython.bint = False, nonempty: cython.bint = False,
                   calling_convention_allowed: cython.bint = False):
    # If empty is true, the declarator must be empty. If nonempty is true,
    # the declarator must be nonempty. Otherwise we don't care.
    # If cmethod_flag is true, then if this declarator declares
    # a function, it's a C method of an extension type.
    pos = s.position()
    if s.sy == '(':
        s.next()
        if s.sy == ')' or looking_at_name(s):
            base = Nodes.CNameDeclaratorNode(pos, name=s.context.intern_ustring(""), cname=None)
            result = p_c_func_declarator(s, pos, ctx, base, cmethod_flag)
        else:
            result = p_c_declarator(s, ctx, empty = empty, is_type = is_type,
                                    cmethod_flag = cmethod_flag,
                                    nonempty = nonempty,
                                    calling_convention_allowed = True)
            s.expect(')')
    else:
        result = p_c_simple_declarator(s, ctx, empty, is_type, cmethod_flag,
                                       assignable, nonempty)
    if not calling_convention_allowed and result.calling_convention and s.sy != '(':
        error(s.position(), "%s on something that is not a function"
            % result.calling_convention)
    while s.sy in ('[', '('):
        pos = s.position()
        if s.sy == '[':
            result = p_c_array_declarator(s, result)
        else:  # sy == '('
            s.next()
            result = p_c_func_declarator(s, pos, ctx, result, cmethod_flag)
        cmethod_flag = 0
    return result


@cython.cfunc
def p_c_array_declarator(s: PyrexScanner, base):
    pos = s.position()
    s.next()  # '['
    if s.sy != ']':
        dim = p_testlist(s)
    else:
        dim = None
    s.expect(']')
    return Nodes.CArrayDeclaratorNode(pos, base = base, dimension = dim)


@cython.cfunc
def p_c_func_declarator(s: PyrexScanner, pos, ctx, base, cmethod_flag: cython.bint):
    # Opening paren has already been skipped
    args = p_c_arg_list(s, ctx, cmethod_flag = cmethod_flag,
                        nonempty_declarators = 0)
    ellipsis = p_optional_ellipsis(s)
    s.expect(')')
    nogil = p_nogil(s)
    exc_val, exc_check, exc_clause = p_exception_value_clause(s, ctx.visibility == 'extern')
    if nogil and exc_clause:
        warning(
            s.position(),
            "The keyword 'nogil' should appear at the end of the "
            "function signature line. Placing it before 'except' "
            "or 'noexcept' will be disallowed in a future version "
            "of Cython.",
            level=2
        )
    nogil = nogil or p_nogil(s)
    with_gil = p_with_gil(s)
    return Nodes.CFuncDeclaratorNode(pos,
        base = base, args = args, has_varargs = ellipsis,
        exception_value = exc_val, exception_check = exc_check,
        nogil = nogil or ctx.nogil or with_gil, with_gil = with_gil, has_explicit_exc_clause=exc_clause)


supported_overloaded_operators = cython.declare(frozenset, frozenset((
    '+', '-', '*', '/', '%',
    '++', '--', '~', '|', '&', '^', '<<', '>>', ',',
    '==', '!=', '>=', '>', '<=', '<',
    '[]', '()', '!', '=',
    'bool',
)))


@cython.cfunc
def p_c_simple_declarator(s: PyrexScanner, ctx,
                          empty: cython.bint, is_type: cython.bint, cmethod_flag: cython.bint,
                          assignable: cython.bint, nonempty: cython.bint):
    pos = s.position()
    calling_convention = p_calling_convention(s)
    if s.sy in ('*', '**'):
        # scanner returns '**' as a single token
        is_ptrptr = s.sy == '**'
        s.next()

        const_pos = s.position()
        is_const = s.systring == 'const' and s.sy == 'IDENT'
        if is_const:
            s.next()

        base = p_c_declarator(s, ctx, empty=empty, is_type=is_type,
                              cmethod_flag=cmethod_flag,
                              assignable=assignable, nonempty=nonempty)
        if is_const:
            base = Nodes.CConstDeclaratorNode(const_pos, base=base)
        if is_ptrptr:
            base = Nodes.CPtrDeclaratorNode(pos, base=base)
        result = Nodes.CPtrDeclaratorNode(pos, base=base)
    elif s.sy == '&' or (s.sy == '&&' and s.context.cpp):
        node_class = Nodes.CppRvalueReferenceDeclaratorNode if s.sy == '&&' else Nodes.CReferenceDeclaratorNode
        s.next()
        base = p_c_declarator(s, ctx, empty=empty, is_type=is_type,
                              cmethod_flag=cmethod_flag,
                              assignable=assignable, nonempty=nonempty)
        result = node_class(pos, base=base)
    else:
        rhs = None
        if s.sy == 'IDENT':
            name = s.systring
            if empty:
                error(s.position(), "Declarator should be empty")
            s.next()
            cname = p_opt_cname(s)
            if name != 'operator' and s.sy == '=' and assignable:
                s.next()
                rhs = p_test(s)
        else:
            if nonempty:
                error(s.position(), "Empty declarator")
            name = ""
            cname = None
        if cname is None and ctx.namespace is not None and nonempty:
            cname = ctx.namespace + "::" + name
        if name == 'operator' and ctx.visibility == 'extern' and nonempty:
            op = s.sy
            if [1 for c in op if c in '+-*/<=>!%&|([^~,']:
                s.next()
                # Handle diphthong operators.
                if op == '(':
                    s.expect(')')
                    op = '()'
                elif op == '[':
                    s.expect(']')
                    op = '[]'
                elif op in ('-', '+', '|', '&') and s.sy == op:
                    op *= 2       # ++, --, ...
                    s.next()
                elif s.sy == '=':
                    op += s.sy    # +=, -=, ...
                    s.next()
                if op not in supported_overloaded_operators:
                    s.error("Overloading operator '%s' not yet supported." % op,
                            fatal=False)
                name += op
            elif op == 'IDENT':
                op = s.systring
                if op not in supported_overloaded_operators:
                    s.error("Overloading operator '%s' not yet supported." % op,
                            fatal=False)
                name = name + ' ' + op
                s.next()
        result = Nodes.CNameDeclaratorNode(pos,
            name = name, cname = cname, default = rhs)
    result.calling_convention = calling_convention
    return result


@cython.cfunc
def p_nogil(s: PyrexScanner) -> cython.bint:
    if s.sy == 'IDENT' and s.systring == 'nogil':
        s.next()
        return True
    else:
        return False


@cython.cfunc
def p_with_gil(s: PyrexScanner) -> cython.bint:
    if s.sy == 'with':
        s.next()
        s.expect_keyword('gil')
        return True
    else:
        return False


@cython.cfunc
def p_exception_value_clause(s: PyrexScanner, is_extern: cython.bint) -> tuple:
    """
    Parse exception value clause.

    Maps clauses to exc_check / exc_value / exc_clause as follows:
     ______________________________________________________________________
    |                             |             |             |            |
    | Clause                      | exc_check   | exc_value   | exc_clause |
    | ___________________________ | ___________ | ___________ | __________ |
    |                             |             |             |            |
    | <nothing> (default func.)   | True        | None        | False      |
    | <nothing> (cdef extern)     | False       | None        | False      |
    | noexcept                    | False       | None        | True       |
    | except <val>                | False       | <val>       | True       |
    | except? <val>               | True        | <val>       | True       |
    | except *                    | True        | None        | True       |
    | except +                    | '+'         | None        | True       |
    | except +*                   | '+'         | '*'         | True       |
    | except +<PyErr>             | '+'         | <PyErr>     | True       |
    | ___________________________ | ___________ | ___________ | __________ |

    Note that the only reason we need `exc_clause` is to raise a
    warning when `'except'` or `'noexcept'` is placed after the
    `'nogil'` keyword.
    """
    exc_clause: cython.bint = False
    exc_val = None
    exc_check = False if is_extern else True

    if s.sy == 'IDENT' and s.systring == 'noexcept':
        exc_clause = True
        s.next()
        exc_check = False
    elif s.sy == 'except':
        exc_clause = True
        s.next()
        if s.sy == '*':
            exc_check = True
            s.next()
        elif s.sy == '+':
            exc_check = '+'
            plus_char_pos = s.position()[2]
            s.next()
            if s.sy == 'IDENT':
                name = s.systring
                if name == 'nogil':
                    if s.position()[2] == plus_char_pos + 1:
                        error(s.position(),
                              "'except +nogil' defines an exception handling function. Use 'except + nogil' for the 'nogil' modifier.")
                    # 'except + nogil' is parsed outside
                else:
                    exc_val = p_name(s, name)
                    s.next()
            elif s.sy == '*':
                exc_val = ExprNodes.CharNode(s.position(), value='*')
                s.next()
        else:
            if s.sy == '?':
                exc_check = True
                s.next()
            else:
                exc_check = False
            # exc_val can be non-None even if exc_check is False, c.f. "except -1"
            exc_val = p_test(s)

    return exc_val, exc_check, exc_clause


c_arg_list_terminators = cython.declare(frozenset, frozenset((
    '*', '**', '...', ')', ':', '/')))


@cython.ccall
def p_c_arg_list(s: PyrexScanner, ctx = Ctx(),
                 in_pyfunc: cython.bint = False, cmethod_flag: cython.bint = False,
                 nonempty_declarators: cython.bint = False, kw_only: cython.bint = False,
                 annotated: cython.bint = True) -> list:
    #  Comma-separated list of C argument declarations, possibly empty.
    #  May have a trailing comma.
    args = []
    is_self_arg = cmethod_flag
    while s.sy not in c_arg_list_terminators:
        args.append(p_c_arg_decl(s, ctx, in_pyfunc, is_self_arg,
            nonempty = nonempty_declarators, kw_only = kw_only,
            annotated = annotated))
        if s.sy != ',':
            break
        s.next()
        is_self_arg = 0
    return args


@cython.cfunc
def p_optional_ellipsis(s: PyrexScanner) -> cython.bint:
    if s.sy == '...':
        expect_ellipsis(s)
        return True
    else:
        return False


@cython.cfunc
def p_c_arg_decl(s: PyrexScanner, ctx, in_pyfunc: cython.bint, cmethod_flag: cython.bint = False,
                 nonempty: cython.bint = False,
                 kw_only: cython.bint = False, annotated: cython.bint = True):
    pos = s.position()
    not_none = or_none = False
    default = None
    annotation = None
    if s.in_python_file:
        # empty type declaration
        base_type = Nodes.CSimpleBaseTypeNode(pos,
            name = None, module_path = [],
            is_basic_c_type = False, signed = 0,
            complex = False, longness = 0,
            is_self_arg = cmethod_flag, templates = None)
    else:
        base_type = p_c_base_type(s, nonempty=nonempty)
    declarator = p_c_declarator(s, ctx, nonempty = nonempty)
    if s.sy in ('not', 'or') and not s.in_python_file:
        kind = s.sy
        s.next()
        if s.sy == 'IDENT' and s.systring == 'None':
            s.next()
        else:
            s.error("Expected 'None'")
        if not in_pyfunc:
            error(pos, "'%s None' only allowed in Python functions" % kind)
        or_none = kind == 'or'
        not_none = kind == 'not'
    if annotated and s.sy == ':':
        s.next()
        annotation = p_annotation(s)
    if s.sy == '=':
        s.next()
        if 'pxd' in ctx.level:
            if s.sy in ['*', '?']:
                # TODO(github/1736): Make this an error for inline declarations.
                default = ExprNodes.NoneNode(pos)
                s.next()
            elif 'inline' in ctx.modifiers:
                default = p_test(s)
            else:
                error(pos, "default values cannot be specified in pxd files, use ? or *")
        else:
            default = p_test(s)
    return Nodes.CArgDeclNode(pos,
        base_type = base_type,
        declarator = declarator,
        not_none = not_none,
        or_none = or_none,
        default = default,
        annotation = annotation,
        kw_only = kw_only)


@cython.cfunc
def p_annotation(s: PyrexScanner):
    """An annotation just has the "test" syntax, but also stores the string it came from

    Note that the string is *allowed* to be changed/processed (although isn't here)
    so may not exactly match the string generated by Python, and if it doesn't
    then it is not a bug.
    """
    pos = s.position()
    expr = p_test(s)
    return ExprNodes.AnnotationNode(pos, expr=expr)


@cython.cfunc
def p_api(s: PyrexScanner) -> cython.bint:
    if s.sy == 'IDENT' and s.systring == 'api':
        s.next()
        return True
    else:
        return False


@cython.cfunc
def p_cdef_statement(s: PyrexScanner, pos, ctx):
    ctx.visibility = p_visibility(s, ctx.visibility)
    ctx.api = ctx.api or p_api(s)
    if ctx.api:
        if ctx.visibility not in ('private', 'public'):
            error(pos, "Cannot combine 'api' with '%s'" % ctx.visibility)
    if (ctx.visibility == 'extern') and s.sy == 'from':
        return p_cdef_extern_block(s, pos, ctx)
    elif s.sy == 'import':
        s.next()
        return p_cdef_extern_block(s, pos, ctx)
    elif p_nogil(s):
        ctx.nogil = True
        if ctx.overridable:
            error(pos, "cdef blocks cannot be declared cpdef")
        return p_cdef_block(s, ctx)
    elif s.sy == ':':
        if ctx.overridable:
            error(pos, "cdef blocks cannot be declared cpdef")
        return p_cdef_block(s, ctx)
    elif s.sy == 'class':
        if ctx.level not in ('module', 'module_pxd'):
            error(pos, "Extension type definition not allowed here")
        if ctx.overridable:
            error(pos, "Extension types cannot be declared cpdef")
        return p_c_class_definition(s, pos, ctx)
    elif s.sy == 'IDENT' and s.systring == 'cppclass':
        return p_cpp_class_definition(s, pos, ctx)
    elif s.sy == 'IDENT' and s.systring in struct_enum_union:
        if ctx.level not in ('module', 'module_pxd'):
            error(pos, "C struct/union/enum definition not allowed here")
        if ctx.overridable:
            if s.systring != 'enum':
                error(pos, "C struct/union cannot be declared cpdef")
        return p_struct_enum(s, pos, ctx)
    elif s.sy == 'IDENT' and s.systring == 'fused':
        return p_fused_definition(s, pos, ctx)
    else:
        return p_c_func_or_var_declaration(s, pos, ctx)


@cython.cfunc
def p_cdef_block(s: PyrexScanner, ctx):
    return p_suite(s, ctx(cdef_flag = True))


@cython.cfunc
def p_cdef_extern_block(s: PyrexScanner, pos, ctx):
    if ctx.overridable:
        error(pos, "cdef extern blocks cannot be declared cpdef")
    include_file = None
    s.expect('from')
    if s.sy == '*':
        s.next()
    else:
        include_file = p_string_literal(s, 'u')[2]
    ctx = ctx(cdef_flag = True, visibility = 'extern')
    if s.systring == "namespace":
        s.next()
        ctx.namespace = p_string_literal(s, 'u')[2]
    if p_nogil(s):
        ctx.nogil = True

    # Use "docstring" as verbatim string to include
    verbatim_include, body = p_suite_with_docstring(s, ctx, True)

    return Nodes.CDefExternNode(pos,
        include_file = include_file,
        verbatim_include = verbatim_include,
        body = body,
        namespace = ctx.namespace)


@cython.cfunc
def p_c_enum_definition(s: PyrexScanner, pos, ctx):
    # s.sy == ident 'enum'
    s.next()

    scoped = False
    if s.context.cpp and (s.sy == 'class' or (s.sy == 'IDENT' and s.systring == 'struct')):
        scoped = True
        s.next()

    if s.sy == 'IDENT':
        name = s.systring
        s.next()
        cname = p_opt_cname(s)
        if cname is None and ctx.namespace is not None:
            cname = ctx.namespace + "::" + name
    else:
        name = cname = None
        if scoped:
            s.error("Unnamed scoped enum not allowed")

    if scoped and s.sy == '(':
        s.next()
        underlying_type = p_c_base_type(s)
        s.expect(')')
    else:
        underlying_type = Nodes.CSimpleBaseTypeNode(
            pos,
            name="int",
            module_path = [],
            is_basic_c_type = True,
            signed = 1,
            complex = False,
            longness = 0
        )

    s.expect(':')
    items = []

    doc = None
    if s.sy != 'NEWLINE':
        p_c_enum_line(s, ctx, items)
    else:
        s.next()  # 'NEWLINE'
        s.expect_indent()
        doc = p_doc_string(s)

        while s.sy not in ('DEDENT', 'EOF'):
            p_c_enum_line(s, ctx, items)

        s.expect_dedent()

    if not items and ctx.visibility != "extern":
        error(pos, "Empty enum definition not allowed outside a 'cdef extern from' block")

    return Nodes.CEnumDefNode(
        pos, name=name, cname=cname,
        scoped=scoped, items=items,
        underlying_type=underlying_type,
        typedef_flag=ctx.typedef_flag, visibility=ctx.visibility,
        create_wrapper=ctx.overridable,
        api=ctx.api, in_pxd=ctx.level == 'module_pxd', doc=doc)


@cython.cfunc
def p_c_enum_line(s: PyrexScanner, ctx, items: list):
    if s.sy != 'pass':
        p_c_enum_item(s, ctx, items)
        while s.sy == ',':
            s.next()
            if s.sy in ('NEWLINE', 'EOF'):
                break
            p_c_enum_item(s, ctx, items)
    else:
        s.next()
    s.expect_newline("Syntax error in enum item list")


@cython.cfunc
def p_c_enum_item(s: PyrexScanner, ctx, items: list):
    pos = s.position()
    name = p_ident(s)
    cname = p_opt_cname(s)
    if cname is None and ctx.namespace is not None:
        cname = ctx.namespace + "::" + name
    value = None
    if s.sy == '=':
        s.next()
        value = p_test(s)
    items.append(Nodes.CEnumDefItemNode(pos,
        name = name, cname = cname, value = value))


@cython.cfunc
def p_c_struct_or_union_definition(s: PyrexScanner, pos, ctx):
    packed = False
    if s.systring == 'packed':
        packed = True
        s.next()
        if s.sy != 'IDENT' or s.systring != 'struct':
            s.expected('struct')
    # s.sy == ident 'struct' or 'union'
    kind = s.systring
    s.next()
    name = p_ident(s)
    cname = p_opt_cname(s)
    if cname is None and ctx.namespace is not None:
        cname = ctx.namespace + "::" + name
    attributes = None
    if s.sy == ':':
        s.next()
        attributes = []
        if s.sy == 'pass':
            s.next()
            s.expect_newline("Expected a newline", ignore_semicolon=True)
        else:
            s.expect('NEWLINE')
            s.expect_indent()
            body_ctx = Ctx(visibility=ctx.visibility)
            while s.sy != 'DEDENT':
                if s.sy != 'pass':
                    attributes.append(
                        p_c_func_or_var_declaration(s, s.position(), body_ctx))
                else:
                    s.next()
                    s.expect_newline("Expected a newline")
            s.expect_dedent()

        if not attributes and ctx.visibility != "extern":
            error(pos, "Empty struct or union definition not allowed outside a 'cdef extern from' block")
    else:
        s.expect_newline("Syntax error in struct or union definition")

    return Nodes.CStructOrUnionDefNode(pos,
        name = name, cname = cname, kind = kind, attributes = attributes,
        typedef_flag = ctx.typedef_flag, visibility = ctx.visibility,
        api = ctx.api, in_pxd = ctx.level == 'module_pxd', packed = packed)


@cython.cfunc
def p_fused_definition(s: PyrexScanner, pos, ctx):
    """
    c(type)def fused my_fused_type:
        ...
    """
    # s.systring == 'fused'

    if ctx.level not in ('module', 'module_pxd'):
        error(pos, "Fused type definition not allowed here")

    s.next()
    name = p_ident(s)

    s.expect(":")
    s.expect_newline()
    s.expect_indent()

    types = []
    while s.sy != 'DEDENT':
        if s.sy != 'pass':
            #types.append(p_c_declarator(s))
            types.append(p_c_base_type(s))  #, nonempty=1))
        else:
            s.next()

        s.expect_newline()

    s.expect_dedent()

    if not types:
        error(pos, "Need at least one type")

    return Nodes.FusedTypeNode(pos, name=name, types=types)


@cython.cfunc
def p_struct_enum(s: PyrexScanner, pos, ctx):
    if s.systring == 'enum':
        return p_c_enum_definition(s, pos, ctx)
    else:
        return p_c_struct_or_union_definition(s, pos, ctx)


@cython.cfunc
def p_visibility(s: PyrexScanner, prev_visibility):
    visibility = prev_visibility
    if s.sy == 'IDENT' and s.systring in ('extern', 'public', 'readonly'):
        visibility = s.systring
        if prev_visibility != 'private' and visibility != prev_visibility:
            s.error("Conflicting visibility options '%s' and '%s'"
                % (prev_visibility, visibility), fatal=False)
        s.next()
    return visibility


@cython.cfunc
def p_c_modifiers(s: PyrexScanner) -> list:
    if s.sy == 'IDENT' and s.systring in ('inline',):
        modifier = s.systring
        s.next()
        return [modifier] + p_c_modifiers(s)
    return []


@cython.cfunc
def p_c_func_or_var_declaration(s: PyrexScanner, pos, ctx):
    cmethod_flag: cython.bint = ctx.level in ('c_class', 'c_class_pxd')
    modifiers = p_c_modifiers(s)
    base_type = p_c_base_type(s, nonempty=True, templates = ctx.templates)
    declarator = p_c_declarator(s, ctx(modifiers=modifiers), cmethod_flag = cmethod_flag,
                                assignable=True, nonempty =True)
    declarator.overridable = ctx.overridable

    if s.sy == 'IDENT' and s.systring == 'const' and ctx.level == 'cpp_class':
        s.next()
        is_const_method = True
    else:
        is_const_method = False

    if s.sy == '->':
        # Special enough to give a better error message and keep going.
        s.error(
            "Return type annotation is not allowed in cdef/cpdef signatures. "
            "Please define it before the function name, as in C signatures.",
            fatal=False)
        s.next()
        p_test(s)  # Keep going, but ignore result.

    if s.sy == ':':
        if ctx.level not in ('module', 'c_class', 'module_pxd', 'c_class_pxd', 'cpp_class') and not ctx.templates:
            s.error("C function definition not allowed here")
        doc, suite = p_suite_with_docstring(s, Ctx(level='function'))
        result = Nodes.CFuncDefNode(pos,
            visibility = ctx.visibility,
            base_type = base_type,
            declarator = declarator,
            body = suite,
            doc = doc,
            modifiers = modifiers,
            api = ctx.api,
            overridable = ctx.overridable,
            is_const_method = is_const_method)
    else:
        #if api:
        #    s.error("'api' not allowed with variable declaration")
        if is_const_method:
            declarator.is_const_method = is_const_method
        declarators = [declarator]
        while s.sy == ',':
            s.next()
            if s.sy == 'NEWLINE':
                break
            declarator = p_c_declarator(s, ctx, cmethod_flag = cmethod_flag,
                                        assignable=True, nonempty=True)
            declarators.append(declarator)
        doc_line = s.start_line + 1
        s.expect_newline("Syntax error in C variable declaration", ignore_semicolon=True)
        if ctx.level in ('c_class', 'c_class_pxd') and s.start_line == doc_line:
            doc = p_doc_string(s)
        else:
            doc = None
        result = Nodes.CVarDefNode(pos,
            visibility = ctx.visibility,
            base_type = base_type,
            declarators = declarators,
            in_pxd = ctx.level in ('module_pxd', 'c_class_pxd'),
            doc = doc,
            api = ctx.api,
            modifiers = modifiers,
            overridable = ctx.overridable)
    return result


@cython.cfunc
def p_ctypedef_statement(s: PyrexScanner, ctx):
    # s.sy == 'ctypedef'
    pos = s.position()
    s.next()
    visibility = p_visibility(s, ctx.visibility)
    api = p_api(s)
    ctx = ctx(typedef_flag=True, visibility = visibility)
    if api:
        ctx.api = True
    if s.sy == 'class':
        return p_c_class_definition(s, pos, ctx)
    elif s.sy == 'IDENT' and s.systring in struct_enum_union:
        return p_struct_enum(s, pos, ctx)
    elif s.sy == 'IDENT' and s.systring == 'fused':
        return p_fused_definition(s, pos, ctx)
    else:
        base_type = p_c_base_type(s, nonempty=True)
        declarator = p_c_declarator(s, ctx, is_type=True, nonempty=True)
        s.expect_newline("Syntax error in ctypedef statement", ignore_semicolon=True)
        return Nodes.CTypeDefNode(
            pos, base_type = base_type,
            declarator = declarator,
            visibility = visibility, api = api,
            in_pxd = ctx.level == 'module_pxd')


@cython.cfunc
def p_decorators(s: PyrexScanner) -> list:
    decorators = []
    while s.sy == '@':
        pos = s.position()
        s.next()
        decorator = p_namedexpr_test(s)
        decorators.append(Nodes.DecoratorNode(pos, decorator=decorator))
        s.expect_newline("Expected a newline after decorator")
    return decorators


@cython.cfunc
def _reject_cdef_modifier_in_py(s: PyrexScanner, name):
    """Step over incorrectly placed cdef modifiers (@see _CDEF_MODIFIERS) to provide a good error message for them.
    """
    if s.sy == 'IDENT' and name in _CDEF_MODIFIERS:
        # Special enough to provide a good error message.
        s.error("Cannot use cdef modifier '%s' in Python function signature. Use a decorator instead." % name, fatal=False)
        return p_ident(s)  # Keep going, in case there are other errors.
    return name


@cython.cfunc
def p_def_statement(s: PyrexScanner, decorators: list = None, is_async_def: cython.bint = False):
    # s.sy == 'def'
    pos = decorators[0].pos if decorators else s.position()
    # PEP 492 switches the async/await keywords on in "async def" functions
    if is_async_def:
        s.enter_async()
    s.next()
    name = _reject_cdef_modifier_in_py(s, p_ident(s))
    s.expect(
        '(',
        "Expected '(', found '%s'. Did you use cdef syntax in a Python declaration? "
        "Use decorators and Python type annotations instead." % (
            s.systring if s.sy == 'IDENT' else s.sy))
    args, star_arg, starstar_arg = p_varargslist(s, terminator=')')
    s.expect(')')
    _reject_cdef_modifier_in_py(s, s.systring)
    return_type_annotation = None
    if s.sy == '->':
        s.next()
        return_type_annotation = p_annotation(s)
        _reject_cdef_modifier_in_py(s, s.systring)

    doc, body = p_suite_with_docstring(s, Ctx(level='function'))
    if is_async_def:
        s.exit_async()

    return Nodes.DefNode(
        pos, name=name, args=args, star_arg=star_arg, starstar_arg=starstar_arg,
        doc=doc, body=body, decorators=decorators, is_async_def=is_async_def,
        return_type_annotation=return_type_annotation)


@cython.cfunc
def p_varargslist(s: PyrexScanner, terminator: cython.Py_UCS4 = ')', annotated: cython.bint = True) -> tuple:
    args = p_c_arg_list(s, in_pyfunc=True, nonempty_declarators=True,
                        annotated = annotated)
    star_arg = None
    starstar_arg = None
    if s.sy == '/':
        if len(args) == 0:
            s.error("Got zero positional-only arguments despite presence of "
                    "positional-only specifier '/'")
        s.next()
        # Mark all args to the left as pos only
        for arg in args:
            arg.pos_only = 1
        if s.sy == ',':
            s.next()
            args.extend(p_c_arg_list(
                s, in_pyfunc=True, nonempty_declarators=True, annotated = annotated))
        elif s.sy != terminator:
            s.error("Syntax error in Python function argument list")
    if s.sy == '*':
        s.next()
        if s.sy == 'IDENT':
            star_arg = p_py_arg_decl(s, annotated=annotated)
        if s.sy == ',':
            s.next()
            args.extend(p_c_arg_list(
                s, in_pyfunc =True, nonempty_declarators=True, kw_only=True, annotated = annotated))
        elif s.sy != terminator:
            s.error("Syntax error in Python function argument list")
    if s.sy == '**':
        s.next()
        starstar_arg = p_py_arg_decl(s, annotated=annotated)
    if s.sy == ',':
        s.next()
    return (args, star_arg, starstar_arg)


@cython.cfunc
def p_py_arg_decl(s: PyrexScanner, annotated: cython.bint = True):
    pos = s.position()
    name = p_ident(s)
    annotation = None
    if annotated and s.sy == ':':
        s.next()
        annotation = p_annotation(s)
    return Nodes.PyArgDeclNode(pos, name = name, annotation = annotation)


@cython.cfunc
def p_class_statement(s: PyrexScanner, decorators):
    # s.sy == 'class'
    pos = s.position()
    s.next()
    class_name = EncodedString(p_ident(s))
    class_name.encoding = s.source_encoding  # FIXME: why is this needed?
    arg_tuple = None
    keyword_dict = None
    if s.sy == '(':
        positional_args, keyword_args = p_call_parse_args(s, allow_genexp=False)
        arg_tuple, keyword_dict = p_call_build_packed_args(pos, positional_args, keyword_args)
    if arg_tuple is None:
        # XXX: empty arg_tuple
        arg_tuple = ExprNodes.TupleNode(pos, args=[])
    doc, body = p_suite_with_docstring(s, Ctx(level='class'))
    return Nodes.PyClassDefNode(
        pos, name=class_name,
        bases=arg_tuple,
        keyword_args=keyword_dict,
        doc=doc, body=body, decorators=decorators,
        force_py3_semantics=s.context.language_level >= 3)


@cython.cfunc
def p_c_class_definition(s: PyrexScanner, pos,  ctx):
    # s.sy == 'class'
    s.next()
    module_path = []
    class_name = p_ident(s)
    while s.sy == '.':
        s.next()
        module_path.append(class_name)
        class_name = p_ident(s)
    if module_path and ctx.visibility != 'extern':
        error(pos, "Qualified class name only allowed for 'extern' C class")
    if module_path and s.sy == 'IDENT' and s.systring == 'as':
        s.next()
        as_name = p_ident(s)
    else:
        as_name = class_name
    objstruct_name = None
    typeobj_name = None
    bases = None
    check_size = None
    if s.sy == '(':
        positional_args, keyword_args = p_call_parse_args(s, allow_genexp=False)
        if keyword_args:
            s.error("C classes cannot take keyword bases.")
        bases, _ = p_call_build_packed_args(pos, positional_args, keyword_args)
    if bases is None:
        bases = ExprNodes.TupleNode(pos, args=[])

    if s.sy == '[':
        if ctx.visibility not in ('public', 'extern') and not ctx.api:
            error(s.position(), "Name options only allowed for 'public', 'api', or 'extern' C class")
        objstruct_name, typeobj_name, check_size = p_c_class_options(s)
    if s.sy == ':':
        if ctx.level == 'module_pxd':
            body_level = 'c_class_pxd'
        else:
            body_level = 'c_class'
        doc, body = p_suite_with_docstring(s, Ctx(level=body_level))
    else:
        s.expect_newline("Syntax error in C class definition")
        doc = None
        body = None
    if ctx.visibility == 'extern':
        if not module_path:
            error(pos, "Module name required for 'extern' C class")
        if typeobj_name:
            error(pos, "Type object name specification not allowed for 'extern' C class")
    elif ctx.visibility == 'public':
        if not objstruct_name:
            error(pos, "Object struct name specification required for 'public' C class")
        if not typeobj_name:
            error(pos, "Type object name specification required for 'public' C class")
    elif ctx.visibility == 'private':
        if ctx.api:
            if not objstruct_name:
                error(pos, "Object struct name specification required for 'api' C class")
            if not typeobj_name:
                error(pos, "Type object name specification required for 'api' C class")
    else:
        error(pos, "Invalid class visibility '%s'" % ctx.visibility)
    return Nodes.CClassDefNode(pos,
        visibility = ctx.visibility,
        typedef_flag = ctx.typedef_flag,
        api = ctx.api,
        module_name = ".".join(module_path),
        class_name = class_name,
        as_name = as_name,
        bases = bases,
        objstruct_name = objstruct_name,
        typeobj_name = typeobj_name,
        check_size = check_size,
        in_pxd = ctx.level == 'module_pxd',
        doc = doc,
        body = body)


@cython.cfunc
def p_c_class_options(s: PyrexScanner) -> tuple:
    objstruct_name = None
    typeobj_name = None
    check_size = None
    s.expect('[')
    while 1:
        if s.sy != 'IDENT':
            break
        if s.systring == 'object':
            s.next()
            objstruct_name = p_ident(s)
        elif s.systring == 'type':
            s.next()
            typeobj_name = p_ident(s)
        elif s.systring == 'check_size':
            s.next()
            check_size = p_ident(s)
            if check_size not in ('ignore', 'warn', 'error'):
                s.error("Expected one of ignore, warn or error, found %r" % check_size)
        if s.sy != ',':
            break
        s.next()
    s.expect(']', "Expected 'object', 'type' or 'check_size'")
    return objstruct_name, typeobj_name, check_size


@cython.cfunc
def p_property_decl(s: PyrexScanner):
    pos = s.position()
    s.next()  # 'property'
    name = p_ident(s)
    doc, body = p_suite_with_docstring(
        s, Ctx(level='property'), with_doc_only=True)
    return Nodes.PropertyNode(pos, name=name, doc=doc, body=body)


@cython.cfunc
def p_ignorable_statement(s: PyrexScanner):
    """
    Parses any kind of ignorable statement that is allowed in .pxd files.
    """
    if s.sy == 'BEGIN_STRING':
        pos = s.position()
        string_node = p_atom(s)
        s.expect_newline("Syntax error in string", ignore_semicolon=True)
        return Nodes.ExprStatNode(pos, expr=string_node)
    return None


@cython.cfunc
def p_doc_string(s: PyrexScanner):
    if s.sy == 'BEGIN_STRING':
        pos = s.position()
        kind, bytes_result, unicode_result = p_cat_string_literal(s)
        s.expect_newline("Syntax error in doc string", ignore_semicolon=True)
        if kind in ('u', ''):
            return unicode_result
        warning(pos, "Python 3 requires docstrings to be unicode strings")
        return bytes_result
    else:
        return None


@cython.cfunc
def _extract_docstring(node) -> tuple:
    """
    Extract a docstring from a statement or from the first statement
    in a list.  Remove the statement if found.  Return a tuple
    (plain-docstring or None, node).
    """
    doc_node = None
    if node is None:
        pass
    elif isinstance(node, Nodes.ExprStatNode):
        if node.expr.is_string_literal:
            doc_node = node.expr
            node = Nodes.StatListNode(node.pos, stats=[])
    elif isinstance(node, Nodes.StatListNode) and node.stats:
        stats = node.stats
        if isinstance(stats[0], Nodes.ExprStatNode):
            if stats[0].expr.is_string_literal:
                doc_node = stats[0].expr
                del stats[0]

    if doc_node is None:
        doc = None
    elif isinstance(doc_node, ExprNodes.BytesNode):
        warning(node.pos,
                "Python 3 requires docstrings to be unicode strings")
        doc = doc_node.value
    else:
        doc = doc_node.value
    return doc, node


@cython.ccall
def p_code(s: PyrexScanner, level=None, ctx=Ctx):
    body = p_statement_list(s, ctx(level = level), first_statement=True)
    if s.sy != 'EOF':
        s.error("Syntax error in statement [%s,%s]" % (
            repr(s.sy), repr(s.systring)))
    return body


_match_compiler_directive_comment = cython.declare(object, re.compile(
    r"^#\s*cython\s*:\s*((\w|[.])+\s*=.*)$").match)


@cython.cfunc
def p_compiler_directive_comments(s: PyrexScanner) -> dict:
    result = {}
    while s.sy == 'commentline':
        pos = s.position()
        m = _match_compiler_directive_comment(s.systring)
        if m:
            directives_string = m.group(1).strip()
            try:
                new_directives = Options.parse_directive_list(directives_string, ignore_unknown=True)
            except ValueError as e:
                s.error(e.args[0], fatal=False)
                s.next()
                continue

            for name in new_directives:
                if name not in result:
                    pass
                elif Options.directive_types.get(name) is list:
                    result[name] += new_directives[name]
                    new_directives[name] = result[name]
                elif new_directives[name] == result[name]:
                    warning(pos, "Duplicate directive found: %s" % (name,))
                else:
                    s.error("Conflicting settings found for top-level directive %s: %r and %r" % (
                        name, result[name], new_directives[name]), pos=pos)

            if 'language_level' in new_directives:
                # Make sure we apply the language level already to the first token that follows the comments.
                s.context.set_language_level(new_directives['language_level'])
            if 'legacy_implicit_noexcept' in new_directives:
                s.context.legacy_implicit_noexcept = new_directives['legacy_implicit_noexcept']


            result.update(new_directives)

        s.next()
    return result


@cython.ccall
def p_module(s: PyrexScanner, pxd, full_module_name, ctx=Ctx):
    pos = s.position()

    directive_comments = p_compiler_directive_comments(s)
    s.parse_comments = False

    if s.context.language_level is None:
        s.context.set_language_level('3')

    level = 'module_pxd' if pxd else 'module'
    doc = p_doc_string(s)
    body = p_statement_list(s, ctx(level=level), first_statement=True)
    if s.sy != 'EOF':
        s.error("Syntax error in statement [%s,%s]" % (
            repr(s.sy), repr(s.systring)))
    return ModuleNode(pos, doc = doc, body = body,
                      full_module_name = full_module_name,
                      directive_comments = directive_comments)


@cython.cfunc
def p_template_definition(s: PyrexScanner) -> tuple:
    name = p_ident(s)
    if s.sy == '=':
        s.expect('=')
        s.expect('*')
        required = False
    else:
        required = True
    return name, required


@cython.cfunc
def p_cpp_class_definition(s: PyrexScanner, pos,  ctx):
    # s.sy == 'cppclass'
    s.next()
    class_name = p_ident(s)
    cname = p_opt_cname(s)
    if cname is None and ctx.namespace is not None:
        cname = ctx.namespace + "::" + class_name
    if s.sy == '.':
        error(pos, "Qualified class name not allowed C++ class")
    if s.sy == '[':
        s.next()
        templates = [p_template_definition(s)]
        while s.sy == ',':
            s.next()
            templates.append(p_template_definition(s))
        s.expect(']')
        template_names = [name for name, required in templates]
    else:
        templates = None
        template_names = None
    if s.sy == '(':
        s.next()
        base_classes = [p_c_base_type(s, templates = template_names)]
        while s.sy == ',':
            s.next()
            base_classes.append(p_c_base_type(s, templates = template_names))
        s.expect(')')
    else:
        base_classes = []
    if s.sy == '[':
        error(s.position(), "Name options not allowed for C++ class")
    nogil = p_nogil(s)
    if s.sy == ':':
        s.next()
        s.expect('NEWLINE')
        s.expect_indent()
        # Allow a cppclass to have docstrings. It will be discarded as comment.
        # The goal of this is consistency: we can make docstrings inside cppclass methods,
        # so why not on the cppclass itself ?
        p_doc_string(s)
        attributes = []
        body_ctx = Ctx(visibility = ctx.visibility, level='cpp_class', nogil=nogil or ctx.nogil)
        body_ctx.templates = template_names
        while s.sy != 'DEDENT':
            if s.sy != 'pass':
                attributes.append(p_cpp_class_attribute(s, body_ctx))
            else:
                s.next()
                s.expect_newline("Expected a newline")
        s.expect_dedent()
    else:
        attributes = None
        s.expect_newline("Syntax error in C++ class definition")
    return Nodes.CppClassNode(pos,
        name = class_name,
        cname = cname,
        base_classes = base_classes,
        visibility = ctx.visibility,
        in_pxd = ctx.level == 'module_pxd',
        attributes = attributes,
        templates = templates)


@cython.cfunc
def p_cpp_class_attribute(s: PyrexScanner, ctx):
    pos = s.position()
    decorators = None
    if s.sy == '@':
        decorators = p_decorators(s)
    if s.systring == 'cppclass':
        return p_cpp_class_definition(s, pos, ctx)
    elif s.systring == 'ctypedef':
        return p_ctypedef_statement(s, ctx)
    elif s.sy == 'IDENT' and s.systring in struct_enum_union:
        if s.systring != 'enum':
            return p_cpp_class_definition(s, pos, ctx)
        else:
            return p_struct_enum(s, pos, ctx)
    else:
        node = p_c_func_or_var_declaration(s, pos, ctx)
        if decorators is not None:
            tup = Nodes.CFuncDefNode, Nodes.CVarDefNode, Nodes.CClassDefNode
            if ctx.allow_struct_enum_decorator:
                tup += Nodes.CStructOrUnionDefNode, Nodes.CEnumDefNode
            if not isinstance(node, tup):
                s.error("Decorators can only be followed by functions or classes")
            node.decorators = decorators
        return node


@cython.cfunc
def p_match_statement(s: PyrexScanner, ctx):
    assert s.sy == "IDENT" and s.systring == "match"
    pos = s.position()
    with tentatively_scan(s) as errors:
        s.next()
        subject = p_namedexpr_test(s)
        subjects = None
        if s.sy == ",":
            subjects = [subject]
        while s.sy == ",":
            s.next()
            if s.sy == ":":
                break
            subjects.append(p_test(s))
        if subjects is not None:
            subject = ExprNodes.TupleNode(pos, args=subjects)
        s.expect(":")
    if errors:
        return None

    # at this stage we are committed to it being a match block so continue
    # outside "with tentatively_scan"
    # (I think this deviates from the PEG parser slightly, and it'd
    # backtrack on the whole thing)
    s.expect_newline()
    s.expect_indent()
    cases = []
    while s.sy != "DEDENT":
        cases.append(p_case_block(s, ctx))
    s.expect_dedent()
    return MatchCaseNodes.MatchNode(pos, subject=subject, cases=cases)


@cython.cfunc
def p_case_block(s: PyrexScanner, ctx):
    if not (s.sy == "IDENT" and s.systring == "case"):
        s.expected("case")
    s.next()
    pos = s.position()
    pattern = p_patterns(s)
    guard = None
    if s.sy == 'if':
        s.next()
        guard = p_test(s)
    body = p_suite(s, ctx)

    return MatchCaseNodes.MatchCaseNode(pos, pattern=pattern, body=body, guard=guard)


@cython.cfunc
def p_patterns(s: PyrexScanner):
    # note - in slight contrast to the name (which comes from the Python grammar),
    # returns a single pattern
    patterns = []
    seq = False
    pos = s.position()
    while True:
        with tentatively_scan(s) as errors:
            pattern = p_maybe_star_pattern(s)
        if errors:
            if patterns:
                break  # all is good provided we have at least 1 pattern
            else:
                e = errors[0]
                s.error(e.args[1], pos=e.args[0])
        patterns.append(pattern)

        if s.sy == ",":
            seq = True
            s.next()
            if s.sy in [":", "if"]:
                break  # common reasons to break
        else:
            break

    if seq:
        return MatchCaseNodes.MatchSequencePatternNode(pos, patterns=patterns)
    else:
        return patterns[0]


@cython.cfunc
def p_maybe_star_pattern(s: PyrexScanner):
    # For match case. Either star_pattern or pattern
    if s.sy == "*":
        # star pattern
        s.next()
        target = None
        if s.systring != "_":  # for match-case '_' is treated as a special wildcard
            target = p_pattern_capture_target(s)
        else:
            s.next()
        pattern = MatchCaseNodes.MatchAndAssignPatternNode(
            s.position(), target=target, is_star=True
        )
        return pattern
    else:
        pattern = p_pattern(s)
        return pattern


@cython.cfunc
def p_pattern(s: PyrexScanner):
    # try "as_pattern" then "or_pattern"
    # (but practically "as_pattern" starts with "or_pattern" too)
    patterns = []
    pos = s.position()
    while True:
        patterns.append(p_closed_pattern(s))
        if s.sy != "|":
            break
        s.next()

    if len(patterns) > 1:
        pattern = MatchCaseNodes.OrPatternNode(
            pos,
            alternatives=patterns
        )
    else:
        pattern = patterns[0]

    if s.sy == 'IDENT' and s.systring == 'as':
        s.next()
        with tentatively_scan(s) as errors:
            pattern.as_targets.append(p_pattern_capture_target(s))
        if errors and s.sy == "_":
            s.next()
            # make this a specific error
            return Nodes.ErrorNode(errors[0].args[0], what=errors[0].args[1])
        elif errors:
            with tentatively_scan(s):
                expr = p_test(s)
                return Nodes.ErrorNode(expr.pos, what="Invalid pattern target")
            s.error(errors[0])
    return pattern


@cython.cfunc
def p_closed_pattern(s: PyrexScanner):
    """
    The PEG parser specifies it as
    | literal_pattern
    | capture_pattern
    | wildcard_pattern
    | value_pattern
    | group_pattern
    | sequence_pattern
    | mapping_pattern
    | class_pattern

    For the sake avoiding too much backtracking, we know:
    * starts with "{" is a sequence_pattern
    * starts with "[" is a mapping_pattern
    * starts with "(" is a group_pattern or sequence_pattern
    * wildcard pattern is just identifier=='_'
    The rest are then tried in order with backtracking
    """
    if s.sy == 'IDENT' and s.systring == '_':
        pos = s.position()
        s.next()
        return MatchCaseNodes.MatchAndAssignPatternNode(pos)
    elif s.sy == '{':
        return p_mapping_pattern(s)
    elif s.sy == '[':
        return p_sequence_pattern(s)
    elif s.sy == '(':
        with tentatively_scan(s) as errors:
            result = p_group_pattern(s)
            if not errors:
                return result
        return p_sequence_pattern(s)

    with tentatively_scan(s) as errors:
        result = p_literal_pattern(s)
        if not errors:
            return result
    with tentatively_scan(s) as errors:
        result = p_capture_pattern(s)
        if not errors:
            return result
    with tentatively_scan(s) as errors:
        result = p_value_pattern(s)
        if not errors:
            return result
    return p_class_pattern(s)


@cython.cfunc
def p_literal_pattern(s: PyrexScanner):
    # a lot of duplication in this function with "p_atom"
    next_must_be_a_number = False
    sign = ''
    if s.sy == '-':
        sign = s.sy
        sign_pos = s.position()
        s.next()
        next_must_be_a_number = True

    sy = s.sy
    pos = s.position()

    res = None
    if sy == 'INT':
        res = p_int_literal(s)
    elif sy == 'FLOAT':
        value = s.systring
        s.next()
        res = ExprNodes.FloatNode(pos, value=value)

    if res is not None and sign == "-":
        res = ExprNodes.UnaryMinusNode(sign_pos, operand=res)

    if res is not None and s.sy in ['+', '-']:
        sign = s.sy
        s.next()
        if s.sy != 'IMAG':
            s.error("Expected imaginary number")
        else:
            add_pos = s.position()
            value = s.systring[:-1]
            s.next()
            res = ExprNodes.binop_node(
                add_pos,
                sign,
                operand1=res,
                operand2=ExprNodes.ImagNode(s.position(), value=value)
            )

    if res is None and sy == 'IMAG':
        value = s.systring[:-1]
        s.next()
        res = ExprNodes.ImagNode(pos, value=sign+value)
        if sign == "-":
            res = ExprNodes.UnaryMinusNode(sign_pos, operand=res)

    if res is not None:
        return MatchCaseNodes.MatchValuePatternNode(pos, value=res)

    if next_must_be_a_number:
        s.error("Expected a number")
    if sy == 'BEGIN_STRING':
        res = p_atom_string(s)
        # Whether f-strings are suitable is validated in PostParse.
        return MatchCaseNodes.MatchValuePatternNode(pos, value=res)
    elif sy == 'IDENT':
        # Note that p_atom_ident_constants includes NULL.
        # This is a deliberate Cython addition to the pattern matching specification
        result = p_atom_ident_constants(s)
        if result:
            return MatchCaseNodes.MatchValuePatternNode(pos, value=result, is_is_check=True)

    s.error("Failed to match literal")


@cython.cfunc
def p_capture_pattern(s: PyrexScanner):
    return MatchCaseNodes.MatchAndAssignPatternNode(
        s.position(),
        target=p_pattern_capture_target(s)
    )


@cython.cfunc
def p_value_pattern(s: PyrexScanner):
    if s.sy != "IDENT":
        s.error("Expected identifier")
    pos = s.position()
    res = p_name(s, s.systring)
    s.next()
    if s.sy != '.':
        s.error(".")
    while s.sy == '.':
        attr_pos = s.position()
        s.next()
        attr = p_ident(s)
        res = ExprNodes.AttributeNode(attr_pos, obj=res, attribute=attr)
    if s.sy in ['(', '=']:
        s.error("Unexpected symbol '%s'" % s.sy)
    return MatchCaseNodes.MatchValuePatternNode(pos, value=res)


@cython.cfunc
def p_group_pattern(s: PyrexScanner):
    s.expect("(")
    pattern = p_pattern(s)
    s.expect(")")
    return pattern


@cython.cfunc
def p_sequence_pattern(s: PyrexScanner):
    opener = s.sy
    pos = s.position()
    if opener in ['[', '(']:
        closer = ']' if opener == '[' else ')'
        s.next()
        # maybe_sequence_pattern and open_sequence_pattern
        patterns = []
        while s.sy != closer:
            patterns.append(p_maybe_star_pattern(s))
            if s.sy == ",":
                s.next()
            else:
                if opener == '(' and len(patterns) == 1:
                    s.error("tuple-like pattern of length 1 must finish with ','")
                break
        s.expect(closer)
        return MatchCaseNodes.MatchSequencePatternNode(pos, patterns=patterns)
    else:
        s.error("Expected '[' or '('")


@cython.cfunc
def p_mapping_pattern(s: PyrexScanner):
    pos = s.position()
    s.expect('{')
    if s.sy == '}':
        # trivial empty mapping
        s.next()
        return MatchCaseNodes.MatchMappingPatternNode(pos)

    double_star_capture_target = None
    items_patterns = []
    star_star_arg_pos = None
    while s.sy != '}':
        if double_star_capture_target and not star_star_arg_pos:
            star_star_arg_pos = s.position()
        if s.sy == '**':
            s.next()
            double_star_capture_target = p_pattern_capture_target(s)
        else:
            # key=(literal_expr | attr)
            with tentatively_scan(s) as errors:
                pattern = p_literal_pattern(s)
                key = pattern.value
            if errors:
                pattern = p_value_pattern(s)
                key = pattern.value
            s.expect(':')
            value = p_pattern(s)
            items_patterns.append((key, value))
        if s.sy != ',':
            break
        s.next()
    s.expect('}')

    if star_star_arg_pos is not None:
        return Nodes.ErrorNode(
            star_star_arg_pos,
            what = "** pattern must be the final part of a mapping pattern."
        )
    return MatchCaseNodes.MatchMappingPatternNode(
        pos,
        keys = [kv[0] for kv in items_patterns],
        value_patterns = [kv[1] for kv in items_patterns],
        double_star_capture_target = double_star_capture_target
    )


@cython.cfunc
def p_class_pattern(s: PyrexScanner):
    # start by parsing the class as name_or_attr
    pos = s.position()
    res = p_name(s, s.systring)
    s.next()
    while s.sy == '.':
        attr_pos = s.position()
        s.next()
        attr = p_ident(s)
        res = ExprNodes.AttributeNode(attr_pos, obj=res, attribute=attr)
    class_ = res

    s.expect("(")
    if s.sy == ")":
        # trivial case with no arguments matched
        s.next()
        return MatchCaseNodes.ClassPatternNode(pos, class_=class_)

    # parse the arguments
    positional_patterns = []
    keyword_patterns = []
    keyword_patterns_error = None
    while s.sy != ')':
        with tentatively_scan(s) as errors:
            positional_patterns.append(p_pattern(s))
        if not errors:
            if keyword_patterns:
                keyword_patterns_error = s.position()
        else:
            with tentatively_scan(s) as errors:
                keyword_patterns.append(p_keyword_pattern(s))
        if s.sy != ",":
            break
        s.next()
    s.expect(")")

    if keyword_patterns_error is not None:
        return Nodes.ErrorNode(
            keyword_patterns_error,
            what="Positional patterns follow keyword patterns"
        )
    return MatchCaseNodes.ClassPatternNode(
        pos, class_ = class_,
        positional_patterns = positional_patterns,
        keyword_pattern_names = [kv[0] for kv in keyword_patterns],
        keyword_pattern_patterns = [kv[1] for kv in keyword_patterns],
    )


@cython.cfunc
def p_keyword_pattern(s: PyrexScanner) -> tuple:
    if s.sy != "IDENT":
        s.error("Expected identifier")
    arg = p_name(s, s.systring)
    s.next()
    s.expect("=")
    value = p_pattern(s)
    return arg, value


@cython.cfunc
def p_pattern_capture_target(s: PyrexScanner):
    # any name but '_', and with some constraints on what follows
    if s.sy != 'IDENT':
        s.error("Expected identifier")
    if s.systring == '_':
        s.error("Pattern capture target cannot be '_'")
    target = p_name(s, s.systring)
    s.next()
    if s.sy in ['.', '(', '=']:
        s.error("Illegal next symbol '%s'" % s.sy)
    return target



#----------------------------------------------
#
#   Debugging
#
#----------------------------------------------

@cython.ccall
def print_parse_tree(f, node, level: cython.long, key = None):
    ind: str = "  " * level
    f.write(ind)
    if key:
        f.write(f"{key}: ")
    if not node:
        f.write("None\n")
    elif type(node) is tuple:
        f.write(f"({node[0]} @ {node[1]}\n")
        for item in node[2:]:
            print_parse_tree(f, item, level+1)
        f.write(f"{ind})\n")
    elif isinstance(node, Nodes.Node):
        try:
            tag = node.tag
        except AttributeError:
            tag = node.__class__.__name__
        f.write(f"{tag} @ {node.pos}\n")
        for name, value in sorted(node.__dict__.items()):
            if name != 'tag' and name != 'pos':
                print_parse_tree(f, value, level+1, name)
    elif type(node) is list:
        f.write("[\n")
        for item in node:
            print_parse_tree(f, item, level+1)
        f.write(f"{ind}]\n")
    else:
        f.write(f"{ind}{node}\n")
