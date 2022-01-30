from typing import Optional, Union, Any, Tuple, Iterable
from typing_extensions import Final

from mypy.nodes import (
    Expression, ComparisonExpr, OpExpr, MemberExpr, UnaryExpr, StarExpr, IndexExpr, LITERAL_YES,
    LITERAL_NO, NameExpr, LITERAL_TYPE, IntExpr, FloatExpr, ComplexExpr, StrExpr, BytesExpr,
    UnicodeExpr, ListExpr, TupleExpr, SetExpr, DictExpr, CallExpr, SliceExpr, CastExpr,
    ConditionalExpr, EllipsisExpr, YieldFromExpr, YieldExpr, RevealExpr, SuperExpr,
    TypeApplication, LambdaExpr, ListComprehension, SetComprehension, DictionaryComprehension,
    GeneratorExpr, BackquoteExpr, TypeVarExpr, TypeAliasExpr, NamedTupleExpr, EnumCallExpr,
    TypedDictExpr, NewTypeExpr, PromoteExpr, AwaitExpr, TempNode, AssignmentExpr, ParamSpecExpr
)
from mypy.visitor import ExpressionVisitor

# [Note Literals and literal_hash]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Mypy uses the term "literal" to refer to any expression built out of
# the following:
#
# * Plain literal expressions, like `1` (integer, float, string, etc.)
#
# * Compound literal expressions, like `(lit1, lit2)` (list, dict,
#   set, or tuple)
#
# * Operator expressions, like `lit1 + lit2`
#
# * Variable references, like `x`
#
# * Member references, like `lit.m`
#
# * Index expressions, like `lit[0]`
#
# A typical "literal" looks like `x[(i,j+1)].m`.
#
# An expression that is a literal has a `literal_hash`, with the
# following properties.
#
# * `literal_hash` is a Key: a tuple containing basic data types and
#   possibly other Keys. So it can be used as a key in a dictionary
#   that will be compared by value (as opposed to the Node itself,
#   which is compared by identity).
#
# * Two expressions have equal `literal_hash`es if and only if they
#   are syntactically equal expressions. (NB: Actually, we also
#   identify as equal expressions like `3` and `3.0`; is this a good
#   idea?)
#
# * The elements of `literal_hash` that are tuples are exactly the
#   subexpressions of the original expression (e.g. the base and index
#   of an index expression, or the operands of an operator expression).


def literal(e: Expression) -> int:
    if isinstance(e, ComparisonExpr):
        return min(literal(o) for o in e.operands)

    elif isinstance(e, OpExpr):
        return min(literal(e.left), literal(e.right))

    elif isinstance(e, (MemberExpr, UnaryExpr, StarExpr)):
        return literal(e.expr)

    elif isinstance(e, AssignmentExpr):
        return literal(e.target)

    elif isinstance(e, IndexExpr):
        if literal(e.index) == LITERAL_YES:
            return literal(e.base)
        else:
            return LITERAL_NO

    elif isinstance(e, NameExpr):
        return LITERAL_TYPE

    if isinstance(e, (IntExpr, FloatExpr, ComplexExpr, StrExpr, BytesExpr, UnicodeExpr)):
        return LITERAL_YES

    if literal_hash(e):
        return LITERAL_YES

    return LITERAL_NO


Key = Tuple[Any, ...]


def subkeys(key: Key) -> Iterable[Key]:
    return [elt for elt in key if isinstance(elt, tuple)]


def literal_hash(e: Expression) -> Optional[Key]:
    return e.accept(_hasher)


class _Hasher(ExpressionVisitor[Optional[Key]]):
    def visit_int_expr(self, e: IntExpr) -> Key:
        return ('Literal', e.value)

    def visit_str_expr(self, e: StrExpr) -> Key:
        return ('Literal', e.value, e.from_python_3)

    def visit_bytes_expr(self, e: BytesExpr) -> Key:
        return ('Literal', e.value)

    def visit_unicode_expr(self, e: UnicodeExpr) -> Key:
        return ('Literal', e.value)

    def visit_float_expr(self, e: FloatExpr) -> Key:
        return ('Literal', e.value)

    def visit_complex_expr(self, e: ComplexExpr) -> Key:
        return ('Literal', e.value)

    def visit_star_expr(self, e: StarExpr) -> Key:
        return ('Star', literal_hash(e.expr))

    def visit_name_expr(self, e: NameExpr) -> Key:
        # N.B: We use the node itself as the key, and not the name,
        # because using the name causes issues when there is shadowing
        # (for example, in list comprehensions).
        return ('Var', e.node)

    def visit_member_expr(self, e: MemberExpr) -> Key:
        return ('Member', literal_hash(e.expr), e.name)

    def visit_op_expr(self, e: OpExpr) -> Key:
        return ('Binary', e.op, literal_hash(e.left), literal_hash(e.right))

    def visit_comparison_expr(self, e: ComparisonExpr) -> Key:
        rest: Any = tuple(e.operators)
        rest += tuple(literal_hash(o) for o in e.operands)
        return ('Comparison',) + rest

    def visit_unary_expr(self, e: UnaryExpr) -> Key:
        return ('Unary', e.op, literal_hash(e.expr))

    def seq_expr(self, e: Union[ListExpr, TupleExpr, SetExpr], name: str) -> Optional[Key]:
        if all(literal(x) == LITERAL_YES for x in e.items):
            rest: Any = tuple(literal_hash(x) for x in e.items)
            return (name,) + rest
        return None

    def visit_list_expr(self, e: ListExpr) -> Optional[Key]:
        return self.seq_expr(e, 'List')

    def visit_dict_expr(self, e: DictExpr) -> Optional[Key]:
        if all(a and literal(a) == literal(b) == LITERAL_YES for a, b in e.items):
            rest: Any = tuple(
                (literal_hash(a) if a else None, literal_hash(b)) for a, b in e.items
            )
            return ("Dict",) + rest
        return None

    def visit_tuple_expr(self, e: TupleExpr) -> Optional[Key]:
        return self.seq_expr(e, 'Tuple')

    def visit_set_expr(self, e: SetExpr) -> Optional[Key]:
        return self.seq_expr(e, 'Set')

    def visit_index_expr(self, e: IndexExpr) -> Optional[Key]:
        if literal(e.index) == LITERAL_YES:
            return ('Index', literal_hash(e.base), literal_hash(e.index))
        return None

    def visit_assignment_expr(self, e: AssignmentExpr) -> Optional[Key]:
        return literal_hash(e.target)

    def visit_call_expr(self, e: CallExpr) -> None:
        return None

    def visit_slice_expr(self, e: SliceExpr) -> None:
        return None

    def visit_cast_expr(self, e: CastExpr) -> None:
        return None

    def visit_conditional_expr(self, e: ConditionalExpr) -> None:
        return None

    def visit_ellipsis(self, e: EllipsisExpr) -> None:
        return None

    def visit_yield_from_expr(self, e: YieldFromExpr) -> None:
        return None

    def visit_yield_expr(self, e: YieldExpr) -> None:
        return None

    def visit_reveal_expr(self, e: RevealExpr) -> None:
        return None

    def visit_super_expr(self, e: SuperExpr) -> None:
        return None

    def visit_type_application(self, e: TypeApplication) -> None:
        return None

    def visit_lambda_expr(self, e: LambdaExpr) -> None:
        return None

    def visit_list_comprehension(self, e: ListComprehension) -> None:
        return None

    def visit_set_comprehension(self, e: SetComprehension) -> None:
        return None

    def visit_dictionary_comprehension(self, e: DictionaryComprehension) -> None:
        return None

    def visit_generator_expr(self, e: GeneratorExpr) -> None:
        return None

    def visit_backquote_expr(self, e: BackquoteExpr) -> None:
        return None

    def visit_type_var_expr(self, e: TypeVarExpr) -> None:
        return None

    def visit_paramspec_expr(self, e: ParamSpecExpr) -> None:
        return None

    def visit_type_alias_expr(self, e: TypeAliasExpr) -> None:
        return None

    def visit_namedtuple_expr(self, e: NamedTupleExpr) -> None:
        return None

    def visit_enum_call_expr(self, e: EnumCallExpr) -> None:
        return None

    def visit_typeddict_expr(self, e: TypedDictExpr) -> None:
        return None

    def visit_newtype_expr(self, e: NewTypeExpr) -> None:
        return None

    def visit__promote_expr(self, e: PromoteExpr) -> None:
        return None

    def visit_await_expr(self, e: AwaitExpr) -> None:
        return None

    def visit_temp_node(self, e: TempNode) -> None:
        return None


_hasher: Final = _Hasher()
