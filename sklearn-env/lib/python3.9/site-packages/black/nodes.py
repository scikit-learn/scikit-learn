"""
blib2to3 Node/Leaf transformation-related utility functions.
"""

import sys
from typing import (
    Collection,
    Generic,
    Iterator,
    List,
    Optional,
    Set,
    Tuple,
    TypeVar,
    Union,
)

if sys.version_info < (3, 8):
    from typing_extensions import Final
else:
    from typing import Final

# lib2to3 fork
from blib2to3.pytree import Node, Leaf, type_repr
from blib2to3 import pygram
from blib2to3.pgen2 import token

from black.cache import CACHE_DIR
from black.strings import has_triple_quotes


pygram.initialize(CACHE_DIR)
syms = pygram.python_symbols


# types
T = TypeVar("T")
LN = Union[Leaf, Node]
LeafID = int
NodeType = int


WHITESPACE: Final = {token.DEDENT, token.INDENT, token.NEWLINE}
STATEMENT: Final = {
    syms.if_stmt,
    syms.while_stmt,
    syms.for_stmt,
    syms.try_stmt,
    syms.except_clause,
    syms.with_stmt,
    syms.funcdef,
    syms.classdef,
}
STANDALONE_COMMENT: Final = 153
token.tok_name[STANDALONE_COMMENT] = "STANDALONE_COMMENT"
LOGIC_OPERATORS: Final = {"and", "or"}
COMPARATORS: Final = {
    token.LESS,
    token.GREATER,
    token.EQEQUAL,
    token.NOTEQUAL,
    token.LESSEQUAL,
    token.GREATEREQUAL,
}
MATH_OPERATORS: Final = {
    token.VBAR,
    token.CIRCUMFLEX,
    token.AMPER,
    token.LEFTSHIFT,
    token.RIGHTSHIFT,
    token.PLUS,
    token.MINUS,
    token.STAR,
    token.SLASH,
    token.DOUBLESLASH,
    token.PERCENT,
    token.AT,
    token.TILDE,
    token.DOUBLESTAR,
}
STARS: Final = {token.STAR, token.DOUBLESTAR}
VARARGS_SPECIALS: Final = STARS | {token.SLASH}
VARARGS_PARENTS: Final = {
    syms.arglist,
    syms.argument,  # double star in arglist
    syms.trailer,  # single argument to call
    syms.typedargslist,
    syms.varargslist,  # lambdas
}
UNPACKING_PARENTS: Final = {
    syms.atom,  # single element of a list or set literal
    syms.dictsetmaker,
    syms.listmaker,
    syms.testlist_gexp,
    syms.testlist_star_expr,
}
TEST_DESCENDANTS: Final = {
    syms.test,
    syms.lambdef,
    syms.or_test,
    syms.and_test,
    syms.not_test,
    syms.comparison,
    syms.star_expr,
    syms.expr,
    syms.xor_expr,
    syms.and_expr,
    syms.shift_expr,
    syms.arith_expr,
    syms.trailer,
    syms.term,
    syms.power,
}
ASSIGNMENTS: Final = {
    "=",
    "+=",
    "-=",
    "*=",
    "@=",
    "/=",
    "%=",
    "&=",
    "|=",
    "^=",
    "<<=",
    ">>=",
    "**=",
    "//=",
}

IMPLICIT_TUPLE = {syms.testlist, syms.testlist_star_expr, syms.exprlist}
BRACKET = {token.LPAR: token.RPAR, token.LSQB: token.RSQB, token.LBRACE: token.RBRACE}
OPENING_BRACKETS = set(BRACKET.keys())
CLOSING_BRACKETS = set(BRACKET.values())
BRACKETS = OPENING_BRACKETS | CLOSING_BRACKETS
ALWAYS_NO_SPACE = CLOSING_BRACKETS | {token.COMMA, STANDALONE_COMMENT}


class Visitor(Generic[T]):
    """Basic lib2to3 visitor that yields things of type `T` on `visit()`."""

    def visit(self, node: LN) -> Iterator[T]:
        """Main method to visit `node` and its children.

        It tries to find a `visit_*()` method for the given `node.type`, like
        `visit_simple_stmt` for Node objects or `visit_INDENT` for Leaf objects.
        If no dedicated `visit_*()` method is found, chooses `visit_default()`
        instead.

        Then yields objects of type `T` from the selected visitor.
        """
        if node.type < 256:
            name = token.tok_name[node.type]
        else:
            name = str(type_repr(node.type))
        # We explicitly branch on whether a visitor exists (instead of
        # using self.visit_default as the default arg to getattr) in order
        # to save needing to create a bound method object and so mypyc can
        # generate a native call to visit_default.
        visitf = getattr(self, f"visit_{name}", None)
        if visitf:
            yield from visitf(node)
        else:
            yield from self.visit_default(node)

    def visit_default(self, node: LN) -> Iterator[T]:
        """Default `visit_*()` implementation. Recurses to children of `node`."""
        if isinstance(node, Node):
            for child in node.children:
                yield from self.visit(child)


def whitespace(leaf: Leaf, *, complex_subscript: bool) -> str:  # noqa: C901
    """Return whitespace prefix if needed for the given `leaf`.

    `complex_subscript` signals whether the given leaf is part of a subscription
    which has non-trivial arguments, like arithmetic expressions or function calls.
    """
    NO = ""
    SPACE = " "
    DOUBLESPACE = "  "
    t = leaf.type
    p = leaf.parent
    v = leaf.value
    if t in ALWAYS_NO_SPACE:
        return NO

    if t == token.COMMENT:
        return DOUBLESPACE

    assert p is not None, f"INTERNAL ERROR: hand-made leaf without parent: {leaf!r}"
    if t == token.COLON and p.type not in {
        syms.subscript,
        syms.subscriptlist,
        syms.sliceop,
    }:
        return NO

    prev = leaf.prev_sibling
    if not prev:
        prevp = preceding_leaf(p)
        if not prevp or prevp.type in OPENING_BRACKETS:
            return NO

        if t == token.COLON:
            if prevp.type == token.COLON:
                return NO

            elif prevp.type != token.COMMA and not complex_subscript:
                return NO

            return SPACE

        if prevp.type == token.EQUAL:
            if prevp.parent:
                if prevp.parent.type in {
                    syms.arglist,
                    syms.argument,
                    syms.parameters,
                    syms.varargslist,
                }:
                    return NO

                elif prevp.parent.type == syms.typedargslist:
                    # A bit hacky: if the equal sign has whitespace, it means we
                    # previously found it's a typed argument.  So, we're using
                    # that, too.
                    return prevp.prefix

        elif prevp.type in VARARGS_SPECIALS:
            if is_vararg(prevp, within=VARARGS_PARENTS | UNPACKING_PARENTS):
                return NO

        elif prevp.type == token.COLON:
            if prevp.parent and prevp.parent.type in {syms.subscript, syms.sliceop}:
                return SPACE if complex_subscript else NO

        elif (
            prevp.parent
            and prevp.parent.type == syms.factor
            and prevp.type in MATH_OPERATORS
        ):
            return NO

        elif (
            prevp.type == token.RIGHTSHIFT
            and prevp.parent
            and prevp.parent.type == syms.shift_expr
            and prevp.prev_sibling
            and prevp.prev_sibling.type == token.NAME
            and prevp.prev_sibling.value == "print"  # type: ignore
        ):
            # Python 2 print chevron
            return NO
        elif prevp.type == token.AT and p.parent and p.parent.type == syms.decorator:
            # no space in decorators
            return NO

    elif prev.type in OPENING_BRACKETS:
        return NO

    if p.type in {syms.parameters, syms.arglist}:
        # untyped function signatures or calls
        if not prev or prev.type != token.COMMA:
            return NO

    elif p.type == syms.varargslist:
        # lambdas
        if prev and prev.type != token.COMMA:
            return NO

    elif p.type == syms.typedargslist:
        # typed function signatures
        if not prev:
            return NO

        if t == token.EQUAL:
            if prev.type != syms.tname:
                return NO

        elif prev.type == token.EQUAL:
            # A bit hacky: if the equal sign has whitespace, it means we
            # previously found it's a typed argument.  So, we're using that, too.
            return prev.prefix

        elif prev.type != token.COMMA:
            return NO

    elif p.type == syms.tname:
        # type names
        if not prev:
            prevp = preceding_leaf(p)
            if not prevp or prevp.type != token.COMMA:
                return NO

    elif p.type == syms.trailer:
        # attributes and calls
        if t == token.LPAR or t == token.RPAR:
            return NO

        if not prev:
            if t == token.DOT:
                prevp = preceding_leaf(p)
                if not prevp or prevp.type != token.NUMBER:
                    return NO

            elif t == token.LSQB:
                return NO

        elif prev.type != token.COMMA:
            return NO

    elif p.type == syms.argument:
        # single argument
        if t == token.EQUAL:
            return NO

        if not prev:
            prevp = preceding_leaf(p)
            if not prevp or prevp.type == token.LPAR:
                return NO

        elif prev.type in {token.EQUAL} | VARARGS_SPECIALS:
            return NO

    elif p.type == syms.decorator:
        # decorators
        return NO

    elif p.type == syms.dotted_name:
        if prev:
            return NO

        prevp = preceding_leaf(p)
        if not prevp or prevp.type == token.AT or prevp.type == token.DOT:
            return NO

    elif p.type == syms.classdef:
        if t == token.LPAR:
            return NO

        if prev and prev.type == token.LPAR:
            return NO

    elif p.type in {syms.subscript, syms.sliceop}:
        # indexing
        if not prev:
            assert p.parent is not None, "subscripts are always parented"
            if p.parent.type == syms.subscriptlist:
                return SPACE

            return NO

        elif not complex_subscript:
            return NO

    elif p.type == syms.atom:
        if prev and t == token.DOT:
            # dots, but not the first one.
            return NO

    elif p.type == syms.dictsetmaker:
        # dict unpacking
        if prev and prev.type == token.DOUBLESTAR:
            return NO

    elif p.type in {syms.factor, syms.star_expr}:
        # unary ops
        if not prev:
            prevp = preceding_leaf(p)
            if not prevp or prevp.type in OPENING_BRACKETS:
                return NO

            prevp_parent = prevp.parent
            assert prevp_parent is not None
            if prevp.type == token.COLON and prevp_parent.type in {
                syms.subscript,
                syms.sliceop,
            }:
                return NO

            elif prevp.type == token.EQUAL and prevp_parent.type == syms.argument:
                return NO

        elif t in {token.NAME, token.NUMBER, token.STRING}:
            return NO

    elif p.type == syms.import_from:
        if t == token.DOT:
            if prev and prev.type == token.DOT:
                return NO

        elif t == token.NAME:
            if v == "import":
                return SPACE

            if prev and prev.type == token.DOT:
                return NO

    elif p.type == syms.sliceop:
        return NO

    return SPACE


def preceding_leaf(node: Optional[LN]) -> Optional[Leaf]:
    """Return the first leaf that precedes `node`, if any."""
    while node:
        res = node.prev_sibling
        if res:
            if isinstance(res, Leaf):
                return res

            try:
                return list(res.leaves())[-1]

            except IndexError:
                return None

        node = node.parent
    return None


def prev_siblings_are(node: Optional[LN], tokens: List[Optional[NodeType]]) -> bool:
    """Return if the `node` and its previous siblings match types against the provided
    list of tokens; the provided `node`has its type matched against the last element in
    the list.  `None` can be used as the first element to declare that the start of the
    list is anchored at the start of its parent's children."""
    if not tokens:
        return True
    if tokens[-1] is None:
        return node is None
    if not node:
        return False
    if node.type != tokens[-1]:
        return False
    return prev_siblings_are(node.prev_sibling, tokens[:-1])


def last_two_except(leaves: List[Leaf], omit: Collection[LeafID]) -> Tuple[Leaf, Leaf]:
    """Return (penultimate, last) leaves skipping brackets in `omit` and contents."""
    stop_after = None
    last = None
    for leaf in reversed(leaves):
        if stop_after:
            if leaf is stop_after:
                stop_after = None
            continue

        if last:
            return leaf, last

        if id(leaf) in omit:
            stop_after = leaf.opening_bracket
        else:
            last = leaf
    else:
        raise LookupError("Last two leaves were also skipped")


def parent_type(node: Optional[LN]) -> Optional[NodeType]:
    """
    Returns:
        @node.parent.type, if @node is not None and has a parent.
            OR
        None, otherwise.
    """
    if node is None or node.parent is None:
        return None

    return node.parent.type


def child_towards(ancestor: Node, descendant: LN) -> Optional[LN]:
    """Return the child of `ancestor` that contains `descendant`."""
    node: Optional[LN] = descendant
    while node and node.parent != ancestor:
        node = node.parent
    return node


def replace_child(old_child: LN, new_child: LN) -> None:
    """
    Side Effects:
        * If @old_child.parent is set, replace @old_child with @new_child in
        @old_child's underlying Node structure.
            OR
        * Otherwise, this function does nothing.
    """
    parent = old_child.parent
    if not parent:
        return

    child_idx = old_child.remove()
    if child_idx is not None:
        parent.insert_child(child_idx, new_child)


def container_of(leaf: Leaf) -> LN:
    """Return `leaf` or one of its ancestors that is the topmost container of it.

    By "container" we mean a node where `leaf` is the very first child.
    """
    same_prefix = leaf.prefix
    container: LN = leaf
    while container:
        parent = container.parent
        if parent is None:
            break

        if parent.children[0].prefix != same_prefix:
            break

        if parent.type == syms.file_input:
            break

        if parent.prev_sibling is not None and parent.prev_sibling.type in BRACKETS:
            break

        container = parent
    return container


def first_leaf_column(node: Node) -> Optional[int]:
    """Returns the column of the first leaf child of a node."""
    for child in node.children:
        if isinstance(child, Leaf):
            return child.column
    return None


def first_child_is_arith(node: Node) -> bool:
    """Whether first child is an arithmetic or a binary arithmetic expression"""
    expr_types = {
        syms.arith_expr,
        syms.shift_expr,
        syms.xor_expr,
        syms.and_expr,
    }
    return bool(node.children and node.children[0].type in expr_types)


def is_docstring(leaf: Leaf) -> bool:
    if prev_siblings_are(
        leaf.parent, [None, token.NEWLINE, token.INDENT, syms.simple_stmt]
    ):
        return True

    # Multiline docstring on the same line as the `def`.
    if prev_siblings_are(leaf.parent, [syms.parameters, token.COLON, syms.simple_stmt]):
        # `syms.parameters` is only used in funcdefs and async_funcdefs in the Python
        # grammar. We're safe to return True without further checks.
        return True

    return False


def is_empty_tuple(node: LN) -> bool:
    """Return True if `node` holds an empty tuple."""
    return (
        node.type == syms.atom
        and len(node.children) == 2
        and node.children[0].type == token.LPAR
        and node.children[1].type == token.RPAR
    )


def is_one_tuple(node: LN) -> bool:
    """Return True if `node` holds a tuple with one element, with or without parens."""
    if node.type == syms.atom:
        gexp = unwrap_singleton_parenthesis(node)
        if gexp is None or gexp.type != syms.testlist_gexp:
            return False

        return len(gexp.children) == 2 and gexp.children[1].type == token.COMMA

    return (
        node.type in IMPLICIT_TUPLE
        and len(node.children) == 2
        and node.children[1].type == token.COMMA
    )


def is_one_tuple_between(opening: Leaf, closing: Leaf, leaves: List[Leaf]) -> bool:
    """Return True if content between `opening` and `closing` looks like a one-tuple."""
    if opening.type != token.LPAR and closing.type != token.RPAR:
        return False

    depth = closing.bracket_depth + 1
    for _opening_index, leaf in enumerate(leaves):
        if leaf is opening:
            break

    else:
        raise LookupError("Opening paren not found in `leaves`")

    commas = 0
    _opening_index += 1
    for leaf in leaves[_opening_index:]:
        if leaf is closing:
            break

        bracket_depth = leaf.bracket_depth
        if bracket_depth == depth and leaf.type == token.COMMA:
            commas += 1
            if leaf.parent and leaf.parent.type in {
                syms.arglist,
                syms.typedargslist,
            }:
                commas += 1
                break

    return commas < 2


def is_walrus_assignment(node: LN) -> bool:
    """Return True iff `node` is of the shape ( test := test )"""
    inner = unwrap_singleton_parenthesis(node)
    return inner is not None and inner.type == syms.namedexpr_test


def is_simple_decorator_trailer(node: LN, last: bool = False) -> bool:
    """Return True iff `node` is a trailer valid in a simple decorator"""
    return node.type == syms.trailer and (
        (
            len(node.children) == 2
            and node.children[0].type == token.DOT
            and node.children[1].type == token.NAME
        )
        # last trailer can be an argument-less parentheses pair
        or (
            last
            and len(node.children) == 2
            and node.children[0].type == token.LPAR
            and node.children[1].type == token.RPAR
        )
        # last trailer can be arguments
        or (
            last
            and len(node.children) == 3
            and node.children[0].type == token.LPAR
            # and node.children[1].type == syms.argument
            and node.children[2].type == token.RPAR
        )
    )


def is_simple_decorator_expression(node: LN) -> bool:
    """Return True iff `node` could be a 'dotted name' decorator

    This function takes the node of the 'namedexpr_test' of the new decorator
    grammar and test if it would be valid under the old decorator grammar.

    The old grammar was: decorator: @ dotted_name [arguments] NEWLINE
    The new grammar is : decorator: @ namedexpr_test NEWLINE
    """
    if node.type == token.NAME:
        return True
    if node.type == syms.power:
        if node.children:
            return (
                node.children[0].type == token.NAME
                and all(map(is_simple_decorator_trailer, node.children[1:-1]))
                and (
                    len(node.children) < 2
                    or is_simple_decorator_trailer(node.children[-1], last=True)
                )
            )
    return False


def is_yield(node: LN) -> bool:
    """Return True if `node` holds a `yield` or `yield from` expression."""
    if node.type == syms.yield_expr:
        return True

    if node.type == token.NAME and node.value == "yield":  # type: ignore
        return True

    if node.type != syms.atom:
        return False

    if len(node.children) != 3:
        return False

    lpar, expr, rpar = node.children
    if lpar.type == token.LPAR and rpar.type == token.RPAR:
        return is_yield(expr)

    return False


def is_vararg(leaf: Leaf, within: Set[NodeType]) -> bool:
    """Return True if `leaf` is a star or double star in a vararg or kwarg.

    If `within` includes VARARGS_PARENTS, this applies to function signatures.
    If `within` includes UNPACKING_PARENTS, it applies to right hand-side
    extended iterable unpacking (PEP 3132) and additional unpacking
    generalizations (PEP 448).
    """
    if leaf.type not in VARARGS_SPECIALS or not leaf.parent:
        return False

    p = leaf.parent
    if p.type == syms.star_expr:
        # Star expressions are also used as assignment targets in extended
        # iterable unpacking (PEP 3132).  See what its parent is instead.
        if not p.parent:
            return False

        p = p.parent

    return p.type in within


def is_multiline_string(leaf: Leaf) -> bool:
    """Return True if `leaf` is a multiline string that actually spans many lines."""
    return has_triple_quotes(leaf.value) and "\n" in leaf.value


def is_stub_suite(node: Node) -> bool:
    """Return True if `node` is a suite with a stub body."""
    if (
        len(node.children) != 4
        or node.children[0].type != token.NEWLINE
        or node.children[1].type != token.INDENT
        or node.children[3].type != token.DEDENT
    ):
        return False

    return is_stub_body(node.children[2])


def is_stub_body(node: LN) -> bool:
    """Return True if `node` is a simple statement containing an ellipsis."""
    if not isinstance(node, Node) or node.type != syms.simple_stmt:
        return False

    if len(node.children) != 2:
        return False

    child = node.children[0]
    return (
        child.type == syms.atom
        and len(child.children) == 3
        and all(leaf == Leaf(token.DOT, ".") for leaf in child.children)
    )


def is_atom_with_invisible_parens(node: LN) -> bool:
    """Given a `LN`, determines whether it's an atom `node` with invisible
    parens. Useful in dedupe-ing and normalizing parens.
    """
    if isinstance(node, Leaf) or node.type != syms.atom:
        return False

    first, last = node.children[0], node.children[-1]
    return (
        isinstance(first, Leaf)
        and first.type == token.LPAR
        and first.value == ""
        and isinstance(last, Leaf)
        and last.type == token.RPAR
        and last.value == ""
    )


def is_empty_par(leaf: Leaf) -> bool:
    return is_empty_lpar(leaf) or is_empty_rpar(leaf)


def is_empty_lpar(leaf: Leaf) -> bool:
    return leaf.type == token.LPAR and leaf.value == ""


def is_empty_rpar(leaf: Leaf) -> bool:
    return leaf.type == token.RPAR and leaf.value == ""


def is_import(leaf: Leaf) -> bool:
    """Return True if the given leaf starts an import statement."""
    p = leaf.parent
    t = leaf.type
    v = leaf.value
    return bool(
        t == token.NAME
        and (
            (v == "import" and p and p.type == syms.import_name)
            or (v == "from" and p and p.type == syms.import_from)
        )
    )


def is_type_comment(leaf: Leaf, suffix: str = "") -> bool:
    """Return True if the given leaf is a special comment.
    Only returns true for type comments for now."""
    t = leaf.type
    v = leaf.value
    return t in {token.COMMENT, STANDALONE_COMMENT} and v.startswith("# type:" + suffix)


def wrap_in_parentheses(parent: Node, child: LN, *, visible: bool = True) -> None:
    """Wrap `child` in parentheses.

    This replaces `child` with an atom holding the parentheses and the old
    child.  That requires moving the prefix.

    If `visible` is False, the leaves will be valueless (and thus invisible).
    """
    lpar = Leaf(token.LPAR, "(" if visible else "")
    rpar = Leaf(token.RPAR, ")" if visible else "")
    prefix = child.prefix
    child.prefix = ""
    index = child.remove() or 0
    new_child = Node(syms.atom, [lpar, child, rpar])
    new_child.prefix = prefix
    parent.insert_child(index, new_child)


def unwrap_singleton_parenthesis(node: LN) -> Optional[LN]:
    """Returns `wrapped` if `node` is of the shape ( wrapped ).

    Parenthesis can be optional. Returns None otherwise"""
    if len(node.children) != 3:
        return None

    lpar, wrapped, rpar = node.children
    if not (lpar.type == token.LPAR and rpar.type == token.RPAR):
        return None

    return wrapped


def ensure_visible(leaf: Leaf) -> None:
    """Make sure parentheses are visible.

    They could be invisible as part of some statements (see
    :func:`normalize_invisible_parens` and :func:`visit_import_from`).
    """
    if leaf.type == token.LPAR:
        leaf.value = "("
    elif leaf.type == token.RPAR:
        leaf.value = ")"
