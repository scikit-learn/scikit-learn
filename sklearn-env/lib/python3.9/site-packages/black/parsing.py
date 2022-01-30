"""
Parse Python code and perform AST validation.
"""
import ast
import sys
from typing import Iterable, Iterator, List, Set, Union

# lib2to3 fork
from blib2to3.pytree import Node, Leaf
from blib2to3 import pygram, pytree
from blib2to3.pgen2 import driver
from blib2to3.pgen2.grammar import Grammar
from blib2to3.pgen2.parse import ParseError

from black.mode import TargetVersion, Feature, supports_feature
from black.nodes import syms

try:
    from typed_ast import ast3, ast27
except ImportError:
    if sys.version_info < (3, 8):
        print(
            "The typed_ast package is required but not installed.\n"
            "You can upgrade to Python 3.8+ or install typed_ast with\n"
            "`python3 -m pip install typed-ast`.",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        ast3 = ast27 = ast


class InvalidInput(ValueError):
    """Raised when input source code fails all parse attempts."""


def get_grammars(target_versions: Set[TargetVersion]) -> List[Grammar]:
    if not target_versions:
        # No target_version specified, so try all grammars.
        return [
            # Python 3.7+
            pygram.python_grammar_no_print_statement_no_exec_statement_async_keywords,
            # Python 3.0-3.6
            pygram.python_grammar_no_print_statement_no_exec_statement,
            # Python 2.7 with future print_function import
            pygram.python_grammar_no_print_statement,
            # Python 2.7
            pygram.python_grammar,
        ]

    if all(version.is_python2() for version in target_versions):
        # Python 2-only code, so try Python 2 grammars.
        return [
            # Python 2.7 with future print_function import
            pygram.python_grammar_no_print_statement,
            # Python 2.7
            pygram.python_grammar,
        ]

    # Python 3-compatible code, so only try Python 3 grammar.
    grammars = []
    # If we have to parse both, try to parse async as a keyword first
    if not supports_feature(target_versions, Feature.ASYNC_IDENTIFIERS):
        # Python 3.7+
        grammars.append(
            pygram.python_grammar_no_print_statement_no_exec_statement_async_keywords
        )
    if not supports_feature(target_versions, Feature.ASYNC_KEYWORDS):
        # Python 3.0-3.6
        grammars.append(pygram.python_grammar_no_print_statement_no_exec_statement)
    # At least one of the above branches must have been taken, because every Python
    # version has exactly one of the two 'ASYNC_*' flags
    return grammars


def lib2to3_parse(src_txt: str, target_versions: Iterable[TargetVersion] = ()) -> Node:
    """Given a string with source, return the lib2to3 Node."""
    if not src_txt.endswith("\n"):
        src_txt += "\n"

    for grammar in get_grammars(set(target_versions)):
        drv = driver.Driver(grammar, pytree.convert)
        try:
            result = drv.parse_string(src_txt, True)
            break

        except ParseError as pe:
            lineno, column = pe.context[1]
            lines = src_txt.splitlines()
            try:
                faulty_line = lines[lineno - 1]
            except IndexError:
                faulty_line = "<line number missing in source>"
            exc = InvalidInput(f"Cannot parse: {lineno}:{column}: {faulty_line}")
    else:
        raise exc from None

    if isinstance(result, Leaf):
        result = Node(syms.file_input, [result])
    return result


def lib2to3_unparse(node: Node) -> str:
    """Given a lib2to3 node, return its string representation."""
    code = str(node)
    return code


def parse_ast(src: str) -> Union[ast.AST, ast3.AST, ast27.AST]:
    filename = "<unknown>"
    if sys.version_info >= (3, 8):
        # TODO: support Python 4+ ;)
        for minor_version in range(sys.version_info[1], 4, -1):
            try:
                return ast.parse(src, filename, feature_version=(3, minor_version))
            except SyntaxError:
                continue
    else:
        for feature_version in (7, 6):
            try:
                return ast3.parse(src, filename, feature_version=feature_version)
            except SyntaxError:
                continue
    if ast27.__name__ == "ast":
        raise SyntaxError(
            "The requested source code has invalid Python 3 syntax.\n"
            "If you are trying to format Python 2 files please reinstall Black"
            " with the 'python2' extra: `python3 -m pip install black[python2]`."
        )
    return ast27.parse(src)


def stringify_ast(
    node: Union[ast.AST, ast3.AST, ast27.AST], depth: int = 0
) -> Iterator[str]:
    """Simple visitor generating strings to compare ASTs by content."""

    node = fixup_ast_constants(node)

    yield f"{'  ' * depth}{node.__class__.__name__}("

    for field in sorted(node._fields):  # noqa: F402
        # TypeIgnore has only one field 'lineno' which breaks this comparison
        type_ignore_classes = (ast3.TypeIgnore, ast27.TypeIgnore)
        if sys.version_info >= (3, 8):
            type_ignore_classes += (ast.TypeIgnore,)
        if isinstance(node, type_ignore_classes):
            break

        try:
            value = getattr(node, field)
        except AttributeError:
            continue

        yield f"{'  ' * (depth+1)}{field}="

        if isinstance(value, list):
            for item in value:
                # Ignore nested tuples within del statements, because we may insert
                # parentheses and they change the AST.
                if (
                    field == "targets"
                    and isinstance(node, (ast.Delete, ast3.Delete, ast27.Delete))
                    and isinstance(item, (ast.Tuple, ast3.Tuple, ast27.Tuple))
                ):
                    for item in item.elts:
                        yield from stringify_ast(item, depth + 2)

                elif isinstance(item, (ast.AST, ast3.AST, ast27.AST)):
                    yield from stringify_ast(item, depth + 2)

        elif isinstance(value, (ast.AST, ast3.AST, ast27.AST)):
            yield from stringify_ast(value, depth + 2)

        else:
            # Constant strings may be indented across newlines, if they are
            # docstrings; fold spaces after newlines when comparing. Similarly,
            # trailing and leading space may be removed.
            # Note that when formatting Python 2 code, at least with Windows
            # line-endings, docstrings can end up here as bytes instead of
            # str so make sure that we handle both cases.
            if (
                isinstance(node, ast.Constant)
                and field == "value"
                and isinstance(value, (str, bytes))
            ):
                lineend = "\n" if isinstance(value, str) else b"\n"
                # To normalize, we strip any leading and trailing space from
                # each line...
                stripped = [line.strip() for line in value.splitlines()]
                normalized = lineend.join(stripped)  # type: ignore[attr-defined]
                # ...and remove any blank lines at the beginning and end of
                # the whole string
                normalized = normalized.strip()
            else:
                normalized = value
            yield f"{'  ' * (depth+2)}{normalized!r},  # {value.__class__.__name__}"

    yield f"{'  ' * depth})  # /{node.__class__.__name__}"


def fixup_ast_constants(
    node: Union[ast.AST, ast3.AST, ast27.AST]
) -> Union[ast.AST, ast3.AST, ast27.AST]:
    """Map ast nodes deprecated in 3.8 to Constant."""
    if isinstance(node, (ast.Str, ast3.Str, ast27.Str, ast.Bytes, ast3.Bytes)):
        return ast.Constant(value=node.s)

    if isinstance(node, (ast.Num, ast3.Num, ast27.Num)):
        return ast.Constant(value=node.n)

    if isinstance(node, (ast.NameConstant, ast3.NameConstant)):
        return ast.Constant(value=node.value)

    return node
