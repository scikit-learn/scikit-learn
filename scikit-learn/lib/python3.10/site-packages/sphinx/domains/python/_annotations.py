from __future__ import annotations

import ast
import functools
import operator
import token
from collections import deque
from inspect import Parameter
from typing import TYPE_CHECKING, Any

from docutils import nodes

from sphinx import addnodes
from sphinx.addnodes import desc_signature, pending_xref, pending_xref_condition
from sphinx.pycode.parser import Token, TokenProcessor
from sphinx.util.inspect import signature_from_str

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator

    from docutils.nodes import Element, Node

    from sphinx.environment import BuildEnvironment


def parse_reftarget(reftarget: str, suppress_prefix: bool = False,
                    ) -> tuple[str, str, str, bool]:
    """Parse a type string and return (reftype, reftarget, title, refspecific flag)"""
    refspecific = False
    if reftarget.startswith('.'):
        reftarget = reftarget[1:]
        title = reftarget
        refspecific = True
    elif reftarget.startswith('~'):
        reftarget = reftarget[1:]
        title = reftarget.split('.')[-1]
    elif suppress_prefix:
        title = reftarget.split('.')[-1]
    elif reftarget.startswith('typing.'):
        title = reftarget[7:]
    else:
        title = reftarget

    if reftarget == 'None' or reftarget.startswith('typing.'):
        # typing module provides non-class types.  Obj reference is good to refer them.
        reftype = 'obj'
    else:
        reftype = 'class'

    return reftype, reftarget, title, refspecific


def type_to_xref(target: str, env: BuildEnvironment, *,
                 suppress_prefix: bool = False) -> addnodes.pending_xref:
    """Convert a type string to a cross reference node."""
    if env:
        kwargs = {'py:module': env.ref_context.get('py:module'),
                  'py:class': env.ref_context.get('py:class')}
    else:
        kwargs = {}

    reftype, target, title, refspecific = parse_reftarget(target, suppress_prefix)

    if env.config.python_use_unqualified_type_names:
        # Note: It would be better to use qualname to describe the object to support support
        # nested classes.  But python domain can't access the real python object because this
        # module should work not-dynamically.
        shortname = title.split('.')[-1]
        contnodes: list[Node] = [pending_xref_condition('', shortname, condition='resolved'),
                                 pending_xref_condition('', title, condition='*')]
    else:
        contnodes = [nodes.Text(title)]

    return pending_xref('', *contnodes,
                        refdomain='py', reftype=reftype, reftarget=target,
                        refspecific=refspecific, **kwargs)


def _parse_annotation(annotation: str, env: BuildEnvironment) -> list[Node]:
    """Parse type annotation."""
    short_literals = env.config.python_display_short_literal_types

    def unparse(node: ast.AST) -> list[Node]:
        if isinstance(node, ast.Attribute):
            return [nodes.Text(f"{unparse(node.value)[0]}.{node.attr}")]
        if isinstance(node, ast.BinOp):
            result: list[Node] = unparse(node.left)
            result.extend(unparse(node.op))
            result.extend(unparse(node.right))
            return result
        if isinstance(node, ast.BitOr):
            return [addnodes.desc_sig_space(),
                    addnodes.desc_sig_punctuation('', '|'),
                    addnodes.desc_sig_space()]
        if isinstance(node, ast.Constant):
            if node.value is Ellipsis:
                return [addnodes.desc_sig_punctuation('', "...")]
            if isinstance(node.value, bool):
                return [addnodes.desc_sig_keyword('', repr(node.value))]
            if isinstance(node.value, int):
                return [addnodes.desc_sig_literal_number('', repr(node.value))]
            if isinstance(node.value, str):
                return [addnodes.desc_sig_literal_string('', repr(node.value))]
            else:
                # handles None, which is further handled by type_to_xref later
                # and fallback for other types that should be converted
                return [nodes.Text(repr(node.value))]
        if isinstance(node, ast.Expr):
            return unparse(node.value)
        if isinstance(node, ast.Invert):
            return [addnodes.desc_sig_punctuation('', '~')]
        if isinstance(node, ast.USub):
            return [addnodes.desc_sig_punctuation('', '-')]
        if isinstance(node, ast.List):
            result = [addnodes.desc_sig_punctuation('', '[')]
            if node.elts:
                # check if there are elements in node.elts to only pop the
                # last element of result if the for-loop was run at least
                # once
                for elem in node.elts:
                    result.extend(unparse(elem))
                    result.append(addnodes.desc_sig_punctuation('', ','))
                    result.append(addnodes.desc_sig_space())
                result.pop()
                result.pop()
            result.append(addnodes.desc_sig_punctuation('', ']'))
            return result
        if isinstance(node, ast.Module):
            return functools.reduce(operator.iadd, (unparse(e) for e in node.body), [])
        if isinstance(node, ast.Name):
            return [nodes.Text(node.id)]
        if isinstance(node, ast.Subscript):
            if getattr(node.value, 'id', '') in {'Optional', 'Union'}:
                return _unparse_pep_604_annotation(node)
            if short_literals and getattr(node.value, 'id', '') == 'Literal':
                return _unparse_pep_604_annotation(node)
            result = unparse(node.value)
            result.append(addnodes.desc_sig_punctuation('', '['))
            result.extend(unparse(node.slice))
            result.append(addnodes.desc_sig_punctuation('', ']'))

            # Wrap the Text nodes inside brackets by literal node if the subscript is a Literal
            if result[0] in ('Literal', 'typing.Literal'):
                for i, subnode in enumerate(result[1:], start=1):
                    if isinstance(subnode, nodes.Text):
                        result[i] = nodes.literal('', '', subnode)
            return result
        if isinstance(node, ast.UnaryOp):
            return unparse(node.op) + unparse(node.operand)
        if isinstance(node, ast.Tuple):
            if node.elts:
                result = []
                for elem in node.elts:
                    result.extend(unparse(elem))
                    result.append(addnodes.desc_sig_punctuation('', ','))
                    result.append(addnodes.desc_sig_space())
                result.pop()
                result.pop()
            else:
                result = [addnodes.desc_sig_punctuation('', '('),
                          addnodes.desc_sig_punctuation('', ')')]

            return result
        if isinstance(node, ast.Call):
            # Call nodes can be used in Annotated type metadata,
            # for example Annotated[str, ArbitraryTypeValidator(str, len=10)]
            args = []
            for arg in node.args:
                args += unparse(arg)
                args.append(addnodes.desc_sig_punctuation('', ','))
                args.append(addnodes.desc_sig_space())
            for kwd in node.keywords:
                args.append(addnodes.desc_sig_name(kwd.arg, kwd.arg))  # type: ignore[arg-type]
                args.append(addnodes.desc_sig_operator('', '='))
                args += unparse(kwd.value)
                args.append(addnodes.desc_sig_punctuation('', ','))
                args.append(addnodes.desc_sig_space())
            result = [
                *unparse(node.func),
                addnodes.desc_sig_punctuation('', '('),
                *args[:-2],  # skip the final comma and space
                addnodes.desc_sig_punctuation('', ')'),
            ]
            return result
        msg = f'unsupported syntax: {node}'
        raise SyntaxError(msg)  # unsupported syntax

    def _unparse_pep_604_annotation(node: ast.Subscript) -> list[Node]:
        subscript = node.slice

        flattened: list[Node] = []
        if isinstance(subscript, ast.Tuple):
            flattened.extend(unparse(subscript.elts[0]))
            for elt in subscript.elts[1:]:
                flattened.extend(unparse(ast.BitOr()))
                flattened.extend(unparse(elt))
        else:
            # e.g. a Union[] inside an Optional[]
            flattened.extend(unparse(subscript))

        if getattr(node.value, 'id', '') == 'Optional':
            flattened.extend(unparse(ast.BitOr()))
            flattened.append(nodes.Text('None'))

        return flattened

    try:
        tree = ast.parse(annotation, type_comments=True)
        result: list[Node] = []
        for node in unparse(tree):
            if isinstance(node, nodes.literal):
                result.append(node[0])
            elif isinstance(node, nodes.Text) and node.strip():
                if (result and isinstance(result[-1], addnodes.desc_sig_punctuation) and
                        result[-1].astext() == '~'):
                    result.pop()
                    result.append(type_to_xref(str(node), env, suppress_prefix=True))
                else:
                    result.append(type_to_xref(str(node), env))
            else:
                result.append(node)
        return result
    except SyntaxError:
        return [type_to_xref(annotation, env)]


class _TypeParameterListParser(TokenProcessor):
    def __init__(self, sig: str) -> None:
        signature = sig.replace('\n', '').strip()
        super().__init__([signature])
        # Each item is a tuple (name, kind, default, annotation) mimicking
        # ``inspect.Parameter`` to allow default values on VAR_POSITIONAL
        # or VAR_KEYWORD parameters.
        self.type_params: list[tuple[str, int, Any, Any]] = []

    def fetch_type_param_spec(self) -> list[Token]:
        tokens = []
        while current := self.fetch_token():
            tokens.append(current)
            for ldelim, rdelim in ('(', ')'), ('{', '}'), ('[', ']'):
                if current == [token.OP, ldelim]:
                    tokens += self.fetch_until([token.OP, rdelim])
                    break
            else:
                if current == token.INDENT:
                    tokens += self.fetch_until(token.DEDENT)
                elif current.match(
                        [token.OP, ':'], [token.OP, '='], [token.OP, ',']):
                    tokens.pop()
                    break
        return tokens

    def parse(self) -> None:
        while current := self.fetch_token():
            if current == token.NAME:
                tp_name = current.value.strip()
                if self.previous and self.previous.match([token.OP, '*'], [token.OP, '**']):
                    if self.previous == [token.OP, '*']:
                        tp_kind = Parameter.VAR_POSITIONAL
                    else:
                        tp_kind = Parameter.VAR_KEYWORD  # type: ignore[assignment]
                else:
                    tp_kind = Parameter.POSITIONAL_OR_KEYWORD  # type: ignore[assignment]

                tp_ann: Any = Parameter.empty
                tp_default: Any = Parameter.empty

                current = self.fetch_token()
                if current and current.match([token.OP, ':'], [token.OP, '=']):
                    if current == [token.OP, ':']:
                        tokens = self.fetch_type_param_spec()
                        tp_ann = self._build_identifier(tokens)

                    if self.current and self.current == [token.OP, '=']:
                        tokens = self.fetch_type_param_spec()
                        tp_default = self._build_identifier(tokens)

                if tp_kind != Parameter.POSITIONAL_OR_KEYWORD and tp_ann != Parameter.empty:
                    msg = ('type parameter bound or constraint is not allowed '
                           f'for {tp_kind.description} parameters')
                    raise SyntaxError(msg)

                type_param = (tp_name, tp_kind, tp_default, tp_ann)
                self.type_params.append(type_param)

    def _build_identifier(self, tokens: list[Token]) -> str:
        from itertools import chain, islice

        def triplewise(iterable: Iterable[Token]) -> Iterator[tuple[Token, ...]]:
            # sliding_window('ABCDEFG', 4) --> ABCD BCDE CDEF DEFG
            it = iter(iterable)
            window = deque(islice(it, 3), maxlen=3)
            if len(window) == 3:
                yield tuple(window)
            for x in it:
                window.append(x)
                yield tuple(window)

        idents: list[str] = []
        tokens: Iterable[Token] = iter(tokens)  # type: ignore[no-redef]
        # do not format opening brackets
        for tok in tokens:
            if not tok.match([token.OP, '('], [token.OP, '['], [token.OP, '{']):
                # check if the first non-delimiter character is an unpack operator
                is_unpack_operator = tok.match([token.OP, '*'], [token.OP, ['**']])
                idents.append(self._pformat_token(tok, native=is_unpack_operator))
                break
            idents.append(tok.value)

        # check the remaining tokens
        stop = Token(token.ENDMARKER, '', (-1, -1), (-1, -1), '<sentinel>')
        is_unpack_operator = False
        for tok, op, after in triplewise(chain(tokens, [stop, stop])):
            ident = self._pformat_token(tok, native=is_unpack_operator)
            idents.append(ident)
            # determine if the next token is an unpack operator depending
            # on the left and right hand side of the operator symbol
            is_unpack_operator = (
                op.match([token.OP, '*'], [token.OP, '**']) and not (
                    tok.match(token.NAME, token.NUMBER, token.STRING,
                              [token.OP, ')'], [token.OP, ']'], [token.OP, '}'])
                    and after.match(token.NAME, token.NUMBER, token.STRING,
                                    [token.OP, '('], [token.OP, '['], [token.OP, '{'])
                )
            )

        return ''.join(idents).strip()

    def _pformat_token(self, tok: Token, native: bool = False) -> str:
        if native:
            return tok.value

        if tok.match(token.NEWLINE, token.ENDMARKER):
            return ''

        if tok.match([token.OP, ':'], [token.OP, ','], [token.OP, '#']):
            return f'{tok.value} '

        # Arithmetic operators are allowed because PEP 695 specifies the
        # default type parameter to be *any* expression (so "T1 << T2" is
        # allowed if it makes sense). The caller is responsible to ensure
        # that a multiplication operator ("*") is not to be confused with
        # an unpack operator (which will not be surrounded by spaces).
        #
        # The operators are ordered according to how likely they are to
        # be used and for (possible) future implementations (e.g., "&" for
        # an intersection type).
        if tok.match(
            # Most likely operators to appear
            [token.OP, '='], [token.OP, '|'],
            # Type composition (future compatibility)
            [token.OP, '&'], [token.OP, '^'], [token.OP, '<'], [token.OP, '>'],
            # Unlikely type composition
            [token.OP, '+'], [token.OP, '-'], [token.OP, '*'], [token.OP, '**'],
            # Unlikely operators but included for completeness
            [token.OP, '@'], [token.OP, '/'], [token.OP, '//'], [token.OP, '%'],
            [token.OP, '<<'], [token.OP, '>>'], [token.OP, '>>>'],
            [token.OP, '<='], [token.OP, '>='], [token.OP, '=='], [token.OP, '!='],
        ):
            return f' {tok.value} '

        return tok.value


def _parse_type_list(
    tp_list: str, env: BuildEnvironment,
    multi_line_parameter_list: bool = False,
) -> addnodes.desc_type_parameter_list:
    """Parse a list of type parameters according to PEP 695."""
    type_params = addnodes.desc_type_parameter_list(tp_list)
    type_params['multi_line_parameter_list'] = multi_line_parameter_list
    # formal parameter names are interpreted as type parameter names and
    # type annotations are interpreted as type parameter bound or constraints
    parser = _TypeParameterListParser(tp_list)
    parser.parse()
    for (tp_name, tp_kind, tp_default, tp_ann) in parser.type_params:
        # no positional-only or keyword-only allowed in a type parameters list
        if tp_kind in {Parameter.POSITIONAL_ONLY, Parameter.KEYWORD_ONLY}:
            msg = ('positional-only or keyword-only parameters '
                   'are prohibited in type parameter lists')
            raise SyntaxError(msg)

        node = addnodes.desc_type_parameter()
        if tp_kind == Parameter.VAR_POSITIONAL:
            node += addnodes.desc_sig_operator('', '*')
        elif tp_kind == Parameter.VAR_KEYWORD:
            node += addnodes.desc_sig_operator('', '**')
        node += addnodes.desc_sig_name('', tp_name)

        if tp_ann is not Parameter.empty:
            annotation = _parse_annotation(tp_ann, env)
            if not annotation:
                continue

            node += addnodes.desc_sig_punctuation('', ':')
            node += addnodes.desc_sig_space()

            type_ann_expr = addnodes.desc_sig_name('', '',
                                                   *annotation)  # type: ignore[arg-type]
            # a type bound is ``T: U`` whereas type constraints
            # must be enclosed with parentheses. ``T: (U, V)``
            if tp_ann.startswith('(') and tp_ann.endswith(')'):
                type_ann_text = type_ann_expr.astext()
                if type_ann_text.startswith('(') and type_ann_text.endswith(')'):
                    node += type_ann_expr
                else:
                    # surrounding braces are lost when using _parse_annotation()
                    node += addnodes.desc_sig_punctuation('', '(')
                    node += type_ann_expr  # type constraint
                    node += addnodes.desc_sig_punctuation('', ')')
            else:
                node += type_ann_expr  # type bound

        if tp_default is not Parameter.empty:
            # Always surround '=' with spaces, even if there is no annotation
            node += addnodes.desc_sig_space()
            node += addnodes.desc_sig_operator('', '=')
            node += addnodes.desc_sig_space()
            node += nodes.inline('', tp_default,
                                 classes=['default_value'],
                                 support_smartquotes=False)

        type_params += node
    return type_params


def _parse_arglist(
    arglist: str, env: BuildEnvironment, multi_line_parameter_list: bool = False,
) -> addnodes.desc_parameterlist:
    """Parse a list of arguments using AST parser"""
    params = addnodes.desc_parameterlist(arglist)
    params['multi_line_parameter_list'] = multi_line_parameter_list
    sig = signature_from_str('(%s)' % arglist)
    last_kind = None
    for param in sig.parameters.values():
        if param.kind != param.POSITIONAL_ONLY and last_kind == param.POSITIONAL_ONLY:
            # PEP-570: Separator for Positional Only Parameter: /
            params += addnodes.desc_parameter('', '', addnodes.desc_sig_operator('', '/'))
        if param.kind == param.KEYWORD_ONLY and last_kind in (param.POSITIONAL_OR_KEYWORD,
                                                              param.POSITIONAL_ONLY,
                                                              None):
            # PEP-3102: Separator for Keyword Only Parameter: *
            params += addnodes.desc_parameter('', '', addnodes.desc_sig_operator('', '*'))

        node = addnodes.desc_parameter()
        if param.kind == param.VAR_POSITIONAL:
            node += addnodes.desc_sig_operator('', '*')
            node += addnodes.desc_sig_name('', param.name)
        elif param.kind == param.VAR_KEYWORD:
            node += addnodes.desc_sig_operator('', '**')
            node += addnodes.desc_sig_name('', param.name)
        else:
            node += addnodes.desc_sig_name('', param.name)

        if param.annotation is not param.empty:
            children = _parse_annotation(param.annotation, env)
            node += addnodes.desc_sig_punctuation('', ':')
            node += addnodes.desc_sig_space()
            node += addnodes.desc_sig_name('', '', *children)  # type: ignore[arg-type]
        if param.default is not param.empty:
            if param.annotation is not param.empty:
                node += addnodes.desc_sig_space()
                node += addnodes.desc_sig_operator('', '=')
                node += addnodes.desc_sig_space()
            else:
                node += addnodes.desc_sig_operator('', '=')
            node += nodes.inline('', param.default, classes=['default_value'],
                                 support_smartquotes=False)

        params += node
        last_kind = param.kind

    if last_kind == Parameter.POSITIONAL_ONLY:
        # PEP-570: Separator for Positional Only Parameter: /
        params += addnodes.desc_parameter('', '', addnodes.desc_sig_operator('', '/'))

    return params


def _pseudo_parse_arglist(
    signode: desc_signature, arglist: str, multi_line_parameter_list: bool = False,
) -> None:
    """"Parse" a list of arguments separated by commas.

    Arguments can have "optional" annotations given by enclosing them in
    brackets.  Currently, this will split at any comma, even if it's inside a
    string literal (e.g. default argument value).
    """
    paramlist = addnodes.desc_parameterlist()
    paramlist['multi_line_parameter_list'] = multi_line_parameter_list
    stack: list[Element] = [paramlist]
    try:
        for argument in arglist.split(','):
            argument = argument.strip()
            ends_open = ends_close = 0
            while argument.startswith('['):
                stack.append(addnodes.desc_optional())
                stack[-2] += stack[-1]
                argument = argument[1:].strip()
            while argument.startswith(']'):
                stack.pop()
                argument = argument[1:].strip()
            while argument.endswith(']') and not argument.endswith('[]'):
                ends_close += 1
                argument = argument[:-1].strip()
            while argument.endswith('['):
                ends_open += 1
                argument = argument[:-1].strip()
            if argument:
                stack[-1] += addnodes.desc_parameter(
                    '', '', addnodes.desc_sig_name(argument, argument))
            while ends_open:
                stack.append(addnodes.desc_optional())
                stack[-2] += stack[-1]
                ends_open -= 1
            while ends_close:
                stack.pop()
                ends_close -= 1
        if len(stack) != 1:
            raise IndexError
    except IndexError:
        # if there are too few or too many elements on the stack, just give up
        # and treat the whole argument list as one argument, discarding the
        # already partially populated paramlist node
        paramlist = addnodes.desc_parameterlist()
        paramlist += addnodes.desc_parameter(arglist, arglist)
        signode += paramlist
    else:
        signode += paramlist
