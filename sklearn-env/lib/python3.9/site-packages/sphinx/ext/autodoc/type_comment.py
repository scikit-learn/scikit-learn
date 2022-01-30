"""
    sphinx.ext.autodoc.type_comment
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Update annotations info of living objects using type_comments.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

from inspect import Parameter, Signature, getsource
from typing import Any, Dict, List, cast

import sphinx
from sphinx.application import Sphinx
from sphinx.locale import __
from sphinx.pycode.ast import ast
from sphinx.pycode.ast import parse as ast_parse
from sphinx.pycode.ast import unparse as ast_unparse
from sphinx.util import inspect, logging

logger = logging.getLogger(__name__)


def not_suppressed(argtypes: List[ast.AST] = []) -> bool:
    """Check given *argtypes* is suppressed type_comment or not."""
    if len(argtypes) == 0:  # no argtypees
        return False
    elif len(argtypes) == 1 and ast_unparse(argtypes[0]) == "...":  # suppressed
        # Note: To support multiple versions of python, this uses ``ast_unparse()`` for
        # comparison with Ellipsis.  Since 3.8, ast.Constant has been used to represent
        # Ellipsis node instead of ast.Ellipsis.
        return False
    else:  # not suppressed
        return True


def signature_from_ast(node: ast.FunctionDef, bound_method: bool,
                       type_comment: ast.FunctionDef) -> Signature:
    """Return a Signature object for the given *node*.

    :param bound_method: Specify *node* is a bound method or not
    """
    params = []
    if hasattr(node.args, "posonlyargs"):  # for py38+
        for arg in node.args.posonlyargs:  # type: ignore
            param = Parameter(arg.arg, Parameter.POSITIONAL_ONLY, annotation=arg.type_comment)
            params.append(param)

    for arg in node.args.args:
        param = Parameter(arg.arg, Parameter.POSITIONAL_OR_KEYWORD,
                          annotation=arg.type_comment or Parameter.empty)
        params.append(param)

    if node.args.vararg:
        param = Parameter(node.args.vararg.arg, Parameter.VAR_POSITIONAL,
                          annotation=node.args.vararg.type_comment or Parameter.empty)
        params.append(param)

    for arg in node.args.kwonlyargs:
        param = Parameter(arg.arg, Parameter.KEYWORD_ONLY,
                          annotation=arg.type_comment or Parameter.empty)
        params.append(param)

    if node.args.kwarg:
        param = Parameter(node.args.kwarg.arg, Parameter.VAR_KEYWORD,
                          annotation=node.args.kwarg.type_comment or Parameter.empty)
        params.append(param)

    # Remove first parameter when *obj* is bound_method
    if bound_method and params:
        params.pop(0)

    # merge type_comment into signature
    if not_suppressed(type_comment.argtypes):  # type: ignore
        for i, param in enumerate(params):
            params[i] = param.replace(annotation=type_comment.argtypes[i])  # type: ignore

    if node.returns:
        return Signature(params, return_annotation=node.returns)
    elif type_comment.returns:
        return Signature(params, return_annotation=ast_unparse(type_comment.returns))
    else:
        return Signature(params)


def get_type_comment(obj: Any, bound_method: bool = False) -> Signature:
    """Get type_comment'ed FunctionDef object from living object.

    This tries to parse original code for living object and returns
    Signature for given *obj*.  It requires py38+ or typed_ast module.
    """
    try:
        source = getsource(obj)
        if source.startswith((' ', r'\t')):
            # subject is placed inside class or block.  To read its docstring,
            # this adds if-block before the declaration.
            module = ast_parse('if True:\n' + source)
            subject = cast(ast.FunctionDef, module.body[0].body[0])  # type: ignore
        else:
            module = ast_parse(source)
            subject = cast(ast.FunctionDef, module.body[0])  # type: ignore

        if getattr(subject, "type_comment", None):
            function = ast_parse(subject.type_comment, mode='func_type')
            return signature_from_ast(subject, bound_method, function)  # type: ignore
        else:
            return None
    except (OSError, TypeError):  # failed to load source code
        return None
    except SyntaxError:  # failed to parse type_comments
        return None


def update_annotations_using_type_comments(app: Sphinx, obj: Any, bound_method: bool) -> None:
    """Update annotations info of *obj* using type_comments."""
    try:
        type_sig = get_type_comment(obj, bound_method)
        if type_sig:
            sig = inspect.signature(obj, bound_method)
            for param in sig.parameters.values():
                if param.name not in obj.__annotations__:
                    annotation = type_sig.parameters[param.name].annotation
                    if annotation is not Parameter.empty:
                        obj.__annotations__[param.name] = ast_unparse(annotation)

            if 'return' not in obj.__annotations__:
                obj.__annotations__['return'] = type_sig.return_annotation
    except KeyError as exc:
        logger.warning(__("Failed to update signature for %r: parameter not found: %s"),
                       obj, exc)
    except NotImplementedError as exc:  # failed to ast.unparse()
        logger.warning(__("Failed to parse type_comment for %r: %s"), obj, exc)


def setup(app: Sphinx) -> Dict[str, Any]:
    app.connect('autodoc-before-process-signature', update_annotations_using_type_comments)

    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}
