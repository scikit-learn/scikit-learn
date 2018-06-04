"""
PEP 0484 ( https://www.python.org/dev/peps/pep-0484/ ) describes type hints
through function annotations. There is a strong suggestion in this document
that only the type of type hinting defined in PEP0484 should be allowed
as annotations in future python versions.

The (initial / probably incomplete) implementation todo list for pep-0484:
v Function parameter annotations with builtin/custom type classes
v Function returntype annotations with builtin/custom type classes
v Function parameter annotations with strings (forward reference)
v Function return type annotations with strings (forward reference)
v Local variable type hints
v Assigned types: `Url = str\ndef get(url:Url) -> str:`
v Type hints in `with` statements
x Stub files support
x support `@no_type_check` and `@no_type_check_decorator`
x support for typing.cast() operator
x support for type hint comments for functions, `# type: (int, str) -> int`.
    See comment from Guido https://github.com/davidhalter/jedi/issues/662
"""

import os
import re

from parso import ParserSyntaxError, parse, split_lines
from parso.python import tree

from jedi._compatibility import unicode, force_unicode
from jedi.evaluate.cache import evaluator_method_cache
from jedi.evaluate import compiled
from jedi.evaluate.base_context import NO_CONTEXTS, ContextSet
from jedi.evaluate.lazy_context import LazyTreeContext
from jedi.evaluate.context import ModuleContext
from jedi.evaluate.helpers import is_string
from jedi import debug
from jedi import parser_utils


def _evaluate_for_annotation(context, annotation, index=None):
    """
    Evaluates a string-node, looking for an annotation
    If index is not None, the annotation is expected to be a tuple
    and we're interested in that index
    """
    context_set = context.eval_node(_fix_forward_reference(context, annotation))
    return context_set.execute_evaluated()


def _evaluate_annotation_string(context, string, index=None):
    node = _get_forward_reference_node(context, string)
    if node is None:
        return NO_CONTEXTS

    context_set = context.eval_node(node)
    if index is not None:
        context_set = context_set.filter(
            lambda context: context.array_type == u'tuple'
                            and len(list(context.py__iter__())) >= index
        ).py__getitem__(index)
    return context_set.execute_evaluated()


def _fix_forward_reference(context, node):
    evaled_nodes = context.eval_node(node)
    if len(evaled_nodes) != 1:
        debug.warning("Eval'ed typing index %s should lead to 1 object, "
                      " not %s" % (node, evaled_nodes))
        return node

    evaled_context = list(evaled_nodes)[0]
    if is_string(evaled_context):
        result = _get_forward_reference_node(context, evaled_context.get_safe_value())
        if result is not None:
            return result

    return node


def _get_forward_reference_node(context, string):
    try:
        new_node = context.evaluator.grammar.parse(
            force_unicode(string),
            start_symbol='eval_input',
            error_recovery=False
        )
    except ParserSyntaxError:
        debug.warning('Annotation not parsed: %s' % string)
        return None
    else:
        module = context.tree_node.get_root_node()
        parser_utils.move(new_node, module.end_pos[0])
        new_node.parent = context.tree_node
        return new_node


def _split_comment_param_declaration(decl_text):
    """
    Split decl_text on commas, but group generic expressions
    together.

    For example, given "foo, Bar[baz, biz]" we return
    ['foo', 'Bar[baz, biz]'].

    """
    try:
        node = parse(decl_text, error_recovery=False).children[0]
    except ParserSyntaxError:
        debug.warning('Comment annotation is not valid Python: %s' % decl_text)
        return []

    if node.type == 'name':
        return [node.get_code().strip()]

    params = []
    try:
        children = node.children
    except AttributeError:
        return []
    else:
        for child in children:
            if child.type in ['name', 'atom_expr', 'power']:
                params.append(child.get_code().strip())

    return params


@evaluator_method_cache()
def infer_param(execution_context, param):
    """
    Infers the type of a function parameter, using type annotations.
    """
    annotation = param.annotation
    if annotation is None:
        # If no Python 3-style annotation, look for a Python 2-style comment
        # annotation.
        # Identify parameters to function in the same sequence as they would
        # appear in a type comment.
        all_params = [child for child in param.parent.children
                      if child.type == 'param']

        node = param.parent.parent
        comment = parser_utils.get_following_comment_same_line(node)
        if comment is None:
            return NO_CONTEXTS

        match = re.match(r"^#\s*type:\s*\(([^#]*)\)\s*->", comment)
        if not match:
            return NO_CONTEXTS
        params_comments = _split_comment_param_declaration(match.group(1))

        # Find the specific param being investigated
        index = all_params.index(param)
        # If the number of parameters doesn't match length of type comment,
        # ignore first parameter (assume it's self).
        if len(params_comments) != len(all_params):
            debug.warning(
                "Comments length != Params length %s %s",
                params_comments, all_params
            )
        from jedi.evaluate.context.instance import BaseInstanceFunctionExecution
        if isinstance(execution_context, BaseInstanceFunctionExecution):
            if index == 0:
                # Assume it's self, which is already handled
                return NO_CONTEXTS
            index -= 1
        if index >= len(params_comments):
            return NO_CONTEXTS

        param_comment = params_comments[index]
        return _evaluate_annotation_string(
            execution_context.get_root_context(),
            param_comment
        )
    module_context = execution_context.get_root_context()
    return _evaluate_for_annotation(module_context, annotation)


def py__annotations__(funcdef):
    return_annotation = funcdef.annotation
    if return_annotation:
        dct = {'return': return_annotation}
    else:
        dct = {}
    for function_param in funcdef.get_params():
        param_annotation = function_param.annotation
        if param_annotation is not None:
            dct[function_param.name.value] = param_annotation
    return dct


@evaluator_method_cache()
def infer_return_types(function_context):
    """
    Infers the type of a function's return value,
    according to type annotations.
    """
    annotation = py__annotations__(function_context.tree_node).get("return", None)
    if annotation is None:
        # If there is no Python 3-type annotation, look for a Python 2-type annotation
        node = function_context.tree_node
        comment = parser_utils.get_following_comment_same_line(node)
        if comment is None:
            return NO_CONTEXTS

        match = re.match(r"^#\s*type:\s*\([^#]*\)\s*->\s*([^#]*)", comment)
        if not match:
            return NO_CONTEXTS

        return _evaluate_annotation_string(
            function_context.get_root_context(),
            match.group(1).strip()
        )

    module_context = function_context.get_root_context()
    return _evaluate_for_annotation(module_context, annotation)


_typing_module = None
_typing_module_code_lines = None


def _get_typing_replacement_module(grammar):
    """
    The idea is to return our jedi replacement for the PEP-0484 typing module
    as discussed at https://github.com/davidhalter/jedi/issues/663
    """
    global _typing_module, _typing_module_code_lines
    if _typing_module is None:
        typing_path = \
            os.path.abspath(os.path.join(__file__, "../jedi_typing.py"))
        with open(typing_path) as f:
            code = unicode(f.read())
        _typing_module = grammar.parse(code)
        _typing_module_code_lines = split_lines(code, keepends=True)
    return _typing_module, _typing_module_code_lines


def py__getitem__(context, typ, node):
    if not typ.get_root_context().name.string_name == "typing":
        return None
    # we assume that any class using [] in a module called
    # "typing" with a name for which we have a replacement
    # should be replaced by that class. This is not 100%
    # airtight but I don't have a better idea to check that it's
    # actually the PEP-0484 typing module and not some other
    if node.type == "subscriptlist":
        nodes = node.children[::2]  # skip the commas
    else:
        nodes = [node]
    del node

    nodes = [_fix_forward_reference(context, node) for node in nodes]
    type_name = typ.name.string_name

    # hacked in Union and Optional, since it's hard to do nicely in parsed code
    if type_name in ("Union", '_Union'):
        # In Python 3.6 it's still called typing.Union but it's an instance
        # called _Union.
        return ContextSet.from_sets(context.eval_node(node) for node in nodes)
    if type_name in ("Optional", '_Optional'):
        # Here we have the same issue like in Union. Therefore we also need to
        # check for the instance typing._Optional (Python 3.6).
        return context.eval_node(nodes[0])

    module_node, code_lines = _get_typing_replacement_module(context.evaluator.latest_grammar)
    typing = ModuleContext(
        context.evaluator,
        module_node=module_node,
        path=None,
        code_lines=code_lines,
    )
    factories = typing.py__getattribute__("factory")
    assert len(factories) == 1
    factory = list(factories)[0]
    assert factory
    function_body_nodes = factory.tree_node.children[4].children
    valid_classnames = set(child.name.value
                           for child in function_body_nodes
                           if isinstance(child, tree.Class))
    if type_name not in valid_classnames:
        return None
    compiled_classname = compiled.create_simple_object(context.evaluator, type_name)

    from jedi.evaluate.context.iterable import FakeSequence
    args = FakeSequence(
        context.evaluator,
        u'tuple',
        [LazyTreeContext(context, n) for n in nodes]
    )

    result = factory.execute_evaluated(compiled_classname, args)
    return result


def find_type_from_comment_hint_for(context, node, name):
    return _find_type_from_comment_hint(context, node, node.children[1], name)


def find_type_from_comment_hint_with(context, node, name):
    assert len(node.children[1].children) == 3, \
        "Can only be here when children[1] is 'foo() as f'"
    varlist = node.children[1].children[2]
    return _find_type_from_comment_hint(context, node, varlist, name)


def find_type_from_comment_hint_assign(context, node, name):
    return _find_type_from_comment_hint(context, node, node.children[0], name)


def _find_type_from_comment_hint(context, node, varlist, name):
    index = None
    if varlist.type in ("testlist_star_expr", "exprlist", "testlist"):
        # something like "a, b = 1, 2"
        index = 0
        for child in varlist.children:
            if child == name:
                break
            if child.type == "operator":
                continue
            index += 1
        else:
            return []

    comment = parser_utils.get_following_comment_same_line(node)
    if comment is None:
        return []
    match = re.match(r"^#\s*type:\s*([^#]*)", comment)
    if match is None:
        return []
    return _evaluate_annotation_string(context, match.group(1).strip(), index)
