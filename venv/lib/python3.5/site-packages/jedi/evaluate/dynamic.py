"""
One of the really important features of |jedi| is to have an option to
understand code like this::

    def foo(bar):
        bar. # completion here
    foo(1)

There's no doubt wheter bar is an ``int`` or not, but if there's also a call
like ``foo('str')``, what would happen? Well, we'll just show both. Because
that's what a human would expect.

It works as follows:

- |Jedi| sees a param
- search for function calls named ``foo``
- execute these calls and check the input.
"""

from parso.python import tree
from jedi import settings
from jedi import debug
from jedi.evaluate.cache import evaluator_function_cache
from jedi.evaluate import imports
from jedi.evaluate.arguments import TreeArguments
from jedi.evaluate.param import create_default_params
from jedi.evaluate.helpers import is_stdlib_path
from jedi.evaluate.utils import to_list
from jedi.parser_utils import get_parent_scope
from jedi.evaluate.context import ModuleContext, instance
from jedi.evaluate.base_context import ContextSet



MAX_PARAM_SEARCHES = 20


class MergedExecutedParams(object):
    """
    Simulates being a parameter while actually just being multiple params.
    """
    def __init__(self, executed_params):
        self._executed_params = executed_params

    def infer(self):
        return ContextSet.from_sets(p.infer() for p in self._executed_params)


@debug.increase_indent
def search_params(evaluator, execution_context, funcdef):
    """
    A dynamic search for param values. If you try to complete a type:

    >>> def func(foo):
    ...     foo
    >>> func(1)
    >>> func("")

    It is not known what the type ``foo`` without analysing the whole code. You
    have to look for all calls to ``func`` to find out what ``foo`` possibly
    is.
    """
    if not settings.dynamic_params:
        return create_default_params(execution_context, funcdef)

    evaluator.dynamic_params_depth += 1
    try:
        path = execution_context.get_root_context().py__file__()
        if path is not None and is_stdlib_path(path):
            # We don't want to search for usages in the stdlib. Usually people
            # don't work with it (except if you are a core maintainer, sorry).
            # This makes everything slower. Just disable it and run the tests,
            # you will see the slowdown, especially in 3.6.
            return create_default_params(execution_context, funcdef)

        if funcdef.type == 'lambdef':
            string_name = _get_lambda_name(funcdef)
            if string_name is None:
                return create_default_params(execution_context, funcdef)
        else:
            string_name = funcdef.name.value
        debug.dbg('Dynamic param search in %s.', string_name, color='MAGENTA')

        try:
            module_context = execution_context.get_root_context()
            function_executions = _search_function_executions(
                evaluator,
                module_context,
                funcdef,
                string_name=string_name,
            )
            if function_executions:
                zipped_params = zip(*list(
                    function_execution.get_params()
                    for function_execution in function_executions
                ))
                params = [MergedExecutedParams(executed_params) for executed_params in zipped_params]
                # Evaluate the ExecutedParams to types.
            else:
                return create_default_params(execution_context, funcdef)
        finally:
            debug.dbg('Dynamic param result finished', color='MAGENTA')
        return params
    finally:
        evaluator.dynamic_params_depth -= 1


@evaluator_function_cache(default=None)
@to_list
def _search_function_executions(evaluator, module_context, funcdef, string_name):
    """
    Returns a list of param names.
    """
    compare_node = funcdef
    if string_name == '__init__':
        cls = get_parent_scope(funcdef)
        if isinstance(cls, tree.Class):
            string_name = cls.name.value
            compare_node = cls

    found_executions = False
    i = 0
    for for_mod_context in imports.get_modules_containing_name(
            evaluator, [module_context], string_name):
        if not isinstance(module_context, ModuleContext):
            return
        for name, trailer in _get_possible_nodes(for_mod_context, string_name):
            i += 1

            # This is a simple way to stop Jedi's dynamic param recursion
            # from going wild: The deeper Jedi's in the recursion, the less
            # code should be evaluated.
            if i * evaluator.dynamic_params_depth > MAX_PARAM_SEARCHES:
                return

            random_context = evaluator.create_context(for_mod_context, name)
            for function_execution in _check_name_for_execution(
                    evaluator, random_context, compare_node, name, trailer):
                found_executions = True
                yield function_execution

        # If there are results after processing a module, we're probably
        # good to process. This is a speed optimization.
        if found_executions:
            return


def _get_lambda_name(node):
    stmt = node.parent
    if stmt.type == 'expr_stmt':
        first_operator = next(stmt.yield_operators(), None)
        if first_operator == '=':
            first = stmt.children[0]
            if first.type == 'name':
                return first.value

    return None


def _get_possible_nodes(module_context, func_string_name):
    try:
        names = module_context.tree_node.get_used_names()[func_string_name]
    except KeyError:
        return

    for name in names:
        bracket = name.get_next_leaf()
        trailer = bracket.parent
        if trailer.type == 'trailer' and bracket == '(':
            yield name, trailer


def _check_name_for_execution(evaluator, context, compare_node, name, trailer):
    from jedi.evaluate.context.function import FunctionExecutionContext

    def create_func_excs():
        arglist = trailer.children[1]
        if arglist == ')':
            arglist = None
        args = TreeArguments(evaluator, context, arglist, trailer)
        if value_node.type == 'classdef':
            created_instance = instance.TreeInstance(
                evaluator,
                value.parent_context,
                value,
                args
            )
            for execution in created_instance.create_init_executions():
                yield execution
        else:
            yield value.get_function_execution(args)

    for value in evaluator.goto_definitions(context, name):
        value_node = value.tree_node
        if compare_node == value_node:
            for func_execution in create_func_excs():
                yield func_execution
        elif isinstance(value.parent_context, FunctionExecutionContext) and \
                compare_node.type == 'funcdef':
            # Here we're trying to find decorators by checking the first
            # parameter. It's not very generic though. Should find a better
            # solution that also applies to nested decorators.
            params = value.parent_context.get_params()
            if len(params) != 1:
                continue
            values = params[0].infer()
            nodes = [v.tree_node for v in values]
            if nodes == [compare_node]:
                # Found a decorator.
                module_context = context.get_root_context()
                execution_context = next(create_func_excs())
                for name, trailer in _get_possible_nodes(module_context, params[0].string_name):
                    if value_node.start_pos < name.start_pos < value_node.end_pos:
                        random_context = evaluator.create_context(execution_context, name)
                        iterator = _check_name_for_execution(
                            evaluator,
                            random_context,
                            compare_node,
                            name,
                            trailer
                        )
                        for function_execution in iterator:
                            yield function_execution
