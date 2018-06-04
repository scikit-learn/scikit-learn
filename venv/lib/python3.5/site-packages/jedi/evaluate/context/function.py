from parso.python import tree

from jedi._compatibility import use_metaclass
from jedi import debug
from jedi.evaluate.cache import evaluator_method_cache, CachedMetaClass
from jedi.evaluate import compiled
from jedi.evaluate import recursion
from jedi.evaluate import docstrings
from jedi.evaluate import pep0484
from jedi.evaluate import flow_analysis
from jedi.evaluate import helpers
from jedi.evaluate.arguments import AnonymousArguments
from jedi.evaluate.filters import ParserTreeFilter, FunctionExecutionFilter, \
    ContextName, AbstractNameDefinition, ParamName
from jedi.evaluate.base_context import ContextualizedNode, NO_CONTEXTS, \
    ContextSet, TreeContext
from jedi.evaluate.lazy_context import LazyKnownContexts, LazyKnownContext, \
    LazyTreeContext
from jedi.evaluate.context import iterable
from jedi.evaluate.context import asynchronous
from jedi import parser_utils
from jedi.evaluate.parser_cache import get_yield_exprs


class LambdaName(AbstractNameDefinition):
    string_name = '<lambda>'
    api_type = u'function'

    def __init__(self, lambda_context):
        self._lambda_context = lambda_context
        self.parent_context = lambda_context.parent_context

    @property
    def start_pos(self):
        return self._lambda_context.tree_node.start_pos

    def infer(self):
        return ContextSet(self._lambda_context)


class FunctionContext(use_metaclass(CachedMetaClass, TreeContext)):
    """
    Needed because of decorators. Decorators are evaluated here.
    """
    api_type = u'function'

    def __init__(self, evaluator, parent_context, funcdef):
        """ This should not be called directly """
        super(FunctionContext, self).__init__(evaluator, parent_context)
        self.tree_node = funcdef

    def get_filters(self, search_global, until_position=None, origin_scope=None):
        if search_global:
            yield ParserTreeFilter(
                self.evaluator,
                context=self,
                until_position=until_position,
                origin_scope=origin_scope
            )
        else:
            scope = self.py__class__()
            for filter in scope.get_filters(search_global=False, origin_scope=origin_scope):
                yield filter

    def infer_function_execution(self, function_execution):
        """
        Created to be used by inheritance.
        """
        is_coroutine = self.tree_node.parent.type == 'async_stmt'
        is_generator = bool(get_yield_exprs(self.evaluator, self.tree_node))

        if is_coroutine:
            if is_generator:
                if self.evaluator.environment.version_info < (3, 6):
                    return NO_CONTEXTS
                return ContextSet(asynchronous.AsyncGenerator(self.evaluator, function_execution))
            else:
                if self.evaluator.environment.version_info < (3, 5):
                    return NO_CONTEXTS
                return ContextSet(asynchronous.Coroutine(self.evaluator, function_execution))
        else:
            if is_generator:
                return ContextSet(iterable.Generator(self.evaluator, function_execution))
            else:
                return function_execution.get_return_values()

    def get_function_execution(self, arguments=None):
        if arguments is None:
            arguments = AnonymousArguments()

        return FunctionExecutionContext(self.evaluator, self.parent_context, self, arguments)

    def py__call__(self, arguments):
        function_execution = self.get_function_execution(arguments)
        return self.infer_function_execution(function_execution)

    def py__class__(self):
        # This differentiation is only necessary for Python2. Python3 does not
        # use a different method class.
        if isinstance(parser_utils.get_parent_scope(self.tree_node), tree.Class):
            name = u'METHOD_CLASS'
        else:
            name = u'FUNCTION_CLASS'
        return compiled.get_special_object(self.evaluator, name)

    @property
    def name(self):
        if self.tree_node.type == 'lambdef':
            return LambdaName(self)
        return ContextName(self, self.tree_node.name)

    def get_param_names(self):
        function_execution = self.get_function_execution()
        return [ParamName(function_execution, param.name)
                for param in self.tree_node.get_params()]


class FunctionExecutionContext(TreeContext):
    """
    This class is used to evaluate functions and their returns.

    This is the most complicated class, because it contains the logic to
    transfer parameters. It is even more complicated, because there may be
    multiple calls to functions and recursion has to be avoided. But this is
    responsibility of the decorators.
    """
    function_execution_filter = FunctionExecutionFilter

    def __init__(self, evaluator, parent_context, function_context, var_args):
        super(FunctionExecutionContext, self).__init__(evaluator, parent_context)
        self.function_context = function_context
        self.tree_node = function_context.tree_node
        self.var_args = var_args

    @evaluator_method_cache(default=NO_CONTEXTS)
    @recursion.execution_recursion_decorator()
    def get_return_values(self, check_yields=False):
        funcdef = self.tree_node
        if funcdef.type == 'lambdef':
            return self.eval_node(funcdef.children[-1])

        if check_yields:
            context_set = NO_CONTEXTS
            returns = get_yield_exprs(self.evaluator, funcdef)
        else:
            returns = funcdef.iter_return_stmts()
            context_set = docstrings.infer_return_types(self.function_context)
            context_set |= pep0484.infer_return_types(self.function_context)

        for r in returns:
            check = flow_analysis.reachability_check(self, funcdef, r)
            if check is flow_analysis.UNREACHABLE:
                debug.dbg('Return unreachable: %s', r)
            else:
                if check_yields:
                    context_set |= ContextSet.from_sets(
                        lazy_context.infer()
                        for lazy_context in self._get_yield_lazy_context(r)
                    )
                else:
                    try:
                        children = r.children
                    except AttributeError:
                        ctx = compiled.builtin_from_name(self.evaluator, u'None')
                        context_set |= ContextSet(ctx)
                    else:
                        context_set |= self.eval_node(children[1])
            if check is flow_analysis.REACHABLE:
                debug.dbg('Return reachable: %s', r)
                break
        return context_set

    def _get_yield_lazy_context(self, yield_expr):
        if yield_expr.type == 'keyword':
            # `yield` just yields None.
            ctx = compiled.builtin_from_name(self.evaluator, u'None')
            yield LazyKnownContext(ctx)
            return

        node = yield_expr.children[1]
        if node.type == 'yield_arg':  # It must be a yield from.
            cn = ContextualizedNode(self, node.children[1])
            for lazy_context in cn.infer().iterate(cn):
                yield lazy_context
        else:
            yield LazyTreeContext(self, node)

    @recursion.execution_recursion_decorator(default=iter([]))
    def get_yield_lazy_contexts(self, is_async=False):
        # TODO: if is_async, wrap yield statements in Awaitable/async_generator_asend
        for_parents = [(y, tree.search_ancestor(y, 'for_stmt', 'funcdef',
                                                'while_stmt', 'if_stmt'))
                       for y in get_yield_exprs(self.evaluator, self.tree_node)]

        # Calculate if the yields are placed within the same for loop.
        yields_order = []
        last_for_stmt = None
        for yield_, for_stmt in for_parents:
            # For really simple for loops we can predict the order. Otherwise
            # we just ignore it.
            parent = for_stmt.parent
            if parent.type == 'suite':
                parent = parent.parent
            if for_stmt.type == 'for_stmt' and parent == self.tree_node \
                    and parser_utils.for_stmt_defines_one_name(for_stmt):  # Simplicity for now.
                if for_stmt == last_for_stmt:
                    yields_order[-1][1].append(yield_)
                else:
                    yields_order.append((for_stmt, [yield_]))
            elif for_stmt == self.tree_node:
                yields_order.append((None, [yield_]))
            else:
                types = self.get_return_values(check_yields=True)
                if types:
                    yield LazyKnownContexts(types)
                return
            last_for_stmt = for_stmt

        for for_stmt, yields in yields_order:
            if for_stmt is None:
                # No for_stmt, just normal yields.
                for yield_ in yields:
                    for result in self._get_yield_lazy_context(yield_):
                        yield result
            else:
                input_node = for_stmt.get_testlist()
                cn = ContextualizedNode(self, input_node)
                ordered = cn.infer().iterate(cn)
                ordered = list(ordered)
                for lazy_context in ordered:
                    dct = {str(for_stmt.children[1].value): lazy_context.infer()}
                    with helpers.predefine_names(self, for_stmt, dct):
                        for yield_in_same_for_stmt in yields:
                            for result in self._get_yield_lazy_context(yield_in_same_for_stmt):
                                yield result

    def get_filters(self, search_global, until_position=None, origin_scope=None):
        yield self.function_execution_filter(self.evaluator, self,
                                             until_position=until_position,
                                             origin_scope=origin_scope)

    @evaluator_method_cache()
    def get_params(self):
        return self.var_args.get_params(self)
