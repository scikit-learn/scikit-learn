from parso.python import tree

from jedi._compatibility import zip_longest
from jedi import debug
from jedi.evaluate import analysis
from jedi.evaluate.lazy_context import LazyKnownContext, LazyKnownContexts, \
    LazyTreeContext, get_merged_lazy_context
from jedi.evaluate.filters import ParamName
from jedi.evaluate.base_context import NO_CONTEXTS
from jedi.evaluate.context import iterable
from jedi.evaluate.param import get_params, ExecutedParam


def try_iter_content(types, depth=0):
    """Helper method for static analysis."""
    if depth > 10:
        # It's possible that a loop has references on itself (especially with
        # CompiledObject). Therefore don't loop infinitely.
        return

    for typ in types:
        try:
            f = typ.py__iter__
        except AttributeError:
            pass
        else:
            for lazy_context in f():
                try_iter_content(lazy_context.infer(), depth + 1)


class AbstractArguments(object):
    context = None
    argument_node = None
    trailer = None

    def eval_argument_clinic(self, parameters):
        """Uses a list with argument clinic information (see PEP 436)."""
        iterator = self.unpack()
        for i, (name, optional, allow_kwargs) in enumerate(parameters):
            key, argument = next(iterator, (None, None))
            if key is not None:
                raise NotImplementedError
            if argument is None and not optional:
                debug.warning('TypeError: %s expected at least %s arguments, got %s',
                              name, len(parameters), i)
                raise ValueError
            values = NO_CONTEXTS if argument is None else argument.infer()

            if not values and not optional:
                # For the stdlib we always want values. If we don't get them,
                # that's ok, maybe something is too hard to resolve, however,
                # we will not proceed with the evaluation of that function.
                debug.warning('argument_clinic "%s" not resolvable.', name)
                raise ValueError
            yield values

    def eval_all(self, funcdef=None):
        """
        Evaluates all arguments as a support for static analysis
        (normally Jedi).
        """
        for key, lazy_context in self.unpack():
            types = lazy_context.infer()
            try_iter_content(types)

    def get_calling_nodes(self):
        raise NotImplementedError

    def unpack(self, funcdef=None):
        raise NotImplementedError

    def get_params(self, execution_context):
        return get_params(execution_context, self)


class AnonymousArguments(AbstractArguments):
    def get_params(self, execution_context):
        from jedi.evaluate.dynamic import search_params
        return search_params(
            execution_context.evaluator,
            execution_context,
            execution_context.tree_node
        )


class TreeArguments(AbstractArguments):
    def __init__(self, evaluator, context, argument_node, trailer=None):
        """
        The argument_node is either a parser node or a list of evaluated
        objects. Those evaluated objects may be lists of evaluated objects
        themselves (one list for the first argument, one for the second, etc).

        :param argument_node: May be an argument_node or a list of nodes.
        """
        self.argument_node = argument_node
        self.context = context
        self._evaluator = evaluator
        self.trailer = trailer  # Can be None, e.g. in a class definition.

    def _split(self):
        if self.argument_node is None:
            return

        # Allow testlist here as well for Python2's class inheritance
        # definitions.
        if not (self.argument_node.type in ('arglist', 'testlist') or (
                # in python 3.5 **arg is an argument, not arglist
                (self.argument_node.type == 'argument') and
                 self.argument_node.children[0] in ('*', '**'))):
            yield 0, self.argument_node
            return

        iterator = iter(self.argument_node.children)
        for child in iterator:
            if child == ',':
                continue
            elif child in ('*', '**'):
                yield len(child.value), next(iterator)
            elif child.type == 'argument' and \
                    child.children[0] in ('*', '**'):
                assert len(child.children) == 2
                yield len(child.children[0].value), child.children[1]
            else:
                yield 0, child

    def unpack(self, funcdef=None):
        named_args = []
        for star_count, el in self._split():
            if star_count == 1:
                arrays = self.context.eval_node(el)
                iterators = [_iterate_star_args(self.context, a, el, funcdef)
                             for a in arrays]
                for values in list(zip_longest(*iterators)):
                    # TODO zip_longest yields None, that means this would raise
                    # an exception?
                    yield None, get_merged_lazy_context(
                        [v for v in values if v is not None]
                    )
            elif star_count == 2:
                arrays = self.context.eval_node(el)
                for dct in arrays:
                    for key, values in _star_star_dict(self.context, dct, el, funcdef):
                        yield key, values
            else:
                if el.type == 'argument':
                    c = el.children
                    if len(c) == 3:  # Keyword argument.
                        named_args.append((c[0].value, LazyTreeContext(self.context, c[2]),))
                    else:  # Generator comprehension.
                        # Include the brackets with the parent.
                        comp = iterable.GeneratorComprehension(
                            self._evaluator, self.context, self.argument_node.parent)
                        yield None, LazyKnownContext(comp)
                else:
                    yield None, LazyTreeContext(self.context, el)

        # Reordering var_args is necessary, because star args sometimes appear
        # after named argument, but in the actual order it's prepended.
        for named_arg in named_args:
            yield named_arg

    def as_tree_tuple_objects(self):
        for star_count, argument in self._split():
            if argument.type == 'argument':
                argument, default = argument.children[::2]
            else:
                default = None
            yield argument, default, star_count

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.argument_node)

    def get_calling_nodes(self):
        from jedi.evaluate.dynamic import MergedExecutedParams
        old_arguments_list = []
        arguments = self

        while arguments not in old_arguments_list:
            if not isinstance(arguments, TreeArguments):
                break

            old_arguments_list.append(arguments)
            for name, default, star_count in reversed(list(arguments.as_tree_tuple_objects())):
                if not star_count or not isinstance(name, tree.Name):
                    continue

                names = self._evaluator.goto(arguments.context, name)
                if len(names) != 1:
                    break
                if not isinstance(names[0], ParamName):
                    break
                param = names[0].get_param()
                if isinstance(param, MergedExecutedParams):
                    # For dynamic searches we don't even want to see errors.
                    return []
                if not isinstance(param, ExecutedParam):
                    break
                if param.var_args is None:
                    break
                arguments = param.var_args
                break

        if arguments.argument_node is not None:
            return [arguments.argument_node]
        if arguments.trailer is not None:
            return [arguments.trailer]
        return []


class ValuesArguments(AbstractArguments):
    def __init__(self, values_list):
        self._values_list = values_list

    def unpack(self, funcdef=None):
        for values in self._values_list:
            yield None, LazyKnownContexts(values)

    def get_calling_nodes(self):
        return []

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self._values_list)


def _iterate_star_args(context, array, input_node, funcdef=None):
    try:
        iter_ = array.py__iter__
    except AttributeError:
        if funcdef is not None:
            # TODO this funcdef should not be needed.
            m = "TypeError: %s() argument after * must be a sequence, not %s" \
                % (funcdef.name.value, array)
            analysis.add(context, 'type-error-star', input_node, message=m)
    else:
        for lazy_context in iter_():
            yield lazy_context


def _star_star_dict(context, array, input_node, funcdef):
    from jedi.evaluate.context.instance import CompiledInstance
    if isinstance(array, CompiledInstance) and array.name.string_name == 'dict':
        # For now ignore this case. In the future add proper iterators and just
        # make one call without crazy isinstance checks.
        return {}
    elif isinstance(array, iterable.Sequence) and array.array_type == 'dict':
        return array.exact_key_items()
    else:
        if funcdef is not None:
            m = "TypeError: %s argument after ** must be a mapping, not %s" \
                % (funcdef.name.value, array)
            analysis.add(context, 'type-error-star-star', input_node, message=m)
        return {}
