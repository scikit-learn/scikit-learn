"""
Filters are objects that you can use to filter names in different scopes. They
are needed for name resolution.
"""
from abc import abstractmethod

from parso.tree import search_ancestor

from jedi._compatibility import use_metaclass, Parameter
from jedi.cache import memoize_method
from jedi.evaluate import flow_analysis
from jedi.evaluate.base_context import ContextSet, Context
from jedi.parser_utils import get_parent_scope
from jedi.evaluate.utils import to_list


class AbstractNameDefinition(object):
    start_pos = None
    string_name = None
    parent_context = None
    tree_name = None

    @abstractmethod
    def infer(self):
        raise NotImplementedError

    @abstractmethod
    def goto(self):
        # Typically names are already definitions and therefore a goto on that
        # name will always result on itself.
        return {self}

    def get_root_context(self):
        return self.parent_context.get_root_context()

    def __repr__(self):
        if self.start_pos is None:
            return '<%s: %s>' % (self.__class__.__name__, self.string_name)
        return '<%s: %s@%s>' % (self.__class__.__name__, self.string_name, self.start_pos)

    def execute(self, arguments):
        return self.infer().execute(arguments)

    def execute_evaluated(self, *args, **kwargs):
        return self.infer().execute_evaluated(*args, **kwargs)

    def is_import(self):
        return False

    @property
    def api_type(self):
        return self.parent_context.api_type


class AbstractTreeName(AbstractNameDefinition):
    def __init__(self, parent_context, tree_name):
        self.parent_context = parent_context
        self.tree_name = tree_name

    def goto(self):
        return self.parent_context.evaluator.goto(self.parent_context, self.tree_name)

    def is_import(self):
        imp = search_ancestor(self.tree_name, 'import_from', 'import_name')
        return imp is not None

    @property
    def string_name(self):
        return self.tree_name.value

    @property
    def start_pos(self):
        return self.tree_name.start_pos


class ContextNameMixin(object):
    def infer(self):
        return ContextSet(self._context)

    def get_root_context(self):
        if self.parent_context is None:
            return self._context
        return super(ContextNameMixin, self).get_root_context()

    @property
    def api_type(self):
        return self._context.api_type


class ContextName(ContextNameMixin, AbstractTreeName):
    def __init__(self, context, tree_name):
        super(ContextName, self).__init__(context.parent_context, tree_name)
        self._context = context


class TreeNameDefinition(AbstractTreeName):
    _API_TYPES = dict(
        import_name='module',
        import_from='module',
        funcdef='function',
        param='param',
        classdef='class',
    )

    def infer(self):
        # Refactor this, should probably be here.
        from jedi.evaluate.syntax_tree import tree_name_to_contexts
        return tree_name_to_contexts(self.parent_context.evaluator, self.parent_context, self.tree_name)

    @property
    def api_type(self):
        definition = self.tree_name.get_definition(import_name_always=True)
        if definition is None:
            return 'statement'
        return self._API_TYPES.get(definition.type, 'statement')


class ParamName(AbstractTreeName):
    api_type = u'param'

    def __init__(self, parent_context, tree_name):
        self.parent_context = parent_context
        self.tree_name = tree_name

    def get_kind(self):
        tree_param = search_ancestor(self.tree_name, 'param')
        if tree_param.star_count == 1:  # *args
            return Parameter.VAR_POSITIONAL
        if tree_param.star_count == 2:  # **kwargs
            return Parameter.VAR_KEYWORD

        parent = tree_param.parent
        for p in parent.children:
            if p.type == 'param':
                if p.star_count:
                    return Parameter.KEYWORD_ONLY
                if p == tree_param:
                    break
        return Parameter.POSITIONAL_OR_KEYWORD

    def infer(self):
        return self.get_param().infer()

    def get_param(self):
        params = self.parent_context.get_params()
        param_node = search_ancestor(self.tree_name, 'param')
        return params[param_node.position_index]


class AnonymousInstanceParamName(ParamName):
    def infer(self):
        param_node = search_ancestor(self.tree_name, 'param')
        # TODO I think this should not belong here. It's not even really true,
        #      because classmethod and other descriptors can change it.
        if param_node.position_index == 0:
            # This is a speed optimization, to return the self param (because
            # it's known). This only affects anonymous instances.
            return ContextSet(self.parent_context.instance)
        else:
            return self.get_param().infer()


class AbstractFilter(object):
    _until_position = None

    def _filter(self, names):
        if self._until_position is not None:
            return [n for n in names if n.start_pos < self._until_position]
        return names

    @abstractmethod
    def get(self, name):
        raise NotImplementedError

    @abstractmethod
    def values(self):
        raise NotImplementedError


class AbstractUsedNamesFilter(AbstractFilter):
    name_class = TreeNameDefinition

    def __init__(self, context, parser_scope):
        self._parser_scope = parser_scope
        self._used_names = self._parser_scope.get_root_node().get_used_names()
        self.context = context

    def get(self, name):
        try:
            names = self._used_names[name]
        except KeyError:
            return []

        return self._convert_names(self._filter(names))

    def _convert_names(self, names):
        return [self.name_class(self.context, name) for name in names]

    def values(self):
        return self._convert_names(name for name_list in self._used_names.values()
                                   for name in self._filter(name_list))

    def __repr__(self):
        return '<%s: %s>' % (self.__class__.__name__, self.context)


class ParserTreeFilter(AbstractUsedNamesFilter):
    def __init__(self, evaluator, context, node_context=None, until_position=None,
                 origin_scope=None):
        """
        node_context is an option to specify a second context for use cases
        like the class mro where the parent class of a new name would be the
        context, but for some type inference it's important to have a local
        context of the other classes.
        """
        if node_context is None:
            node_context = context
        super(ParserTreeFilter, self).__init__(context, node_context.tree_node)
        self._node_context = node_context
        self._origin_scope = origin_scope
        self._until_position = until_position

    def _filter(self, names):
        names = super(ParserTreeFilter, self)._filter(names)
        names = [n for n in names if self._is_name_reachable(n)]
        return list(self._check_flows(names))

    def _is_name_reachable(self, name):
        if not name.is_definition():
            return False
        parent = name.parent
        if parent.type == 'trailer':
            return False
        base_node = parent if parent.type in ('classdef', 'funcdef') else name
        return get_parent_scope(base_node) == self._parser_scope

    def _check_flows(self, names):
        for name in sorted(names, key=lambda name: name.start_pos, reverse=True):
            check = flow_analysis.reachability_check(
                context=self._node_context,
                context_scope=self._parser_scope,
                node=name,
                origin_scope=self._origin_scope
            )
            if check is not flow_analysis.UNREACHABLE:
                yield name

            if check is flow_analysis.REACHABLE:
                break


class FunctionExecutionFilter(ParserTreeFilter):
    param_name = ParamName

    def __init__(self, evaluator, context, node_context=None,
                 until_position=None, origin_scope=None):
        super(FunctionExecutionFilter, self).__init__(
            evaluator,
            context,
            node_context,
            until_position,
            origin_scope
        )

    @to_list
    def _convert_names(self, names):
        for name in names:
            param = search_ancestor(name, 'param')
            if param:
                yield self.param_name(self.context, name)
            else:
                yield TreeNameDefinition(self.context, name)


class AnonymousInstanceFunctionExecutionFilter(FunctionExecutionFilter):
    param_name = AnonymousInstanceParamName


class GlobalNameFilter(AbstractUsedNamesFilter):
    def __init__(self, context, parser_scope):
        super(GlobalNameFilter, self).__init__(context, parser_scope)

    @to_list
    def _filter(self, names):
        for name in names:
            if name.parent.type == 'global_stmt':
                yield name


class DictFilter(AbstractFilter):
    def __init__(self, dct):
        self._dct = dct

    def get(self, name):
        try:
            value = self._convert(name, self._dct[name])
        except KeyError:
            return []
        else:
            return list(self._filter([value]))

    def values(self):
        def yielder():
            for item in self._dct.items():
                try:
                    yield self._convert(*item)
                except KeyError:
                    pass
        return self._filter(yielder())

    def _convert(self, name, value):
        return value


class MergedFilter(object):
    def __init__(self, *filters):
        self._filters = filters

    def get(self, name):
        return [n for filter in self._filters for n in filter.get(name)]

    def values(self):
        return [n for filter in self._filters for n in filter.values()]

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, ', '.join(str(f) for f in self._filters))


class _BuiltinMappedMethod(Context):
    """``Generator.__next__`` ``dict.values`` methods and so on."""
    api_type = u'function'

    def __init__(self, builtin_context, method, builtin_func):
        super(_BuiltinMappedMethod, self).__init__(
            builtin_context.evaluator,
            parent_context=builtin_context
        )
        self._method = method
        self._builtin_func = builtin_func

    def py__call__(self, params):
        # TODO add TypeError if params are given/or not correct.
        return self._method(self.parent_context)

    def __getattr__(self, name):
        return getattr(self._builtin_func, name)


class SpecialMethodFilter(DictFilter):
    """
    A filter for methods that are defined in this module on the corresponding
    classes like Generator (for __next__, etc).
    """
    class SpecialMethodName(AbstractNameDefinition):
        api_type = u'function'

        def __init__(self, parent_context, string_name, value, builtin_context):
            callable_, python_version = value
            if python_version is not None and \
                    python_version != parent_context.evaluator.environment.version_info.major:
                raise KeyError

            self.parent_context = parent_context
            self.string_name = string_name
            self._callable = callable_
            self._builtin_context = builtin_context

        def infer(self):
            for filter in self._builtin_context.get_filters():
                # We can take the first index, because on builtin methods there's
                # always only going to be one name. The same is true for the
                # inferred values.
                for name in filter.get(self.string_name):
                    builtin_func = next(iter(name.infer()))
                    break
                else:
                    continue
                break
            return ContextSet(
                _BuiltinMappedMethod(self.parent_context, self._callable, builtin_func)
            )

    def __init__(self, context, dct, builtin_context):
        super(SpecialMethodFilter, self).__init__(dct)
        self.context = context
        self._builtin_context = builtin_context
        """
        This context is what will be used to introspect the name, where as the
        other context will be used to execute the function.

        We distinguish, because we have to.
        """

    def _convert(self, name, value):
        return self.SpecialMethodName(self.context, name, value, self._builtin_context)


class _OverwriteMeta(type):
    def __init__(cls, name, bases, dct):
        super(_OverwriteMeta, cls).__init__(name, bases, dct)

        base_dct = {}
        for base_cls in reversed(cls.__bases__):
            try:
                base_dct.update(base_cls.overwritten_methods)
            except AttributeError:
                pass

        for func in cls.__dict__.values():
            try:
                base_dct.update(func.registered_overwritten_methods)
            except AttributeError:
                pass
        cls.overwritten_methods = base_dct


class AbstractObjectOverwrite(use_metaclass(_OverwriteMeta, object)):
    def get_object(self):
        raise NotImplementedError

    def get_filters(self, search_global, *args, **kwargs):
        yield SpecialMethodFilter(self, self.overwritten_methods, self.get_object())

        for filter in self.get_object().get_filters(search_global):
            yield filter


class BuiltinOverwrite(Context, AbstractObjectOverwrite):
    special_object_identifier = None

    def __init__(self, evaluator):
        super(BuiltinOverwrite, self).__init__(evaluator, evaluator.builtins_module)

    @memoize_method
    def get_object(self):
        from jedi.evaluate import compiled
        assert self.special_object_identifier
        return compiled.get_special_object(self.evaluator, self.special_object_identifier)

    def py__class__(self):
        return self.get_object().py__class__()


def publish_method(method_name, python_version_match=None):
    def decorator(func):
        dct = func.__dict__.setdefault('registered_overwritten_methods', {})
        dct[method_name] = func, python_version_match
        return func
    return decorator


def get_global_filters(evaluator, context, until_position, origin_scope):
    """
    Returns all filters in order of priority for name resolution.

    For global name lookups. The filters will handle name resolution
    themselves, but here we gather possible filters downwards.

    >>> from jedi._compatibility import u, no_unicode_pprint
    >>> from jedi import Script
    >>> script = Script(u('''
    ... x = ['a', 'b', 'c']
    ... def func():
    ...     y = None
    ... '''))
    >>> module_node = script._module_node
    >>> scope = next(module_node.iter_funcdefs())
    >>> scope
    <Function: func@3-5>
    >>> context = script._get_module().create_context(scope)
    >>> filters = list(get_global_filters(context.evaluator, context, (4, 0), None))

    First we get the names from the function scope.

    >>> no_unicode_pprint(filters[0])                    #doctest: +ELLIPSIS
    MergedFilter(<ParserTreeFilter: ...>, <GlobalNameFilter: ...>)
    >>> sorted(str(n) for n in filters[0].values())
    ['<TreeNameDefinition: func@(3, 4)>', '<TreeNameDefinition: x@(2, 0)>']
    >>> filters[0]._filters[0]._until_position
    (4, 0)
    >>> filters[0]._filters[1]._until_position

    Then it yields the names from one level "lower". In this example, this is
    the module scope (including globals).
    As a side note, you can see, that the position in the filter is None on the
    globals filter, because there the whole module is searched.

    >>> list(filters[1].values())  # package modules -> Also empty.
    []
    >>> sorted(name.string_name for name in filters[2].values())  # Module attributes
    ['__doc__', '__file__', '__name__', '__package__']

    Finally, it yields the builtin filter, if `include_builtin` is
    true (default).

    >>> filters[3].values()                              #doctest: +ELLIPSIS
    [<CompiledName: ...>, ...]
    """
    from jedi.evaluate.context.function import FunctionExecutionContext
    while context is not None:
        # Names in methods cannot be resolved within the class.
        for filter in context.get_filters(
                search_global=True,
                until_position=until_position,
                origin_scope=origin_scope):
            yield filter
        if isinstance(context, FunctionExecutionContext):
            # The position should be reset if the current scope is a function.
            until_position = None

        context = context.parent_context

    # Add builtins to the global scope.
    for filter in evaluator.builtins_module.get_filters(search_global=True):
        yield filter
