"""
Contains all classes and functions to deal with lists, dicts, generators and
iterators in general.

Array modifications
*******************

If the content of an array (``set``/``list``) is requested somewhere, the
current module will be checked for appearances of ``arr.append``,
``arr.insert``, etc.  If the ``arr`` name points to an actual array, the
content will be added

This can be really cpu intensive, as you can imagine. Because |jedi| has to
follow **every** ``append`` and check wheter it's the right array. However this
works pretty good, because in *slow* cases, the recursion detector and other
settings will stop this process.

It is important to note that:

1. Array modfications work only in the current module.
2. Jedi only checks Array additions; ``list.pop``, etc are ignored.
"""
from jedi import debug
from jedi import settings
from jedi._compatibility import force_unicode, is_py3
from jedi.cache import memoize_method
from jedi.evaluate import compiled
from jedi.evaluate import analysis
from jedi.evaluate import recursion
from jedi.evaluate.lazy_context import LazyKnownContext, LazyKnownContexts, \
    LazyTreeContext
from jedi.evaluate.helpers import get_int_or_none, is_string, \
    predefine_names, evaluate_call_of_leaf
from jedi.evaluate.utils import safe_property
from jedi.evaluate.utils import to_list
from jedi.evaluate.cache import evaluator_method_cache
from jedi.evaluate.filters import ParserTreeFilter, BuiltinOverwrite, \
    publish_method
from jedi.evaluate.base_context import ContextSet, NO_CONTEXTS, Context, \
    TreeContext, ContextualizedNode
from jedi.parser_utils import get_comp_fors


class IterableMixin(object):
    def py__stop_iteration_returns(self):
        return ContextSet(compiled.builtin_from_name(self.evaluator, u'None'))


class GeneratorBase(BuiltinOverwrite, IterableMixin):
    array_type = None
    special_object_identifier = u'GENERATOR_OBJECT'

    @publish_method('send')
    @publish_method('next', python_version_match=2)
    @publish_method('__next__', python_version_match=3)
    def py__next__(self):
        return ContextSet.from_sets(lazy_context.infer() for lazy_context in self.py__iter__())

    @property
    def name(self):
        return compiled.CompiledContextName(self, 'generator')


class Generator(GeneratorBase):
    """Handling of `yield` functions."""
    def __init__(self, evaluator, func_execution_context):
        super(Generator, self).__init__(evaluator)
        self._func_execution_context = func_execution_context

    def py__iter__(self):
        return self._func_execution_context.get_yield_lazy_contexts()

    def py__stop_iteration_returns(self):
        return self._func_execution_context.get_return_values()

    def __repr__(self):
        return "<%s of %s>" % (type(self).__name__, self._func_execution_context)


class CompForContext(TreeContext):
    @classmethod
    def from_comp_for(cls, parent_context, comp_for):
        return cls(parent_context.evaluator, parent_context, comp_for)

    def __init__(self, evaluator, parent_context, comp_for):
        super(CompForContext, self).__init__(evaluator, parent_context)
        self.tree_node = comp_for

    def get_node(self):
        return self.tree_node

    def get_filters(self, search_global, until_position=None, origin_scope=None):
        yield ParserTreeFilter(self.evaluator, self)


def comprehension_from_atom(evaluator, context, atom):
    bracket = atom.children[0]
    if bracket == '{':
        if atom.children[1].children[1] == ':':
            cls = DictComprehension
        else:
            cls = SetComprehension
    elif bracket == '(':
        cls = GeneratorComprehension
    elif bracket == '[':
        cls = ListComprehension
    return cls(evaluator, context, atom)


class ComprehensionMixin(object):
    def __init__(self, evaluator, defining_context, atom):
        super(ComprehensionMixin, self).__init__(evaluator)
        self._defining_context = defining_context
        self._atom = atom

    def _get_comprehension(self):
        "return 'a for a in b'"
        # The atom contains a testlist_comp
        return self._atom.children[1]

    def _get_comp_for(self):
        "return CompFor('for a in b')"
        return self._get_comprehension().children[1]

    def _eval_node(self, index=0):
        """
        The first part `x + 1` of the list comprehension:

            [x + 1 for x in foo]
        """
        return self._get_comprehension().children[index]

    @evaluator_method_cache()
    def _get_comp_for_context(self, parent_context, comp_for):
        # TODO shouldn't this be part of create_context?
        return CompForContext.from_comp_for(parent_context, comp_for)

    def _nested(self, comp_fors, parent_context=None):
        comp_for = comp_fors[0]

        is_async = 'async' == comp_for.children[comp_for.children.index('for') - 1]

        input_node = comp_for.children[comp_for.children.index('in') + 1]
        parent_context = parent_context or self._defining_context
        input_types = parent_context.eval_node(input_node)
        # TODO: simulate await if self.is_async

        cn = ContextualizedNode(parent_context, input_node)
        iterated = input_types.iterate(cn, is_async=is_async)
        exprlist = comp_for.children[comp_for.children.index('for') + 1]
        for i, lazy_context in enumerate(iterated):
            types = lazy_context.infer()
            dct = unpack_tuple_to_dict(parent_context, types, exprlist)
            context_ = self._get_comp_for_context(
                parent_context,
                comp_for,
            )
            with predefine_names(context_, comp_for, dct):
                try:
                    for result in self._nested(comp_fors[1:], context_):
                        yield result
                except IndexError:
                    iterated = context_.eval_node(self._eval_node())
                    if self.array_type == 'dict':
                        yield iterated, context_.eval_node(self._eval_node(2))
                    else:
                        yield iterated

    @evaluator_method_cache(default=[])
    @to_list
    def _iterate(self):
        comp_fors = tuple(get_comp_fors(self._get_comp_for()))
        for result in self._nested(comp_fors):
            yield result

    def py__iter__(self):
        for set_ in self._iterate():
            yield LazyKnownContexts(set_)

    def __repr__(self):
        return "<%s of %s>" % (type(self).__name__, self._atom)


class Sequence(BuiltinOverwrite, IterableMixin):
    api_type = u'instance'

    @property
    def name(self):
        return compiled.CompiledContextName(self, self.array_type)

    @memoize_method
    def get_object(self):
        compiled_obj = compiled.builtin_from_name(self.evaluator, self.array_type)
        only_obj, = compiled_obj.execute_evaluated(self)
        return only_obj

    def py__bool__(self):
        return None  # We don't know the length, because of appends.

    def py__class__(self):
        return compiled.builtin_from_name(self.evaluator, self.array_type)

    @safe_property
    def parent(self):
        return self.evaluator.builtins_module

    def dict_values(self):
        return ContextSet.from_sets(
            self._defining_context.eval_node(v)
            for k, v in self._items()
        )


class ListComprehension(ComprehensionMixin, Sequence):
    array_type = u'list'

    def py__getitem__(self, index):
        if isinstance(index, slice):
            return ContextSet(self)

        all_types = list(self.py__iter__())
        return all_types[index].infer()


class SetComprehension(ComprehensionMixin, Sequence):
    array_type = u'set'


class DictComprehension(ComprehensionMixin, Sequence):
    array_type = u'dict'

    def _get_comp_for(self):
        return self._get_comprehension().children[3]

    def py__iter__(self):
        for keys, values in self._iterate():
            yield LazyKnownContexts(keys)

    def py__getitem__(self, index):
        for keys, values in self._iterate():
            for k in keys:
                if isinstance(k, compiled.CompiledObject):
                    if k.get_safe_value(default=object()) == index:
                        return values
        return self.dict_values()

    def dict_values(self):
        return ContextSet.from_sets(values for keys, values in self._iterate())

    @publish_method('values')
    def _imitate_values(self):
        lazy_context = LazyKnownContexts(self.dict_values())
        return ContextSet(FakeSequence(self.evaluator, u'list', [lazy_context]))

    @publish_method('items')
    def _imitate_items(self):
        items = ContextSet.from_iterable(
            FakeSequence(
                self.evaluator, u'tuple'
                (LazyKnownContexts(keys), LazyKnownContexts(values))
            ) for keys, values in self._iterate()
        )

        return create_evaluated_sequence_set(self.evaluator, items, sequence_type=u'list')


class GeneratorComprehension(ComprehensionMixin, GeneratorBase):
    pass


class SequenceLiteralContext(Sequence):
    mapping = {'(': u'tuple',
               '[': u'list',
               '{': u'set'}

    def __init__(self, evaluator, defining_context, atom):
        super(SequenceLiteralContext, self).__init__(evaluator)
        self.atom = atom
        self._defining_context = defining_context

        if self.atom.type in ('testlist_star_expr', 'testlist'):
            self.array_type = u'tuple'
        else:
            self.array_type = SequenceLiteralContext.mapping[atom.children[0]]
            """The builtin name of the array (list, set, tuple or dict)."""

    def py__getitem__(self, index):
        """Here the index is an int/str. Raises IndexError/KeyError."""
        if self.array_type == u'dict':
            compiled_obj_index = compiled.create_simple_object(self.evaluator, index)
            for key, value in self._items():
                for k in self._defining_context.eval_node(key):
                    if isinstance(k, compiled.CompiledObject) \
                            and k.execute_operation(compiled_obj_index, u'==').get_safe_value():
                        return self._defining_context.eval_node(value)
            raise KeyError('No key found in dictionary %s.' % self)

        # Can raise an IndexError
        if isinstance(index, slice):
            return ContextSet(self)
        else:
            return self._defining_context.eval_node(self._items()[index])

    def py__iter__(self):
        """
        While values returns the possible values for any array field, this
        function returns the value for a certain index.
        """
        if self.array_type == u'dict':
            # Get keys.
            types = ContextSet()
            for k, _ in self._items():
                types |= self._defining_context.eval_node(k)
            # We don't know which dict index comes first, therefore always
            # yield all the types.
            for _ in types:
                yield LazyKnownContexts(types)
        else:
            for node in self._items():
                yield LazyTreeContext(self._defining_context, node)

            for addition in check_array_additions(self._defining_context, self):
                yield addition

    def _values(self):
        """Returns a list of a list of node."""
        if self.array_type == u'dict':
            return ContextSet.from_sets(v for k, v in self._items())
        else:
            return self._items()

    def _items(self):
        c = self.atom.children

        if self.atom.type in ('testlist_star_expr', 'testlist'):
            return c[::2]

        array_node = c[1]
        if array_node in (']', '}', ')'):
            return []  # Direct closing bracket, doesn't contain items.

        if array_node.type == 'testlist_comp':
            return array_node.children[::2]
        elif array_node.type == 'dictorsetmaker':
            kv = []
            iterator = iter(array_node.children)
            for key in iterator:
                op = next(iterator, None)
                if op is None or op == ',':
                    kv.append(key)  # A set.
                else:
                    assert op == ':'  # A dict.
                    kv.append((key, next(iterator)))
                    next(iterator, None)  # Possible comma.
            return kv
        else:
            return [array_node]

    def exact_key_items(self):
        """
        Returns a generator of tuples like dict.items(), where the key is
        resolved (as a string) and the values are still lazy contexts.
        """
        for key_node, value in self._items():
            for key in self._defining_context.eval_node(key_node):
                if is_string(key):
                    yield key.get_safe_value(), LazyTreeContext(self._defining_context, value)

    def __repr__(self):
        return "<%s of %s>" % (self.__class__.__name__, self.atom)


class DictLiteralContext(SequenceLiteralContext):
    array_type = u'dict'

    def __init__(self, evaluator, defining_context, atom):
        super(SequenceLiteralContext, self).__init__(evaluator)
        self._defining_context = defining_context
        self.atom = atom

    @publish_method('values')
    def _imitate_values(self):
        lazy_context = LazyKnownContexts(self.dict_values())
        return ContextSet(FakeSequence(self.evaluator, u'list', [lazy_context]))

    @publish_method('items')
    def _imitate_items(self):
        lazy_contexts = [
            LazyKnownContext(FakeSequence(
                self.evaluator, u'tuple',
                (LazyTreeContext(self._defining_context, key_node),
                 LazyTreeContext(self._defining_context, value_node))
            )) for key_node, value_node in self._items()
        ]

        return ContextSet(FakeSequence(self.evaluator, u'list', lazy_contexts))


class _FakeArray(SequenceLiteralContext):
    def __init__(self, evaluator, container, type):
        super(SequenceLiteralContext, self).__init__(evaluator)
        self.array_type = type
        self.atom = container
        # TODO is this class really needed?


class FakeSequence(_FakeArray):
    def __init__(self, evaluator, array_type, lazy_context_list):
        """
        type should be one of "tuple", "list"
        """
        super(FakeSequence, self).__init__(evaluator, None, array_type)
        self._lazy_context_list = lazy_context_list

    def py__getitem__(self, index):
        return self._lazy_context_list[index].infer()

    def py__iter__(self):
        return self._lazy_context_list

    def py__bool__(self):
        return bool(len(self._lazy_context_list))

    def __repr__(self):
        return "<%s of %s>" % (type(self).__name__, self._lazy_context_list)


class FakeDict(_FakeArray):
    def __init__(self, evaluator, dct):
        super(FakeDict, self).__init__(evaluator, dct, u'dict')
        self._dct = dct

    def py__iter__(self):
        for key in self._dct:
            yield LazyKnownContext(compiled.create_simple_object(self.evaluator, key))

    def py__getitem__(self, index):
        if is_py3 and self.evaluator.environment.version_info.major == 2:
            # In Python 2 bytes and unicode compare.
            if isinstance(index, bytes):
                index_unicode = force_unicode(index)
                try:
                    return self._dct[index_unicode].infer()
                except KeyError:
                    pass
            elif isinstance(index, str):
                index_bytes = index.encode('utf-8')
                try:
                    return self._dct[index_bytes].infer()
                except KeyError:
                    pass

        return self._dct[index].infer()

    @publish_method('values')
    def _values(self):
        return ContextSet(FakeSequence(
            self.evaluator, u'tuple',
            [LazyKnownContexts(self.dict_values())]
        ))

    def dict_values(self):
        return ContextSet.from_sets(lazy_context.infer() for lazy_context in self._dct.values())

    def exact_key_items(self):
        return self._dct.items()


class MergedArray(_FakeArray):
    def __init__(self, evaluator, arrays):
        super(MergedArray, self).__init__(evaluator, arrays, arrays[-1].array_type)
        self._arrays = arrays

    def py__iter__(self):
        for array in self._arrays:
            for lazy_context in array.py__iter__():
                yield lazy_context

    def py__getitem__(self, index):
        return ContextSet.from_sets(lazy_context.infer() for lazy_context in self.py__iter__())

    def _items(self):
        for array in self._arrays:
            for a in array._items():
                yield a

    def __len__(self):
        return sum(len(a) for a in self._arrays)


def unpack_tuple_to_dict(context, types, exprlist):
    """
    Unpacking tuple assignments in for statements and expr_stmts.
    """
    if exprlist.type == 'name':
        return {exprlist.value: types}
    elif exprlist.type == 'atom' and exprlist.children[0] in '([':
        return unpack_tuple_to_dict(context, types, exprlist.children[1])
    elif exprlist.type in ('testlist', 'testlist_comp', 'exprlist',
                           'testlist_star_expr'):
        dct = {}
        parts = iter(exprlist.children[::2])
        n = 0
        for lazy_context in types.iterate(exprlist):
            n += 1
            try:
                part = next(parts)
            except StopIteration:
                # TODO this context is probably not right.
                analysis.add(context, 'value-error-too-many-values', part,
                             message="ValueError: too many values to unpack (expected %s)" % n)
            else:
                dct.update(unpack_tuple_to_dict(context, lazy_context.infer(), part))
        has_parts = next(parts, None)
        if types and has_parts is not None:
            # TODO this context is probably not right.
            analysis.add(context, 'value-error-too-few-values', has_parts,
                         message="ValueError: need more than %s values to unpack" % n)
        return dct
    elif exprlist.type == 'power' or exprlist.type == 'atom_expr':
        # Something like ``arr[x], var = ...``.
        # This is something that is not yet supported, would also be difficult
        # to write into a dict.
        return {}
    elif exprlist.type == 'star_expr':  # `a, *b, c = x` type unpackings
        # Currently we're not supporting them.
        return {}
    raise NotImplementedError


def check_array_additions(context, sequence):
    """ Just a mapper function for the internal _check_array_additions """
    if sequence.array_type not in ('list', 'set'):
        # TODO also check for dict updates
        return NO_CONTEXTS

    return _check_array_additions(context, sequence)


@evaluator_method_cache(default=NO_CONTEXTS)
@debug.increase_indent
def _check_array_additions(context, sequence):
    """
    Checks if a `Array` has "add" (append, insert, extend) statements:

    >>> a = [""]
    >>> a.append(1)
    """
    from jedi.evaluate import arguments

    debug.dbg('Dynamic array search for %s' % sequence, color='MAGENTA')
    module_context = context.get_root_context()
    if not settings.dynamic_array_additions or isinstance(module_context, compiled.CompiledObject):
        debug.dbg('Dynamic array search aborted.', color='MAGENTA')
        return ContextSet()

    def find_additions(context, arglist, add_name):
        params = list(arguments.TreeArguments(context.evaluator, context, arglist).unpack())
        result = set()
        if add_name in ['insert']:
            params = params[1:]
        if add_name in ['append', 'add', 'insert']:
            for key, whatever in params:
                result.add(whatever)
        elif add_name in ['extend', 'update']:
            for key, lazy_context in params:
                result |= set(lazy_context.infer().iterate())
        return result

    temp_param_add, settings.dynamic_params_for_other_modules = \
        settings.dynamic_params_for_other_modules, False

    is_list = sequence.name.string_name == 'list'
    search_names = (['append', 'extend', 'insert'] if is_list else ['add', 'update'])

    added_types = set()
    for add_name in search_names:
        try:
            possible_names = module_context.tree_node.get_used_names()[add_name]
        except KeyError:
            continue
        else:
            for name in possible_names:
                context_node = context.tree_node
                if not (context_node.start_pos < name.start_pos < context_node.end_pos):
                    continue
                trailer = name.parent
                power = trailer.parent
                trailer_pos = power.children.index(trailer)
                try:
                    execution_trailer = power.children[trailer_pos + 1]
                except IndexError:
                    continue
                else:
                    if execution_trailer.type != 'trailer' \
                            or execution_trailer.children[0] != '(' \
                            or execution_trailer.children[1] == ')':
                        continue

                random_context = context.create_context(name)

                with recursion.execution_allowed(context.evaluator, power) as allowed:
                    if allowed:
                        found = evaluate_call_of_leaf(
                            random_context,
                            name,
                            cut_own_trailer=True
                        )
                        if sequence in found:
                            # The arrays match. Now add the results
                            added_types |= find_additions(
                                random_context,
                                execution_trailer.children[1],
                                add_name
                            )

    # reset settings
    settings.dynamic_params_for_other_modules = temp_param_add
    debug.dbg('Dynamic array result %s' % added_types, color='MAGENTA')
    return added_types


def get_dynamic_array_instance(instance):
    """Used for set() and list() instances."""
    if not settings.dynamic_array_additions:
        return instance.var_args

    ai = _ArrayInstance(instance)
    from jedi.evaluate import arguments
    return arguments.ValuesArguments([ContextSet(ai)])


class _ArrayInstance(object):
    """
    Used for the usage of set() and list().
    This is definitely a hack, but a good one :-)
    It makes it possible to use set/list conversions.

    In contrast to Array, ListComprehension and all other iterable types, this
    is something that is only used inside `evaluate/compiled/fake/builtins.py`
    and therefore doesn't need filters, `py__bool__` and so on, because
    we don't use these operations in `builtins.py`.
    """
    def __init__(self, instance):
        self.instance = instance
        self.var_args = instance.var_args

    def py__iter__(self):
        var_args = self.var_args
        try:
            _, lazy_context = next(var_args.unpack())
        except StopIteration:
            pass
        else:
            for lazy in lazy_context.infer().iterate():
                yield lazy

        from jedi.evaluate import arguments
        if isinstance(var_args, arguments.TreeArguments):
            additions = _check_array_additions(var_args.context, self.instance)
            for addition in additions:
                yield addition

    def iterate(self, contextualized_node=None, is_async=False):
        return self.py__iter__()


class Slice(Context):
    def __init__(self, context, start, stop, step):
        super(Slice, self).__init__(
            context.evaluator,
            parent_context=context.evaluator.builtins_module
        )
        self._context = context
        # all of them are either a Precedence or None.
        self._start = start
        self._stop = stop
        self._step = step

    @property
    def obj(self):
        """
        Imitate CompiledObject.obj behavior and return a ``builtin.slice()``
        object.
        """
        def get(element):
            if element is None:
                return None

            result = self._context.eval_node(element)
            if len(result) != 1:
                # For simplicity, we want slices to be clear defined with just
                # one type.  Otherwise we will return an empty slice object.
                raise IndexError

            context, = result
            return get_int_or_none(context)

        try:
            return slice(get(self._start), get(self._stop), get(self._step))
        except IndexError:
            return slice(None, None, None)
