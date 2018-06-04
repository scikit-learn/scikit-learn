"""
Implementations of standard library functions, because it's not possible to
understand them with Jedi.

To add a new implementation, create a function and add it to the
``_implemented`` dict at the bottom of this module.

Note that this module exists only to implement very specific functionality in
the standard library. The usual way to understand the standard library is the
compiled module that returns the types for C-builtins.
"""
import re

import parso

from jedi._compatibility import force_unicode
from jedi import debug
from jedi.evaluate.arguments import ValuesArguments
from jedi.evaluate import analysis
from jedi.evaluate import compiled
from jedi.evaluate.context.instance import InstanceFunctionExecution, \
    AbstractInstanceContext, CompiledInstance, BoundMethod, \
    AnonymousInstanceFunctionExecution
from jedi.evaluate.base_context import ContextualizedNode, \
    NO_CONTEXTS, ContextSet
from jedi.evaluate.context import ClassContext, ModuleContext
from jedi.evaluate.context import iterable
from jedi.evaluate.lazy_context import LazyTreeContext
from jedi.evaluate.syntax_tree import is_string

# Now this is all part of fake tuples in Jedi. However super doesn't work on
# __init__ and __new__ doesn't work at all. So adding this to nametuples is
# just the easiest way.
_NAMEDTUPLE_INIT = """
    def __init__(_cls, {arg_list}):
        'A helper function for namedtuple.'
        self.__iterable = ({arg_list})

    def __iter__(self):
        for i in self.__iterable:
            yield i

    def __getitem__(self, y):
        return self.__iterable[y]

"""


class NotInStdLib(LookupError):
    pass


def execute(evaluator, obj, arguments):
    if isinstance(obj, BoundMethod):
        raise NotInStdLib()

    try:
        obj_name = obj.name.string_name
    except AttributeError:
        pass
    else:
        if obj.parent_context == evaluator.builtins_module:
            module_name = 'builtins'
        elif isinstance(obj.parent_context, ModuleContext):
            module_name = obj.parent_context.name.string_name
        else:
            module_name = ''

        # for now we just support builtin functions.
        try:
            func = _implemented[module_name][obj_name]
        except KeyError:
            pass
        else:
            return func(evaluator, obj, arguments)
    raise NotInStdLib()


def _follow_param(evaluator, arguments, index):
    try:
        key, lazy_context = list(arguments.unpack())[index]
    except IndexError:
        return NO_CONTEXTS
    else:
        return lazy_context.infer()


def argument_clinic(string, want_obj=False, want_context=False, want_arguments=False):
    """
    Works like Argument Clinic (PEP 436), to validate function params.
    """
    clinic_args = []
    allow_kwargs = False
    optional = False
    while string:
        # Optional arguments have to begin with a bracket. And should always be
        # at the end of the arguments. This is therefore not a proper argument
        # clinic implementation. `range()` for exmple allows an optional start
        # value at the beginning.
        match = re.match('(?:(?:(\[),? ?|, ?|)(\w+)|, ?/)\]*', string)
        string = string[len(match.group(0)):]
        if not match.group(2):  # A slash -> allow named arguments
            allow_kwargs = True
            continue
        optional = optional or bool(match.group(1))
        word = match.group(2)
        clinic_args.append((word, optional, allow_kwargs))

    def f(func):
        def wrapper(evaluator, obj, arguments):
            debug.dbg('builtin start %s' % obj, color='MAGENTA')
            result = NO_CONTEXTS
            try:
                lst = list(arguments.eval_argument_clinic(clinic_args))
            except ValueError:
                pass
            else:
                kwargs = {}
                if want_context:
                    kwargs['context'] = arguments.context
                if want_obj:
                    kwargs['obj'] = obj
                if want_arguments:
                    kwargs['arguments'] = arguments
                result = func(evaluator, *lst, **kwargs)
            finally:
                debug.dbg('builtin end: %s', result, color='MAGENTA')
            return result

        return wrapper
    return f


@argument_clinic('iterator[, default], /')
def builtins_next(evaluator, iterators, defaults):
    """
    TODO this function is currently not used. It's a stab at implementing next
    in a different way than fake objects. This would be a bit more flexible.
    """
    if evaluator.environment.version_info.major == 2:
        name = 'next'
    else:
        name = '__next__'

    context_set = NO_CONTEXTS
    for iterator in iterators:
        if isinstance(iterator, AbstractInstanceContext):
            context_set = ContextSet.from_sets(
                n.infer()
                for filter in iterator.get_filters(include_self_names=True)
                for n in filter.get(name)
            ).execute_evaluated()
    if context_set:
        return context_set
    return defaults


@argument_clinic('object, name[, default], /')
def builtins_getattr(evaluator, objects, names, defaults=None):
    # follow the first param
    for obj in objects:
        for name in names:
            if is_string(name):
                return obj.py__getattribute__(force_unicode(name.get_safe_value()))
            else:
                debug.warning('getattr called without str')
                continue
    return NO_CONTEXTS


@argument_clinic('object[, bases, dict], /')
def builtins_type(evaluator, objects, bases, dicts):
    if bases or dicts:
        # It's a type creation... maybe someday...
        return NO_CONTEXTS
    else:
        return objects.py__class__()


class SuperInstance(AbstractInstanceContext):
    """To be used like the object ``super`` returns."""
    def __init__(self, evaluator, cls):
        su = cls.py_mro()[1]
        super().__init__(evaluator, su and su[0] or self)


@argument_clinic('[type[, obj]], /', want_context=True)
def builtins_super(evaluator, types, objects, context):
    # TODO make this able to detect multiple inheritance super
    if isinstance(context, (InstanceFunctionExecution,
                            AnonymousInstanceFunctionExecution)):
        su = context.instance.py__class__().py__bases__()
        return su[0].infer().execute_evaluated()
    return NO_CONTEXTS


@argument_clinic('sequence, /', want_obj=True, want_arguments=True)
def builtins_reversed(evaluator, sequences, obj, arguments):
    # While we could do without this variable (just by using sequences), we
    # want static analysis to work well. Therefore we need to generated the
    # values again.
    key, lazy_context = next(arguments.unpack())
    cn = None
    if isinstance(lazy_context, LazyTreeContext):
        # TODO access private
        cn = ContextualizedNode(lazy_context._context, lazy_context.data)
    ordered = list(sequences.iterate(cn))

    rev = list(reversed(ordered))
    # Repack iterator values and then run it the normal way. This is
    # necessary, because `reversed` is a function and autocompletion
    # would fail in certain cases like `reversed(x).__iter__` if we
    # just returned the result directly.
    seq = iterable.FakeSequence(evaluator, u'list', rev)
    arguments = ValuesArguments([ContextSet(seq)])
    return ContextSet(CompiledInstance(evaluator, evaluator.builtins_module, obj, arguments))


@argument_clinic('obj, type, /', want_arguments=True)
def builtins_isinstance(evaluator, objects, types, arguments):
    bool_results = set()
    for o in objects:
        cls = o.py__class__()
        try:
            mro_func = cls.py__mro__
        except AttributeError:
            # This is temporary. Everything should have a class attribute in
            # Python?! Maybe we'll leave it here, because some numpy objects or
            # whatever might not.
            bool_results = set([True, False])
            break

        mro = mro_func()

        for cls_or_tup in types:
            if cls_or_tup.is_class():
                bool_results.add(cls_or_tup in mro)
            elif cls_or_tup.name.string_name == 'tuple' \
                    and cls_or_tup.get_root_context() == evaluator.builtins_module:
                # Check for tuples.
                classes = ContextSet.from_sets(
                    lazy_context.infer()
                    for lazy_context in cls_or_tup.iterate()
                )
                bool_results.add(any(cls in mro for cls in classes))
            else:
                _, lazy_context = list(arguments.unpack())[1]
                if isinstance(lazy_context, LazyTreeContext):
                    node = lazy_context.data
                    message = 'TypeError: isinstance() arg 2 must be a ' \
                              'class, type, or tuple of classes and types, ' \
                              'not %s.' % cls_or_tup
                    analysis.add(lazy_context._context, 'type-error-isinstance', node, message)

    return ContextSet.from_iterable(
        compiled.builtin_from_name(evaluator, force_unicode(str(b)))
        for b in bool_results
    )


def collections_namedtuple(evaluator, obj, arguments):
    """
    Implementation of the namedtuple function.

    This has to be done by processing the namedtuple class template and
    evaluating the result.

    """
    collections_context = obj.parent_context
    _class_template_set = collections_context.py__getattribute__(u'_class_template')
    if not _class_template_set:
        # Namedtuples are not supported on Python 2.6, early 2.7, because the
        # _class_template variable is not defined, there.
        return NO_CONTEXTS

    # Process arguments
    # TODO here we only use one of the types, we should use all.
    # TODO this is buggy, doesn't need to be a string
    name = list(_follow_param(evaluator, arguments, 0))[0].get_safe_value()
    _fields = list(_follow_param(evaluator, arguments, 1))[0]
    if isinstance(_fields, compiled.CompiledObject):
        fields = _fields.get_safe_value().replace(',', ' ').split()
    elif isinstance(_fields, iterable.Sequence):
        fields = [
            v.get_safe_value()
            for lazy_context in _fields.py__iter__()
            for v in lazy_context.infer() if is_string(v)
        ]
    else:
        return NO_CONTEXTS

    def get_var(name):
        x, = collections_context.py__getattribute__(name)
        return x.get_safe_value()

    base = next(iter(_class_template_set)).get_safe_value()
    base += _NAMEDTUPLE_INIT
    # Build source code
    code = base.format(
        typename=name,
        field_names=tuple(fields),
        num_fields=len(fields),
        arg_list=repr(tuple(fields)).replace("u'", "").replace("'", "")[1:-1],
        repr_fmt=', '.join(get_var(u'_repr_template').format(name=name) for name in fields),
        field_defs='\n'.join(get_var(u'_field_template').format(index=index, name=name)
                             for index, name in enumerate(fields))
    )

    # Parse source code
    module = evaluator.grammar.parse(code)
    generated_class = next(module.iter_classdefs())
    parent_context = ModuleContext(
        evaluator, module, None,
        code_lines=parso.split_lines(code, keepends=True),
    )
    return ContextSet(ClassContext(evaluator, parent_context, generated_class))


@argument_clinic('first, /')
def _return_first_param(evaluator, firsts):
    return firsts


_implemented = {
    'builtins': {
        'getattr': builtins_getattr,
        'type': builtins_type,
        'super': builtins_super,
        'reversed': builtins_reversed,
        'isinstance': builtins_isinstance,
    },
    'copy': {
        'copy': _return_first_param,
        'deepcopy': _return_first_param,
    },
    'json': {
        'load': lambda *args: NO_CONTEXTS,
        'loads': lambda *args: NO_CONTEXTS,
    },
    'collections': {
        'namedtuple': collections_namedtuple,
    },
}
