import inspect
import types
import sys
from textwrap import dedent
import operator as op
from collections import namedtuple

from jedi import debug
from jedi._compatibility import unicode, is_py3, is_py34, builtins, \
    py_version, force_unicode, print_to_stderr
from jedi.evaluate.compiled.getattr_static import getattr_static
from jedi.evaluate.utils import dotted_from_fs_path


MethodDescriptorType = type(str.replace)
# These are not considered classes and access is granted even though they have
# a __class__ attribute.
NOT_CLASS_TYPES = (
    types.BuiltinFunctionType,
    types.CodeType,
    types.FrameType,
    types.FunctionType,
    types.GeneratorType,
    types.GetSetDescriptorType,
    types.LambdaType,
    types.MemberDescriptorType,
    types.MethodType,
    types.ModuleType,
    types.TracebackType,
    MethodDescriptorType
)

if is_py3:
    NOT_CLASS_TYPES += (
        types.MappingProxyType,
        types.SimpleNamespace
    )
    if is_py34:
        NOT_CLASS_TYPES += (types.DynamicClassAttribute,)


# Those types don't exist in typing.
MethodDescriptorType = type(str.replace)
WrapperDescriptorType = type(set.__iter__)
# `object.__subclasshook__` is an already executed descriptor.
object_class_dict = type.__dict__["__dict__"].__get__(object)
ClassMethodDescriptorType = type(object_class_dict['__subclasshook__'])

def _a_generator(foo):
    """Used to have an object to return for generators."""
    yield 42
    yield foo


_sentinel = object()

# Maps Python syntax to the operator module.
COMPARISON_OPERATORS = {
    '==': op.eq,
    '!=': op.ne,
    'is': op.is_,
    'is not': op.is_not,
    '<': op.lt,
    '<=': op.le,
    '>': op.gt,
    '>=': op.ge,
}

_OPERATORS = {
    '+': op.add,
    '-': op.sub,
}
_OPERATORS.update(COMPARISON_OPERATORS)

ALLOWED_DESCRIPTOR_ACCESS = (
    types.FunctionType,
    types.GetSetDescriptorType,
    types.MemberDescriptorType,
    MethodDescriptorType,
    WrapperDescriptorType,
    ClassMethodDescriptorType,
    staticmethod,
    classmethod,
)


def safe_getattr(obj, name, default=_sentinel):
    try:
        attr, is_get_descriptor = getattr_static(obj, name)
    except AttributeError:
        if default is _sentinel:
            raise
        return default
    else:
        if type(attr) in ALLOWED_DESCRIPTOR_ACCESS:
            # In case of descriptors that have get methods we cannot return
            # it's value, because that would mean code execution.
            return getattr(obj, name)
    return attr


SignatureParam = namedtuple(
    'SignatureParam',
    'name has_default default has_annotation annotation kind_name'
)


def compiled_objects_cache(attribute_name):
    def decorator(func):
        """
        This decorator caches just the ids, oopposed to caching the object itself.
        Caching the id has the advantage that an object doesn't need to be
        hashable.
        """
        def wrapper(evaluator, obj, parent_context=None):
            cache = getattr(evaluator, attribute_name)
            # Do a very cheap form of caching here.
            key = id(obj)
            try:
                cache[key]
                return cache[key][0]
            except KeyError:
                # TODO wuaaaarrghhhhhhhh
                if attribute_name == 'mixed_cache':
                    result = func(evaluator, obj, parent_context)
                else:
                    result = func(evaluator, obj)
                # Need to cache all of them, otherwise the id could be overwritten.
                cache[key] = result, obj, parent_context
                return result
        return wrapper

    return decorator


def create_access(evaluator, obj):
    return evaluator.compiled_subprocess.get_or_create_access_handle(obj)


def load_module(evaluator, path=None, name=None, sys_path=None):
    if sys_path is None:
        sys_path = list(evaluator.get_sys_path())
    if path is not None:
        dotted_path = dotted_from_fs_path(path, sys_path=sys_path)
    else:
        dotted_path = name

    temp, sys.path = sys.path, sys_path
    try:
        __import__(dotted_path)
    except ImportError:
        # If a module is "corrupt" or not really a Python module or whatever.
        debug.warning('Module %s not importable in path %s.', dotted_path, path)
        return None
    except Exception:
        # Since __import__ pretty much makes code execution possible, just
        # catch any error here and print it.
        import traceback
        print_to_stderr("Cannot import:\n%s" % traceback.format_exc())
        return None
    finally:
        sys.path = temp

    # Just access the cache after import, because of #59 as well as the very
    # complicated import structure of Python.
    module = sys.modules[dotted_path]
    return create_access_path(evaluator, module)


class AccessPath(object):
    def __init__(self, accesses):
        self.accesses = accesses

    # Writing both of these methods here looks a bit ridiculous. However with
    # the differences of Python 2/3 it's actually necessary, because we will
    # otherwise have a accesses attribute that is bytes instead of unicode.
    def __getstate__(self):
        return self.accesses

    def __setstate__(self, value):
        self.accesses = value


def create_access_path(evaluator, obj):
    access = create_access(evaluator, obj)
    return AccessPath(access.get_access_path_tuples())


def _force_unicode_decorator(func):
    return lambda *args, **kwargs: force_unicode(func(*args, **kwargs))


class DirectObjectAccess(object):
    def __init__(self, evaluator, obj):
        self._evaluator = evaluator
        self._obj = obj

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__, self.get_repr())

    def _create_access(self, obj):
        return create_access(self._evaluator, obj)

    def _create_access_path(self, obj):
        return create_access_path(self._evaluator, obj)

    def py__bool__(self):
        return bool(self._obj)

    def py__file__(self):
        try:
            return self._obj.__file__
        except AttributeError:
            return None

    def py__doc__(self, include_call_signature=False):
        return force_unicode(inspect.getdoc(self._obj)) or u''

    def py__name__(self):
        if not _is_class_instance(self._obj) or \
                inspect.ismethoddescriptor(self._obj):  # slots
            cls = self._obj
        else:
            try:
                cls = self._obj.__class__
            except AttributeError:
                # happens with numpy.core.umath._UFUNC_API (you get it
                # automatically by doing `import numpy`.
                return None

        try:
            return force_unicode(cls.__name__)
        except AttributeError:
            return None

    def py__mro__accesses(self):
        return tuple(self._create_access_path(cls) for cls in self._obj.__mro__[1:])

    def py__getitem__(self, index):
        if type(self._obj) not in (str, list, tuple, unicode, bytes, bytearray, dict):
            # Get rid of side effects, we won't call custom `__getitem__`s.
            return None

        return self._create_access_path(self._obj[index])

    def py__iter__list(self):
        if type(self._obj) not in (str, list, tuple, unicode, bytes, bytearray, dict):
            # Get rid of side effects, we won't call custom `__getitem__`s.
            return []

        lst = []
        for i, part in enumerate(self._obj):
            if i > 20:
                # Should not go crazy with large iterators
                break
            lst.append(self._create_access_path(part))
        return lst

    def py__class__(self):
        return self._create_access_path(self._obj.__class__)

    def py__bases__(self):
        return [self._create_access_path(base) for base in self._obj.__bases__]

    @_force_unicode_decorator
    def get_repr(self):
        builtins = 'builtins', '__builtin__'

        if inspect.ismodule(self._obj):
            return repr(self._obj)
        # Try to avoid execution of the property.
        if safe_getattr(self._obj, '__module__', default='') in builtins:
            return repr(self._obj)

        type_ = type(self._obj)
        if type_ == type:
            return type.__repr__(self._obj)

        if safe_getattr(type_, '__module__', default='') in builtins:
            # Allow direct execution of repr for builtins.
            return repr(self._obj)
        return object.__repr__(self._obj)

    def is_class(self):
        return inspect.isclass(self._obj)

    def ismethoddescriptor(self):
        return inspect.ismethoddescriptor(self._obj)

    def dir(self):
        return list(map(force_unicode, dir(self._obj)))

    def has_iter(self):
        try:
            iter(self._obj)
            return True
        except TypeError:
            return False

    def is_allowed_getattr(self, name):
        # TODO this API is ugly.
        try:
            attr, is_get_descriptor = getattr_static(self._obj, name)
        except AttributeError:
            return False, False
        else:
            if is_get_descriptor and type(attr) not in ALLOWED_DESCRIPTOR_ACCESS:
                # In case of descriptors that have get methods we cannot return
                # it's value, because that would mean code execution.
                return True, True
        return True, False

    def getattr(self, name, default=_sentinel):
        try:
            return self._create_access(getattr(self._obj, name))
        except AttributeError:
            # Happens e.g. in properties of
            # PyQt4.QtGui.QStyleOptionComboBox.currentText
            # -> just set it to None
            if default is _sentinel:
                raise
            return self._create_access(default)

    def get_safe_value(self):
        if type(self._obj) in (bool, bytes, float, int, str, unicode, slice):
            return self._obj
        raise ValueError("Object is type %s and not simple" % type(self._obj))

    def get_api_type(self):
        obj = self._obj
        if self.is_class():
            return u'class'
        elif inspect.ismodule(obj):
            return u'module'
        elif inspect.isbuiltin(obj) or inspect.ismethod(obj) \
                or inspect.ismethoddescriptor(obj) or inspect.isfunction(obj):
            return u'function'
        # Everything else...
        return u'instance'

    def get_access_path_tuples(self):
        accesses = [create_access(self._evaluator, o) for o in self._get_objects_path()]
        return [(access.py__name__(), access) for access in accesses]

    def _get_objects_path(self):
        def get():
            obj = self._obj
            yield obj
            try:
                obj = obj.__objclass__
            except AttributeError:
                pass
            else:
                yield obj

            try:
                # Returns a dotted string path.
                imp_plz = obj.__module__
            except AttributeError:
                # Unfortunately in some cases like `int` there's no __module__
                if not inspect.ismodule(obj):
                    yield builtins
            else:
                if imp_plz is None:
                    # Happens for example in `(_ for _ in []).send.__module__`.
                    yield builtins
                else:
                    try:
                        # TODO use sys.modules, __module__ can be faked.
                        yield sys.modules[imp_plz]
                    except KeyError:
                        # __module__ can be something arbitrary that doesn't exist.
                        yield builtins

        return list(reversed(list(get())))

    def execute_operation(self, other_access_handle, operator):
        other_access = other_access_handle.access
        op = _OPERATORS[operator]
        return self._create_access_path(op(self._obj, other_access._obj))

    def needs_type_completions(self):
        return inspect.isclass(self._obj) and self._obj != type

    def get_signature_params(self):
        obj = self._obj
        if py_version < 33:
            raise ValueError("inspect.signature was introduced in 3.3")
        if py_version == 34:
            # In 3.4 inspect.signature are wrong for str and int. This has
            # been fixed in 3.5. The signature of object is returned,
            # because no signature was found for str. Here we imitate 3.5
            # logic and just ignore the signature if the magic methods
            # don't match object.
            # 3.3 doesn't even have the logic and returns nothing for str
            # and classes that inherit from object.
            user_def = inspect._signature_get_user_defined_method
            if (inspect.isclass(obj)
                    and not user_def(type(obj), '__init__')
                    and not user_def(type(obj), '__new__')
                    and (obj.__init__ != object.__init__
                         or obj.__new__ != object.__new__)):
                raise ValueError

        try:
            signature = inspect.signature(obj)
        except (RuntimeError, TypeError):
            # Reading the code of the function in Python 3.6 implies there are
            # at least these errors that might occur if something is wrong with
            # the signature. In that case we just want a simple escape for now.
            raise ValueError
        return [
            SignatureParam(
                name=p.name,
                has_default=p.default is not p.empty,
                default=self._create_access_path(p.default),
                has_annotation=p.annotation is not p.empty,
                annotation=self._create_access_path(p.annotation),
                kind_name=str(p.kind)
            ) for p in signature.parameters.values()
        ]

    def negate(self):
        return self._create_access_path(-self._obj)

    def dict_values(self):
        return [self._create_access_path(v) for v in self._obj.values()]

    def is_super_class(self, exception):
        return issubclass(exception, self._obj)

    def get_dir_infos(self):
        """
        Used to return a couple of infos that are needed when accessing the sub
        objects of an objects
        """
        # TODO is_allowed_getattr might raise an AttributeError
        tuples = dict(
            (force_unicode(name), self.is_allowed_getattr(name))
            for name in self.dir()
        )
        return self.needs_type_completions(), tuples


def _is_class_instance(obj):
    """Like inspect.* methods."""
    try:
        cls = obj.__class__
    except AttributeError:
        return False
    else:
        return cls != type and not issubclass(cls, NOT_CLASS_TYPES)


if py_version >= 35:
    exec(compile(dedent("""
        async def _coroutine(): pass
        _coroutine = _coroutine()
        CoroutineType = type(_coroutine)
        _coroutine.close()  # Prevent ResourceWarning
    """), 'blub', 'exec'))
    _coroutine_wrapper = _coroutine.__await__()
else:
    _coroutine = None
    _coroutine_wrapper = None

if py_version >= 36:
    exec(compile(dedent("""
        async def _async_generator():
            yield
        _async_generator = _async_generator()
        AsyncGeneratorType = type(_async_generator)
    """), 'blub', 'exec'))
else:
    _async_generator = None

class _SPECIAL_OBJECTS(object):
    FUNCTION_CLASS = types.FunctionType
    METHOD_CLASS = type(DirectObjectAccess.py__bool__)
    MODULE_CLASS = types.ModuleType
    GENERATOR_OBJECT = _a_generator(1.0)
    BUILTINS = builtins
    COROUTINE = _coroutine
    COROUTINE_WRAPPER = _coroutine_wrapper
    ASYNC_GENERATOR = _async_generator


def get_special_object(evaluator, identifier):
    obj = getattr(_SPECIAL_OBJECTS, identifier)
    return create_access_path(evaluator, obj)
