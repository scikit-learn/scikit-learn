import functools
import inspect
import itertools
import operator
import toolz
from toolz.functoolz import (curry, is_valid_args, is_partial_args, is_arity,
                             num_required_args, has_varargs, has_keywords)
from toolz._signatures import builtins
import toolz._signatures as _sigs
from toolz.compatibility import PY3, PY33
from toolz.utils import raises


def make_func(param_string, raise_if_called=True):
    if not param_string.startswith('('):
        param_string = '(%s)' % param_string
    if raise_if_called:
        body = 'raise ValueError("function should not be called")'
    else:
        body = 'return True'
    d = {}
    exec('def func%s:\n    %s' % (param_string, body), globals(), d)
    return d['func']


def test_make_func():
    f = make_func('')
    assert raises(ValueError, lambda: f())
    assert raises(TypeError, lambda: f(1))

    f = make_func('', raise_if_called=False)
    assert f()
    assert raises(TypeError, lambda: f(1))

    f = make_func('x, y=1', raise_if_called=False)
    assert f(1)
    assert f(x=1)
    assert f(1, 2)
    assert f(x=1, y=2)
    assert raises(TypeError, lambda: f(1, 2, 3))

    f = make_func('(x, y=1)', raise_if_called=False)
    assert f(1)
    assert f(x=1)
    assert f(1, 2)
    assert f(x=1, y=2)
    assert raises(TypeError, lambda: f(1, 2, 3))


def test_is_valid(check_valid=is_valid_args, incomplete=False):
    orig_check_valid = check_valid
    check_valid = lambda func, *args, **kwargs: orig_check_valid(func, args, kwargs)

    f = make_func('')
    assert check_valid(f)
    assert check_valid(f, 1) is False
    assert check_valid(f, x=1) is False

    f = make_func('x')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1)
    assert check_valid(f, x=1)
    assert check_valid(f, 1, x=2) is False
    assert check_valid(f, 1, y=2) is False
    assert check_valid(f, 1, 2) is False
    assert check_valid(f, x=1, y=2) is False

    f = make_func('x=1')
    assert check_valid(f)
    assert check_valid(f, 1)
    assert check_valid(f, x=1)
    assert check_valid(f, 1, x=2) is False
    assert check_valid(f, 1, y=2) is False
    assert check_valid(f, 1, 2) is False
    assert check_valid(f, x=1, y=2) is False

    f = make_func('*args')
    assert check_valid(f)
    assert check_valid(f, 1)
    assert check_valid(f, 1, 2)
    assert check_valid(f, x=1) is False

    f = make_func('**kwargs')
    assert check_valid(f)
    assert check_valid(f, x=1)
    assert check_valid(f, x=1, y=2)
    assert check_valid(f, 1) is False

    f = make_func('x, *args')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1)
    assert check_valid(f, 1, 2)
    assert check_valid(f, x=1)
    assert check_valid(f, 1, x=1) is False
    assert check_valid(f, 1, y=1) is False

    f = make_func('x, y=1, **kwargs')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1)
    assert check_valid(f, x=1)
    assert check_valid(f, 1, 2)
    assert check_valid(f, x=1, y=2, z=3)
    assert check_valid(f, 1, 2, y=3) is False

    f = make_func('a, b, c=3, d=4')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1) is incomplete
    assert check_valid(f, 1, 2)
    assert check_valid(f, 1, c=3) is incomplete
    assert check_valid(f, 1, e=3) is False
    assert check_valid(f, 1, 2, e=3) is False
    assert check_valid(f, 1, 2, b=3) is False

    assert check_valid(1) is False


def test_is_valid_py3(check_valid=is_valid_args, incomplete=False):
    if not PY3:
        return
    orig_check_valid = check_valid
    check_valid = lambda func, *args, **kwargs: orig_check_valid(func, args, kwargs)

    f = make_func('x, *, y=1')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1)
    assert check_valid(f, x=1)
    assert check_valid(f, 1, y=2)
    assert check_valid(f, 1, 2) is False
    assert check_valid(f, 1, z=2) is False

    f = make_func('x, *args, y=1')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1)
    assert check_valid(f, x=1)
    assert check_valid(f, 1, y=2)
    assert check_valid(f, 1, 2, y=2)
    assert check_valid(f, 1, 2)
    assert check_valid(f, 1, z=2) is False

    f = make_func('*, y=1')
    assert check_valid(f)
    assert check_valid(f, 1) is False
    assert check_valid(f, y=1)
    assert check_valid(f, z=1) is False

    f = make_func('x, *, y')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1) is incomplete
    assert check_valid(f, x=1) is incomplete
    assert check_valid(f, 1, y=2)
    assert check_valid(f, x=1, y=2)
    assert check_valid(f, 1, 2) is False
    assert check_valid(f, 1, z=2) is False
    assert check_valid(f, 1, y=1, z=2) is False

    f = make_func('x=1, *, y, z=3')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1, z=3) is incomplete
    assert check_valid(f, y=2)
    assert check_valid(f, 1, y=2)
    assert check_valid(f, x=1, y=2)
    assert check_valid(f, x=1, y=2, z=3)
    assert check_valid(f, 1, x=1, y=2) is False
    assert check_valid(f, 1, 3, y=2) is False

    f = make_func('w, x=2, *args, y, z=4')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1) is incomplete
    assert check_valid(f, 1, y=3)

    f = make_func('a, b, c=3, d=4, *args, e=5, f=6, g, h')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1) is incomplete
    assert check_valid(f, 1, 2) is incomplete
    assert check_valid(f, 1, 2, g=7) is incomplete
    assert check_valid(f, 1, 2, g=7, h=8)
    assert check_valid(f, 1, 2, 3, 4, 5, 6, 7, 8, 9) is incomplete

    f = make_func('a: int, b: float')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1) is incomplete
    assert check_valid(f, b=1) is incomplete
    assert check_valid(f, 1, 2)

    f = make_func('(a: int, b: float) -> float')
    assert check_valid(f) is incomplete
    assert check_valid(f, 1) is incomplete
    assert check_valid(f, b=1) is incomplete
    assert check_valid(f, 1, 2)

    f.__signature__ = 34
    assert check_valid(f) is False

    class RaisesValueError(object):
        def __call__(self):
            pass
        @property
        def __signature__(self):
            raise ValueError('Testing Python 3.4')

    f = RaisesValueError()
    assert check_valid(f) is None


def test_is_partial():
    test_is_valid(check_valid=is_partial_args, incomplete=True)
    test_is_valid_py3(check_valid=is_partial_args, incomplete=True)


def test_is_valid_curry():
    def check_curry(func, args, kwargs, incomplete=True):
        try:
            curry(func)(*args, **kwargs)
            curry(func, *args)(**kwargs)
            curry(func, **kwargs)(*args)
            curry(func, *args, **kwargs)()
            if not isinstance(func, type(lambda: None)):
                return None
            return incomplete
        except ValueError:
            return True
        except TypeError:
            return False

    check_valid = functools.partial(check_curry, incomplete=True)
    test_is_valid(check_valid=check_valid, incomplete=True)
    test_is_valid_py3(check_valid=check_valid, incomplete=True)

    check_valid = functools.partial(check_curry, incomplete=False)
    test_is_valid(check_valid=check_valid, incomplete=False)
    test_is_valid_py3(check_valid=check_valid, incomplete=False)


def test_func_keyword():
    def f(func=None):
        pass
    assert is_valid_args(f, (), {})
    assert is_valid_args(f, (None,), {})
    assert is_valid_args(f, (), {'func': None})
    assert is_valid_args(f, (None,), {'func': None}) is False
    assert is_partial_args(f, (), {})
    assert is_partial_args(f, (None,), {})
    assert is_partial_args(f, (), {'func': None})
    assert is_partial_args(f, (None,), {'func': None}) is False


def test_has_unknown_args():
    assert has_varargs(1) is False
    assert has_varargs(map)
    assert has_varargs(make_func('')) is False
    assert has_varargs(make_func('x, y, z')) is False
    assert has_varargs(make_func('*args'))
    assert has_varargs(make_func('**kwargs')) is False
    assert has_varargs(make_func('x, y, *args, **kwargs'))
    assert has_varargs(make_func('x, y, z=1')) is False
    assert has_varargs(make_func('x, y, z=1, **kwargs')) is False

    if PY3:
        f = make_func('*args')
        f.__signature__ = 34
        assert has_varargs(f) is False

        class RaisesValueError(object):
            def __call__(self):
                pass
            @property
            def __signature__(self):
                raise ValueError('Testing Python 3.4')

        f = RaisesValueError()
        assert has_varargs(f) is None


def test_num_required_args():
    assert num_required_args(lambda: None) == 0
    assert num_required_args(lambda x: None) == 1
    assert num_required_args(lambda x, *args: None) == 1
    assert num_required_args(lambda x, **kwargs: None) == 1
    assert num_required_args(lambda x, y, *args, **kwargs: None) == 2
    assert num_required_args(map) == 2
    assert num_required_args(dict) is None


def test_has_keywords():
    assert has_keywords(lambda: None) is False
    assert has_keywords(lambda x: None) is False
    assert has_keywords(lambda x=1: None)
    assert has_keywords(lambda **kwargs: None)
    assert has_keywords(int)
    assert has_keywords(sorted)
    assert has_keywords(max)
    assert has_keywords(map) is False
    assert has_keywords(bytearray) is None


def test_has_varargs():
    assert has_varargs(lambda: None) is False
    assert has_varargs(lambda *args: None)
    assert has_varargs(lambda **kwargs: None) is False
    assert has_varargs(map)
    if PY3:
        assert has_varargs(max) is None


def test_is_arity():
    assert is_arity(0, lambda: None)
    assert is_arity(1, lambda: None) is False
    assert is_arity(1, lambda x: None)
    assert is_arity(3, lambda x, y, z: None)
    assert is_arity(1, lambda x, *args: None) is False
    assert is_arity(1, lambda x, **kwargs: None) is False
    assert is_arity(1, all)
    assert is_arity(2, map) is False
    assert is_arity(2, range) is None


def test_introspect_curry_valid_py3(check_valid=is_valid_args, incomplete=False):
    if not PY3:
        return
    orig_check_valid = check_valid
    check_valid = lambda _func, *args, **kwargs: orig_check_valid(_func, args, kwargs)

    f = toolz.curry(make_func('x, y, z=0'))
    assert check_valid(f)
    assert check_valid(f, 1)
    assert check_valid(f, 1, 2)
    assert check_valid(f, 1, 2, 3)
    assert check_valid(f, 1, 2, 3, 4) is False
    assert check_valid(f, invalid_keyword=True) is False
    assert check_valid(f(1))
    assert check_valid(f(1), 2)
    assert check_valid(f(1), 2, 3)
    assert check_valid(f(1), 2, 3, 4) is False
    assert check_valid(f(1), x=2) is False
    assert check_valid(f(1), y=2)
    assert check_valid(f(x=1), 2) is False
    assert check_valid(f(x=1), y=2)
    assert check_valid(f(y=2), 1)
    assert check_valid(f(y=2), 1, z=3)
    assert check_valid(f(y=2), 1, 3) is False

    f = toolz.curry(make_func('x, y, z=0'), 1, x=1)
    assert check_valid(f) is False
    assert check_valid(f, z=3) is False

    f = toolz.curry(make_func('x, y, *args, z'))
    assert check_valid(f)
    assert check_valid(f, 0)
    assert check_valid(f(1), 0)
    assert check_valid(f(1, 2), 0)
    assert check_valid(f(1, 2, 3), 0)
    assert check_valid(f(1, 2, 3, 4), 0)
    assert check_valid(f(1, 2, 3, 4), z=4)
    assert check_valid(f(x=1))
    assert check_valid(f(x=1), 1) is False
    assert check_valid(f(x=1), y=2)


def test_introspect_curry_partial_py3():
    test_introspect_curry_valid_py3(check_valid=is_partial_args, incomplete=True)


def test_introspect_curry_py3():
    if not PY3:
        return
    f = toolz.curry(make_func(''))
    assert num_required_args(f) == 0
    assert is_arity(0, f)
    assert has_varargs(f) is False
    assert has_keywords(f) is False

    f = toolz.curry(make_func('x'))
    assert num_required_args(f) == 0
    assert is_arity(0, f) is False
    assert is_arity(1, f) is False
    assert has_varargs(f) is False
    assert has_keywords(f)  # A side-effect of being curried

    f = toolz.curry(make_func('x, y, z=0'))
    assert num_required_args(f) == 0
    assert is_arity(0, f) is False
    assert is_arity(1, f) is False
    assert is_arity(2, f) is False
    assert is_arity(3, f) is False
    assert has_varargs(f) is False
    assert has_keywords(f)

    f = toolz.curry(make_func('*args, **kwargs'))
    assert num_required_args(f) == 0
    assert has_varargs(f)
    assert has_keywords(f)


def test_introspect_builtin_modules():
    mods = [builtins, functools, itertools, operator, toolz,
            toolz.functoolz, toolz.itertoolz, toolz.dicttoolz, toolz.recipes]

    blacklist = set()

    def add_blacklist(mod, attr):
        if hasattr(mod, attr):
            blacklist.add(getattr(mod, attr))

    add_blacklist(builtins, 'basestring')
    add_blacklist(builtins, 'breakpoint')
    add_blacklist(builtins, 'NoneType')
    add_blacklist(builtins, '__metaclass__')
    add_blacklist(builtins, 'sequenceiterator')

    def is_missing(modname, name, func):
        if name.startswith('_') and not name.startswith('__'):
            return False
        if name.startswith('__pyx_unpickle_') or name.endswith('_cython__'):
            return False
        try:
            if issubclass(func, BaseException):
                return False
        except TypeError:
            pass
        try:
            return (callable(func)
                    and func.__module__ is not None
                    and modname in func.__module__
                    and is_partial_args(func, (), {}) is not True
                    and func not in blacklist)
        except AttributeError:
            return False

    missing = {}
    for mod in mods:
        modname = mod.__name__
        for name, func in vars(mod).items():
            if is_missing(modname, name, func):
                if modname not in missing:
                    missing[modname] = []
                missing[modname].append(name)
    if missing:
        messages = []
        for modname, names in sorted(missing.items()):
            msg = '{0}:\n    {1}'.format(modname, '\n    '.join(sorted(names)))
            messages.append(msg)
        message = 'Missing introspection for the following callables:\n\n'
        raise AssertionError(message + '\n\n'.join(messages))


def test_inspect_signature_property():
    if not PY3:
        return

    # By adding AddX to our signature registry, we can inspect the class
    # itself and objects of the class.  `inspect.signature` doesn't like
    # it when `obj.__signature__` is a property.
    class AddX(object):
        def __init__(self, func):
            self.func = func

        def __call__(self, addx, *args, **kwargs):
            return addx + self.func(*args, **kwargs)

        @property
        def __signature__(self):
            sig = inspect.signature(self.func)
            params = list(sig.parameters.values())
            kind = inspect.Parameter.POSITIONAL_OR_KEYWORD
            newparam = inspect.Parameter('addx', kind)
            params = [newparam] + params
            return sig.replace(parameters=params)

    addx = AddX(lambda x: x)
    sig = inspect.signature(addx)
    assert sig == inspect.Signature(parameters=[
        inspect.Parameter('addx', inspect.Parameter.POSITIONAL_OR_KEYWORD),
        inspect.Parameter('x', inspect.Parameter.POSITIONAL_OR_KEYWORD)])

    assert num_required_args(AddX) is False
    _sigs.signatures[AddX] = (_sigs.expand_sig((0, lambda func: None)),)
    assert num_required_args(AddX) == 1
    del _sigs.signatures[AddX]


def test_inspect_wrapped_property():
    class Wrapped(object):
        def __init__(self, func):
            self.func = func

        def __call__(self, *args, **kwargs):
            return self.func(*args, **kwargs)

        @property
        def __wrapped__(self):
            return self.func

    func = lambda x: x
    wrapped = Wrapped(func)
    if PY3:
        assert inspect.signature(func) == inspect.signature(wrapped)

    assert num_required_args(Wrapped) == (False if PY33 else None)
    _sigs.signatures[Wrapped] = (_sigs.expand_sig((0, lambda func: None)),)
    assert num_required_args(Wrapped) == 1
