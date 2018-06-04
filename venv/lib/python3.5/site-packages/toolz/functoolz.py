from functools import reduce, partial
import inspect
import operator
from operator import attrgetter
from textwrap import dedent

from .compatibility import PY3, PY33, PY34, PYPY, import_module
from .utils import no_default


__all__ = ('identity', 'thread_first', 'thread_last', 'memoize', 'compose',
           'pipe', 'complement', 'juxt', 'do', 'curry', 'flip', 'excepts')


def identity(x):
    """ Identity function. Return x

    >>> identity(3)
    3
    """
    return x


def thread_first(val, *forms):
    """ Thread value through a sequence of functions/forms

    >>> def double(x): return 2*x
    >>> def inc(x):    return x + 1
    >>> thread_first(1, inc, double)
    4

    If the function expects more than one input you can specify those inputs
    in a tuple.  The value is used as the first input.

    >>> def add(x, y): return x + y
    >>> def pow(x, y): return x**y
    >>> thread_first(1, (add, 4), (pow, 2))  # pow(add(1, 4), 2)
    25

    So in general
        thread_first(x, f, (g, y, z))
    expands to
        g(f(x), y, z)

    See Also:
        thread_last
    """
    def evalform_front(val, form):
        if callable(form):
            return form(val)
        if isinstance(form, tuple):
            func, args = form[0], form[1:]
            args = (val,) + args
            return func(*args)
    return reduce(evalform_front, forms, val)


def thread_last(val, *forms):
    """ Thread value through a sequence of functions/forms

    >>> def double(x): return 2*x
    >>> def inc(x):    return x + 1
    >>> thread_last(1, inc, double)
    4

    If the function expects more than one input you can specify those inputs
    in a tuple.  The value is used as the last input.

    >>> def add(x, y): return x + y
    >>> def pow(x, y): return x**y
    >>> thread_last(1, (add, 4), (pow, 2))  # pow(2, add(4, 1))
    32

    So in general
        thread_last(x, f, (g, y, z))
    expands to
        g(y, z, f(x))

    >>> def iseven(x):
    ...     return x % 2 == 0
    >>> list(thread_last([1, 2, 3], (map, inc), (filter, iseven)))
    [2, 4]

    See Also:
        thread_first
    """
    def evalform_back(val, form):
        if callable(form):
            return form(val)
        if isinstance(form, tuple):
            func, args = form[0], form[1:]
            args = args + (val,)
            return func(*args)
    return reduce(evalform_back, forms, val)


def instanceproperty(fget=None, fset=None, fdel=None, doc=None, classval=None):
    """ Like @property, but returns ``classval`` when used as a class attribute

    >>> class MyClass(object):
    ...     '''The class docstring'''
    ...     @instanceproperty(classval=__doc__)
    ...     def __doc__(self):
    ...         return 'An object docstring'
    ...     @instanceproperty
    ...     def val(self):
    ...         return 42
    ...
    >>> MyClass.__doc__
    'The class docstring'
    >>> MyClass.val is None
    True
    >>> obj = MyClass()
    >>> obj.__doc__
    'An object docstring'
    >>> obj.val
    42
    """
    if fget is None:
        return partial(instanceproperty, fset=fset, fdel=fdel, doc=doc,
                       classval=classval)
    return InstanceProperty(fget=fget, fset=fset, fdel=fdel, doc=doc,
                            classval=classval)


class InstanceProperty(property):
    """ Like @property, but returns ``classval`` when used as a class attribute

    Should not be used directly.  Use ``instanceproperty`` instead.
    """
    def __init__(self, fget=None, fset=None, fdel=None, doc=None,
                 classval=None):
        self.classval = classval
        property.__init__(self, fget=fget, fset=fset, fdel=fdel, doc=doc)

    def __get__(self, obj, type=None):
        if obj is None:
            return self.classval
        return property.__get__(self, obj, type)

    def __reduce__(self):
        state = (self.fget, self.fset, self.fdel, self.__doc__, self.classval)
        return InstanceProperty, state


class curry(object):
    """ Curry a callable function

    Enables partial application of arguments through calling a function with an
    incomplete set of arguments.

    >>> def mul(x, y):
    ...     return x * y
    >>> mul = curry(mul)

    >>> double = mul(2)
    >>> double(10)
    20

    Also supports keyword arguments

    >>> @curry                  # Can use curry as a decorator
    ... def f(x, y, a=10):
    ...     return a * (x + y)

    >>> add = f(a=1)
    >>> add(2, 3)
    5

    See Also:
        toolz.curried - namespace of curried functions
                        https://toolz.readthedocs.io/en/latest/curry.html
    """
    def __init__(self, *args, **kwargs):
        if not args:
            raise TypeError('__init__() takes at least 2 arguments (1 given)')
        func, args = args[0], args[1:]
        if not callable(func):
            raise TypeError("Input must be callable")

        # curry- or functools.partial-like object?  Unpack and merge arguments
        if (
            hasattr(func, 'func')
            and hasattr(func, 'args')
            and hasattr(func, 'keywords')
            and isinstance(func.args, tuple)
        ):
            _kwargs = {}
            if func.keywords:
                _kwargs.update(func.keywords)
            _kwargs.update(kwargs)
            kwargs = _kwargs
            args = func.args + args
            func = func.func

        if kwargs:
            self._partial = partial(func, *args, **kwargs)
        else:
            self._partial = partial(func, *args)

        self.__doc__ = getattr(func, '__doc__', None)
        self.__name__ = getattr(func, '__name__', '<curry>')
        self.__module__ = getattr(func, '__module__', None)
        self.__qualname__ = getattr(func, '__qualname__', None)
        self._sigspec = None
        self._has_unknown_args = None

    @instanceproperty
    def func(self):
        return self._partial.func

    if PY3:  # pragma: py2 no cover
        @instanceproperty
        def __signature__(self):
            sig = inspect.signature(self.func)
            args = self.args or ()
            keywords = self.keywords or {}
            if is_partial_args(self.func, args, keywords, sig) is False:
                raise TypeError('curry object has incorrect arguments')

            params = list(sig.parameters.values())
            skip = 0
            for param in params[:len(args)]:
                if param.kind == param.VAR_POSITIONAL:
                    break
                skip += 1

            kwonly = False
            newparams = []
            for param in params[skip:]:
                kind = param.kind
                default = param.default
                if kind == param.VAR_KEYWORD:
                    pass
                elif kind == param.VAR_POSITIONAL:
                    if kwonly:
                        continue
                elif param.name in keywords:
                    default = keywords[param.name]
                    kind = param.KEYWORD_ONLY
                    kwonly = True
                else:
                    if kwonly:
                        kind = param.KEYWORD_ONLY
                    if default is param.empty:
                        default = no_default
                newparams.append(param.replace(default=default, kind=kind))

            return sig.replace(parameters=newparams)

    @instanceproperty
    def args(self):
        return self._partial.args

    @instanceproperty
    def keywords(self):
        return self._partial.keywords

    @instanceproperty
    def func_name(self):
        return self.__name__

    def __str__(self):
        return str(self.func)

    def __repr__(self):
        return repr(self.func)

    def __hash__(self):
        return hash((self.func, self.args,
                     frozenset(self.keywords.items()) if self.keywords
                     else None))

    def __eq__(self, other):
        return (isinstance(other, curry) and self.func == other.func and
                self.args == other.args and self.keywords == other.keywords)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __call__(self, *args, **kwargs):
        try:
            return self._partial(*args, **kwargs)
        except TypeError as exc:
            if self._should_curry(args, kwargs, exc):
                return self.bind(*args, **kwargs)
            raise

    def _should_curry(self, args, kwargs, exc=None):
        func = self.func
        args = self.args + args
        if self.keywords:
            kwargs = dict(self.keywords, **kwargs)
        if self._sigspec is None:
            sigspec = self._sigspec = _sigs.signature_or_spec(func)
            self._has_unknown_args = has_varargs(func, sigspec) is not False
        else:
            sigspec = self._sigspec

        if is_partial_args(func, args, kwargs, sigspec) is False:
            # Nothing can make the call valid
            return False
        elif self._has_unknown_args:
            # The call may be valid and raised a TypeError, but we curry
            # anyway because the function may have `*args`.  This is useful
            # for decorators with signature `func(*args, **kwargs)`.
            return True
        elif not is_valid_args(func, args, kwargs, sigspec):
            # Adding more arguments may make the call valid
            return True
        else:
            # There was a genuine TypeError
            return False

    def bind(self, *args, **kwargs):
        return type(self)(self, *args, **kwargs)

    def call(self, *args, **kwargs):
        return self._partial(*args, **kwargs)

    def __get__(self, instance, owner):
        if instance is None:
            return self
        return curry(self, instance)

    def __reduce__(self):
        func = self.func
        modname = getattr(func, '__module__', None)
        qualname = getattr(func, '__qualname__', None)
        if qualname is None:  # pragma: py3 no cover
            qualname = getattr(func, '__name__', None)
        is_decorated = None
        if modname and qualname:
            attrs = []
            obj = import_module(modname)
            for attr in qualname.split('.'):
                if isinstance(obj, curry):  # pragma: py2 no cover
                    attrs.append('func')
                    obj = obj.func
                obj = getattr(obj, attr, None)
                if obj is None:
                    break
                attrs.append(attr)
            if isinstance(obj, curry) and obj.func is func:
                is_decorated = obj is self
                qualname = '.'.join(attrs)
                func = '%s:%s' % (modname, qualname)

        # functools.partial objects can't be pickled
        userdict = tuple((k, v) for k, v in self.__dict__.items()
                         if k not in ('_partial', '_sigspec'))
        state = (type(self), func, self.args, self.keywords, userdict,
                 is_decorated)
        return (_restore_curry, state)


def _restore_curry(cls, func, args, kwargs, userdict, is_decorated):
    if isinstance(func, str):
        modname, qualname = func.rsplit(':', 1)
        obj = import_module(modname)
        for attr in qualname.split('.'):
            obj = getattr(obj, attr)
        if is_decorated:
            return obj
        func = obj.func
    obj = cls(func, *args, **(kwargs or {}))
    obj.__dict__.update(userdict)
    return obj


@curry
def memoize(func, cache=None, key=None):
    """ Cache a function's result for speedy future evaluation

    Considerations:
        Trades memory for speed.
        Only use on pure functions.

    >>> def add(x, y):  return x + y
    >>> add = memoize(add)

    Or use as a decorator

    >>> @memoize
    ... def add(x, y):
    ...     return x + y

    Use the ``cache`` keyword to provide a dict-like object as an initial cache

    >>> @memoize(cache={(1, 2): 3})
    ... def add(x, y):
    ...     return x + y

    Note that the above works as a decorator because ``memoize`` is curried.

    It is also possible to provide a ``key(args, kwargs)`` function that
    calculates keys used for the cache, which receives an ``args`` tuple and
    ``kwargs`` dict as input, and must return a hashable value.  However,
    the default key function should be sufficient most of the time.

    >>> # Use key function that ignores extraneous keyword arguments
    >>> @memoize(key=lambda args, kwargs: args)
    ... def add(x, y, verbose=False):
    ...     if verbose:
    ...         print('Calculating %s + %s' % (x, y))
    ...     return x + y
    """
    if cache is None:
        cache = {}

    try:
        may_have_kwargs = has_keywords(func) is not False
        # Is unary function (single arg, no variadic argument or keywords)?
        is_unary = is_arity(1, func)
    except TypeError:  # pragma: no cover
        may_have_kwargs = True
        is_unary = False

    if key is None:
        if is_unary:
            def key(args, kwargs):
                return args[0]
        elif may_have_kwargs:
            def key(args, kwargs):
                return (
                    args or None,
                    frozenset(kwargs.items()) if kwargs else None,
                )
        else:
            def key(args, kwargs):
                return args

    def memof(*args, **kwargs):
        k = key(args, kwargs)
        try:
            return cache[k]
        except TypeError:
            raise TypeError("Arguments to memoized function must be hashable")
        except KeyError:
            cache[k] = result = func(*args, **kwargs)
            return result

    try:
        memof.__name__ = func.__name__
    except AttributeError:
        pass
    memof.__doc__ = func.__doc__
    memof.__wrapped__ = func
    return memof


class Compose(object):
    """ A composition of functions

    See Also:
        compose
    """
    __slots__ = 'first', 'funcs'

    def __init__(self, funcs):
        funcs = tuple(reversed(funcs))
        self.first = funcs[0]
        self.funcs = funcs[1:]

    def __call__(self, *args, **kwargs):
        ret = self.first(*args, **kwargs)
        for f in self.funcs:
            ret = f(ret)
        return ret

    def __getstate__(self):
        return self.first, self.funcs

    def __setstate__(self, state):
        self.first, self.funcs = state

    @instanceproperty(classval=__doc__)
    def __doc__(self):
        def composed_doc(*fs):
            """Generate a docstring for the composition of fs.
            """
            if not fs:
                # Argument name for the docstring.
                return '*args, **kwargs'

            return '{f}({g})'.format(f=fs[0].__name__, g=composed_doc(*fs[1:]))

        try:
            return (
                'lambda *args, **kwargs: ' +
                composed_doc(*reversed((self.first,) + self.funcs))
            )
        except AttributeError:
            # One of our callables does not have a `__name__`, whatever.
            return 'A composition of functions'

    @property
    def __name__(self):
        try:
            return '_of_'.join(
                f.__name__ for f in reversed((self.first,) + self.funcs)
            )
        except AttributeError:
            return type(self).__name__


def compose(*funcs):
    """ Compose functions to operate in series.

    Returns a function that applies other functions in sequence.

    Functions are applied from right to left so that
    ``compose(f, g, h)(x, y)`` is the same as ``f(g(h(x, y)))``.

    If no arguments are provided, the identity function (f(x) = x) is returned.

    >>> inc = lambda i: i + 1
    >>> compose(str, inc)(3)
    '4'

    See Also:
        pipe
    """
    if not funcs:
        return identity
    if len(funcs) == 1:
        return funcs[0]
    else:
        return Compose(funcs)


def pipe(data, *funcs):
    """ Pipe a value through a sequence of functions

    I.e. ``pipe(data, f, g, h)`` is equivalent to ``h(g(f(data)))``

    We think of the value as progressing through a pipe of several
    transformations, much like pipes in UNIX

    ``$ cat data | f | g | h``

    >>> double = lambda i: 2 * i
    >>> pipe(3, double, str)
    '6'

    See Also:
        compose
        thread_first
        thread_last
    """
    for func in funcs:
        data = func(data)
    return data


def complement(func):
    """ Convert a predicate function to its logical complement.

    In other words, return a function that, for inputs that normally
    yield True, yields False, and vice-versa.

    >>> def iseven(n): return n % 2 == 0
    >>> isodd = complement(iseven)
    >>> iseven(2)
    True
    >>> isodd(2)
    False
    """
    return compose(operator.not_, func)


class juxt(object):
    """ Creates a function that calls several functions with the same arguments

    Takes several functions and returns a function that applies its arguments
    to each of those functions then returns a tuple of the results.

    Name comes from juxtaposition: the fact of two things being seen or placed
    close together with contrasting effect.

    >>> inc = lambda x: x + 1
    >>> double = lambda x: x * 2
    >>> juxt(inc, double)(10)
    (11, 20)
    >>> juxt([inc, double])(10)
    (11, 20)
    """
    __slots__ = ['funcs']

    def __init__(self, *funcs):
        if len(funcs) == 1 and not callable(funcs[0]):
            funcs = funcs[0]
        self.funcs = tuple(funcs)

    def __call__(self, *args, **kwargs):
        return tuple(func(*args, **kwargs) for func in self.funcs)

    def __getstate__(self):
        return self.funcs

    def __setstate__(self, state):
        self.funcs = state


def do(func, x):
    """ Runs ``func`` on ``x``, returns ``x``

    Because the results of ``func`` are not returned, only the side
    effects of ``func`` are relevant.

    Logging functions can be made by composing ``do`` with a storage function
    like ``list.append`` or ``file.write``

    >>> from toolz import compose
    >>> from toolz.curried import do

    >>> log = []
    >>> inc = lambda x: x + 1
    >>> inc = compose(inc, do(log.append))
    >>> inc(1)
    2
    >>> inc(11)
    12
    >>> log
    [1, 11]
    """
    func(x)
    return x


@curry
def flip(func, a, b):
    """ Call the function call with the arguments flipped

    This function is curried.

    >>> def div(a, b):
    ...     return a // b
    ...
    >>> flip(div, 2, 6)
    3
    >>> div_by_two = flip(div, 2)
    >>> div_by_two(4)
    2

    This is particularly useful for built in functions and functions defined
    in C extensions that accept positional only arguments. For example:
    isinstance, issubclass.

    >>> data = [1, 'a', 'b', 2, 1.5, object(), 3]
    >>> only_ints = list(filter(flip(isinstance, int), data))
    >>> only_ints
    [1, 2, 3]
    """
    return func(b, a)


def return_none(exc):
    """ Returns None.
    """
    return None


class excepts(object):
    """A wrapper around a function to catch exceptions and
    dispatch to a handler.

    This is like a functional try/except block, in the same way that
    ifexprs are functional if/else blocks.

    Examples
    --------
    >>> excepting = excepts(
    ...     ValueError,
    ...     lambda a: [1, 2].index(a),
    ...     lambda _: -1,
    ... )
    >>> excepting(1)
    0
    >>> excepting(3)
    -1

    Multiple exceptions and default except clause.
    >>> excepting = excepts((IndexError, KeyError), lambda a: a[0])
    >>> excepting([])
    >>> excepting([1])
    1
    >>> excepting({})
    >>> excepting({0: 1})
    1
    """
    def __init__(self, exc, func, handler=return_none):
        self.exc = exc
        self.func = func
        self.handler = handler

    def __call__(self, *args, **kwargs):
        try:
            return self.func(*args, **kwargs)
        except self.exc as e:
            return self.handler(e)

    @instanceproperty(classval=__doc__)
    def __doc__(self):
        exc = self.exc
        try:
            if isinstance(exc, tuple):
                exc_name = '(%s)' % ', '.join(
                    map(attrgetter('__name__'), exc),
                )
            else:
                exc_name = exc.__name__

            return dedent(
                """\
                A wrapper around {inst.func.__name__!r} that will except:
                {exc}
                and handle any exceptions with {inst.handler.__name__!r}.

                Docs for {inst.func.__name__!r}:
                {inst.func.__doc__}

                Docs for {inst.handler.__name__!r}:
                {inst.handler.__doc__}
                """
            ).format(
                inst=self,
                exc=exc_name,
            )
        except AttributeError:
            return type(self).__doc__

    @property
    def __name__(self):
        exc = self.exc
        try:
            if isinstance(exc, tuple):
                exc_name = '_or_'.join(map(attrgetter('__name__'), exc))
            else:
                exc_name = exc.__name__
            return '%s_excepting_%s' % (self.func.__name__, exc_name)
        except AttributeError:
            return 'excepting'


if PY3:  # pragma: py2 no cover
    def _check_sigspec(sigspec, func, builtin_func, *builtin_args):
        if sigspec is None:
            try:
                sigspec = inspect.signature(func)
            except (ValueError, TypeError) as e:
                sigspec = e
        if isinstance(sigspec, ValueError):
            return None, builtin_func(*builtin_args)
        elif not isinstance(sigspec, inspect.Signature):
            if (
                func in _sigs.signatures
                and ((
                    hasattr(func, '__signature__')
                    and hasattr(func.__signature__, '__get__')
                ) or (
                    PY33
                    and hasattr(func, '__wrapped__')
                    and hasattr(func.__wrapped__, '__get__')
                    and not callable(func.__wrapped__)
                ))
            ):  # pragma: no cover (not covered in Python 3.4)
                val = builtin_func(*builtin_args)
                return None, val
            return None, False
        return sigspec, None

else:  # pragma: py3 no cover
    def _check_sigspec(sigspec, func, builtin_func, *builtin_args):
        if sigspec is None:
            try:
                sigspec = inspect.getargspec(func)
            except TypeError as e:
                sigspec = e
        if isinstance(sigspec, TypeError):
            if not callable(func):
                return None, False
            return None, builtin_func(*builtin_args)
        return sigspec, None


if PY34 or PYPY:  # pragma: no cover
    _check_sigspec_orig = _check_sigspec

    def _check_sigspec(sigspec, func, builtin_func, *builtin_args):
        # Python 3.4 and PyPy may lie, so use our registry for builtins instead
        if func in _sigs.signatures:
            val = builtin_func(*builtin_args)
            return None, val
        return _check_sigspec_orig(sigspec, func, builtin_func, *builtin_args)

_check_sigspec.__doc__ = """ \
Private function to aid in introspection compatibly across Python versions.

If a callable doesn't have a signature (Python 3) or an argspec (Python 2),
the signature registry in toolz._signatures is used.
"""

if PY3:  # pragma: py2 no cover
    def num_required_args(func, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._num_required_args,
                                     func)
        if sigspec is None:
            return rv
        return sum(1 for p in sigspec.parameters.values()
                   if p.default is p.empty
                   and p.kind in (p.POSITIONAL_OR_KEYWORD, p.POSITIONAL_ONLY))

    def has_varargs(func, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._has_varargs, func)
        if sigspec is None:
            return rv
        return any(p.kind == p.VAR_POSITIONAL
                   for p in sigspec.parameters.values())

    def has_keywords(func, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._has_keywords, func)
        if sigspec is None:
            return rv
        return any(p.default is not p.empty
                   or p.kind in (p.KEYWORD_ONLY, p.VAR_KEYWORD)
                   for p in sigspec.parameters.values())

    def is_valid_args(func, args, kwargs, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._is_valid_args,
                                     func, args, kwargs)
        if sigspec is None:
            return rv
        try:
            sigspec.bind(*args, **kwargs)
        except TypeError:
            return False
        return True

    def is_partial_args(func, args, kwargs, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._is_partial_args,
                                     func, args, kwargs)
        if sigspec is None:
            return rv
        try:
            sigspec.bind_partial(*args, **kwargs)
        except TypeError:
            return False
        return True

else:  # pragma: py3 no cover
    def num_required_args(func, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._num_required_args,
                                     func)
        if sigspec is None:
            return rv
        num_defaults = len(sigspec.defaults) if sigspec.defaults else 0
        return len(sigspec.args) - num_defaults

    def has_varargs(func, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._has_varargs, func)
        if sigspec is None:
            return rv
        return sigspec.varargs is not None

    def has_keywords(func, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._has_keywords, func)
        if sigspec is None:
            return rv
        return sigspec.defaults is not None or sigspec.keywords is not None

    def is_valid_args(func, args, kwargs, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._is_valid_args,
                                     func, args, kwargs)
        if sigspec is None:
            return rv
        spec = sigspec
        defaults = spec.defaults or ()
        num_pos = len(spec.args) - len(defaults)
        missing_pos = spec.args[len(args):num_pos]
        if any(arg not in kwargs for arg in missing_pos):
            return False

        if spec.varargs is None:
            num_extra_pos = max(0, len(args) - num_pos)
        else:
            num_extra_pos = 0

        kwargs = dict(kwargs)

        # Add missing keyword arguments (unless already included in `args`)
        missing_kwargs = spec.args[num_pos + num_extra_pos:]
        kwargs.update(zip(missing_kwargs, defaults[num_extra_pos:]))

        # Convert call to use positional arguments
        args = args + tuple(kwargs.pop(key) for key in spec.args[len(args):])

        if (
            not spec.keywords and kwargs
            or not spec.varargs and len(args) > len(spec.args)
            or set(spec.args[:len(args)]) & set(kwargs)
        ):
            return False
        else:
            return True

    def is_partial_args(func, args, kwargs, sigspec=None):
        sigspec, rv = _check_sigspec(sigspec, func, _sigs._is_partial_args,
                                     func, args, kwargs)
        if sigspec is None:
            return rv
        spec = sigspec
        defaults = spec.defaults or ()
        num_pos = len(spec.args) - len(defaults)
        if spec.varargs is None:
            num_extra_pos = max(0, len(args) - num_pos)
        else:
            num_extra_pos = 0

        kwargs = dict(kwargs)

        # Add missing keyword arguments (unless already included in `args`)
        missing_kwargs = spec.args[num_pos + num_extra_pos:]
        kwargs.update(zip(missing_kwargs, defaults[num_extra_pos:]))

        # Add missing position arguments as keywords (may already be in kwargs)
        missing_args = spec.args[len(args):num_pos + num_extra_pos]
        kwargs.update((x, None) for x in missing_args)

        # Convert call to use positional arguments
        args = args + tuple(kwargs.pop(key) for key in spec.args[len(args):])

        if (
            not spec.keywords and kwargs
            or not spec.varargs and len(args) > len(spec.args)
            or set(spec.args[:len(args)]) & set(kwargs)
        ):
            return False
        else:
            return True


def is_arity(n, func, sigspec=None):
    """ Does a function have only n positional arguments?

    This function relies on introspection and does not call the function.
    Returns None if validity can't be determined.

    >>> def f(x):
    ...     return x
    >>> is_arity(1, f)
    True
    >>> def g(x, y=1):
    ...     return x + y
    >>> is_arity(1, g)
    False
    """
    sigspec, rv = _check_sigspec(sigspec, func, _sigs._is_arity, n, func)
    if sigspec is None:
        return rv
    num = num_required_args(func, sigspec)
    if num is not None:
        num = num == n
        if not num:
            return False
    varargs = has_varargs(func, sigspec)
    if varargs:
        return False
    keywords = has_keywords(func, sigspec)
    if keywords:
        return False
    if num is None or varargs is None or keywords is None:  # pragma: no cover
        return None
    return True


num_required_args.__doc__ = """ \
Number of required positional arguments

    This function relies on introspection and does not call the function.
    Returns None if validity can't be determined.

    >>> def f(x, y, z=3):
    ...     return x + y + z
    >>> num_required_args(f)
    2
    >>> def g(*args, **kwargs):
    ...     pass
    >>> num_required_args(g)
    0
    """

has_varargs.__doc__ = """ \
Does a function have variadic positional arguments?

    This function relies on introspection and does not call the function.
    Returns None if validity can't be determined.

    >>> def f(*args):
    ...    return args
    >>> has_varargs(f)
    True
    >>> def g(**kwargs):
    ...    return kwargs
    >>> has_varargs(g)
    False
    """

has_keywords.__doc__ = """ \
Does a function have keyword arguments?

    This function relies on introspection and does not call the function.
    Returns None if validity can't be determined.

    >>> def f(x, y=0):
    ...     return x + y

    >>> has_keywords(f)
    True
    """

is_valid_args.__doc__ = """ \
Is ``func(*args, **kwargs)`` a valid function call?

    This function relies on introspection and does not call the function.
    Returns None if validity can't be determined.

    >>> def add(x, y):
    ...     return x + y

    >>> is_valid_args(add, (1,), {})
    False
    >>> is_valid_args(add, (1, 2), {})
    True
    >>> is_valid_args(map, (), {})
    False

    **Implementation notes**
    Python 2 relies on ``inspect.getargspec``, which only works for
    user-defined functions.  Python 3 uses ``inspect.signature``, which
    works for many more types of callables.

    Many builtins in the standard library are also supported.
    """

is_partial_args.__doc__ = """ \
Can partial(func, *args, **kwargs)(*args2, **kwargs2) be a valid call?

    Returns True *only* if the call is valid or if it is possible for the
    call to become valid by adding more positional or keyword arguments.

    This function relies on introspection and does not call the function.
    Returns None if validity can't be determined.

    >>> def add(x, y):
    ...     return x + y

    >>> is_partial_args(add, (1,), {})
    True
    >>> is_partial_args(add, (1, 2), {})
    True
    >>> is_partial_args(add, (1, 2, 3), {})
    False
    >>> is_partial_args(map, (), {})
    True

    **Implementation notes**
    Python 2 relies on ``inspect.getargspec``, which only works for
    user-defined functions.  Python 3 uses ``inspect.signature``, which
    works for many more types of callables.

    Many builtins in the standard library are also supported.
    """

from . import _signatures as _sigs
