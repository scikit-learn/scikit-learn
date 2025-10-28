# cython.* namespace for pure mode.

# Possible version formats: "3.1.0", "3.1.0a1", "3.1.0a1.dev0"
__version__ = "3.1.6"


# BEGIN shameless copy from Cython/minivect/minitypes.py

class _ArrayType:

    is_array = True
    subtypes = ['dtype']

    def __init__(self, dtype, ndim, is_c_contig=False, is_f_contig=False,
                 inner_contig=False, broadcasting=None):
        self.dtype = dtype
        self.ndim = ndim
        self.is_c_contig = is_c_contig
        self.is_f_contig = is_f_contig
        self.inner_contig = inner_contig or is_c_contig or is_f_contig
        self.broadcasting = broadcasting

    def __repr__(self):
        axes = [":"] * self.ndim
        if self.is_c_contig:
            axes[-1] = "::1"
        elif self.is_f_contig:
            axes[0] = "::1"

        return "%s[%s]" % (self.dtype, ", ".join(axes))


def index_type(base_type, item):
    """
    Support array type creation by slicing, e.g. double[:, :] specifies
    a 2D strided array of doubles. The syntax is the same as for
    Cython memoryviews.
    """
    class InvalidTypeSpecification(Exception):
        pass

    def verify_slice(s):
        if s.start or s.stop or s.step not in (None, 1):
            raise InvalidTypeSpecification(
                "Only a step of 1 may be provided to indicate C or "
                "Fortran contiguity")

    if isinstance(item, tuple):
        step_idx = None
        for idx, s in enumerate(item):
            verify_slice(s)
            if s.step and (step_idx or idx not in (0, len(item) - 1)):
                raise InvalidTypeSpecification(
                    "Step may only be provided once, and only in the "
                    "first or last dimension.")

            if s.step == 1:
                step_idx = idx

        return _ArrayType(base_type, len(item),
                          is_c_contig=step_idx == len(item) - 1,
                          is_f_contig=step_idx == 0)
    elif isinstance(item, slice):
        verify_slice(item)
        return _ArrayType(base_type, 1, is_c_contig=bool(item.step))
    else:
        # int[8] etc.
        assert int(item) == item  # array size must be a plain integer
        return array(base_type, item)

# END shameless copy


compiled = False

_Unspecified = object()

# Function decorators

def _empty_decorator(x):
    return x

def locals(**arg_types):
    return _empty_decorator

def test_assert_path_exists(*paths):
    return _empty_decorator

def test_fail_if_path_exists(*paths):
    return _empty_decorator

class _EmptyDecoratorAndManager:
    def __call__(self, x):
        return x
    def __enter__(self):
        pass
    def __exit__(self, exc_type, exc_value, traceback):
        pass

class _Optimization:
    pass

cclass = ccall = cfunc = _EmptyDecoratorAndManager()

annotation_typing = returns = wraparound = boundscheck = initializedcheck = \
    nonecheck = embedsignature = cdivision = cdivision_warnings = \
    always_allow_keywords = profile = linetrace = infer_types = \
    unraisable_tracebacks = freelist = auto_pickle = cpow = trashcan = \
    auto_cpdef = c_api_binop_methods = \
    allow_none_for_extension_args = callspec = show_performance_hints = \
    cpp_locals = py2_import = iterable_coroutine = remove_unreachable = \
    overflowcheck = \
        lambda _: _EmptyDecoratorAndManager()

# Note that fast_getattr is untested and undocumented!
fast_getattr = lambda _: _EmptyDecoratorAndManager()
# c_compile_guard is largely for internal use
c_compile_guard = lambda _:_EmptyDecoratorAndManager()

exceptval = lambda _=None, check=True: _EmptyDecoratorAndManager()

optimize = _Optimization()


embedsignature.format = overflowcheck.fold = optimize.use_switch = \
    optimize.unpack_method_calls = lambda arg: _EmptyDecoratorAndManager()

final = internal = type_version_tag = no_gc_clear = no_gc = total_ordering = \
    ufunc = _empty_decorator

binding = lambda _: _empty_decorator

class warn:
    undeclared = unreachable = maybe_uninitialized = unused = \
        unused_arg = unused_result = \
            lambda _: _EmptyDecoratorAndManager()


_cython_inline = None
def inline(f, *args, **kwds):
    if isinstance(f, str):
        global _cython_inline
        if _cython_inline is None:
            from Cython.Build.Inline import cython_inline as _cython_inline
        return _cython_inline(f, *args, **kwds)
    else:
        assert len(args) == len(kwds) == 0
        return f


def compile(f):
    from Cython.Build.Inline import RuntimeCompiledFunction
    return RuntimeCompiledFunction(f)


# Special functions

def cdiv(a, b):
    if a < 0:
        a = -a
        b = -b
    if b < 0:
        return (a + b + 1) // b
    return a // b

def cmod(a, b):
    r = a % b
    if (a * b) < 0 and r:
        r -= b
    return r


# Emulated language constructs

def cast(t, *args, **kwargs):
    kwargs.pop('typecheck', None)
    assert not kwargs

    if isinstance(t, typedef):
        return t(*args)
    elif isinstance(t, type):  # Doesn't work with old-style classes of Python 2.x
        if len(args) != 1 or not (args[0] is None or isinstance(args[0], t)):
            return t(*args)

    return args[0]

def sizeof(arg):
    return 1

def typeof(arg):
    return arg.__class__.__name__
    # return type(arg)

def address(arg):
    return pointer(type(arg))([arg])

def _is_value_type(t):
    if isinstance(t, typedef):
        return _is_value_type(t._basetype)

    return isinstance(t, type) and issubclass(t, (StructType, UnionType, ArrayType))

def declare(t=None, value=_Unspecified, **kwds):
    if value is not _Unspecified:
        return cast(t, value)
    elif _is_value_type(t):
        return t()
    else:
        return None

class _nogil:
    """Support for 'with nogil' statement and @nogil decorator.
    """
    def __call__(self, x):
        if callable(x):
            # Used as function decorator => return the function unchanged.
            return x
        # Used as conditional context manager or to create an "@nogil(True/False)" decorator => keep going.
        return self

    def __enter__(self):
        pass
    def __exit__(self, exc_class, exc, tb):
        return exc_class is None

nogil = _nogil()
gil = _nogil()
with_gil = _nogil()  # Actually not a context manager, but compilation will give the right error.
del _nogil


class critical_section:
    def __init__(self, arg0, arg1=None):
        # It's ambiguous if this is being used as a decorator or context manager
        # even with a callable arg.
        self.arg0 = arg0
    def __call__(self, *args, **kwds):
        return self.arg0(*args, **kwds)
    def __enter__(self):
        pass
    def __exit__(self, exc_class, exc, tb):
        return False


# Emulated types

class CythonMetaType(type):

    def __getitem__(type, ix):
        return array(type, ix)

CythonTypeObject = CythonMetaType('CythonTypeObject', (object,), {})

class CythonType(CythonTypeObject):

    def _pointer(self, n=1):
        for i in range(n):
            self = pointer(self)
        return self

class PointerType(CythonType):

    def __init__(self, value=None):
        if isinstance(value, (ArrayType, PointerType)):
            self._items = [cast(self._basetype, a) for a in value._items]
        elif isinstance(value, list):
            self._items = [cast(self._basetype, a) for a in value]
        elif value is None or value == 0:
            self._items = []
        else:
            raise ValueError

    def __getitem__(self, ix):
        if ix < 0:
            raise IndexError("negative indexing not allowed in C")
        return self._items[ix]

    def __setitem__(self, ix, value):
        if ix < 0:
            raise IndexError("negative indexing not allowed in C")
        self._items[ix] = cast(self._basetype, value)

    def __eq__(self, value):
        if value is None and not self._items:
            return True
        elif type(self) != type(value):
            return False
        else:
            return not self._items and not value._items

    def __repr__(self):
        return f"{self._basetype} *"


class ArrayType(PointerType):

    def __init__(self, value=None):
        if value is None:
            self._items = [None] * self._n
        else:
            super().__init__(value)


class StructType(CythonType):

    def __init__(self, *posargs, **data):
        if not (posargs or data):
            return
        if posargs and data:
            raise ValueError('Cannot accept both positional and keyword arguments.')

        # Allow 'cast_from' as single positional or keyword argument.
        if data and len(data) == 1 and 'cast_from' in data:
            cast_from = data.pop('cast_from')
        elif len(posargs) == 1 and type(posargs[0]) is type(self):
            cast_from, posargs = posargs[0], ()
        elif posargs:
            for key, arg in zip(self._members, posargs):
                setattr(self, key, arg)
            return
        else:
            for key, value in data.items():
                if key not in self._members:
                    raise ValueError("Invalid struct attribute for %s: %s" % (
                        self.__class__.__name__, key))
                setattr(self, key, value)
            return

        # do cast
        if data:
            raise ValueError('Cannot accept keyword arguments when casting.')
        if type(cast_from) is not type(self):
            raise ValueError('Cannot cast from %s' % cast_from)
        for key, value in cast_from.__dict__.items():
            setattr(self, key, value)

    def __setattr__(self, key, value):
        if key in self._members:
            self.__dict__[key] = cast(self._members[key], value)
        else:
            raise AttributeError("Struct has no member '%s'" % key)


class UnionType(CythonType):

    def __init__(self, cast_from=_Unspecified, **data):
        if cast_from is not _Unspecified:
            # do type cast
            if len(data) > 0:
                raise ValueError('Cannot accept keyword arguments when casting.')
            if isinstance(cast_from, dict):
                datadict = cast_from
            elif type(cast_from) is type(self):
                datadict = cast_from.__dict__
            else:
                raise ValueError('Cannot cast from %s' % cast_from)
        else:
            datadict = data
        if len(datadict) > 1:
            raise AttributeError("Union can only store one field at a time.")
        for key, value in datadict.items():
            setattr(self, key, value)

    def __setattr__(self, key, value):
        if key == '__dict__':
            CythonType.__setattr__(self, key, value)
        elif key in self._members:
            self.__dict__ = {key: cast(self._members[key], value)}
        else:
            raise AttributeError("Union has no member '%s'" % key)


class pointer(PointerType):
    # Implemented as class to support both 'pointer(int)' and 'pointer[int]'.
    def __new__(cls, basetype):
        class PointerInstance(PointerType):
            _basetype = basetype
        return PointerInstance

    def __class_getitem__(cls, basetype):
        return cls(basetype)


class array(ArrayType):
    # Implemented as class to support both 'array(int, 5)' and 'array[int, 5]'.
    def __new__(cls, basetype, n):
        class ArrayInstance(ArrayType):
            _basetype = basetype
            _n = n
        return ArrayInstance

    def __class_getitem__(cls, item):
        basetype, n = item
        return cls(basetype, item)


def struct(**members):
    class StructInstance(StructType):
        _members = members
    for key in members:
        setattr(StructInstance, key, None)
    return StructInstance

def union(**members):
    class UnionInstance(UnionType):
        _members = members
    for key in members:
        setattr(UnionInstance, key, None)
    return UnionInstance


class typedef(CythonType):

    def __init__(self, type, name=None):
        self._basetype = type
        self.name = name

    def __call__(self, *arg):
        value = cast(self._basetype, *arg)
        return value

    def __repr__(self):
        return self.name or str(self._basetype)

    __getitem__ = index_type


class const(typedef):
    def __init__(self, type, name=None):
        name = f"const {name or repr(type)}"
        super().__init__(type, name)

    def __class_getitem__(cls, base_type):
        return const(base_type)


class volatile(typedef):
    def __init__(self, type, name=None):
        name = f"volatile {name or repr(type)}"
        super().__init__(type, name)

    def __class_getitem__(cls, base_type):
        return volatile(base_type)


class _FusedType(CythonType):
    __getitem__ = index_type


def fused_type(*args):
    if not args:
        raise TypeError("Expected at least one type as argument")

    # Find the numeric type with biggest rank if all types are numeric
    rank = -1
    for type in args:
        if type not in (py_int, py_long, py_float, py_complex):
            break

        if type_ordering.index(type) > rank:
            result_type = type
    else:
        return result_type

    # Not a simple numeric type, return a fused type instance. The result
    # isn't really meant to be used, as we can't keep track of the context in
    # pure-mode. Casting won't do anything in this case.
    return _FusedType()


def _specialized_from_args(signatures, args, kwargs):
    "Perhaps this should be implemented in a TreeFragment in Cython code"
    raise Exception("yet to be implemented")


py_int = typedef(int, "int")
py_long = typedef(int, "long")  # for legacy Py2 code only
py_float = typedef(float, "float")
py_complex = typedef(complex, "double complex")


# Predefined types

int_types = [
    'char',
    'short',
    'Py_UNICODE',
    'int',
    'Py_UCS4',
    'long',
    'longlong',
    'Py_hash_t',
    'Py_ssize_t',
    'size_t',
    'ssize_t',
    'ptrdiff_t',
]
float_types = [
    'longdouble',
    'double',
    'float',
]
complex_types = [
    'longdoublecomplex',
    'doublecomplex',
    'floatcomplex',
    'complex',
]
other_types = [
    'bint',
    'void',
    'Py_tss_t',
]

to_repr = {
    'longlong': 'long long',
    'longdouble': 'long double',
    'longdoublecomplex': 'long double complex',
    'doublecomplex': 'double complex',
    'floatcomplex': 'float complex',
}.get

gs = globals()

gs['unicode'] = typedef(str, 'unicode')

for name in int_types:
    reprname = to_repr(name, name)
    gs[name] = typedef(py_int, reprname)
    if name not in ('Py_UNICODE', 'Py_UCS4', 'Py_hash_t', 'ptrdiff_t') and not name.endswith('size_t'):
        gs['u'+name] = typedef(py_int, "unsigned " + reprname)
        gs['s'+name] = typedef(py_int, "signed " + reprname)

for name in float_types:
    gs[name] = typedef(py_float, to_repr(name, name))

for name in complex_types:
    gs[name] = typedef(py_complex, to_repr(name, name))

del name, reprname

bint = typedef(bool, "bint")
void = typedef(None, "void")
Py_tss_t = typedef(None, "Py_tss_t")

# Generate const types.
for t in int_types + float_types + complex_types + other_types:
    for t in (t, f'u{t}', f's{t}'):
        if t in gs:
            gs[f"const_{t}"] = const(gs[t], t)

# Generate pointer types: p_int, p_const_char, etc.
for i in range(1, 4):
    for const_ in ('', 'const_'):
        for t in int_types:
            for t in (t, f'u{t}', f's{t}'):
                if t in gs:
                    gs[f"{'p'*i}_{const_}{t}"] = pointer(gs[f"{'p'*(i-1)}{'_' if i > 1 else ''}{const_}{t}"])

        for t in float_types + complex_types:
            gs[f"{'p'*i}_{const_}{t}"] = pointer(gs[f"{'p'*(i-1)}{'_' if i > 1 else ''}{const_}{t}"])

    gs[f"{'p'*i}_const_bint"] = pointer(gs[f"{'p'*(i-1)}{'_' if i > 1 else ''}const_bint"])
    for t in other_types:
        gs[f"{'p'*i}_{t}"] = pointer(gs[f"{'p'*(i-1)}{'_' if i > 1 else ''}{t}"])

del t, const_, i

NULL = gs['p_void'](0)

del gs


def __getattr__(name):
    # looks like 'gs' has some users out there by now...
    if name == 'gs':
        import warnings
        warnings.warn(
            "'gs' is not a publicly exposed name in cython.*. Use vars() or globals() instead.",
            DeprecationWarning)
        return globals()
    raise AttributeError(f"'cython' has no attribute {name!r}")


integral = floating = numeric = _FusedType()

type_ordering = [py_int, py_long, py_float, py_complex]

class CythonDotParallel:
    """
    The cython.parallel module.
    """

    __all__ = ['parallel', 'prange', 'threadid']

    def parallel(self, num_threads=None):
        return nogil

    def prange(self, start=0, stop=None, step=1, nogil=False, schedule=None, chunksize=None, num_threads=None):
        if stop is None:
            stop = start
            start = 0
        return range(start, stop, step)

    def threadid(self):
        return 0

    # def threadsavailable(self):
        # return 1

class CythonDotImportedFromElsewhere:
    """
    cython.dataclasses just shadows the standard library modules of the same name
    """
    def __init__(self, module):
        self.__path__ = []
        self.__file__ = None
        self.__name__ = module
        self.__package__ = module

    def __getattr__(self, attr):
        # we typically only expect this to be called once
        from importlib import import_module
        import sys
        try:
            mod = import_module(self.__name__)
        except ImportError:
            # but if they don't exist (Python is not sufficiently up-to-date) then
            # you can't use them
            raise AttributeError("%s: the standard library module %s is not available" %
                                 (attr, self.__name__))
        sys.modules['cython.%s' % self.__name__] = mod
        return getattr(mod, attr)

class CythonCImports:
    """
    Simplistic module mock to make cimports sort-of work in Python code.
    """
    def __init__(self, module, **attributes):
        self.__path__ = []
        self.__file__ = None
        self.__name__ = module
        self.__package__ = module
        if attributes:
            self.__dict__.update(attributes)

    def __getattr__(self, item):
        if item.startswith('__') and item.endswith('__'):
            raise AttributeError(item)

        package = self.__package__[len('cython.cimports.'):]

        from importlib import import_module
        try:
            return import_module(item, package or None)
        except ImportError:
            ex = AttributeError(item)
            ex.__cause__ = None
            raise ex


import math, sys
sys.modules['cython.parallel'] = CythonDotParallel()
sys.modules['cython.cimports.libc.math'] = math
sys.modules['cython.cimports.libc'] = CythonCImports('cython.cimports.libc', math=math)
sys.modules['cython.cimports'] = CythonCImports('cython.cimports', libc=sys.modules['cython.cimports.libc'])

# In pure Python mode @cython.dataclasses.dataclass and dataclass field should just
# shadow the standard library ones (if they are available)
dataclasses = sys.modules['cython.dataclasses'] = CythonDotImportedFromElsewhere('dataclasses')
del math, sys

class pymutex:
    def __init__(self):
        import threading
        self._l = threading.Lock()

    def acquire(self):
        return self._l.acquire()

    def release(self):
        return self._l.release()

    def __enter__(self):
        return self._l.__enter__()

    def __exit__(self, exc_type, exc_value, traceback):
        return self._l.__exit__(exc_type, exc_value, traceback)

pythread_type_lock = pymutex
