import dataclasses as dataclasses
from builtins import (int as py_int, float as py_float,
                      bool as py_bool, str as py_str, complex as py_complex)
from types import TracebackType
from typing import (
    Any, Iterable, Sequence, Optional, Type, TypeVar, Generic, Callable, overload,
    TypeAlias, Annotated,
)

# Type checkers assume typing_extensions is always present
from typing_extensions import Literal, ParamSpec, overload, final as final

# Normally all imports aren't exported in stub files but we can explicitly expose
# imports by using import ... as ... (with the same name) which was done for
# dataclasses and the final decorator.

__version__: str

# Predefined types

Py_UCS4 = py_int | str
Py_UNICODE = py_int | str

bint = py_bool
void = Type[None]
basestring = py_str
unicode = py_str

compiled: bool


_T = TypeVar('_T')
_P = ParamSpec('_P')
_C = TypeVar('_C', bound='Callable')
_TypeT = TypeVar('_TypeT', bound='Type')
_Decorator = Callable[[_C], _C]


_func_deco: _Decorator

cfunc = ccall = compile = ufunc = _func_deco

def locals(**kwargs: Any) -> _Decorator: ...

def _class_deco(__cls: _TypeT) -> _TypeT: ...

cclass = internal = c_api_binop_methods = type_version_tag = no_gc_clear = no_gc = total_ordering = _class_deco

# May be a bit hard to read but essentially means:
# > Returns a callable that takes another callable with these parameters and *some*
# > return value, then returns another callable with the same parameters but
# > the return type is the previous 'type' parameter.
# On Python 3.5, the latest version of Mypy available is 0.910 which doesn't understand ParamSpec
def returns(__type: Type[_T]) -> Callable[[Callable[_P, object]], Callable[_P, _T]]: ...  # type: ignore

def exceptval(__val: Any, *, check: bool = False) -> _Decorator: ...

class _EmptyDecoratorAndManager(object):
    @overload
    def __call__(self, __val: bool) -> _Decorator: ...

    @overload
    def __call__(self, __func: _C) -> _C: ...

    def __enter__(self) -> None: ...

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc: Optional[BaseException],
        tb: Optional[TracebackType]
    ) -> None: ...
_empty_decorator_and_manager: _EmptyDecoratorAndManager

@overload
def _compiler_directive(__func: _C) -> _C: ...

@overload
def _compiler_directive(__val: bool = ...) -> _Decorator: ...

# These all come from 'Compiler directives' on Source Files and Compilation.
# The following directives are missing as they need to be global:
# - c_string_type
# - c_string_encoding
# Note that c_api_binop_methods and type_version_tag is defined above.

annotation_typing = returns = wraparound = boundscheck = initializedcheck = \
    nonecheck = cdivision = cdivision_warnings = \
    profile = linetrace = infer_types = \
    freelist = auto_pickle = cpow = trashcan = \
    auto_cpdef = c_api_binop_methods = \
    allow_none_for_extension_args = callspec = show_performance_hints = \
    cpp_locals = py2_import = iterable_coroutine = remove_unreachable = \
    _empty_decorator_and_manager

binding = embedsignature = always_allow_keywords = unraisable_tracebacks = \
    cpp_locals = \
    _compiler_directive

# overflowcheck() has to be specialized because there is also overflowcheck.fold
class _OverflowcheckClass:
    def __call__(self, __val: bool = ...) -> _Decorator: ...

    def fold(self, __val: bool = ...) -> _Decorator: ...

overflowcheck = _OverflowcheckClass()

class optimize:
    @staticmethod
    def use_switch(__val: bool = ...) -> _Decorator: ...

    @staticmethod
    def unpack_method_calls(__val: bool = ...) -> _Decorator: ...

class warn:
    @staticmethod
    def undeclared(__val: bool = ...) -> _Decorator: ...

    @staticmethod
    def unreachable(__val: bool = ...) -> _Decorator: ...

    @staticmethod
    def maybe_uninitialized(__val: bool = ...) -> _Decorator: ...

    @staticmethod
    def unused(__val: bool = ...) -> _Decorator: ...

    @staticmethod
    def unused_argument(__val: bool = ...) -> _Decorator: ...

    @staticmethod
    def multiple_declarators(__val: bool = ...) -> _Decorator: ...

@overload
def inline(__func: _C) -> _C: ...

@overload
def inline(__code: str, *, get_type: Callable[[object, object], str] = ..., lib_dir: str = ...,
           cython_include_dirs: Iterable[str] = ..., cython_compiler_directives: Iterable[str] = ...,
           force: bool = ..., quiet: bool = ..., locals: dict[str, str] = ..., globals: dict[str, str] = ...,
           language_level: str = ...) -> Any: ...

def cdiv(__a: int, __b: int) -> int: ...

def cmod(__a: int, __b: int) -> int: ...

@overload
def cast(__t: Type[_T], __value: Any) -> _T: ...

# On Python 3.5, the latest version of Mypy available is 0.910 which doesn't understand ParamSpec
@overload
def cast(__t: Callable[_P, _T], *args: _P.args, **kwargs: _P.kwargs) -> _T: ...  # type: ignore

def sizeof(__obj: object) -> int: ...

def typeof(__obj: object) -> str: ...

def address(__obj: object) -> PointerType: ...

const: TypeAlias = Annotated[_T, "cython.const"]
volatile: TypeAlias = Annotated[_T, "cython.volatile"]


@overload
def declare(
    t: Optional[Callable[..., _T]] = ...,
    value: Any = ...,
) -> _T:
    ...

# This one is for attributes, they cannot have initializers through cython.declare() currently.
@overload
def declare(
    t: Callable[..., _T],
    *,
    visibility: Literal['public', 'readonly', 'private'] = ...,
) -> _T:
    ...

@overload
def declare(**kwargs: type) -> None: ...


class _nogil:
    @overload
    def __call__(self, __val: bool) -> _Decorator: ...

    @overload
    def __call__(self, __func: _C) -> _C: ...

    @overload
    def __call__(self) -> '_nogil': ...

    def __enter__(self) -> None: ...

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc: Optional[BaseException],
        tb: Optional[TracebackType]
    ) -> None: ...

nogil = gil = _nogil


class _ArrayType(Generic[_T]):
    is_array: bool
    subtypes: Sequence[str]
    dtype: _T
    ndim: int
    is_c_contig: bool
    is_f_contig: bool
    inner_contig: bool
    broadcasting: Any

    # broadcasting is not used, so it's not clear about its type
    def __init__(self, dtype: _T, ndim: int, is_c_contig: bool = ...,
                 is_f_contig: bool = ..., inner_contig: bool = ...,
                 broadcasting: Any = ...) -> None: ...
    def __repr__(self) -> str: ...


class CythonTypeObject(object):
    ...

class CythonType(CythonTypeObject):
    ...

class PointerType(CythonType, Generic[_T]):
    def __init__(
        self,
        value: Optional[ArrayType[_T] | PointerType[_T] | list[_T] | int] = ...
    ) -> None: ...
    def __getitem__(self, ix: int) -> _T: ...
    def __setitem__(self, ix: int, value: _T) -> None: ...
    def __eq__(self, value: object) -> bool: ...
    def __repr__(self) -> str: ...

class ArrayType(PointerType[_T]):
    def __init__(self) -> None: ...

def index_type(
    base_type: _T, item: tuple | slice | int) -> _ArrayType[_T]: ...

class pointer(PointerType[_T]):
    def __new__(cls, basetype: _T) -> Type[PointerType[_T]]: ...
    def __class_getitem__(cls, basetype: _T) -> Type[PointerType[_T]]: ...

class array(ArrayType[_T]):
    def __new__(basetype: _T, n: int) -> Type[ArrayType[_T, int]]: ...
    def __class_getitem__(cls, item: tuple[_T, int]) -> Type[ArrayType[_T, int]]: ...

def struct(**members: type) -> Type[Any]: ...

def union(**members: type) -> Type[Any]: ...

def fused_type(*args: Any) -> Type[Any]: ...


class typedef(CythonType, Generic[_T]):
    name: str

    def __init__(self, type: _T, name: Optional[str] = ...) -> None: ...
    def __call__(self, *arg: Any) -> _T: ...
    def __repr__(self) -> str: ...
    __getitem__ = index_type

NULL: pointer[Any]

##### START: GENERATED LIST OF GENERATED TYPES #####
# Generated by "Tools/cython-generate-shadow-pyi.py" on 2024-09-24 11:45:22.474391

const_bint : TypeAlias = const[bint]
p_const_bint = pointer[const[bint]]
pp_const_bint = pointer[pointer[const[bint]]]
ppp_const_bint = pointer[pointer[pointer[const[bint]]]]
p_bint = pointer[bint]
pp_bint = pointer[pointer[bint]]
ppp_bint = pointer[pointer[pointer[bint]]]
char : TypeAlias = py_int
const_char : TypeAlias = const[py_int]
p_const_char : TypeAlias = pointer[const[py_int]]
pp_const_char = pointer[pointer[const[py_int]]]
ppp_const_char = pointer[pointer[pointer[const[py_int]]]]
p_char = pointer[py_int]
pp_char = pointer[pointer[py_int]]
ppp_char = pointer[pointer[pointer[py_int]]]
complex : TypeAlias = py_complex
const_complex : TypeAlias = const[py_complex]
p_const_complex = pointer[const[py_complex]]
pp_const_complex = pointer[pointer[const[py_complex]]]
ppp_const_complex = pointer[pointer[pointer[const[py_complex]]]]
p_complex = pointer[py_complex]
pp_complex = pointer[pointer[py_complex]]
ppp_complex = pointer[pointer[pointer[py_complex]]]
double : TypeAlias = py_float
const_double : TypeAlias = const[py_float]
p_const_double = pointer[const[py_float]]
pp_const_double = pointer[pointer[const[py_float]]]
ppp_const_double = pointer[pointer[pointer[const[py_float]]]]
p_double = pointer[py_float]
pp_double = pointer[pointer[py_float]]
ppp_double = pointer[pointer[pointer[py_float]]]
doublecomplex : TypeAlias = py_complex
const_doublecomplex : TypeAlias = const[py_complex]
p_const_doublecomplex = pointer[const[py_complex]]
pp_const_doublecomplex = pointer[pointer[const[py_complex]]]
ppp_const_doublecomplex = pointer[pointer[pointer[const[py_complex]]]]
p_doublecomplex = pointer[py_complex]
pp_doublecomplex = pointer[pointer[py_complex]]
ppp_doublecomplex = pointer[pointer[pointer[py_complex]]]
float : TypeAlias = py_float
const_float : TypeAlias = const[py_float]
p_const_float = pointer[const[py_float]]
pp_const_float = pointer[pointer[const[py_float]]]
ppp_const_float = pointer[pointer[pointer[const[py_float]]]]
p_float = pointer[py_float]
pp_float = pointer[pointer[py_float]]
ppp_float = pointer[pointer[pointer[py_float]]]
floatcomplex : TypeAlias = py_complex
const_floatcomplex : TypeAlias = const[py_complex]
p_const_floatcomplex = pointer[const[py_complex]]
pp_const_floatcomplex = pointer[pointer[const[py_complex]]]
ppp_const_floatcomplex = pointer[pointer[pointer[const[py_complex]]]]
p_floatcomplex = pointer[py_complex]
pp_floatcomplex = pointer[pointer[py_complex]]
ppp_floatcomplex = pointer[pointer[pointer[py_complex]]]
int : TypeAlias = py_int
const_int : TypeAlias = const[py_int]
p_const_int = pointer[const[py_int]]
pp_const_int = pointer[pointer[const[py_int]]]
ppp_const_int = pointer[pointer[pointer[const[py_int]]]]
p_int = pointer[py_int]
pp_int = pointer[pointer[py_int]]
ppp_int = pointer[pointer[pointer[py_int]]]
long : TypeAlias = py_int
const_long : TypeAlias = const[py_int]
p_const_long = pointer[const[py_int]]
pp_const_long = pointer[pointer[const[py_int]]]
ppp_const_long = pointer[pointer[pointer[const[py_int]]]]
p_long = pointer[py_int]
pp_long = pointer[pointer[py_int]]
ppp_long = pointer[pointer[pointer[py_int]]]
py_long : TypeAlias = py_int
longdouble : TypeAlias = py_float
const_longdouble : TypeAlias = const[py_float]
p_const_longdouble = pointer[const[py_float]]
pp_const_longdouble = pointer[pointer[const[py_float]]]
ppp_const_longdouble = pointer[pointer[pointer[const[py_float]]]]
p_longdouble = pointer[py_float]
pp_longdouble = pointer[pointer[py_float]]
ppp_longdouble = pointer[pointer[pointer[py_float]]]
longdoublecomplex : TypeAlias = py_complex
const_longdoublecomplex : TypeAlias = const[py_complex]
p_const_longdoublecomplex = pointer[const[py_complex]]
pp_const_longdoublecomplex = pointer[pointer[const[py_complex]]]
ppp_const_longdoublecomplex = pointer[pointer[pointer[const[py_complex]]]]
p_longdoublecomplex = pointer[py_complex]
pp_longdoublecomplex = pointer[pointer[py_complex]]
ppp_longdoublecomplex = pointer[pointer[pointer[py_complex]]]
longlong : TypeAlias = py_int
const_longlong : TypeAlias = const[py_int]
p_const_longlong = pointer[const[py_int]]
pp_const_longlong = pointer[pointer[const[py_int]]]
ppp_const_longlong = pointer[pointer[pointer[const[py_int]]]]
p_longlong = pointer[py_int]
pp_longlong = pointer[pointer[py_int]]
ppp_longlong = pointer[pointer[pointer[py_int]]]
schar : TypeAlias = py_int
const_schar : TypeAlias = const[py_int]
p_const_schar = pointer[const[py_int]]
pp_const_schar = pointer[pointer[const[py_int]]]
ppp_const_schar = pointer[pointer[pointer[const[py_int]]]]
p_schar = pointer[py_int]
pp_schar = pointer[pointer[py_int]]
ppp_schar = pointer[pointer[pointer[py_int]]]
short : TypeAlias = py_int
const_short : TypeAlias = const[py_int]
p_const_short = pointer[const[py_int]]
pp_const_short = pointer[pointer[const[py_int]]]
ppp_const_short = pointer[pointer[pointer[const[py_int]]]]
p_short = pointer[py_int]
pp_short = pointer[pointer[py_int]]
ppp_short = pointer[pointer[pointer[py_int]]]
sint : TypeAlias = py_int
const_sint : TypeAlias = const[py_int]
p_const_sint = pointer[const[py_int]]
pp_const_sint = pointer[pointer[const[py_int]]]
ppp_const_sint = pointer[pointer[pointer[const[py_int]]]]
p_sint = pointer[py_int]
pp_sint = pointer[pointer[py_int]]
ppp_sint = pointer[pointer[pointer[py_int]]]
slong : TypeAlias = py_int
const_slong : TypeAlias = const[py_int]
p_const_slong = pointer[const[py_int]]
pp_const_slong = pointer[pointer[const[py_int]]]
ppp_const_slong = pointer[pointer[pointer[const[py_int]]]]
p_slong = pointer[py_int]
pp_slong = pointer[pointer[py_int]]
ppp_slong = pointer[pointer[pointer[py_int]]]
slonglong : TypeAlias = py_int
const_slonglong : TypeAlias = const[py_int]
p_const_slonglong = pointer[const[py_int]]
pp_const_slonglong = pointer[pointer[const[py_int]]]
ppp_const_slonglong = pointer[pointer[pointer[const[py_int]]]]
p_slonglong = pointer[py_int]
pp_slonglong = pointer[pointer[py_int]]
ppp_slonglong = pointer[pointer[pointer[py_int]]]
sshort : TypeAlias = py_int
const_sshort : TypeAlias = const[py_int]
p_const_sshort = pointer[const[py_int]]
pp_const_sshort = pointer[pointer[const[py_int]]]
ppp_const_sshort = pointer[pointer[pointer[const[py_int]]]]
p_sshort = pointer[py_int]
pp_sshort = pointer[pointer[py_int]]
ppp_sshort = pointer[pointer[pointer[py_int]]]
Py_hash_t : TypeAlias = py_int
const_Py_hash_t : TypeAlias = const[py_int]
p_const_Py_hash_t = pointer[const[py_int]]
pp_const_Py_hash_t = pointer[pointer[const[py_int]]]
ppp_const_Py_hash_t = pointer[pointer[pointer[const[py_int]]]]
p_Py_hash_t = pointer[py_int]
pp_Py_hash_t = pointer[pointer[py_int]]
ppp_Py_hash_t = pointer[pointer[pointer[py_int]]]
ptrdiff_t : TypeAlias = py_int
const_ptrdiff_t : TypeAlias = const[py_int]
p_const_ptrdiff_t = pointer[const[py_int]]
pp_const_ptrdiff_t = pointer[pointer[const[py_int]]]
ppp_const_ptrdiff_t = pointer[pointer[pointer[const[py_int]]]]
p_ptrdiff_t = pointer[py_int]
pp_ptrdiff_t = pointer[pointer[py_int]]
ppp_ptrdiff_t = pointer[pointer[pointer[py_int]]]
size_t : TypeAlias = py_int
const_size_t : TypeAlias = const[py_int]
p_const_size_t = pointer[const[py_int]]
pp_const_size_t = pointer[pointer[const[py_int]]]
ppp_const_size_t = pointer[pointer[pointer[const[py_int]]]]
p_size_t = pointer[py_int]
pp_size_t = pointer[pointer[py_int]]
ppp_size_t = pointer[pointer[pointer[py_int]]]
ssize_t : TypeAlias = py_int
const_ssize_t : TypeAlias = const[py_int]
p_const_ssize_t = pointer[const[py_int]]
pp_const_ssize_t = pointer[pointer[const[py_int]]]
ppp_const_ssize_t = pointer[pointer[pointer[const[py_int]]]]
p_ssize_t = pointer[py_int]
pp_ssize_t = pointer[pointer[py_int]]
ppp_ssize_t = pointer[pointer[pointer[py_int]]]
Py_ssize_t : TypeAlias = py_int
const_Py_ssize_t : TypeAlias = const[py_int]
p_const_Py_ssize_t = pointer[const[py_int]]
pp_const_Py_ssize_t = pointer[pointer[const[py_int]]]
ppp_const_Py_ssize_t = pointer[pointer[pointer[const[py_int]]]]
p_Py_ssize_t = pointer[py_int]
pp_Py_ssize_t = pointer[pointer[py_int]]
ppp_Py_ssize_t = pointer[pointer[pointer[py_int]]]
Py_tss_t : TypeAlias = Any
const_Py_tss_t : TypeAlias = const[Any]
p_Py_tss_t = pointer[Any]
pp_Py_tss_t = pointer[pointer[Any]]
ppp_Py_tss_t = pointer[pointer[pointer[Any]]]
uchar : TypeAlias = py_int
const_uchar : TypeAlias = const[py_int]
p_const_uchar = pointer[const[py_int]]
pp_const_uchar = pointer[pointer[const[py_int]]]
ppp_const_uchar = pointer[pointer[pointer[const[py_int]]]]
p_uchar = pointer[py_int]
pp_uchar = pointer[pointer[py_int]]
ppp_uchar = pointer[pointer[pointer[py_int]]]
const_Py_UCS4 : TypeAlias = const[py_int]
p_const_Py_UCS4 = pointer[const[py_int]]
pp_const_Py_UCS4 = pointer[pointer[const[py_int]]]
ppp_const_Py_UCS4 = pointer[pointer[pointer[const[py_int]]]]
p_Py_UCS4 = pointer[py_int]
pp_Py_UCS4 = pointer[pointer[py_int]]
ppp_Py_UCS4 = pointer[pointer[pointer[py_int]]]
uint : TypeAlias = py_int
const_uint : TypeAlias = const[py_int]
p_const_uint = pointer[const[py_int]]
pp_const_uint = pointer[pointer[const[py_int]]]
ppp_const_uint = pointer[pointer[pointer[const[py_int]]]]
p_uint = pointer[py_int]
pp_uint = pointer[pointer[py_int]]
ppp_uint = pointer[pointer[pointer[py_int]]]
ulong : TypeAlias = py_int
const_ulong : TypeAlias = const[py_int]
p_const_ulong = pointer[const[py_int]]
pp_const_ulong = pointer[pointer[const[py_int]]]
ppp_const_ulong = pointer[pointer[pointer[const[py_int]]]]
p_ulong = pointer[py_int]
pp_ulong = pointer[pointer[py_int]]
ppp_ulong = pointer[pointer[pointer[py_int]]]
ulonglong : TypeAlias = py_int
const_ulonglong : TypeAlias = const[py_int]
p_const_ulonglong = pointer[const[py_int]]
pp_const_ulonglong = pointer[pointer[const[py_int]]]
ppp_const_ulonglong = pointer[pointer[pointer[const[py_int]]]]
p_ulonglong = pointer[py_int]
pp_ulonglong = pointer[pointer[py_int]]
ppp_ulonglong = pointer[pointer[pointer[py_int]]]
const_Py_UNICODE : TypeAlias = const[py_int]
p_const_Py_UNICODE = pointer[const[py_int]]
pp_const_Py_UNICODE = pointer[pointer[const[py_int]]]
ppp_const_Py_UNICODE = pointer[pointer[pointer[const[py_int]]]]
p_Py_UNICODE = pointer[py_int]
pp_Py_UNICODE = pointer[pointer[py_int]]
ppp_Py_UNICODE = pointer[pointer[pointer[py_int]]]
ushort : TypeAlias = py_int
const_ushort : TypeAlias = const[py_int]
p_const_ushort = pointer[const[py_int]]
pp_const_ushort = pointer[pointer[const[py_int]]]
ppp_const_ushort = pointer[pointer[pointer[const[py_int]]]]
p_ushort = pointer[py_int]
pp_ushort = pointer[pointer[py_int]]
ppp_ushort = pointer[pointer[pointer[py_int]]]
const_void : TypeAlias = const[Any]
p_void = pointer[Any]
pp_void = pointer[pointer[Any]]
ppp_void = pointer[pointer[pointer[Any]]]

##### END: GENERATED LIST OF GENERATED TYPES #####
