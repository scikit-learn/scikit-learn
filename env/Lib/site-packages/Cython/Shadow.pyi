from builtins import (int as py_int, float as py_float,
                      bool as py_bool, str as py_str, complex as py_complex)
from typing import (Union, Dict, Any, Sequence, Optional,
                    List, TypeVar, Type, Generic)

int = py_int
long = py_int
longlong = py_int
short = py_int
char = py_int
sint = py_int
slong = py_int
slonglong = py_int
sshort = py_int
schar = py_int
uint = py_int
ulong = py_int
ulonglong = py_int
ushort = py_int
uchar = py_int
size_t = py_int
Py_ssize_t = py_int
Py_UCS4 = Union[py_int, str]
Py_UNICODE = Union[py_int, str]
float = py_float
double = py_float
longdouble = py_float
complex = py_complex
floatcomplex = py_complex
doublecomplex = py_complex
longdoublecomplex = py_complex
bint = py_bool
void = Union[None]
basestring = py_str
unicode = py_str

gs: Dict[str, Any]  # Should match the return type of globals()

_T = TypeVar('_T')

class _ArrayType(object, Generic[_T]):
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
        value: Optional[Union[ArrayType[_T], PointerType[_T], List[_T], int]] = ...
    ) -> None: ...
    def __getitem__(self, ix: int) -> _T: ...
    def __setitem__(self, ix: int, value: _T) -> None: ...
    def __eq__(self, value: object) -> bool: ...
    def __repr__(self) -> str: ...

class ArrayType(PointerType[_T]):
    def __init__(self) -> None: ...

#class StructType(CythonType, Generic[_T]):
#    def __init__(
#        self,
#        value: List[Type[_T]] = ...
#    ) -> None: ...

def index_type(
    base_type: _T, item: Union[tuple, slice, int]) -> _ArrayType[_T]: ...

def pointer(basetype: _T) -> Type[PointerType[_T]]: ...

def array(basetype: _T, n: int) -> Type[ArrayType[_T]]: ...

#def struct(basetype: _T) -> Type[StructType[_T]]: ...

class typedef(CythonType, Generic[_T]):
    name: str

    def __init__(self, type: _T, name: Optional[str] = ...) -> None: ...
    def __call__(self, *arg: Any) -> _T: ...
    def __repr__(self) -> str: ...
    __getitem__ = index_type

#class _FusedType(CythonType, Generic[_T]):
#    def __init__(self) -> None: ...

#def fused_type(*args: Tuple[_T]) -> Type[FusedType[_T]]: ...
