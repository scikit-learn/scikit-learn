from __future__ import annotations

from typing import Any, Sequence, Protocol, Union, TypeVar
from numpy import (
    ndarray,
    dtype,
    generic,
    bool_,
    unsignedinteger,
    integer,
    floating,
    complexfloating,
    number,
    timedelta64,
    datetime64,
    object_,
    void,
    str_,
    bytes_,
)
from ._nested_sequence import _NestedSequence

_T = TypeVar("_T")
_ScalarType = TypeVar("_ScalarType", bound=generic)
_DType = TypeVar("_DType", bound="dtype[Any]")
_DType_co = TypeVar("_DType_co", covariant=True, bound="dtype[Any]")

# The `_SupportsArray` protocol only cares about the default dtype
# (i.e. `dtype=None` or no `dtype` parameter at all) of the to-be returned
# array.
# Concrete implementations of the protocol are responsible for adding
# any and all remaining overloads
class _SupportsArray(Protocol[_DType_co]):
    def __array__(self) -> ndarray[Any, _DType_co]: ...


# TODO: Wait until mypy supports recursive objects in combination with typevars
_FiniteNestedSequence = Union[
    _T,
    Sequence[_T],
    Sequence[Sequence[_T]],
    Sequence[Sequence[Sequence[_T]]],
    Sequence[Sequence[Sequence[Sequence[_T]]]],
]

# A union representing array-like objects; consists of two typevars:
# One representing types that can be parametrized w.r.t. `np.dtype`
# and another one for the rest
_ArrayLike = Union[
    _SupportsArray[_DType],
    _NestedSequence[_SupportsArray[_DType]],
    _T,
    _NestedSequence[_T],
]

# TODO: support buffer protocols once
#
# https://bugs.python.org/issue27501
#
# is resolved. See also the mypy issue:
#
# https://github.com/python/typing/issues/593
ArrayLike = _ArrayLike[
    dtype,
    Union[bool, int, float, complex, str, bytes],
]

# `ArrayLike<X>_co`: array-like objects that can be coerced into `X`
# given the casting rules `same_kind`
_ArrayLikeBool_co = _ArrayLike[
    "dtype[bool_]",
    bool,
]
_ArrayLikeUInt_co = _ArrayLike[
    "dtype[Union[bool_, unsignedinteger[Any]]]",
    bool,
]
_ArrayLikeInt_co = _ArrayLike[
    "dtype[Union[bool_, integer[Any]]]",
    Union[bool, int],
]
_ArrayLikeFloat_co = _ArrayLike[
    "dtype[Union[bool_, integer[Any], floating[Any]]]",
    Union[bool, int, float],
]
_ArrayLikeComplex_co = _ArrayLike[
    "dtype[Union[bool_, integer[Any], floating[Any], complexfloating[Any, Any]]]",
    Union[bool, int, float, complex],
]
_ArrayLikeNumber_co = _ArrayLike[
    "dtype[Union[bool_, number[Any]]]",
    Union[bool, int, float, complex],
]
_ArrayLikeTD64_co = _ArrayLike[
    "dtype[Union[bool_, integer[Any], timedelta64]]",
    Union[bool, int],
]
_ArrayLikeDT64_co = Union[
    _SupportsArray["dtype[datetime64]"],
    _NestedSequence[_SupportsArray["dtype[datetime64]"]],
]
_ArrayLikeObject_co = Union[
    _SupportsArray["dtype[object_]"],
    _NestedSequence[_SupportsArray["dtype[object_]"]],
]

_ArrayLikeVoid_co = Union[
    _SupportsArray["dtype[void]"],
    _NestedSequence[_SupportsArray["dtype[void]"]],
]
_ArrayLikeStr_co = _ArrayLike[
    "dtype[str_]",
    str,
]
_ArrayLikeBytes_co = _ArrayLike[
    "dtype[bytes_]",
    bytes,
]

_ArrayLikeInt = _ArrayLike[
    "dtype[integer[Any]]",
    int,
]
