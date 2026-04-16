""" This modules provides the translation tables from python to c++. """

import gast as ast
from importlib import import_module
import inspect
import logging
import numpy
import sys

from pythran.typing import Dict, Set, List, TypeVar, Union, Optional, NDArray
from pythran.typing import Generator, Fun, Tuple, Iterable, Sized, File

from pythran.conversion import to_ast, ToNotEval
from pythran.intrinsic import Class
from pythran.intrinsic import ClassWithConstConstructor, ExceptionClass
from pythran.intrinsic import ClassWithReadOnceConstructor
from pythran.intrinsic import ConstFunctionIntr, FunctionIntr, UpdateEffect
from pythran.intrinsic import ConstMethodIntr, MethodIntr
from pythran.intrinsic import AttributeIntr, StaticAttributeIntr
from pythran.intrinsic import ReadEffect, ConstantIntr, UFunc
from pythran.intrinsic import ReadOnceMethodIntr
from pythran.intrinsic import ReadOnceFunctionIntr, ConstExceptionIntr
from pythran import interval


logger = logging.getLogger("pythran")

pythran_ward = '__pythran_'

namespace = "pythonic"


cxx_keywords = {
    'and', 'and_eq', 'asm', 'auto', 'bitand', 'bitor',
    'bool', 'break', 'case', 'catch', 'char', 'class',
    'compl', 'const', 'const_cast', 'continue', 'default', 'delete',
    'do', 'double', 'dynamic_cast', 'else', 'enum', 'explicit',
    'export', 'extern', 'false', 'float', 'for', 'friend',
    'goto', 'if', 'inline', 'int', 'long', 'mutable', 'namespace', 'new',
    'not', 'not_eq', 'operator', 'or', 'or_eq', 'private', 'protected',
    'public', 'register', 'reinterpret_cast', 'return', 'short', 'signed',
    'sizeof', 'static', 'static_cast',
    'struct', 'switch', 'template', 'this', 'throw', 'true',
    'try', 'typedef', 'typeid', 'typename', 'union', 'unsigned',
    'using', 'virtual', 'void', 'volatile', 'wchar_t', 'while',
    'xor', 'xor_eq',
    # C++11 additions
    'constexpr', 'decltype', 'thread_local', 'noexcept', 'nullptr', 'static_assert',
    'alignof', 'alignas',
    # C++14 additions (nothing)
    # reserved namespaces
    'std',
}


def make_lazy(exp):
    return '[&] () {{ return {0}; }}'.format(exp)


def make_and(x, y):
    lx, ly = make_lazy(x), make_lazy(y)
    return 'pythonic::builtins::pythran::and_({0}, {1})'.format(lx, ly)


def make_or(x, y):
    lx, ly = make_lazy(x), make_lazy(y)
    return 'pythonic::builtins::pythran::or_({0}, {1})'.format(lx, ly)


operator_to_lambda = {
    # boolop
    ast.And: make_and,
    ast.Or: make_or,
    # operator
    ast.Add: "pythonic::operator_::add({0}, {1})".format,
    ast.Sub: "pythonic::operator_::sub({0}, {1})".format,
    ast.Mult: "pythonic::operator_::mul({0}, {1})".format,
    ast.Div: "pythonic::operator_::div({0}, {1})".format,
    ast.Mod: "pythonic::operator_::mod({0}, {1})".format,
    ast.Pow: "pythonic::builtins::pow({0}, {1})".format,
    ast.LShift: "pythonic::operator_::lshift({0}, {1})".format,
    ast.RShift: "pythonic::operator_::rshift({0}, {1})".format,
    ast.BitOr: "pythonic::operator_::or_({0}, {1})".format,
    ast.BitXor: "pythonic::operator_::xor_({0}, {1})".format,
    ast.BitAnd: "pythonic::operator_::and_({0}, {1})".format,
    ast.MatMult: "pythonic::operator_::functor::matmul()({0}, {1})".format,
    ast.FloorDiv: "pythonic::operator_::functor::floordiv()({0}, {1})".format,
    # unaryop
    ast.Invert: "pythonic::operator_::invert({0})".format,
    ast.Not: "pythonic::operator_::not_({0})".format,
    ast.UAdd: "pythonic::operator_::pos({0})".format,
    ast.USub: "pythonic::operator_::neg({0})".format,
    # cmpop
    ast.Eq: "pythonic::operator_::eq({0}, {1})".format,
    ast.NotEq: "pythonic::operator_::ne({0}, {1})".format,
    ast.Lt: "pythonic::operator_::lt({0}, {1})".format,
    ast.LtE: "pythonic::operator_::le({0}, {1})".format,
    ast.Gt: "pythonic::operator_::gt({0}, {1})".format,
    ast.GtE: "pythonic::operator_::ge({0}, {1})".format,
    ast.Is: "pythonic::operator_::is_({0}, {1})".format,
    ast.IsNot: ("pythonic::operator_::is_not({0}, {1})").format,
    ast.In: "pythonic::operator_::contains({1}, {0})".format,
    ast.NotIn: "(!pythonic::operator_::contains({1}, {0}))".format,
}

update_operator_to_lambda = {
    # operator
    ast.Add: "({0} += {1})".format,
    ast.Sub: "({0} -= {1})".format,
    ast.Mult: "({0} *= {1})".format,
    ast.Div: "(pythonic::operator_::idiv({0}, {1}))".format,
    ast.Mod: "(pythonic::operator_::imod({0}, {1}))".format,
    ast.Pow: "(pythonic::operator_::ipow({0}, {1}))".format,
    ast.LShift: "({0} <<= {1})".format,
    ast.RShift: "({0} >>= {1})".format,
    ast.BitOr: "({0} |= {1})".format,
    ast.BitXor: "({0} ^= {1})".format,
    ast.BitAnd: "({0} &= {1})".format,
    ast.MatMult: "(pythonic::operator_::imatmul({0}, {1}))".format,
    ast.FloorDiv:
        "(pythonic::operator_::functor::ifloordiv{{}}({0}, {1}))".format,
}

T0, T1, T2, T3 = TypeVar('T0'), TypeVar('T1'), TypeVar('T2'), TypeVar('T3')
T4, T5, T6, T7 = TypeVar('T4'), TypeVar('T5'), TypeVar('T6'), TypeVar('T7')


_bool_signature = Union[
    Fun[[], bool],
    Fun[[T0], bool]
]

_int_signature = Union[
    Fun[[], int],
    Fun[[bool], int],
    Fun[[int], int],
    Fun[[float], int],
    Fun[[str], int],
]

_float_signature = Union[
    Fun[[], float],
    Fun[[str], float],
    Fun[[float], float],
]

_complex_signature = Union[
    Fun[[float], complex],
    Fun[[float, float], complex],
]

# workaround changes in numpy interaction with getfullargspec
try:
    inspect.getfullargspec(numpy.asarray)
    # if we have a description, honor it
    extra_numpy_asarray_descr = {}
except TypeError:
    extra_numpy_asarray_descr = {'args':('a', 'dtype'),
                                 'defaults': (None,)}





def update_effects(self, node):
    """
    Combiner when we update the first argument of a function.

    It turn type of first parameter in combination of all others
    parameters types.
    """
    return [self.combine(node.args[0], None, node_args_k)
            for node_args_k in node.args[1:]]


BINARY_UFUNC = {"accumulate": FunctionIntr()}
REDUCED_BINARY_UFUNC = {"accumulate": FunctionIntr(),
                        "reduce": ConstFunctionIntr()}

CLASSES = {
    "dtype": {
        "type": MethodIntr(),
    },
    "array": {
        # array have fixed type, no need for signature
        "append": MethodIntr(),
        "buffer_info": ConstMethodIntr(),
        "byteswap": MethodIntr(),
        "count": ConstMethodIntr(),
        "extend": MethodIntr(),
        "fromfile": MethodIntr(),
        "fromlist": MethodIntr(),
        "frombytes": MethodIntr(),
        "insert": MethodIntr(),
        "pop": MethodIntr(),
        "remove": MethodIntr(),
        "reverse": MethodIntr(),
    },
    "list": {
        "append": MethodIntr(signature=Fun[[List[T0], T0], None]),
        "extend": MethodIntr(update_effects),
        "pop": MethodIntr(
            signature=Union[
                Fun[[List[T0]], T0],
                Fun[[List[T0], int], T0],
            ],
        ),
        "reverse": MethodIntr(signature=Fun[[List[T0]], None]),
        "sort": MethodIntr(
            args=("self", "key",),
        ),
        "count": ConstMethodIntr(signature=Fun[[List[T0], T0], int]),
        "remove": MethodIntr(signature=Fun[[List[T0], T0], None]),
        "insert": MethodIntr(signature=Fun[[List[T0], int, T0], None]),
        "clear": MethodIntr(signature=Fun[[List[T0]], None]),
    },
    "slice": {
        "start": AttributeIntr(signature=Fun[[T0], int]),
        "stop": AttributeIntr(signature=Fun[[T0], int]),
        "step": AttributeIntr(signature=Fun[[T0], int]),
    },
    "str": {
        "__mod__": ConstMethodIntr(
            signature=Union[
                Fun[[str, T0], str],
                Fun[[str, T0, T1], str],
                Fun[[str, T0, T1, T2], str],
                Fun[[str, T0, T1, T2, T3, T4], str],
                Fun[[str, T0, T1, T2, T3, T4, T5], str],
                Fun[[str, T0, T1, T2, T3, T4, T5, T6], str],
            ],
        ),

        "capitalize": ConstMethodIntr(signature=Fun[[str], str]),
        "count": ConstMethodIntr(signature=Union[
            Fun[[str, str], int],
            Fun[[str, str, int], int],
            Fun[[str, str, int, int], int],
        ]),
        "endswith": ConstMethodIntr(
            signature=Union[
                Fun[[str, str], bool],
                Fun[[str, str, Optional[int]], bool],
                Fun[[str, str, Optional[int], Optional[int]], bool],
            ],
        ),
        "startswith": ConstMethodIntr(
            signature=Union[
                Fun[[str, str], bool],
                Fun[[str, str, Optional[int]], bool],
                Fun[[str, str, Optional[int], Optional[int]], bool],
            ],
        ),
        "find": ConstMethodIntr(
            signature=Union[
                Fun[[str, str], int],
                Fun[[str, str, Optional[int]], int],
                Fun[[str, str, Optional[int], Optional[int]], int],
            ],
        ),
        "isalpha": ConstMethodIntr(signature=Fun[[str], bool]),
        "isdigit": ConstMethodIntr(signature=Fun[[str], bool]),
        "join": ConstMethodIntr(signature=Fun[[str, Iterable[str]], str]),
        "lower": ConstMethodIntr(signature=Fun[[str], str]),
        "replace": ConstMethodIntr(
            signature=Union[
                Fun[[str, str, str], str],
                Fun[[str, str, str, int], str],
            ]
        ),
        "split": ConstMethodIntr(
            signature=Union[
                Fun[[str], List[str]],
                Fun[[str, str], List[str]],
                Fun[[str, None], List[str]],
                Fun[[str, str, int], List[str]],
                Fun[[str, None, int], List[str]],
            ]
        ),
        "strip": ConstMethodIntr(
            signature=Union[
                Fun[[str], str],
                Fun[[str, str], str],
            ]
        ),
        "lstrip": ConstMethodIntr(
            signature=Union[
                Fun[[str], str],
                Fun[[str, str], str],
            ]
        ),
        "rstrip": ConstMethodIntr(
            signature=Union[
                Fun[[str], str],
                Fun[[str, str], str],
            ]
        ),
        "upper": ConstMethodIntr(
            signature=Fun[[str], str]
        ),
    },
    "set": {
        "add": MethodIntr(signature=Fun[[Set[T0], T0], None]),
        "clear": MethodIntr(signature=Fun[[Set[T0]], None]),
        "copy": ConstMethodIntr(signature=Fun[[Set[T0]], Iterable[T0]]),
        "discard": MethodIntr(signature=Fun[[Set[T0], T0], None]),
        "remove": MethodIntr(signature=Fun[[Set[T0], T0], None]),
        "isdisjoint": ConstMethodIntr(
            signature=Fun[[Set[T0], Set[T0]], bool]),
        "union": ConstMethodIntr(
            signature=Union[
                Fun[[Set[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0], Iterable[T0]],
                    Set[T0]],
            ]
        ),
        "update": MethodIntr(update_effects),
        "intersection": ConstMethodIntr(
            signature=Union[
                Fun[[Set[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0], Iterable[T0]],
                    Set[T0]],
            ]
        ),
        "intersection_update": MethodIntr(update_effects),
        "difference": ConstMethodIntr(
            signature=Union[
                Fun[[Set[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0], Iterable[T0]],
                    Set[T0]],
            ]
        ),
        "difference_update": MethodIntr(update_effects),
        "symmetric_difference": ConstMethodIntr(
            signature=Union[
                Fun[[Set[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0]], Set[T0]],
                Fun[[Set[T0], Iterable[T0], Iterable[T0], Iterable[T0]],
                    Set[T0]],
            ]
        ),
        "symmetric_difference_update": MethodIntr(update_effects),
        "issuperset": ConstMethodIntr(
            signature=Fun[[Set[T0], Set[T0]], bool]),
        "issubset": ConstMethodIntr(
            signature=Fun[[Set[T0], Set[T0]], bool]),
    },
    "Exception": {
        "args": AttributeIntr(signature=Fun[[T0], str]),
        "errno": AttributeIntr(signature=Fun[[T0], str]),
        "strerror": AttributeIntr(signature=Fun[[T0], str]),
        "filename": AttributeIntr(signature=Fun[[T0], str]),
    },
    "float": {
        "is_integer": ConstMethodIntr(signature=Fun[[float], bool]),
    },
    "complex": {
        "conjugate": ConstMethodIntr(),
        "real": AttributeIntr(
            signature=Union[
                Fun[[complex], float],
                Fun[[NDArray[complex, :]], NDArray[float, :]],
                Fun[[NDArray[complex, :, :]], NDArray[float, :, :]],
                Fun[[NDArray[complex, :, :, :]], NDArray[float, :, :, :]],
                Fun[[NDArray[complex, :, :, :, :]],
                    NDArray[float, :, :, :, :]],
            ]
        ),
        "imag": AttributeIntr(
            signature=Union[
                Fun[[complex], float],
                Fun[[NDArray[complex, :]], NDArray[float, :]],
                Fun[[NDArray[complex, :, :]], NDArray[float, :, :]],
                Fun[[NDArray[complex, :, :, :]], NDArray[float, :, :, :]],
                Fun[[NDArray[complex, :, :, :, :]],
                    NDArray[float, :, :, :, :]],
            ]
        ),
    },
    "dict": {
        "fromkeys": ConstFunctionIntr(
            signature=Union[
                Fun[[Iterable[T0]], Dict[T0, Optional[T1]]],
                Fun[[Iterable[T0], T1], Dict[T0, T1]],
            ],
        ),
        "clear": MethodIntr(signature=Fun[[Dict[T0, T1]], None]),
        "copy": ConstMethodIntr(
            signature=Fun[[Dict[T0, T1]], Dict[T0, T1]]),
        "get": ConstMethodIntr(
            signature=Union[
                Fun[[Dict[T0, T1], T0], Optional[T1]],
                Fun[[Dict[T0, T1], T0, T1], T1],
            ],
        ),
        "items": ConstMethodIntr(
            signature=Fun[[Dict[T0, T1]], List[Tuple[T0, T1]]]),
        "keys": ConstMethodIntr(signature=Fun[[Dict[T0, T1]], List[T0]]),
        "pop": MethodIntr(
            signature=Union[
                Fun[[Dict[T0, T1], T0], T1],
                Fun[[Dict[T0, T1], T0, T1], T1],
            ]
        ),
        "popitem": MethodIntr(
            signature=Fun[[Dict[T0, T1]], Tuple[T0, T1]]),
        "setdefault": MethodIntr(
            signature=Union[
                Fun[[Dict[T0, T1], T0, T1], T1],
                Fun[[Dict[T0, T1], T0], T1]
            ],
            return_alias=lambda args: {
                ast.Subscript(args[0], args[1], ast.Load())
            }.union({args[2]} if len(args) == 3 else set())
        ),
        "update": MethodIntr(update_effects),
        "values": ConstMethodIntr(signature=Fun[[Dict[T0, T1]], List[T1]]),
    },
    "file": {
        # Member variables
        "closed": AttributeIntr(signature=Fun[[File], bool]),
        "mode": AttributeIntr(signature=Fun[[File], str]),
        "name": AttributeIntr(signature=Fun[[File], str]),
        "newlines": AttributeIntr(signature=Fun[[File], str]),
        # Member functions
        "close": MethodIntr(
            signature=Fun[[File], None],
            global_effects=True
        ),
        "flush": MethodIntr(
            signature=Fun[[File], None],
            global_effects=True
        ),
        "fileno": MethodIntr(
            signature=Fun[[File], int],
        ),
        "isatty": MethodIntr(signature=Fun[[File], bool]),
        "next": MethodIntr(global_effects=True),
        "read": MethodIntr(
            signature=Union[
                Fun[[File], str],
                Fun[[File, int], str],
            ],
            global_effects=True
        ),
        "readline": MethodIntr(
            signature=Union[
                Fun[[File], str],
                Fun[[File, int], str],
            ],
            global_effects=True
        ),
        "readlines": MethodIntr(
            signature=Union[
                Fun[[File], List[str]],
                Fun[[File, int], List[str]],
            ],
            global_effects=True
        ),
        "seek": MethodIntr(
            signature=Union[
                Fun[[File, int], None],
                Fun[[File, int, int], None],
            ],
            global_effects=True
        ),
        "tell": MethodIntr(signature=Fun[[File], int]),
        "truncate": MethodIntr(
            signature=Union[
                Fun[[File], None],
                Fun[[File, int], None],
            ],
            global_effects=True
        ),
        "write": MethodIntr(
            signature=Fun[[File, str], None],
            global_effects=True
        ),
        "writelines": MethodIntr(
            signature=Fun[[File, Iterable[str]], None],
            global_effects=True
        ),
    },
    "finfo": {
        "eps": AttributeIntr(signature=float),
    },
    "ndarray": {
        "astype": MethodIntr(
            signature=Union[
                # dtype = bool
                Fun[[NDArray[bool, :], _bool_signature], NDArray[bool, :]],
                Fun[[NDArray[int, :], _bool_signature], NDArray[bool, :]],
                Fun[[NDArray[float, :], _bool_signature], NDArray[bool, :]],
                Fun[[NDArray[complex, :], _bool_signature], NDArray[bool, :]],
                Fun[[NDArray[bool, :, :], _bool_signature],
                    NDArray[bool, :, :]],
                Fun[[NDArray[int, :, :], _bool_signature],
                    NDArray[bool, :, :]],
                Fun[[NDArray[float, :, :], _bool_signature],
                    NDArray[bool, :, :]],
                Fun[[NDArray[complex, :, :], _bool_signature],
                    NDArray[bool, :, :]],
                Fun[[NDArray[bool, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :]],
                Fun[[NDArray[int, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :]],
                Fun[[NDArray[float, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :]],
                Fun[[NDArray[complex, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :]],
                Fun[[NDArray[bool, :, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :, :]],
                Fun[[NDArray[int, :, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :, :]],
                Fun[[NDArray[float, :, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :, :]],
                Fun[[NDArray[complex, :, :, :, :], _bool_signature],
                    NDArray[bool, :, :, :, :]],
                # dtype = int
                Fun[[NDArray[bool, :], _int_signature], NDArray[int, :]],
                Fun[[NDArray[int, :], _int_signature], NDArray[int, :]],
                Fun[[NDArray[float, :], _int_signature], NDArray[int, :]],
                Fun[[NDArray[complex, :], _int_signature], NDArray[int, :]],
                Fun[[NDArray[bool, :, :], _int_signature], NDArray[int, :, :]],
                Fun[[NDArray[int, :, :], _int_signature], NDArray[int, :, :]],
                Fun[[NDArray[float, :, :], _int_signature],
                    NDArray[int, :, :]],
                Fun[[NDArray[complex, :, :], _int_signature],
                    NDArray[int, :, :]],
                Fun[[NDArray[bool, :, :, :], _int_signature],
                    NDArray[int, :, :, :]],
                Fun[[NDArray[int, :, :, :], _int_signature],
                    NDArray[int, :, :, :]],
                Fun[[NDArray[float, :, :, :], _int_signature],
                    NDArray[int, :, :, :]],
                Fun[[NDArray[complex, :, :, :], _int_signature],
                    NDArray[int, :, :, :]],
                Fun[[NDArray[bool, :, :, :, :], _int_signature],
                    NDArray[int, :, :, :, :]],
                Fun[[NDArray[int, :, :, :, :], _int_signature],
                    NDArray[int, :, :, :, :]],
                Fun[[NDArray[float, :, :, :, :], _int_signature],
                    NDArray[int, :, :, :, :]],
                Fun[[NDArray[complex, :, :, :, :], _int_signature],
                    NDArray[int, :, :, :, :]],
                # dtype = float
                Fun[[NDArray[bool, :], _float_signature], NDArray[float, :]],
                Fun[[NDArray[int, :], _float_signature], NDArray[float, :]],
                Fun[[NDArray[float, :], _float_signature], NDArray[float, :]],
                Fun[[NDArray[complex, :], _float_signature],
                    NDArray[float, :]],
                Fun[[NDArray[bool, :, :], _float_signature],
                    NDArray[float, :, :]],
                Fun[[NDArray[int, :, :], _float_signature],
                    NDArray[float, :, :]],
                Fun[[NDArray[float, :, :], _float_signature],
                    NDArray[float, :, :]],
                Fun[[NDArray[complex, :, :], _float_signature],
                    NDArray[float, :, :]],
                Fun[[NDArray[bool, :, :, :], _float_signature],
                    NDArray[float, :, :, :]],
                Fun[[NDArray[int, :, :, :], _float_signature],
                    NDArray[float, :, :, :]],
                Fun[[NDArray[float, :, :, :], _float_signature],
                    NDArray[float, :, :, :]],
                Fun[[NDArray[complex, :, :, :], _float_signature],
                    NDArray[float, :, :, :]],
                Fun[[NDArray[bool, :, :, :, :], _float_signature],
                    NDArray[float, :, :, :, :]],
                Fun[[NDArray[int, :, :, :, :], _float_signature],
                    NDArray[float, :, :, :, :]],
                Fun[[NDArray[float, :, :, :, :], _float_signature],
                    NDArray[float, :, :, :, :]],
                Fun[[NDArray[complex, :, :, :, :], _float_signature],
                    NDArray[float, :, :, :, :]],
                # dtype = complex
                Fun[[NDArray[bool, :], _complex_signature],
                    NDArray[complex, :]],
                Fun[[NDArray[int, :], _complex_signature],
                    NDArray[complex, :]],
                Fun[[NDArray[float, :], _complex_signature],
                    NDArray[complex, :]],
                Fun[[NDArray[complex, :], _complex_signature],
                    NDArray[complex, :]],
                Fun[[NDArray[bool, :, :], _complex_signature],
                    NDArray[complex, :, :]],
                Fun[[NDArray[int, :, :], _complex_signature],
                    NDArray[complex, :, :]],
                Fun[[NDArray[float, :, :], _complex_signature],
                    NDArray[complex, :, :]],
                Fun[[NDArray[complex, :, :], _complex_signature],
                    NDArray[complex, :, :]],
                Fun[[NDArray[bool, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :]],
                Fun[[NDArray[int, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :]],
                Fun[[NDArray[float, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :]],
                Fun[[NDArray[complex, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :]],
                Fun[[NDArray[bool, :, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :, :]],
                Fun[[NDArray[int, :, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :, :]],
                Fun[[NDArray[float, :, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :, :]],
                Fun[[NDArray[complex, :, :, :, :], _complex_signature],
                    NDArray[complex, :, :, :, :]],
            ]
        ),
        "dtype": AttributeIntr(),
        "fill": MethodIntr(
            signature=Union[
                # 1d
                Fun[[NDArray[bool, :], bool], None],
                Fun[[NDArray[int, :], int], None],
                Fun[[NDArray[float, :], float], None],
                Fun[[NDArray[complex, :], complex], None],
                # 2d
                Fun[[NDArray[bool, :, :], bool], None],
                Fun[[NDArray[int, :, :], int], None],
                Fun[[NDArray[float, :, :], float], None],
                Fun[[NDArray[complex, :, :], complex], None],
                # 3d
                Fun[[NDArray[bool, :, :, :], bool], None],
                Fun[[NDArray[int, :, :, :], int], None],
                Fun[[NDArray[float, :, :, :], float], None],
                Fun[[NDArray[complex, :, :, :], complex], None],
                # 4d
                Fun[[NDArray[bool, :, :, :, :], bool], None],
                Fun[[NDArray[int, :, :, :, :], int], None],
                Fun[[NDArray[float, :, :, :, :], float], None],
                Fun[[NDArray[complex, :, :, :, :], complex], None],

            ],
        ),
        "flat": AttributeIntr(
            signature=Union[
                # 1d
                Fun[[NDArray[bool, :]], Generator[bool]],
                Fun[[NDArray[int, :]], Generator[int]],
                Fun[[NDArray[float, :]], Generator[float]],
                Fun[[NDArray[complex, :]], Generator[complex]],
                # 2d
                Fun[[NDArray[bool, :, :]], Generator[bool]],
                Fun[[NDArray[int, :, :]], Generator[int]],
                Fun[[NDArray[float, :, :]], Generator[float]],
                Fun[[NDArray[complex, :, :]], Generator[complex]],
                # 3d
                Fun[[NDArray[bool, :, :, :]], Generator[bool]],
                Fun[[NDArray[int, :, :, :]], Generator[int]],
                Fun[[NDArray[float, :, :, :]], Generator[float]],
                Fun[[NDArray[complex, :, :, :]], Generator[complex]],
                # 4d
                Fun[[NDArray[bool, :, :, :, :]], Generator[bool]],
                Fun[[NDArray[int, :, :, :, :]], Generator[int]],
                Fun[[NDArray[float, :, :, :, :]], Generator[float]],
                Fun[[NDArray[complex, :, :, :, :]], Generator[complex]],
            ]
        ),
        "flatten": MethodIntr(
            signature=Union[
                # 1d
                Fun[[NDArray[bool, :]], NDArray[bool, :]],
                Fun[[NDArray[int, :]], NDArray[int, :]],
                Fun[[NDArray[float, :]], NDArray[float, :]],
                Fun[[NDArray[complex, :]], NDArray[complex, :]],
                # 2d
                Fun[[NDArray[bool, :, :]], NDArray[bool, :]],
                Fun[[NDArray[int, :, :]], NDArray[int, :]],
                Fun[[NDArray[float, :, :]], NDArray[float, :]],
                Fun[[NDArray[complex, :, :]], NDArray[complex, :]],
                # 3d
                Fun[[NDArray[bool, :, :, :]], NDArray[bool, :]],
                Fun[[NDArray[int, :, :, :]], NDArray[int, :]],
                Fun[[NDArray[float, :, :, :]], NDArray[float, :]],
                Fun[[NDArray[complex, :, :, :]], NDArray[complex, :]],
                # 4d
                Fun[[NDArray[bool, :, :, :, :]], NDArray[bool, :]],
                Fun[[NDArray[int, :, :, :, :]], NDArray[int, :]],
                Fun[[NDArray[float, :, :, :, :]], NDArray[float, :]],
                Fun[[NDArray[complex, :, :, :, :]], NDArray[complex, :]],
            ]
        ),
        "item": MethodIntr(
            signature=Union[
                # item = int
                # 1d
                Fun[[NDArray[bool, :], int], bool],
                Fun[[NDArray[int, :], int], int],
                Fun[[NDArray[float, :], int], float],
                Fun[[NDArray[complex, :], int], complex],
                # 2d
                Fun[[NDArray[bool, :, :], int], bool],
                Fun[[NDArray[int, :, :], int], int],
                Fun[[NDArray[float, :, :], int], float],
                Fun[[NDArray[complex, :, :], int], complex],
                # 3d
                Fun[[NDArray[bool, :, :, :], int], bool],
                Fun[[NDArray[int, :, :, :], int], int],
                Fun[[NDArray[float, :, :, :], int], float],
                Fun[[NDArray[complex, :, :, :], int], complex],
                # 4d
                Fun[[NDArray[bool, :, :, :, :], int], bool],
                Fun[[NDArray[int, :, :, :, :], int], int],
                Fun[[NDArray[float, :, :, :, :], int], float],
                Fun[[NDArray[complex, :, :, :, :], int], complex],

                # item = tuple
                # 1d
                Fun[[NDArray[bool, :], Tuple[int]], bool],
                Fun[[NDArray[int, :], Tuple[int]], int],
                Fun[[NDArray[float, :], Tuple[int]], float],
                Fun[[NDArray[complex, :], Tuple[int]], complex],
                # 2d
                Fun[[NDArray[bool, :, :], Tuple[int, int]], bool],
                Fun[[NDArray[int, :, :], Tuple[int, int]], int],
                Fun[[NDArray[float, :, :], Tuple[int, int]], float],
                Fun[[NDArray[complex, :, :], Tuple[int, int]], complex],
                # 3d
                Fun[[NDArray[bool, :, :, :], Tuple[int, int, int]], bool],
                Fun[[NDArray[int, :, :, :], Tuple[int, int, int]], int],
                Fun[[NDArray[float, :, :, :], Tuple[int, int, int]], float],
                Fun[[NDArray[complex, :, :, :], Tuple[int, int, int]],
                    complex],
                # 4d
                Fun[[NDArray[bool, :, :, :, :], Tuple[int, int, int, int]],
                    bool],
                Fun[[NDArray[int, :, :, :, :], Tuple[int, int, int, int]],
                    int],
                Fun[[NDArray[float, :, :, :, :], Tuple[int, int, int, int]],
                    float],
                Fun[[NDArray[complex, :, :, :, :], Tuple[int, int, int, int]],
                    complex],
            ]
        ),
        "itemsize": StaticAttributeIntr(signature=Fun[[NDArray[T0, :]], int],
                                  return_range=interval.positive_values),
        "nbytes": AttributeIntr(
            signature=Fun[[NDArray[T0, :]], int],
            return_range=interval.positive_values
        ),
        "ndim": StaticAttributeIntr(signature=Fun[[NDArray[T0, :]], int],
                              return_range=interval.positive_values),
        "reshape": ConstMethodIntr(
            signature=Union[
                Fun[[NDArray[T0, :], int], NDArray[T1, :]],
                Fun[[NDArray[T0, :], Tuple[int]], NDArray[T1, :]],
                Fun[[NDArray[T0, :], int, int], NDArray[T1, :, :]],
                Fun[[NDArray[T0, :], Tuple[int, int]], NDArray[T1, :, :]],
                Fun[[NDArray[T0, :], int, int, int],
                    NDArray[T1, :, :, :]],
                Fun[[NDArray[T0, :], Tuple[int, int, int]],
                    NDArray[T1, :, :, :]],
                Fun[[NDArray[T0, :], int, int, int, int],
                    NDArray[T1, :, :, :, :]],
                Fun[[NDArray[T0, :], Tuple[int, int, int, int]],
                    NDArray[T1, :, :, :, :]],
            ]
        ),
        "shape": AttributeIntr(
            signature=Union[
                # bool
                Fun[[NDArray[bool, :]], Tuple[int]],
                Fun[[NDArray[bool, :, :]], Tuple[int, int]],
                Fun[[NDArray[bool, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[bool, :, :, :, :]], Tuple[int, int, int, int]],
                # int
                Fun[[NDArray[int, :]], Tuple[int]],
                Fun[[NDArray[int, :, :]], Tuple[int, int]],
                Fun[[NDArray[int, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[int, :, :, :, :]], Tuple[int, int, int, int]],
                # float
                Fun[[NDArray[float, :]], Tuple[int]],
                Fun[[NDArray[float, :, :]], Tuple[int, int]],
                Fun[[NDArray[float, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[float, :, :, :, :]], Tuple[int, int, int, int]],
                # complex
                Fun[[NDArray[complex, :]], Tuple[int]],
                Fun[[NDArray[complex, :, :]], Tuple[int, int]],
                Fun[[NDArray[complex, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[complex, :, :, :, :]], Tuple[int, int, int, int]],
            ],
            return_range_content=interval.positive_values
        ),
        "size": AttributeIntr(signature=Fun[[NDArray[T0, :]], int],
                              return_range=interval.positive_values),
        "sort": MethodIntr(
            args=("self", "axis", "kind"),
            defaults=(-1, None)
        ),
        "strides": AttributeIntr(
            signature=Union[
                # bool
                Fun[[NDArray[bool, :]], Tuple[int]],
                Fun[[NDArray[bool, :, :]], Tuple[int, int]],
                Fun[[NDArray[bool, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[bool, :, :, :, :]], Tuple[int, int, int, int]],
                # int
                Fun[[NDArray[int, :]], Tuple[int]],
                Fun[[NDArray[int, :, :]], Tuple[int, int]],
                Fun[[NDArray[int, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[int, :, :, :, :]], Tuple[int, int, int, int]],
                # float
                Fun[[NDArray[float, :]], Tuple[int]],
                Fun[[NDArray[float, :, :]], Tuple[int, int]],
                Fun[[NDArray[float, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[float, :, :, :, :]], Tuple[int, int, int, int]],
                # complex
                Fun[[NDArray[complex, :]], Tuple[int]],
                Fun[[NDArray[complex, :, :]], Tuple[int, int]],
                Fun[[NDArray[complex, :, :, :]], Tuple[int, int, int]],
                Fun[[NDArray[complex, :, :, :, :]], Tuple[int, int, int, int]],
            ]
        ),
        "T": AttributeIntr(signature=Fun[[NDArray[T0, :]], NDArray[T0, :]]),
        "tofile": ConstMethodIntr(signature=Fun[[NDArray[T0, :]], str, str], global_effects=True),
        "tostring": ConstMethodIntr(signature=Fun[[NDArray[T0, :]], str]),
        "view": MethodIntr(),
    },
}


_numpy_ones_signature = Union[
    # 1d
    Fun[[int], NDArray[float, :]],
    Fun[[int, _bool_signature], NDArray[bool, :]],
    Fun[[int, _int_signature], NDArray[int, :]],
    Fun[[int, _float_signature], NDArray[float, :]],
    Fun[[int, _complex_signature], NDArray[complex, :]],
    # 1D tuple
    Fun[[Tuple[int]], NDArray[float, :]],
    Fun[[Tuple[int], _bool_signature], NDArray[bool, :]],
    Fun[[Tuple[int], _int_signature], NDArray[int, :]],
    Fun[[Tuple[int], _float_signature], NDArray[float, :]],
    Fun[[Tuple[int], _complex_signature], NDArray[complex, :]],
    # 2D tuple
    Fun[[Tuple[int, int]], NDArray[float, :, :]],
    Fun[[Tuple[int, int], _bool_signature], NDArray[bool, :, :]],
    Fun[[Tuple[int, int], _int_signature], NDArray[int, :, :]],
    Fun[[Tuple[int, int], _float_signature], NDArray[float, :, :]],
    Fun[[Tuple[int, int], _complex_signature], NDArray[complex, :, :]],
    # 3D tuple
    Fun[[Tuple[int, int, int]], NDArray[float, :, :, :]],
    Fun[[Tuple[int, int, int], _bool_signature], NDArray[bool, :, :, :]],
    Fun[[Tuple[int, int, int], _int_signature], NDArray[int, :, :, :]],
    Fun[[Tuple[int, int, int], _float_signature], NDArray[float, :, :, :]],
    Fun[[Tuple[int, int, int], _complex_signature], NDArray[complex, :, :, :]],
    # 4D tuple
    Fun[[Tuple[int, int, int, int]], NDArray[float, :, :, :, :]],
    Fun[[Tuple[int, int, int, int], _bool_signature],
        NDArray[bool, :, :, :, :]],
    Fun[[Tuple[int, int, int, int], _int_signature],
        NDArray[int, :, :, :, :]],
    Fun[[Tuple[int, int, int, int], _float_signature],
        NDArray[float, :, :, :, :]],
    Fun[[Tuple[int, int, int, int], _complex_signature],
        NDArray[complex, :, :, :, :]],
]
_numpy_ones_like_signature = Union[
    # scalar
    Fun[[bool], bool],
    Fun[[int], int],
    Fun[[float], float],
    Fun[[complex], complex],
    # scalar + None
    Fun[[bool, None], bool],
    Fun[[int, None], int],
    Fun[[float, None], float],
    Fun[[complex, None], complex],
    # scalar + dtype
    Fun[[bool, _bool_signature], bool],
    Fun[[bool, _int_signature], int],
    Fun[[bool, _float_signature], float],
    Fun[[bool, _complex_signature], complex],
    Fun[[int, _bool_signature], bool],
    Fun[[int, _int_signature], int],
    Fun[[int, _float_signature], float],
    Fun[[int, _complex_signature], complex],
    Fun[[complex, _bool_signature], bool],
    Fun[[complex, _int_signature], int],
    Fun[[complex, _float_signature], float],
    Fun[[complex, _complex_signature], complex],
    # array 1D
    Fun[[Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex]], NDArray[complex, :]],

    # array 2d
    Fun[[Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]]], NDArray[complex, :, :]],

    # array 3d
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], _float_signature],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],

    # array 4d
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[bool, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
        NDArray[complex, :, :, :, :]],

    # with dtype
]

_numpy_unary_op_signature = Union[
    # 1d
    Fun[[bool], bool],
    Fun[[int], int],
    Fun[[float], float],
    Fun[[complex], complex],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex]], NDArray[complex, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]]], NDArray[complex, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[bool, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
        NDArray[complex, :, :, :, :]],
]

_numpy_float_unary_op_signature = Union[
    # 1d
    Fun[[bool], float],
    Fun[[int], float],
    Fun[[float], float],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],
]

_numpy_int_unary_op_signature = Union[
    # 1d
    Fun[[bool], bool],
    Fun[[int], int],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[int]], NDArray[int, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[int, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[bool, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[int, :, :, :, :]],
]

_numpy_unary_op_angle_signature = Union[
    # no extra option
    # 1d
    Fun[[bool], float],
    Fun[[int], float],
    Fun[[float], float],
    Fun[[complex], float],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex]], NDArray[float, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]]], NDArray[float, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]]], NDArray[float, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
        NDArray[float, :, :, :, :]],

    # extra option
    # 1d
    Fun[[bool, bool], float],
    Fun[[int, bool], float],
    Fun[[float, bool], float],
    Fun[[complex, bool], float],
    # 1d Iterable
    Fun[[Iterable[bool], bool], NDArray[float, :]],
    Fun[[Iterable[int], bool], NDArray[float, :]],
    Fun[[Iterable[float], bool], NDArray[float, :]],
    Fun[[Iterable[complex], bool], NDArray[float, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], bool], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], bool], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], bool], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]], bool], NDArray[float, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], bool], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], bool], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], bool], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], bool],
        NDArray[float, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], bool],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], bool],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], bool],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], bool],
        NDArray[float, :, :, :, :]],
]

_numpy_array_str_signature = Union[
    tuple(Fun[[NDArray[(dtype,) + slices]], str]
          for dtype in (bool, int, float, complex)
          for slices in [(slice(0),) * i
                         for i in range(1, 5)])
]
_numpy_float_unary_op_float_signature = Union[
    # 1d
    Fun[[bool], float],
    Fun[[int], float],
    Fun[[float], float],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],
]


_numpy_unary_op_float_signature = Union[
    # 1d
    Fun[[bool], float],
    Fun[[int], float],
    Fun[[float], float],
    Fun[[complex], complex],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex]], NDArray[complex, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]]], NDArray[complex, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
        NDArray[complex, :, :, :, :]],
]

_numpy_unary_op_int_signature = Union[
    # 1d
    Fun[[bool], int],
    Fun[[int], int],
    Fun[[float], int],
    Fun[[complex], int],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[int, :]],
    Fun[[Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[float]], NDArray[int, :]],
    Fun[[Iterable[complex]], NDArray[int, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[complex]]], NDArray[int, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]]], NDArray[int, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
        NDArray[int, :, :, :, :]],
]

_numpy_unary_op_axis_signature = Union[
    # no axis
    # 1d
    Fun[[bool], bool],
    Fun[[int], int],
    Fun[[float], float],
    Fun[[complex], complex],
    # 1d Iterable
    Fun[[Iterable[bool]], bool],
    Fun[[Iterable[int]], int],
    Fun[[Iterable[float]], float],
    Fun[[Iterable[complex]], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], bool],
    Fun[[Iterable[Iterable[int]]], int],
    Fun[[Iterable[Iterable[float]]], float],
    Fun[[Iterable[Iterable[complex]]], complex],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], bool],
    Fun[[Iterable[Iterable[Iterable[int]]]], int],
    Fun[[Iterable[Iterable[Iterable[float]]]], float],
    Fun[[Iterable[Iterable[Iterable[complex]]]], complex],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]], bool],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]], int],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]], complex],

    # axis
    # 1d Iterable
    Fun[[Iterable[bool], int], bool],
    Fun[[Iterable[int], int], int],
    Fun[[Iterable[float], int], float],
    Fun[[Iterable[complex], int], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], int], NDArray[bool, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[int, :]],
    Fun[[Iterable[Iterable[float]], int], NDArray[float, :]],
    Fun[[Iterable[Iterable[complex]], int], NDArray[complex, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], int], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], int], NDArray[complex, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], int],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], int],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], int],
        NDArray[complex, :, :, :]],
]

_numpy_unary_op_int_axis_signature = Union[
    # no axis
    # 1d
    Fun[[bool], int],
    Fun[[int], int],
    Fun[[float], int],
    Fun[[complex], int],
    # 1d Iterable
    Fun[[Iterable[bool]], int],
    Fun[[Iterable[int]], int],
    Fun[[Iterable[float]], int],
    Fun[[Iterable[complex]], int],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], int],
    Fun[[Iterable[Iterable[int]]], int],
    Fun[[Iterable[Iterable[float]]], int],
    Fun[[Iterable[Iterable[complex]]], int],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], int],
    Fun[[Iterable[Iterable[Iterable[int]]]], int],
    Fun[[Iterable[Iterable[Iterable[float]]]], int],
    Fun[[Iterable[Iterable[Iterable[complex]]]], int],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]], int],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]], int],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]], int],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]], int],

    # axis
    # 1d Iterable
    Fun[[Iterable[bool], int], int],
    Fun[[Iterable[int], int], int],
    Fun[[Iterable[float], int], int],
    Fun[[Iterable[complex], int], int],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], int], NDArray[int, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[int, :]],
    Fun[[Iterable[Iterable[float]], int], NDArray[int, :]],
    Fun[[Iterable[Iterable[complex]], int], NDArray[int, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], int], NDArray[int, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], int],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], int],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], int],
        NDArray[int, :, :, :]],
]

_numpy_unary_op_sum_axis_signature = Union[
    # no axis
    # 1d
    Fun[[bool], int],
    Fun[[int], int],
    Fun[[float], float],
    Fun[[complex], complex],
    # 1d Iterable
    Fun[[Iterable[bool]], int],
    Fun[[Iterable[int]], int],
    Fun[[Iterable[float]], float],
    Fun[[Iterable[complex]], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], int],
    Fun[[Iterable[Iterable[int]]], int],
    Fun[[Iterable[Iterable[float]]], float],
    Fun[[Iterable[Iterable[complex]]], complex],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], int],
    Fun[[Iterable[Iterable[Iterable[int]]]], int],
    Fun[[Iterable[Iterable[Iterable[float]]]], float],
    Fun[[Iterable[Iterable[Iterable[complex]]]], complex],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]], int],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]], int],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]], complex],

    # axis
    # 1d
    Fun[[bool, int], int],
    Fun[[int, int], int],
    Fun[[float, int], float],
    Fun[[complex, int], complex],
    # 1d Iterable
    Fun[[Iterable[bool], int], int],
    Fun[[Iterable[int], int], int],
    Fun[[Iterable[float], int], float],
    Fun[[Iterable[complex], int], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], int], NDArray[int, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[int, :]],
    Fun[[Iterable[Iterable[float]], int], NDArray[float, :]],
    Fun[[Iterable[Iterable[complex]], int], NDArray[complex, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], int], NDArray[complex, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], int],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], int],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], int],
        NDArray[complex, :, :, :]],
]

_numpy_unary_op_cumsum_axis_signature = Union[
    # no axis
    # 1d
    Fun[[bool], NDArray[int, :]],
    Fun[[int], NDArray[int, :]],
    Fun[[float], NDArray[float, :]],
    Fun[[complex], NDArray[complex, :]],
    # 1d Iterable
    Fun[[Iterable[bool]], NDArray[int, :]],
    Fun[[Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex]], NDArray[complex, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], NDArray[int, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[int, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :]],
    Fun[[Iterable[Iterable[complex]]], NDArray[complex, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[int, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[int, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]]], NDArray[complex, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]], NDArray[int, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]], NDArray[int, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]], NDArray[float, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
        NDArray[complex, :]],

    # axis
    # 1d
    Fun[[bool, int], NDArray[int, :]],
    Fun[[int, int], NDArray[int, :]],
    Fun[[float, int], NDArray[float, :]],
    Fun[[complex, int], NDArray[complex, :]],
    # 1d Iterable
    Fun[[Iterable[bool], int], NDArray[int, :]],
    Fun[[Iterable[int], int], NDArray[int, :]],
    Fun[[Iterable[float], int], NDArray[float, :]],
    Fun[[Iterable[complex], int], NDArray[complex, :]],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[float]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]], int], NDArray[complex, :, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], int], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], int], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], int],
        NDArray[complex, :, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], int],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], int],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], int],
        NDArray[complex, :, :, :, :]],
]
_numpy_unary_op_average_axis_signature = Union[
    # no axis
    # 1d
    Fun[[bool], float],
    Fun[[int], float],
    Fun[[float], float],
    Fun[[complex], complex],
    # 1d Iterable
    Fun[[Iterable[bool]], float],
    Fun[[Iterable[int]], float],
    Fun[[Iterable[float]], float],
    Fun[[Iterable[complex]], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], float],
    Fun[[Iterable[Iterable[int]]], float],
    Fun[[Iterable[Iterable[float]]], float],
    Fun[[Iterable[Iterable[complex]]], complex],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], float],
    Fun[[Iterable[Iterable[Iterable[int]]]], float],
    Fun[[Iterable[Iterable[Iterable[float]]]], float],
    Fun[[Iterable[Iterable[Iterable[complex]]]], complex],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]], complex],

    # axis None
    # 1d
    Fun[[bool, None], float],
    Fun[[int, None], float],
    Fun[[float, None], float],
    Fun[[complex, None], complex],
    # 1d Iterable
    Fun[[Iterable[bool], None], float],
    Fun[[Iterable[int], None], float],
    Fun[[Iterable[float], None], float],
    Fun[[Iterable[complex], None], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], None], float],
    Fun[[Iterable[Iterable[int]], None], float],
    Fun[[Iterable[Iterable[float]], None], float],
    Fun[[Iterable[Iterable[complex]], None], complex],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], None], float],
    Fun[[Iterable[Iterable[Iterable[int]]], None], float],
    Fun[[Iterable[Iterable[Iterable[float]]], None], float],
    Fun[[Iterable[Iterable[Iterable[complex]]], None], complex],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], None], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], None], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], None], float],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], None], complex],

    # axis
    # 1d
    Fun[[bool, int], float],
    Fun[[int, int], float],
    Fun[[float, int], float],
    Fun[[complex, int], complex],
    # 1d Iterable
    Fun[[Iterable[bool], int], float],
    Fun[[Iterable[int], int], float],
    Fun[[Iterable[float], int], float],
    Fun[[Iterable[complex], int], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], int], NDArray[float, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[float, :]],
    Fun[[Iterable[Iterable[float]], int], NDArray[float, :]],
    Fun[[Iterable[Iterable[complex]], int], NDArray[complex, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], int], NDArray[complex, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], int],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], int],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], int],
        NDArray[complex, :, :, :]],

    # axis None + weight
    # 1d
    Fun[[bool, None, float], float],
    Fun[[int, None, float], float],
    Fun[[float, None, float], float],
    Fun[[complex, None, float], complex],
    # 1d Iterable
    Fun[[Iterable[bool], None, Iterable[float]], float],
    Fun[[Iterable[int], None, Iterable[float]], float],
    Fun[[Iterable[float], None, Iterable[float]], float],
    Fun[[Iterable[complex], None, Iterable[float]], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], None, Iterable[float]], float],
    Fun[[Iterable[Iterable[int]], None, Iterable[float]], float],
    Fun[[Iterable[Iterable[float]], None, Iterable[float]], float],
    Fun[[Iterable[Iterable[complex]], None, Iterable[float]], complex],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], None, Iterable[float]], float],
    Fun[[Iterable[Iterable[Iterable[int]]], None, Iterable[float]], float],
    Fun[[Iterable[Iterable[Iterable[float]]], None, Iterable[float]], float],
    Fun[[Iterable[Iterable[Iterable[complex]]], None, Iterable[float]],
        complex],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], None, Iterable[float]],
        float],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], None, Iterable[float]],
        float],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], None, Iterable[float]],
        float],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
         None, Iterable[float]], complex],

    # axis
    # 1d
    Fun[[bool, int, float], float],
    Fun[[int, int, float], float],
    Fun[[float, int, float], float],
    Fun[[complex, int, float], complex],
    # 1d Iterable
    Fun[[Iterable[bool], int, Iterable[float]], float],
    Fun[[Iterable[int], int, Iterable[float]], float],
    Fun[[Iterable[float], int, Iterable[float]], float],
    Fun[[Iterable[complex], int, Iterable[float]], complex],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], int, Iterable[Iterable[bool]]],
        NDArray[float, :]],
    Fun[[Iterable[Iterable[int]], int, Iterable[Iterable[int]]],
        NDArray[float, :]],
    Fun[[Iterable[Iterable[float]], int, Iterable[
        Iterable[float]]], NDArray[float, :]],
    Fun[[Iterable[Iterable[complex]], int, Iterable[
        Iterable[complex]]], NDArray[complex, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], int, Iterable[
        Iterable[Iterable[bool]]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int, Iterable[
        Iterable[Iterable[int]]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], int, Iterable[
        Iterable[Iterable[float]]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], int, Iterable[
        Iterable[Iterable[complex]]]], NDArray[complex, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], int, Iterable[
        Iterable[Iterable[Iterable[bool]]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int, Iterable[
        Iterable[Iterable[Iterable[int]]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], int, Iterable[
        Iterable[Iterable[Iterable[float]]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], int, Iterable[
        Iterable[Iterable[Iterable[complex]]]]], NDArray[complex, :, :, :]],

]

_numpy_unary_op_bool_axis_signature = Union[
    # no axis
    # 1d
    Fun[[bool], bool],
    Fun[[int], bool],
    Fun[[float], bool],
    Fun[[complex], bool],
    # 1d Iterable
    Fun[[Iterable[bool]], bool],
    Fun[[Iterable[int]], bool],
    Fun[[Iterable[float]], bool],
    Fun[[Iterable[complex]], bool],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]]], bool],
    Fun[[Iterable[Iterable[int]]], bool],
    Fun[[Iterable[Iterable[float]]], bool],
    Fun[[Iterable[Iterable[complex]]], bool],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]]], bool],
    Fun[[Iterable[Iterable[Iterable[int]]]], bool],
    Fun[[Iterable[Iterable[Iterable[float]]]], bool],
    Fun[[Iterable[Iterable[Iterable[complex]]]], bool],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]], bool],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]], bool],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]], bool],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]], bool],

    # axis
    # 1d
    Fun[[bool, int], bool],
    Fun[[int, int], bool],
    Fun[[float, int], bool],
    Fun[[complex, int], bool],
    # 1d Iterable
    Fun[[Iterable[bool], int], bool],
    Fun[[Iterable[int], int], bool],
    Fun[[Iterable[float], int], bool],
    Fun[[Iterable[complex], int], bool],
    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], int], NDArray[bool, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[bool, :]],
    Fun[[Iterable[Iterable[float]], int], NDArray[bool, :]],
    Fun[[Iterable[Iterable[complex]], int], NDArray[bool, :]],
    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], int], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], int], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], int], NDArray[bool, :, :]],
    # 4d Iterable
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], int],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], int],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], int],
        NDArray[bool, :, :, :]],
]


_numpy_binary_op_signature = Union[
    # 1d
    Fun[[bool, bool], bool],
    Fun[[int, int], int],
    Fun[[float, float], float],
    Fun[[complex, complex], complex],

    # 1d Iterable
    Fun[[Iterable[bool], Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[bool], bool], NDArray[bool, :]],
    Fun[[bool, Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[int], Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[int], int], NDArray[int, :]],
    Fun[[int, Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[float], Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[float], float], NDArray[float, :]],
    Fun[[float, Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex], Iterable[complex]], NDArray[complex, :]],
    Fun[[Iterable[complex], complex], NDArray[complex, :]],
    Fun[[complex, Iterable[complex]], NDArray[complex, :]],

    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], Iterable[Iterable[bool]]],
        NDArray[bool, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[bool]], NDArray[bool, :, :]],
    Fun[[bool, Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[bool]], bool], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]]],
        NDArray[int, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[int]], NDArray[int, :, :]],
    Fun[[int, Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[Iterable[float]]],
        NDArray[float, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[float]], NDArray[float, :, :]],
    Fun[[float, Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], float], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[
        Iterable[complex]]], NDArray[complex, :, :]],
    Fun[[Iterable[complex], Iterable[Iterable[complex]]],
        NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[complex]],
        NDArray[complex, :, :]],
    Fun[[complex, Iterable[Iterable[complex]]], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], complex], NDArray[complex, :, :]],

    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[Iterable[bool]]]],
        NDArray[bool, :, :, :]],
    Fun[[bool, Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[bool]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[bool]],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], bool], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[
        Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[Iterable[int]]]],
        NDArray[int, :, :, :]],
    Fun[[int, Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[int]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[int]],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[
        Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[Iterable[float]]]],
        NDArray[float, :, :, :]],
    Fun[[float, Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[float]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         Iterable[float]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], float], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
        Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[
        Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],
    Fun[[Iterable[complex], Iterable[Iterable[Iterable[complex]]]],
        NDArray[complex, :, :, :]],
    Fun[[complex, Iterable[Iterable[Iterable[complex]]]],
        NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
        Iterable[complex]]], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         Iterable[complex]], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], complex],
        NDArray[complex, :, :, :]],
]
_numpy_binary_op_bool_signature = Union[
    # 1d
    Fun[[bool, bool], bool],
    Fun[[int, int], bool],
    Fun[[float, float], bool],
    Fun[[complex, complex], bool],

    # 1d Iterable
    Fun[[Iterable[bool], Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[bool], bool], NDArray[bool, :]],
    Fun[[bool, Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[int], Iterable[int]], NDArray[bool, :]],
    Fun[[Iterable[int], int], NDArray[bool, :]],
    Fun[[int, Iterable[int]], NDArray[bool, :]],
    Fun[[Iterable[float], Iterable[float]], NDArray[bool, :]],
    Fun[[Iterable[float], float], NDArray[bool, :]],
    Fun[[float, Iterable[float]], NDArray[bool, :]],
    Fun[[Iterable[complex], Iterable[complex]], NDArray[bool, :]],
    Fun[[Iterable[complex], complex], NDArray[bool, :]],
    Fun[[complex, Iterable[complex]], NDArray[bool, :]],

    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], Iterable[Iterable[bool]]],
        NDArray[bool, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[bool]], NDArray[bool, :, :]],
    Fun[[bool, Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[bool]], bool], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]]],
        NDArray[bool, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[int]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[int]], NDArray[bool, :, :]],
    Fun[[int, Iterable[Iterable[int]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[Iterable[float]]],
        NDArray[bool, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[float]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[float]], NDArray[bool, :, :]],
    Fun[[float, Iterable[Iterable[float]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[float]], float], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[
        Iterable[complex]]], NDArray[bool, :, :]],
    Fun[[Iterable[complex], Iterable[Iterable[complex]]],
        NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[complex]],
        NDArray[bool, :, :]],
    Fun[[complex, Iterable[Iterable[complex]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[complex]], complex], NDArray[bool, :, :]],

    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[Iterable[bool]]]],
        NDArray[bool, :, :, :]],
    Fun[[bool, Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[bool]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[bool]],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], bool], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[Iterable[int]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[
        Iterable[Iterable[int]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[Iterable[int]]]],
        NDArray[bool, :, :, :]],
    Fun[[int, Iterable[Iterable[Iterable[int]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[int]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[int]],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[Iterable[float]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[
        Iterable[Iterable[float]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[Iterable[float]]]],
        NDArray[bool, :, :, :]],
    Fun[[float, Iterable[Iterable[Iterable[float]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[float]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         Iterable[float]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], float], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
        Iterable[Iterable[complex]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[
        Iterable[Iterable[complex]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[complex], Iterable[Iterable[Iterable[complex]]]],
        NDArray[bool, :, :, :]],
    Fun[[complex, Iterable[Iterable[Iterable[complex]]]],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
        Iterable[complex]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         Iterable[complex]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], complex],
        NDArray[bool, :, :, :]],
]
_numpy_binary_op_float_signature = Union[
    # 1d
    Fun[[bool, bool], float],
    Fun[[int, int], float],
    Fun[[float, float], float],
    Fun[[complex, complex], complex],

    # 1d Iterable
    Fun[[Iterable[bool], Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[bool], bool], NDArray[float, :]],
    Fun[[bool, Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[int], Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[int], int], NDArray[float, :]],
    Fun[[int, Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[float], Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[float], float], NDArray[float, :]],
    Fun[[float, Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex], Iterable[complex]], NDArray[complex, :]],
    Fun[[Iterable[complex], complex], NDArray[complex, :]],
    Fun[[complex, Iterable[complex]], NDArray[complex, :]],

    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], Iterable[Iterable[bool]]],
        NDArray[float, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[bool]], NDArray[float, :, :]],
    Fun[[bool, Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[bool]], bool], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]]],
        NDArray[float, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[int]], NDArray[float, :, :]],
    Fun[[int, Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[Iterable[float]]],
        NDArray[float, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[float]], NDArray[float, :, :]],
    Fun[[float, Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], float], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[
        Iterable[complex]]], NDArray[complex, :, :]],
    Fun[[Iterable[complex], Iterable[Iterable[complex]]],
        NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[complex]],
        NDArray[complex, :, :]],
    Fun[[complex, Iterable[Iterable[complex]]], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], complex], NDArray[complex, :, :]],

    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[Iterable[bool]]]],
        NDArray[float, :, :, :]],
    Fun[[bool, Iterable[Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[bool]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]],
         Iterable[bool]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], bool], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[Iterable[int]]]],
        NDArray[float, :, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[Iterable[int]]]],
        NDArray[float, :, :, :]],
    Fun[[int, Iterable[Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[int]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[int]],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[
        Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[Iterable[float]]]],
        NDArray[float, :, :, :]],
    Fun[[float, Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[float]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         Iterable[float]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], float], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
        Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[
        Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],
    Fun[[Iterable[complex], Iterable[Iterable[Iterable[complex]]]],
        NDArray[complex, :, :, :]],
    Fun[[complex, Iterable[Iterable[Iterable[complex]]]],
        NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
        Iterable[complex]]], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         Iterable[complex]], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]], complex],
        NDArray[complex, :, :, :]],
]
_numpy_ternary_op_signature = Union[
    # scalar
    Fun[[int, int, int], int],
    Fun[[float, float, float], float],
    Fun[[complex, complex, complex], complex],

    # 1D
    Fun[[Iterable[int], int, int], NDArray[int, :]],
    Fun[[Iterable[int], Iterable[int], int], NDArray[int, :]],
    Fun[[Iterable[int], int, Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[int], Iterable[int], Iterable[int]], NDArray[int, :]],

    Fun[[Iterable[float], float, float], NDArray[float, :]],
    Fun[[Iterable[float], Iterable[float], float], NDArray[float, :]],
    Fun[[Iterable[float], float, Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[float], Iterable[float], Iterable[float]],
        NDArray[float, :]],

    Fun[[Iterable[complex], complex, complex], NDArray[complex, :]],
    Fun[[Iterable[complex], Iterable[complex], complex], NDArray[complex, :]],
    Fun[[Iterable[complex], complex, Iterable[complex]], NDArray[complex, :]],
    Fun[[Iterable[complex], Iterable[complex],
         Iterable[complex]], NDArray[complex, :]],

    # 2D
    Fun[[Iterable[Iterable[int]], int, int], NDArray[int, :, :]],

    Fun[[Iterable[Iterable[int]], Iterable[int], int], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[int], Iterable[int]],
        NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], int, Iterable[int]], NDArray[int, :, :]],

    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]], int],
        NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], int, Iterable[Iterable[int]]],
        NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]],
         Iterable[int]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[int],
         Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]],
         Iterable[Iterable[int]]], NDArray[int, :, :]],

    Fun[[Iterable[Iterable[float]], float, float], NDArray[float, :, :]],

    Fun[[Iterable[Iterable[float]], Iterable[float], float],
        NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[float],
         Iterable[float]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], float, Iterable[float]],
        NDArray[float, :, :]],

    Fun[[Iterable[Iterable[float]], Iterable[
        Iterable[float]], float], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], float, Iterable[
        Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[Iterable[float]],
         Iterable[float]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[float],
         Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[Iterable[float]],
         Iterable[Iterable[float]]], NDArray[float, :, :]],

    Fun[[Iterable[Iterable[complex]], complex, complex],
        NDArray[complex, :, :]],

    Fun[[Iterable[Iterable[complex]], Iterable[
        complex], complex], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[complex],
         Iterable[complex]], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], complex,
         Iterable[complex]], NDArray[complex, :, :]],

    Fun[[Iterable[Iterable[complex]], Iterable[
        Iterable[complex]], complex], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], complex, Iterable[
        Iterable[complex]]], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[Iterable[complex]],
         Iterable[complex]], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[complex],
         Iterable[Iterable[complex]]], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], Iterable[Iterable[complex]],
         Iterable[Iterable[complex]]], NDArray[complex, :, :]],
]

_numpy_int_binary_op_signature = Union[
    # 1d
    Fun[[bool, bool], bool],
    Fun[[int, int], int],

    # 1d Iterable
    Fun[[Iterable[bool], Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[bool], bool], NDArray[bool, :]],
    Fun[[bool, Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[int], Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[int], int], NDArray[int, :]],
    Fun[[int, Iterable[int]], NDArray[int, :]],

    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], Iterable[Iterable[bool]]],
        NDArray[bool, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[bool]], NDArray[bool, :, :]],
    Fun[[bool, Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[bool]], bool], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]]],
        NDArray[int, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[int]], NDArray[int, :, :]],
    Fun[[int, Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[int, :, :]],

    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[Iterable[bool]]]],
        NDArray[bool, :, :, :]],
    Fun[[bool, Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[bool]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[bool]],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], bool], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[
        Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[Iterable[int]]]],
        NDArray[int, :, :, :]],
    Fun[[int, Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[int]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[int]],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[int, :, :, :]],
]

_numpy_binary_op_float_no_complex_signature = Union[
    # 1d
    Fun[[bool, bool], float],
    Fun[[int, int], float],
    Fun[[float, float], float],

    # 1d Iterable
    Fun[[Iterable[bool], Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[bool], bool], NDArray[float, :]],
    Fun[[bool, Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[int], Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[int], int], NDArray[float, :]],
    Fun[[int, Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[float], Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[float], float], NDArray[float, :]],
    Fun[[float, Iterable[float]], NDArray[float, :]],

    # 2d Iterable
    Fun[[Iterable[Iterable[bool]], Iterable[Iterable[bool]]],
        NDArray[float, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[bool]], NDArray[float, :, :]],
    Fun[[bool, Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[bool]], bool], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[int]]],
        NDArray[float, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[int]], NDArray[float, :, :]],
    Fun[[int, Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], int], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[Iterable[float]]],
        NDArray[float, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[float]], NDArray[float, :, :]],
    Fun[[float, Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], float], NDArray[float, :, :]],

    # 3d Iterable
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[bool]], Iterable[
        Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[bool], Iterable[Iterable[Iterable[bool]]]],
        NDArray[float, :, :, :]],
    Fun[[bool, Iterable[Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
        Iterable[bool]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]],
         Iterable[bool]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[bool]]], bool], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[int]], Iterable[Iterable[Iterable[int]]]],
        NDArray[float, :, :, :]],
    Fun[[Iterable[int], Iterable[Iterable[Iterable[int]]]],
        NDArray[float, :, :, :]],
    Fun[[int, Iterable[Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
        Iterable[int]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], Iterable[int]],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[float]], Iterable[
        Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[float], Iterable[Iterable[Iterable[float]]]],
        NDArray[float, :, :, :]],
    Fun[[float, Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
        Iterable[float]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         Iterable[float]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]], float], NDArray[float, :, :, :]],
]

_numpy_allclose_signature = Union[
    _numpy_binary_op_signature.__args__ +
    tuple([Fun[c.__args__[:-1] + (float,), c.__args__[-1]]
           for c in _numpy_unary_op_signature.__args__] +
          [Fun[c.__args__[:-1] + (float, float), c.__args__[-1]]
           for c in _numpy_unary_op_signature.__args__] +
          [Fun[c.__args__[:-1] + (float, float, bool), c.__args__[-1]]
           for c in _numpy_unary_op_signature.__args__]
          )
]

_numpy_around_signature = Union[
    _numpy_unary_op_float_signature.__args__ +
    tuple([Fun[c.__args__[:-1] + (int,), c.__args__[-1]]
           for c in _numpy_unary_op_float_signature.__args__])
]

_functools_reduce_signature = Union[
    Fun[[Fun[[T0, T0], T0], Iterable[T0]], T0],
    Fun[[Fun[[T0, T1], T0], Iterable[T1], T0], T0],
]


def partialsum(seq):
    s = tuple()
    for elt in seq:
        s += elt,
        yield s


_operator_add_signature = Union[
    _numpy_binary_op_signature.__args__
]

_operator_eq_signature = Union[
    _numpy_binary_op_bool_signature.__args__ +
    (Fun[[str, str], bool],
     Fun[[List[T0], List[T0]], bool],
     Fun[[Set[T0], Set[T0]], bool],
     Fun[[T0, None], bool],
     Fun[[None, T0], bool],
     Fun[[Dict[T0, T1], Dict[T0, T1]], bool],) +
    tuple(Fun[[Tuple[t0], Tuple[t1]], Tuple[t0 + t1]]
          for t0 in partialsum([T0, T1, T2, T3])
          for t1 in partialsum([T4, T5, T6, T7]))
]

_operator_sub_signature = Union[
    _numpy_binary_op_signature.__args__ +
    (Fun[[Set[T0], Set[T0]], Set[T0]],)
]

_operator_mod_signature = Union[
    _numpy_binary_op_signature.__args__ +
    (Fun[[str, T0], str],)
]

_operator_mul_signature = Union[
    _numpy_binary_op_signature.__args__ +
    (Fun[[str, int], str], Fun[[int, str], str],
     Fun[[List[T0], int], List[T0]], Fun[[int, List[T0]], List[T0]])
]

_operator_contains_signature = Fun[[Iterable[T0], T0], bool]

_operator_getitem_signature = Union[
    Fun[[List[T0], int], T0],
    Fun[[List[T0], slice], List[T0]],
    Fun[[Dict[T0, T1], T0], T1],
    Fun[[str, int], str],
    Fun[[str, slice], str],
    # arrays
    Fun[[NDArray[T0, :], T1], T2],
    # large tuple
    Fun[[Iterable[T0], int], T0],
]


_numpy_farray_signature = Union[
    # no dtype
    # scalar
    Fun[[bool], float],
    Fun[[int], float],
    Fun[[float], float],
    # 1D array
    Fun[[Iterable[bool]], NDArray[float, :]],
    Fun[[Iterable[int]], NDArray[float, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],

    # bool dtype
    # scalar
    Fun[[bool, _bool_signature], float],
    Fun[[int, _bool_signature], float],
    Fun[[float, _bool_signature], float],
    Fun[[complex, _bool_signature], float],
    # 1D array
    Fun[[Iterable[bool], _bool_signature], NDArray[float, :]],
    Fun[[Iterable[int], _bool_signature], NDArray[float, :]],
    Fun[[Iterable[float], _bool_signature], NDArray[float, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _bool_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], _bool_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], _bool_signature], NDArray[float, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]],
         _bool_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], _bool_signature],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _bool_signature], NDArray[float, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _bool_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _bool_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _bool_signature], NDArray[float, :, :, :, :]],

    # int dtype
    Fun[[bool, _int_signature], float],
    Fun[[int, _int_signature], float],
    Fun[[float, _int_signature], float],
    Fun[[complex, _int_signature], float],
    # 1D array
    Fun[[Iterable[bool], _int_signature], NDArray[float, :]],
    Fun[[Iterable[int], _int_signature], NDArray[float, :]],
    Fun[[Iterable[float], _int_signature], NDArray[float, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _int_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], _int_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], _int_signature], NDArray[float, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]],
         _int_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], _int_signature],
        NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _int_signature], NDArray[float, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _int_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _int_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _int_signature], NDArray[float, :, :, :, :]],

    # float dtype
    # scalar
    Fun[[bool, _float_signature], float],
    Fun[[int, _float_signature], float],
    Fun[[float, _float_signature], float],
    Fun[[complex, _float_signature], float],
    # 1D array
    Fun[[Iterable[bool], _float_signature], NDArray[float, :]],
    Fun[[Iterable[int], _float_signature], NDArray[float, :]],
    Fun[[Iterable[float], _float_signature], NDArray[float, :]],
    Fun[[Iterable[complex], _float_signature], NDArray[float, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _float_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], _float_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], _float_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]], _float_signature], NDArray[float, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]],
         _float_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]],
         _float_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _float_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         _float_signature], NDArray[float, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _float_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _float_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _float_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
         _float_signature], NDArray[float, :, :, :, :]],

    # complex dtype
    # scalar
    Fun[[bool, _complex_signature], complex],
    Fun[[int, _complex_signature], complex],
    Fun[[float, _complex_signature], complex],
    Fun[[complex, _complex_signature], complex],
    # 1D array
    Fun[[Iterable[bool], _complex_signature], NDArray[complex, :]],
    Fun[[Iterable[int], _complex_signature], NDArray[complex, :]],
    Fun[[Iterable[float], _complex_signature], NDArray[complex, :]],
    Fun[[Iterable[complex], _complex_signature], NDArray[complex, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _complex_signature],
        NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[int]], _complex_signature], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[float]], _complex_signature],
        NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], _complex_signature],
        NDArray[complex, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]], _complex_signature],
        NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], _complex_signature],
        NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _complex_signature], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         _complex_signature], NDArray[complex, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
]

_numpy_array_signature = Union[
    # no dtype
    # scalar
    Fun[[bool], bool],
    Fun[[int], int],
    Fun[[float], float],
    Fun[[complex], complex],
    # 1D array
    Fun[[Iterable[bool]], NDArray[bool, :]],
    Fun[[Iterable[int]], NDArray[int, :]],
    Fun[[Iterable[float]], NDArray[float, :]],
    Fun[[Iterable[complex]], NDArray[complex, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]]], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]]], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]]], NDArray[complex, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]]], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]]], NDArray[complex, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
        NDArray[bool, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
        NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
        NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
        NDArray[complex, :, :, :, :]],

    # bool dtype
    # scalar
    Fun[[bool, _bool_signature], bool],
    Fun[[int, _bool_signature], bool],
    Fun[[float, _bool_signature], bool],
    Fun[[complex, _bool_signature], bool],
    # 1D array
    Fun[[Iterable[bool], _bool_signature], NDArray[bool, :]],
    Fun[[Iterable[int], _bool_signature], NDArray[bool, :]],
    Fun[[Iterable[float], _bool_signature], NDArray[bool, :]],
    Fun[[Iterable[complex], _bool_signature], NDArray[bool, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _bool_signature], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[int]], _bool_signature], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[float]], _bool_signature], NDArray[bool, :, :]],
    Fun[[Iterable[Iterable[complex]], _bool_signature], NDArray[bool, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]],
         _bool_signature], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], _bool_signature],
        NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _bool_signature], NDArray[bool, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         _bool_signature], NDArray[bool, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _bool_signature], NDArray[bool, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _bool_signature], NDArray[bool, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _bool_signature], NDArray[bool, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
         _bool_signature], NDArray[bool, :, :, :, :]],

    # int dtype
    # scalar
    Fun[[bool, _int_signature], int],
    Fun[[int, _int_signature], int],
    Fun[[float, _int_signature], int],
    Fun[[complex, _int_signature], int],
    # 1D array
    Fun[[Iterable[bool], _int_signature], NDArray[int, :]],
    Fun[[Iterable[int], _int_signature], NDArray[int, :]],
    Fun[[Iterable[float], _int_signature], NDArray[int, :]],
    Fun[[Iterable[complex], _int_signature], NDArray[int, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _int_signature], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[int]], _int_signature], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[float]], _int_signature], NDArray[int, :, :]],
    Fun[[Iterable[Iterable[complex]], _int_signature], NDArray[int, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]], _int_signature],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], _int_signature],
        NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _int_signature], NDArray[int, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         _int_signature], NDArray[int, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _int_signature], NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _int_signature], NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _int_signature], NDArray[int, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
         _int_signature], NDArray[int, :, :, :, :]],

    # float dtype
    # scalar
    Fun[[bool, _float_signature], float],
    Fun[[int, _float_signature], float],
    Fun[[float, _float_signature], float],
    Fun[[complex, _float_signature], float],
    # 1D array
    Fun[[Iterable[bool], _float_signature], NDArray[float, :]],
    Fun[[Iterable[int], _float_signature], NDArray[float, :]],
    Fun[[Iterable[float], _float_signature], NDArray[float, :]],
    Fun[[Iterable[complex], _float_signature], NDArray[float, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _float_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[int]], _float_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[float]], _float_signature], NDArray[float, :, :]],
    Fun[[Iterable[Iterable[complex]], _float_signature], NDArray[float, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]],
         _float_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]],
         _float_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _float_signature], NDArray[float, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         _float_signature], NDArray[float, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _float_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _float_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _float_signature], NDArray[float, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
         _float_signature], NDArray[float, :, :, :, :]],

    # complex dtype
    # scalar
    Fun[[bool, _complex_signature], complex],
    Fun[[int, _complex_signature], complex],
    Fun[[float, _complex_signature], complex],
    Fun[[complex, _complex_signature], complex],
    # 1D array
    Fun[[Iterable[bool], _complex_signature], NDArray[complex, :]],
    Fun[[Iterable[int], _complex_signature], NDArray[complex, :]],
    Fun[[Iterable[float], _complex_signature], NDArray[complex, :]],
    Fun[[Iterable[complex], _complex_signature], NDArray[complex, :]],
    # 2D array
    Fun[[Iterable[Iterable[bool]], _complex_signature],
        NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[int]], _complex_signature], NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[float]], _complex_signature],
        NDArray[complex, :, :]],
    Fun[[Iterable[Iterable[complex]], _complex_signature],
        NDArray[complex, :, :]],
    # 3D array
    Fun[[Iterable[Iterable[Iterable[bool]]], _complex_signature],
        NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[int]]], _complex_signature],
        NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[float]]],
         _complex_signature], NDArray[complex, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[complex]]],
         _complex_signature], NDArray[complex, :, :, :]],
    # 4D array
    Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
    Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
         _complex_signature], NDArray[complex, :, :, :, :]],
]


def expand_numpy_2_args(args, defaults=None, force=False):
    if force or numpy.__version__[0] < '2':
        if defaults is not None:
            return {'args': args, 'defaults': defaults}
        else:
            return {'args': args}
    else:
        return {}


# each module consist in a module_name <> set of symbols
MODULES = {
    "builtins": {
        "pythran": {
            "abssqr": ConstFunctionIntr(),
            "static_list": ReadOnceFunctionIntr(
                signature=Fun[[Iterable[T0]], List[T0]],
                return_alias=lambda args: {args[0]}),
            "is_none": ConstFunctionIntr(),
            "kwonly": ConstFunctionIntr(),
            "len_set": ConstFunctionIntr(signature=Fun[[Iterable[T0]], int]),
            "make_shape": ConstFunctionIntr(),
            "restrict_assign": FunctionIntr(),
            "static_if": ConstFunctionIntr(),
            "StaticIfBreak": ConstFunctionIntr(),
            "StaticIfCont": ConstFunctionIntr(),
            "StaticIfNoReturn": ConstFunctionIntr(),
            "StaticIfReturn": ConstFunctionIntr(),
        },
        "abs": ConstFunctionIntr(
            signature=Union[
                Fun[[int], int],
                Fun[[float], float],
                Fun[[complex], float],
                Fun[[NDArray[int, :]], NDArray[int, :]],
                Fun[[NDArray[int, :, :]], NDArray[int, :, :]],
                Fun[[NDArray[int, :, :, :]], NDArray[int, :, :, :]],
                Fun[[NDArray[int, :, :, :, :]], NDArray[int, :, :, :, :]],
                Fun[[NDArray[float, :]], NDArray[float, :]],
                Fun[[NDArray[float, :, :]], NDArray[float, :, :]],
                Fun[[NDArray[float, :, :, :]], NDArray[float, :, :, :]],
                Fun[[NDArray[float, :, :, :, :]], NDArray[float, :, :, :, :]],
                Fun[[NDArray[complex, :]], NDArray[float, :]],
                Fun[[NDArray[complex, :, :]], NDArray[float, :, :]],
                Fun[[NDArray[complex, :, :, :]], NDArray[float, :, :, :]],
                Fun[[NDArray[complex, :, :, :, :]], NDArray[float, :, :, :, :]]
            ],
        ),
        "BaseException": ConstExceptionIntr(),
        "SystemExit": ConstExceptionIntr(),
        "KeyboardInterrupt": ConstExceptionIntr(),
        "GeneratorExit": ConstExceptionIntr(),
        "Exception": ExceptionClass(CLASSES["Exception"]),
        "StopIteration": ConstExceptionIntr(),
        "Warning": ConstExceptionIntr(),
        "BytesWarning": ConstExceptionIntr(),
        "UnicodeWarning": ConstExceptionIntr(),
        "ImportWarning": ConstExceptionIntr(),
        "FutureWarning": ConstExceptionIntr(),
        "UserWarning": ConstExceptionIntr(),
        "SyntaxWarning": ConstExceptionIntr(),
        "RuntimeWarning": ConstExceptionIntr(),
        "PendingDeprecationWarning": ConstExceptionIntr(),
        "DeprecationWarning": ConstExceptionIntr(),
        "BufferError": ConstExceptionIntr(),
        "ArithmeticError": ConstExceptionIntr(),
        "AssertionError": ConstExceptionIntr(),
        "AttributeError": ConstExceptionIntr(),
        "EnvironmentError": ConstExceptionIntr(),
        "EOFError": ConstExceptionIntr(),
        "ImportError": ConstExceptionIntr(),
        "LookupError": ConstExceptionIntr(),
        "MemoryError": ConstExceptionIntr(),
        "NameError": ConstExceptionIntr(),
        "ReferenceError": ConstExceptionIntr(),
        "RuntimeError": ConstExceptionIntr(),
        "SyntaxError": ConstExceptionIntr(),
        "SystemError": ConstExceptionIntr(),
        "TypeError": ConstExceptionIntr(),
        "ValueError": ConstExceptionIntr(),
        "FloatingPointError": ConstExceptionIntr(),
        "OverflowError": ConstExceptionIntr(),
        "ZeroDivisionError": ConstExceptionIntr(),
        "IOError": ConstExceptionIntr(),
        "OSError": ConstExceptionIntr(),
        "IndexError": ConstExceptionIntr(),
        "KeyError": ConstExceptionIntr(),
        "UnboundLocalError": ConstExceptionIntr(),
        "NotImplementedError": ConstExceptionIntr(),
        "IndentationError": ConstExceptionIntr(),
        "TabError": ConstExceptionIntr(),
        "UnicodeError": ConstExceptionIntr(),
        #  "UnicodeDecodeError": ConstExceptionIntr(),
        #  "UnicodeEncodeError": ConstExceptionIntr(),
        #  "UnicodeTranslateError": ConstExceptionIntr(),
        "all": ReadOnceFunctionIntr(signature=Fun[[Iterable[T0]], bool]),
        "any": ReadOnceFunctionIntr(signature=Fun[[Iterable[T0]], bool]),
        "bin": ConstFunctionIntr(signature=Fun[[int], str]),
        "bool": ConstFunctionIntr(signature=_bool_signature),
        "chr": ConstFunctionIntr(signature=Fun[[int], str]),
        "complex": ClassWithConstConstructor(
            CLASSES['complex'],
            signature=_complex_signature
        ),
        "dict": ClassWithReadOnceConstructor(
            CLASSES['dict'],
            signature=Union[
                Fun[[], Dict[T0, T1]],
                Fun[[Iterable[Tuple[T0, T1]]], Dict[T0, T1]],
            ],
        ),
        "divmod": ConstFunctionIntr(
            signature=Union[
                Fun[[int, int], Tuple[int, int]],
                Fun[[float, int], Tuple[float, float]],
                Fun[[int, float], Tuple[float, float]],
                Fun[[float, float], Tuple[float, float]],
            ],
        ),
        "enumerate": ReadOnceFunctionIntr(
            signature=Union[
                Fun[[Iterable[T0]], Generator[Tuple[int, T0]]],
                Fun[[Iterable[T0], int], Generator[Tuple[int, T0]]],
            ],
        ),
        "filter": ReadOnceFunctionIntr(
            signature=Union[
                Fun[[None, Iterable[T0]], List[T0]],
                Fun[[Fun[[T0], bool], Iterable[T0]], List[T0]],
            ],
        ),
        "float": ClassWithConstConstructor(
            CLASSES['float'],
            signature=_float_signature
        ),
        "getattr": ConstFunctionIntr(),
        "hex": ConstFunctionIntr(signature=Fun[[int], str]),
        "id": ConstFunctionIntr(signature=Fun[[T0], int]),
        "int": ConstFunctionIntr(signature=_int_signature),
        "isinstance": ConstFunctionIntr(signature=Fun[[T0, T1], bool]),
        "iter": FunctionIntr(
            signature=Fun[[Iterable[T0]], Generator[T0]]),  # not const
        "len": ConstFunctionIntr(
            signature=Fun[[Sized], int],
            return_range=interval.positive_values
        ),
        "list": ClassWithReadOnceConstructor(
            CLASSES['list'],
            signature=Union[
                Fun[[], List[T0]],
                Fun[[Iterable[T0]], List[T0]]
            ],
        ),
        "map": ReadOnceFunctionIntr(
            signature=Union[
                Fun[[None, Iterable[T0]], List[T0]],
                Fun[[None, Iterable[T0], Iterable[T1]], List[Tuple[T0, T1]]],
                Fun[[None, Iterable[T0], Iterable[T1], Iterable[T2]],
                    List[Tuple[T0, T1, T2]]],
                Fun[[None, Iterable[T0], Iterable[T1], Iterable[T2],
                     Iterable[T3]], List[Tuple[T0, T1, T2, T3]]],
                Fun[[Fun[[T0], T7], Iterable[T0]], List[T7]],
                Fun[[Fun[[T0, T1], T7], Iterable[T0], Iterable[T1]], List[T7]],
                Fun[[Fun[[T0, T1, T2], T7], Iterable[T0], Iterable[T1],
                     Iterable[T2]], List[T7]],
                Fun[[Fun[[T0, T1, T2, T3], T7], Iterable[T0], Iterable[T1],
                     Iterable[T2], Iterable[T3]], List[T7]],
            ]
        ),
        "max": ReadOnceFunctionIntr(
            kwonlyargs=('key',),
            signature=Union[
                Fun[[T0, T0], T0],
                Fun[[T0, T0, T0], T0],
                Fun[[T0, T0, T0, T0], T0],
                Fun[[Iterable[T0]], T0],
            ],
            return_range=interval.max_values
        ),
        "min": ReadOnceFunctionIntr(
            kwonlyargs=('key', 'default'),
            signature=Union[
                Fun[[int, int], int],
                Fun[[float, float], float],
                Fun[[Iterable[T0]], T0],
            ],
            return_range=interval.min_values
        ),
        "next": FunctionIntr(  # not const
            signature=Union[
                Fun[[Iterable[T0]], T0],
                Fun[[Iterable[T0], T0], T0],
            ],
        ),  # not const
        "oct": ConstFunctionIntr(signature=Fun[[int], str]),
        "ord": ConstFunctionIntr(
            signature=Fun[[str], int],
            return_range=interval.ord_values
        ),
        "open": ConstFunctionIntr(
            signature=Union[
                Fun[[str], File],
                Fun[[str, str], File],
            ],
            global_effects=True
        ),
        "print": ConstFunctionIntr(global_effects=True),
        "pow": ConstFunctionIntr(
            signature=Union[
                Fun[[int, int], int],
                Fun[[int, int, int], int],
                Fun[[int, float], float],
                Fun[[int, float, int], float],
                Fun[[float, float], float],
                Fun[[float, float, int], float],
            ],
            immediate_arguments=[1]
        ),
        "range": ConstFunctionIntr(
            signature=Union[
                Fun[[int], List[int]],
                Fun[[int, int], List[int]],
                Fun[[int, int, int], List[int]],
            ],
            return_range_content=interval.range_values
        ),
        "reduce": ReadOnceFunctionIntr(signature=_functools_reduce_signature),
        "reversed": ReadOnceFunctionIntr(
            signature=Fun[[Iterable[T0]], Iterable[T0]]
        ),
        "round": ConstFunctionIntr(
            signature=Union[
                Fun[[float], float],
                Fun[[float, int], float],
            ],
        ),
        "set": ClassWithReadOnceConstructor(
            CLASSES['set'],
            signature=Union[
                Fun[[], Set[T0]],
                Fun[[Iterable[T0]], Set[T0]]
            ],
        ),
        "slice": ClassWithConstConstructor(CLASSES['slice']),
        "sorted": ConstFunctionIntr(signature=Fun[[Iterable[T0]], List[T0]]),
        "str": ClassWithConstConstructor(
            CLASSES['str'],
            signature=Fun[[T0], str],
        ),
        "sum": ReadOnceFunctionIntr(
            signature=Union[
                Fun[[Iterable[int]], int],
                Fun[[Iterable[float]], float],
                Fun[[Iterable[int], int], int],
                Fun[[Iterable[float], float], float],
            ],
        ),
        "tuple": ReadOnceFunctionIntr(
            signature=Union[
                Fun[[], Tuple[()]],
                Fun[[Tuple[T0]], Tuple[T0]],
                Fun[[Tuple[T0, T1]], Tuple[T0, T1]],
                Fun[[Tuple[T0, T1, T2]], Tuple[T0, T1, T2]],
                # FIXME: We accept some type loss here
                Fun[[List[T0]], Iterable[T0]],
            ],
        ),
        "type": ConstFunctionIntr(),
        "zip": ReadOnceFunctionIntr(
            signature=Union[
                Fun[[], List[T0]],
                Fun[[Iterable[T0]], List[Tuple[T0]]],
                Fun[[Iterable[T0], Iterable[T1]], List[Tuple[T0, T1]]],
                Fun[[Iterable[T0], Iterable[T1], Iterable[T2]],
                    List[Tuple[T0, T1, T2]]],
                Fun[[Iterable[T0], Iterable[T1], Iterable[T2], Iterable[T3]],
                    List[Tuple[T0, T1, T2, T3]]],
            ]
        ),
        "False": ConstantIntr(
            signature=bool,
            return_range=lambda _: interval.Range(0, 0)
        ),
        "None": ConstantIntr(signature=None),
        "True": ConstantIntr(
            signature=bool,
            return_range=lambda _: interval.Range(1, 1)
        ),
    },
    "array": {
            "typecodes": ConstantIntr(signature=str),
            "array": ClassWithReadOnceConstructor(
                CLASSES['array'],
                immediate_arguments=[0]),
    },
    "scipy": {
        "special": {
            "binom": ConstFunctionIntr(
                signature=_numpy_binary_op_float_signature
            ),
            "gammaincinv": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "gammaln": ConstFunctionIntr(
                signature=_numpy_unary_op_float_signature
            ),
            "gamma": ConstFunctionIntr(
                signature=_numpy_unary_op_float_signature
            ),
            "hankel1": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "hankel2": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "i0": ConstFunctionIntr(
                signature=_numpy_unary_op_float_signature
            ),
            "i0e": ConstFunctionIntr(
                signature=_numpy_unary_op_float_signature
            ),
            "iv": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "ivp": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "jv": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "jvp": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "kv": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "kvp": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "ndtr": ConstFunctionIntr(
                signature=_numpy_unary_op_float_signature
            ),
            "ndtri": ConstFunctionIntr(
                signature=_numpy_unary_op_float_signature
            ),
            "yv": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "yvp": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "spherical_jn": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
            "spherical_yn": UFunc(
                BINARY_UFUNC,
                signature=_numpy_binary_op_float_signature
            ),
        }
    },
    "numpy": {
        "abs": ConstFunctionIntr(signature=_numpy_unary_op_signature),
        "absolute": ConstFunctionIntr(signature=_numpy_ones_signature),
        "add": UFunc(
            REDUCED_BINARY_UFUNC,
            signature=_numpy_binary_op_signature,
        ),
        "alen": ConstFunctionIntr(
            signature=Union[
                # scalar
                Fun[[bool], int],
                Fun[[int], int],
                Fun[[float], int],
                Fun[[complex], int],
                # Sized
                Fun[[Sized], int],
            ],
            return_range=interval.positive_values
        ),
        "all": ConstMethodIntr(
            signature=_numpy_unary_op_bool_axis_signature,
            return_range=interval.bool_values
        ),
        "allclose": ConstFunctionIntr(
            signature=_numpy_allclose_signature,
            return_range=interval.bool_values
        ),
        "alltrue": ConstFunctionIntr(
            signature=_numpy_unary_op_bool_axis_signature,
            return_range=interval.bool_values,
            args=("a", "axis"),
            defaults=(None,)
        ),
        "amax": ConstFunctionIntr(signature=_numpy_unary_op_axis_signature),
        "amin": ConstFunctionIntr(signature=_numpy_unary_op_axis_signature),
        "angle": ConstFunctionIntr(signature=_numpy_unary_op_angle_signature),
        "any": ConstMethodIntr(
            signature=_numpy_unary_op_bool_axis_signature,
            return_range=interval.bool_values
        ),
        "append": ConstFunctionIntr(
            signature=Union[
                # no axis -> flattened output
                # scalar
                Fun[[bool, bool], NDArray[bool, :]],
                Fun[[int, int], NDArray[int, :]],
                Fun[[float, float], NDArray[float, :]],
                Fun[[complex, float], NDArray[float, :]],
                # 1D Array
                # FIXME: second argument could have a shape larger than first
                Fun[[Iterable[bool], bool], NDArray[bool, :]],
                Fun[[Iterable[int], int], NDArray[int, :]],
                Fun[[Iterable[float], float], NDArray[float, :]],
                Fun[[Iterable[complex], complex], NDArray[complex, :]],

                Fun[[Iterable[bool], Iterable[bool]], NDArray[bool, :]],
                Fun[[Iterable[int], Iterable[int]], NDArray[int, :]],
                Fun[[Iterable[float], Iterable[float]], NDArray[float, :]],
                Fun[[Iterable[complex], Iterable[complex]],
                    NDArray[complex, :]],

                # 2D Array
                Fun[[Iterable[Iterable[bool]], bool], NDArray[bool, :]],
                Fun[[Iterable[Iterable[int]], int], NDArray[int, :]],
                Fun[[Iterable[Iterable[float]], float], NDArray[float, :]],
                Fun[[Iterable[Iterable[complex]], complex],
                    NDArray[complex, :]],

                Fun[[Iterable[Iterable[bool]], Iterable[bool]],
                    NDArray[bool, :]],
                Fun[[Iterable[Iterable[int]], Iterable[int]], NDArray[int, :]],
                Fun[[Iterable[Iterable[float]], Iterable[float]],
                    NDArray[float, :]],
                Fun[[Iterable[Iterable[complex]], Iterable[
                    complex]], NDArray[complex, :]],

                Fun[[Iterable[Iterable[bool]], Iterable[
                    Iterable[bool]]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[int]], Iterable[
                    Iterable[int]]], NDArray[int, :]],
                Fun[[Iterable[Iterable[float]], Iterable[
                    Iterable[float]]], NDArray[float, :]],
                Fun[[Iterable[Iterable[complex]], Iterable[
                    Iterable[complex]]], NDArray[complex, :]],

                Fun[[bool, Iterable[Iterable[bool]]], NDArray[bool, :]],
                Fun[[int, Iterable[Iterable[int]]], NDArray[int, :]],
                Fun[[float, Iterable[Iterable[float]]], NDArray[float, :]],
                Fun[[complex, Iterable[Iterable[complex]]],
                    NDArray[complex, :]],

                Fun[[Iterable[bool], Iterable[Iterable[bool]]],
                    NDArray[bool, :]],
                Fun[[Iterable[int], Iterable[Iterable[int]]], NDArray[int, :]],
                Fun[[Iterable[float], Iterable[Iterable[float]]],
                    NDArray[float, :]],
                Fun[[Iterable[complex], Iterable[Iterable[complex]]],
                    NDArray[complex, :]],

                Fun[[Iterable[Iterable[bool]], Iterable[
                    Iterable[bool]]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[int]], Iterable[
                    Iterable[int]]], NDArray[int, :]],
                Fun[[Iterable[Iterable[float]], Iterable[
                    Iterable[float]]], NDArray[float, :]],
                Fun[[Iterable[Iterable[complex]], Iterable[
                    Iterable[complex]]], NDArray[complex, :]],

                # 3D Array FIXME: same as above
                Fun[[Iterable[Iterable[Iterable[bool]]], bool],
                    NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[int]]], int], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[float]]], float],
                    NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[complex]]],
                     complex], NDArray[complex, :]],

                Fun[[Iterable[Iterable[Iterable[bool]]],
                     Iterable[bool]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[int]]],
                     Iterable[int]], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[float]]],
                     Iterable[float]], NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[complex]]],
                     Iterable[complex]], NDArray[complex, :]],

                Fun[[Iterable[Iterable[Iterable[bool]]],
                     Iterable[Iterable[bool]]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[int]]],
                     Iterable[Iterable[int]]], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
                    Iterable[float]]], NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
                    Iterable[complex]]], NDArray[complex, :]],

                Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
                    Iterable[Iterable[bool]]]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
                    Iterable[Iterable[int]]]], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
                    Iterable[Iterable[float]]]], NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
                    Iterable[Iterable[complex]]]], NDArray[complex, :]],

                # 4D Array FIXME: same as above
                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], bool],
                    NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], int],
                    NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], float],
                    NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
                     complex], NDArray[complex, :]],

                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
                     Iterable[bool]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
                     Iterable[int]], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
                     Iterable[float]], NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
                     Iterable[complex]], NDArray[complex, :]],

                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
                     Iterable[Iterable[bool]]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
                     Iterable[Iterable[int]]], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
                     Iterable[Iterable[float]]], NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
                     Iterable[Iterable[complex]]], NDArray[complex, :]],

                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
                     Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
                     Iterable[Iterable[Iterable[int]]]], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], Iterable[
                    Iterable[Iterable[float]]]], NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], Iterable[
                    Iterable[Iterable[complex]]]], NDArray[complex, :]],

                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]], Iterable[
                    Iterable[Iterable[Iterable[bool]]]]], NDArray[bool, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]], Iterable[
                    Iterable[Iterable[Iterable[int]]]]], NDArray[int, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]], Iterable[
                    Iterable[Iterable[Iterable[float]]]]], NDArray[float, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]], Iterable[
                    Iterable[Iterable[Iterable[complex]]]]],
                    NDArray[complex, :]],

                # FIXME: same as above with None axis

                # axis -> same dims

                # 1D
                Fun[[Iterable[bool], Iterable[bool], int], Iterable[bool]],
                Fun[[Iterable[int], Iterable[int], int], Iterable[int]],
                Fun[[Iterable[float], Iterable[float], int], Iterable[float]],
                Fun[[Iterable[complex], Iterable[complex], int],
                    Iterable[complex]],

                # 2D
                Fun[[Iterable[Iterable[bool]], Iterable[
                    Iterable[bool]], int], Iterable[Iterable[bool]]],
                Fun[[Iterable[Iterable[int]], Iterable[
                    Iterable[int]], int], Iterable[Iterable[int]]],
                Fun[[Iterable[Iterable[float]], Iterable[
                    Iterable[float]], int], Iterable[Iterable[float]]],
                Fun[[Iterable[Iterable[complex]], Iterable[
                    Iterable[complex]], int], Iterable[Iterable[complex]]],

                # 3D
                Fun[[Iterable[Iterable[Iterable[bool]]], Iterable[
                    Iterable[Iterable[bool]]], int],
                    Iterable[Iterable[Iterable[bool]]]],
                Fun[[Iterable[Iterable[Iterable[int]]], Iterable[
                    Iterable[Iterable[int]]], int],
                    Iterable[Iterable[Iterable[int]]]],
                Fun[[Iterable[Iterable[Iterable[float]]], Iterable[
                    Iterable[Iterable[float]]], int],
                    Iterable[Iterable[Iterable[float]]]],
                Fun[[Iterable[Iterable[Iterable[complex]]], Iterable[
                    Iterable[Iterable[complex]]], int],
                    Iterable[Iterable[Iterable[complex]]]],

                # 4D
                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]],
                     Iterable[Iterable[Iterable[Iterable[bool]]]], int],
                    Iterable[Iterable[Iterable[Iterable[bool]]]]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]],
                     Iterable[Iterable[Iterable[Iterable[int]]]], int],
                    Iterable[Iterable[Iterable[Iterable[int]]]]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]],
                     Iterable[Iterable[Iterable[Iterable[float]]]], int],
                    Iterable[Iterable[Iterable[Iterable[float]]]]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]],
                     Iterable[Iterable[Iterable[Iterable[complex]]]],
                     int],
                    Iterable[Iterable[Iterable[Iterable[complex]]]]],
            ]
        ),
        "arange": ConstFunctionIntr(
            signature=Union[
                Fun[[float], NDArray[float, :]],
                Fun[[float, float], NDArray[float, :]],
                Fun[[float, float, float], NDArray[float, :]],
                Fun[[float, float, float, None], NDArray[float, :]],
                Fun[[float, float, float, _bool_signature], NDArray[bool, :]],
                Fun[[float, float, float, _int_signature], NDArray[int, :]],
                Fun[[float, float, float, _float_signature],
                    NDArray[float, :]],
                Fun[[float, float, float, _complex_signature],
                    NDArray[complex, :]],
            ],
            return_range_content=interval.range_values,
            args=('start', 'stop', 'step', 'dtype'),
            defaults=(1, None)
        ),
        "arccos": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "arccosh": ConstFunctionIntr(
            signature=_numpy_unary_op_float_signature),
        "arcsin": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "arcsinh": ConstFunctionIntr(
            signature=_numpy_unary_op_float_signature),
        "arctan": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "arctan2": UFunc(
            BINARY_UFUNC,
            signature=_numpy_binary_op_float_no_complex_signature
        ),
        "arctanh": ConstFunctionIntr(
            signature=_numpy_unary_op_float_signature),
        "argmax": ConstMethodIntr(
            signature=_numpy_unary_op_int_axis_signature,
            return_range=interval.positive_values
        ),
        "argmin": ConstMethodIntr(
            signature=_numpy_unary_op_int_axis_signature,
            return_range=interval.positive_values
        ),
        "argsort": ConstMethodIntr(
            signature=_numpy_unary_op_int_axis_signature,
            return_range=interval.positive_values
        ),
        "argwhere": ConstFunctionIntr(
            signature=_numpy_unary_op_int_signature,
            return_range=interval.positive_values
        ),
        "around": ConstFunctionIntr(signature=_numpy_around_signature),
        "array": FunctionIntr(signature=_numpy_array_signature,
                              args=('object', 'dtype'),
                              defaults=(None,)),
        "array2string": ConstFunctionIntr(
            signature=_numpy_array_str_signature),
        "array_equal": ConstFunctionIntr(signature=Fun[[T0, T1], bool]),
        "array_equiv": ConstFunctionIntr(signature=Fun[[T0, T1], bool]),
        "array_split": ConstFunctionIntr(
            signature=Union[
                # int split
                Fun[[NDArray[T0, :], int], List[NDArray[T0, :]]],
                # array split
                Fun[[NDArray[T0, :], Iterable[int]], List[NDArray[T0, :]]],
            ]
        ),
        "array_str": ConstFunctionIntr(signature=_numpy_array_str_signature),
        "asarray": ReadOnceFunctionIntr(signature=_numpy_array_signature,
                                        **extra_numpy_asarray_descr),
        "asarray_chkfinite": ConstFunctionIntr(
            signature=_numpy_array_signature),
        "ascontiguousarray": ConstFunctionIntr(
            signature=_numpy_array_signature),
        "asfarray": ConstFunctionIntr(signature=_numpy_farray_signature),
        "asscalar": ConstFunctionIntr(
            signature=Union[
                Fun[[NDArray[bool, :]], bool],
                Fun[[NDArray[int, :]], int],
                Fun[[NDArray[float, :]], float],
                Fun[[NDArray[complex, :]], complex],
                Fun[[NDArray[bool, :, :]], bool],
                Fun[[NDArray[int, :, :]], int],
                Fun[[NDArray[float, :, :]], float],
                Fun[[NDArray[complex, :, :]], complex],
                Fun[[NDArray[bool, :, :, :]], bool],
                Fun[[NDArray[int, :, :, :]], int],
                Fun[[NDArray[float, :, :, :]], float],
                Fun[[NDArray[complex, :, :, :]], complex],
            ]
        ),
        "atleast_1d": ConstFunctionIntr(
            signature=Union[
                # scalar
                Fun[[bool], NDArray[bool, :]],
                Fun[[int], NDArray[int, :]],
                Fun[[float], NDArray[float, :]],
                Fun[[complex], NDArray[complex, :]],
                # 1d
                Fun[[Iterable[bool]], NDArray[bool, :]],
                Fun[[Iterable[int]], NDArray[int, :]],
                Fun[[Iterable[float]], NDArray[float, :]],
                Fun[[Iterable[complex]], NDArray[complex, :]],
                # 2d+
                Fun[[NDArray[T0, :, :]], NDArray[T0, :, :]],
                Fun[[Iterable[Iterable[bool]]], NDArray[bool, :, :]],
                Fun[[Iterable[Iterable[int]]], NDArray[int, :, :]],
                Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
                Fun[[Iterable[Iterable[complex]]], NDArray[complex, :, :]],
                # 3d
                Fun[[Iterable[Iterable[Iterable[bool]]]],
                    NDArray[bool, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[int]]]],
                    NDArray[int, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[float]]]],
                    NDArray[float, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :, :]],
                # 4d
                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
                    NDArray[bool, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
                    NDArray[int, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
                    NDArray[float, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
                    NDArray[complex, :, :, :, :]],
            ]
        ),
        "atleast_2d": ConstFunctionIntr(
            signature=Union[
                # scalar
                Fun[[bool], NDArray[bool, :, :]],
                Fun[[int], NDArray[int, :, :]],
                Fun[[float], NDArray[float, :, :]],
                Fun[[complex], NDArray[complex, :, :]],
                # 1d
                Fun[[Iterable[bool]], NDArray[bool, :, :]],
                Fun[[Iterable[int]], NDArray[int, :, :]],
                Fun[[Iterable[float]], NDArray[float, :, :]],
                Fun[[Iterable[complex]], NDArray[complex, :, :]],
                # 2d+
                Fun[[NDArray[T0, :, :]], NDArray[T0, :, :]],
                Fun[[Iterable[Iterable[bool]]], NDArray[bool, :, :]],
                Fun[[Iterable[Iterable[int]]], NDArray[int, :, :]],
                Fun[[Iterable[Iterable[float]]], NDArray[float, :, :]],
                Fun[[Iterable[Iterable[complex]]], NDArray[complex, :, :]],
                # 3d
                Fun[[Iterable[Iterable[Iterable[bool]]]],
                    NDArray[bool, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[int]]]],
                    NDArray[int, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[float]]]],
                    NDArray[float, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :, :]],
                # 4d
                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
                    NDArray[bool, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
                    NDArray[int, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
                    NDArray[float, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
                    NDArray[complex, :, :, :, :]],
            ]
        ),
        "atleast_3d": ConstFunctionIntr(
            signature=Union[
                # scalar
                Fun[[bool], NDArray[bool, :, :, :]],
                Fun[[int], NDArray[int, :, :, :]],
                Fun[[float], NDArray[float, :, :, :]],
                Fun[[complex], NDArray[complex, :, :, :]],
                # 1d
                Fun[[Iterable[bool]], NDArray[bool, :, :, :]],
                Fun[[Iterable[int]], NDArray[int, :, :, :]],
                Fun[[Iterable[float]], NDArray[float, :, :, :]],
                Fun[[Iterable[complex]], NDArray[complex, :, :, :]],
                # 2d+
                Fun[[NDArray[T0, :, :]], NDArray[T0, :, :, :]],
                Fun[[Iterable[Iterable[bool]]], NDArray[bool, :, :, :]],
                Fun[[Iterable[Iterable[int]]], NDArray[int, :, :, :]],
                Fun[[Iterable[Iterable[float]]], NDArray[float, :, :, :]],
                Fun[[Iterable[Iterable[complex]]], NDArray[complex, :, :, :]],
                # 3d
                Fun[[Iterable[Iterable[Iterable[bool]]]],
                    NDArray[bool, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[int]]]],
                    NDArray[int, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[float]]]],
                    NDArray[float, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :, :]],
                # 4d
                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]],
                    NDArray[bool, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]],
                    NDArray[int, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]],
                    NDArray[float, :, :, :, :]],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]],
                    NDArray[complex, :, :, :, :]],
            ]
        ),
        "average": ConstFunctionIntr(
            signature=_numpy_unary_op_average_axis_signature),
        "base_repr": ConstFunctionIntr(
            signature=Union[
                Fun[[bool], str],
                Fun[[bool, int], str],
                Fun[[bool, int, int], str],
                Fun[[int], str],
                Fun[[int, int], str],
                Fun[[int, int, int], str],
            ]
        ),
        "binary_repr": ConstFunctionIntr(
            signature=Union[
                Fun[[bool], str],
                Fun[[bool, int], str],
                Fun[[bool, None], str],
                Fun[[int], str],
                Fun[[int, int], str],
                Fun[[int, None], str],
            ]
        ),
        "bincount": ConstFunctionIntr(
            signature=Union[
                Fun[[Iterable[bool]], NDArray[int, :]],
                Fun[[Iterable[int]], NDArray[int, :]],
                #
                Fun[[Iterable[bool], Iterable[float]], NDArray[int, :]],
                Fun[[Iterable[int], Iterable[float]], NDArray[int, :]],
                #
                Fun[[Iterable[bool], Iterable[float], int], NDArray[int, :]],
                Fun[[Iterable[int], Iterable[float], int], NDArray[int, :]],
            ],
        ),
        "bitwise_and": UFunc(
            REDUCED_BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature
        ),
        "bitwise_not": ConstFunctionIntr(
            signature=_numpy_int_unary_op_signature
        ),
        "bitwise_or": UFunc(
            REDUCED_BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature
        ),
        "bitwise_xor": UFunc(
            REDUCED_BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature
        ),
        "bool_": ConstFunctionIntr(signature=_int_signature),
        "broadcast_to": ConstFunctionIntr(),
        "byte": ConstFunctionIntr(signature=_int_signature),
        "cbrt": ConstFunctionIntr(
            signature=_numpy_unary_op_float_signature
        ),
        "ceil": ConstFunctionIntr(signature=_numpy_float_unary_op_signature),
        "clip": ConstMethodIntr(signature=_numpy_ternary_op_signature),
        "concatenate": ConstFunctionIntr(
            args=('_', 'axis'),
            defaults=(0,),
            signature=Union[
                # 1D
                Fun[[Iterable[Iterable[bool]]], NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool]]], NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool], int]], NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool], Iterable[bool]]], NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool], Iterable[bool], int]],
                    NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool], Iterable[bool],
                           Iterable[bool]]], NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool], Iterable[bool],
                           Iterable[bool], int]], NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool], Iterable[bool], Iterable[
                    bool], Iterable[bool]]], NDArray[bool, :]],
                Fun[[Tuple[Iterable[bool], Iterable[bool], Iterable[
                    bool], Iterable[bool], int]], NDArray[bool, :]],

                Fun[[Iterable[Iterable[int]]], NDArray[int, :]],
                Fun[[Tuple[Iterable[int]]], NDArray[int, :]],
                Fun[[Tuple[Iterable[int], int]], NDArray[int, :]],
                Fun[[Tuple[Iterable[int], Iterable[int]]], NDArray[int, :]],
                Fun[[Tuple[Iterable[int], Iterable[int], int]],
                    NDArray[int, :]],
                Fun[[Tuple[Iterable[int], Iterable[int], Iterable[int]]],
                    NDArray[int, :]],
                Fun[[Tuple[Iterable[int], Iterable[int],
                           Iterable[int], int]], NDArray[int, :]],
                Fun[[Tuple[Iterable[int], Iterable[int], Iterable[
                    int], Iterable[int]]], NDArray[int, :]],
                Fun[[Tuple[Iterable[int], Iterable[int], Iterable[
                    int], Iterable[int], int]], NDArray[int, :]],

                Fun[[Iterable[Iterable[float]]], NDArray[float, :]],
                Fun[[Tuple[Iterable[float]]], NDArray[float, :]],
                Fun[[Tuple[Iterable[float], int]], NDArray[float, :]],
                Fun[[Tuple[Iterable[float], Iterable[float]]],
                    NDArray[float, :]],
                Fun[[Tuple[Iterable[float], Iterable[float], int]],
                    NDArray[float, :]],
                Fun[[Tuple[Iterable[float], Iterable[float],
                           Iterable[float]]], NDArray[float, :]],
                Fun[[Tuple[Iterable[float], Iterable[float],
                           Iterable[float], int]], NDArray[float, :]],
                Fun[[Tuple[Iterable[float], Iterable[float], Iterable[
                    float], Iterable[float]]], NDArray[float, :]],
                Fun[[Tuple[Iterable[float], Iterable[float], Iterable[
                    float], Iterable[float], int]], NDArray[float, :]],

                Fun[[Iterable[Iterable[complex]]], NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex]]], NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex], int]], NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex], Iterable[complex]]],
                    NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex], Iterable[complex], int]],
                    NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex], Iterable[complex],
                           Iterable[complex]]], NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex], Iterable[complex],
                           Iterable[complex], int]], NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex], Iterable[complex], Iterable[
                    complex], Iterable[complex]]], NDArray[complex, :]],
                Fun[[Tuple[Iterable[complex], Iterable[complex], Iterable[
                    complex], Iterable[complex], int]], NDArray[complex, :]],

                # 2D
                Fun[[Iterable[Iterable[Iterable[bool]]]], NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]]]], NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]], int]],
                    NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]], Iterable[
                    Iterable[bool]]]], NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]], Iterable[
                    Iterable[bool]], int]], NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]], Iterable[Iterable[bool]],
                           Iterable[Iterable[bool]]]], NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]], Iterable[Iterable[bool]],
                           Iterable[Iterable[bool]], int]],
                    NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]], Iterable[Iterable[bool]],
                           Iterable[Iterable[bool]],
                           Iterable[Iterable[bool]]]],
                    NDArray[bool, :, :]],
                Fun[[Tuple[Iterable[Iterable[bool]],
                           Iterable[Iterable[bool]],
                           Iterable[Iterable[bool]],
                           Iterable[Iterable[bool]], int]],
                    NDArray[bool, :, :]],

                Fun[[Iterable[Iterable[Iterable[int]]]], NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]]]], NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]], int]], NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]], Iterable[
                    Iterable[int]]]], NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]], Iterable[
                    Iterable[int]], int]], NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]], Iterable[Iterable[int]],
                           Iterable[Iterable[int]]]], NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]], Iterable[Iterable[int]],
                           Iterable[Iterable[int]], int]], NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]], Iterable[Iterable[int]],
                           Iterable[Iterable[int]], Iterable[Iterable[int]]]],
                    NDArray[int, :, :]],
                Fun[[Tuple[Iterable[Iterable[int]], Iterable[Iterable[int]],
                           Iterable[Iterable[int]], Iterable[Iterable[int]],
                           int]],
                    NDArray[int, :, :]],

                Fun[[Iterable[Iterable[Iterable[float]]]],
                    NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]]]], NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]], int]],
                    NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]],
                           Iterable[Iterable[float]]]], NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]],
                           Iterable[Iterable[float]], int]],
                    NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]],
                           Iterable[Iterable[float]],
                           Iterable[Iterable[float]]]],
                    NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]],
                           Iterable[Iterable[float]],
                           Iterable[Iterable[float]], int]],
                    NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]],
                           Iterable[Iterable[float]],
                           Iterable[Iterable[float]],
                           Iterable[Iterable[float]]]], NDArray[float, :, :]],
                Fun[[Tuple[Iterable[Iterable[float]],
                           Iterable[Iterable[float]],
                           Iterable[Iterable[float]],
                           Iterable[Iterable[float]], int]],
                    NDArray[float, :, :]],

                Fun[[Iterable[Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]], int]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]], int]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]], int]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]]]],
                    NDArray[complex, :, :]],
                Fun[[Tuple[Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]],
                           Iterable[Iterable[complex]], int]],
                    NDArray[complex, :, :]],
            ]
        ),
        "complex64": ConstFunctionIntr(signature=_complex_signature),
        "complex128": ConstFunctionIntr(signature=_complex_signature),
        "complex256": ConstFunctionIntr(signature=_complex_signature),
        "conj": ConstMethodIntr(signature=_numpy_unary_op_signature),
        "conjugate": ConstMethodIntr(signature=_numpy_unary_op_signature),
        "convolve": ConstMethodIntr(requires_blas=True),
        "correlate": ConstMethodIntr(requires_blas=True),
        "copy": ConstMethodIntr(signature=_numpy_array_signature),
        "copyto": FunctionIntr(
            argument_effects=[UpdateEffect(), ReadEffect(),
                              ReadEffect(), ReadEffect()],
            signature=Union[
                # 1d
                Fun[[NDArray[bool, :], bool], None],
                Fun[[NDArray[bool, :], Iterable[bool]], None],
                Fun[[NDArray[int, :], int], None],
                Fun[[NDArray[int, :], Iterable[int]], None],
                Fun[[NDArray[float, :], float], None],
                Fun[[NDArray[float, :], Iterable[float]], None],
                Fun[[NDArray[complex, :], complex], None],
                Fun[[NDArray[complex, :], Iterable[complex]], None],

                # 2d
                Fun[[NDArray[bool, :, :], bool], None],
                Fun[[NDArray[bool, :, :], Iterable[bool]], None],
                Fun[[NDArray[bool, :, :], Iterable[Iterable[bool]]], None],
                Fun[[NDArray[int, :, :], int], None],
                Fun[[NDArray[int, :, :], Iterable[int]], None],
                Fun[[NDArray[int, :, :], Iterable[Iterable[int]]], None],
                Fun[[NDArray[float, :, :], float], None],
                Fun[[NDArray[float, :, :], Iterable[float]], None],
                Fun[[NDArray[float, :, :], Iterable[Iterable[float]]], None],
                Fun[[NDArray[complex, :, :], complex], None],
                Fun[[NDArray[complex, :, :], Iterable[complex]], None],
                Fun[[NDArray[complex, :, :],
                     Iterable[Iterable[complex]]], None],

                # 3d
                Fun[[NDArray[bool, :, :, :], bool], None],
                Fun[[NDArray[bool, :, :, :], Iterable[bool]], None],
                Fun[[NDArray[bool, :, :, :], Iterable[Iterable[bool]]], None],
                Fun[[NDArray[bool, :, :, :],
                     Iterable[Iterable[Iterable[bool]]]], None],
                Fun[[NDArray[int, :, :, :], int], None],
                Fun[[NDArray[int, :, :, :], Iterable[int]], None],
                Fun[[NDArray[int, :, :, :], Iterable[Iterable[int]]], None],
                Fun[[NDArray[int, :, :, :],
                     Iterable[Iterable[Iterable[int]]]], None],
                Fun[[NDArray[float, :, :, :], float], None],
                Fun[[NDArray[float, :, :, :],
                     Iterable[float]], None],
                Fun[[NDArray[float, :, :, :],
                     Iterable[Iterable[float]]], None],
                Fun[[NDArray[float, :, :, :],
                     Iterable[Iterable[Iterable[float]]]], None],
                Fun[[NDArray[complex, :, :, :],
                     complex], None],
                Fun[[NDArray[complex, :, :, :], Iterable[complex]], None],
                Fun[[NDArray[complex, :, :, :],
                     Iterable[Iterable[complex]]], None],
                Fun[[NDArray[complex, :, :, :],
                     Iterable[Iterable[Iterable[complex]]]], None],
            ]
        ),
        "copysign": UFunc(BINARY_UFUNC),
        "count_nonzero": ConstFunctionIntr(
            signature=Union[
                # scalar
                Fun[[bool], int],
                Fun[[int], int],
                Fun[[float], int],
                Fun[[complex], int],
                # 1d
                Fun[[Iterable[bool]], int],
                Fun[[Iterable[int]], int],
                Fun[[Iterable[float]], int],
                Fun[[Iterable[complex]], int],
                # 2d
                Fun[[Iterable[Iterable[bool]]], int],
                Fun[[Iterable[Iterable[int]]], int],
                Fun[[Iterable[Iterable[float]]], int],
                Fun[[Iterable[Iterable[complex]]], int],
                # 3d
                Fun[[Iterable[Iterable[Iterable[bool]]]], int],
                Fun[[Iterable[Iterable[Iterable[int]]]], int],
                Fun[[Iterable[Iterable[Iterable[float]]]], int],
                Fun[[Iterable[Iterable[Iterable[complex]]]], int],
                # 4d
                Fun[[Iterable[Iterable[Iterable[Iterable[bool]]]]], int],
                Fun[[Iterable[Iterable[Iterable[Iterable[int]]]]], int],
                Fun[[Iterable[Iterable[Iterable[Iterable[float]]]]], int],
                Fun[[Iterable[Iterable[Iterable[Iterable[complex]]]]], int],
            ]
        ),
        "cos": ConstFunctionIntr(
            signature=_numpy_unary_op_float_signature
        ),
        "cosh": ConstFunctionIntr(
            signature=_numpy_unary_op_float_signature
        ),
        "cross": ConstFunctionIntr(),
        "ctypeslib": {
            "as_array": ConstFunctionIntr()
        },
        "cumprod": ConstMethodIntr(
            signature=_numpy_unary_op_cumsum_axis_signature
        ),
        "cumproduct": ConstFunctionIntr(
            signature=_numpy_unary_op_cumsum_axis_signature,
            args=("a", "axis", "dtype", "out"),
            defaults=(None, None, None),
        ),
        "cumsum": ConstMethodIntr(
            signature=_numpy_unary_op_cumsum_axis_signature
        ),
        "deg2rad": ConstFunctionIntr(
            signature=_numpy_float_unary_op_float_signature
        ),
        "degrees": ConstFunctionIntr(
            signature=_numpy_float_unary_op_float_signature
        ),
        "delete": ConstFunctionIntr(),
        "diag": ConstFunctionIntr(),
        "diagflat": ConstFunctionIntr(),
        "diagonal": ConstMethodIntr(),
        "diff": ConstFunctionIntr(),
        "digitize": ConstFunctionIntr(),
        "divide": UFunc(BINARY_UFUNC),
        "dot": ConstMethodIntr(requires_blas=True),
        "double": ConstFunctionIntr(signature=_float_signature),
        "dtype": ClassWithConstConstructor(CLASSES["dtype"]),
        "e": ConstantIntr(),
        "ediff1d": ConstFunctionIntr(),
        "empty": ConstFunctionIntr(args=('shape', 'dtype'),
                                   defaults=("numpy.float64",),
                                   signature=_numpy_ones_signature,
                                   global_effects=True, # to avoid folding
                                   ),
        "empty_like": ConstFunctionIntr(
            args=('a', 'dtype'),
            defaults=("numpy.float64",),
            signature=_numpy_ones_like_signature
        ),
        "equal": UFunc(BINARY_UFUNC),
        "exp": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "expand_dims": ConstFunctionIntr(),
        "expm1": ConstFunctionIntr(),
        "eye": ConstFunctionIntr(),
        "fabs": ConstFunctionIntr(),
        "fill_diagonal": FunctionIntr(
            argument_effects=[UpdateEffect(), ReadEffect()],
            signature=Union[
                Fun[[NDArray[T0, :], bool], None],
                Fun[[NDArray[T0, :], int], None],
                Fun[[NDArray[T0, :], float], None],
                Fun[[NDArray[T0, :], complex], None],
            ]
        ),
        "finfo": ClassWithConstConstructor(CLASSES['finfo']),
        "fix": ConstFunctionIntr(),
        "flatnonzero": ConstFunctionIntr(),
        "fliplr": ConstFunctionIntr(),
        "flip": ConstFunctionIntr(),
        "flipud": ConstFunctionIntr(),
        "float32": ConstFunctionIntr(signature=_float_signature),
        "float64": ConstFunctionIntr(signature=_float_signature),
        "float128": ConstFunctionIntr(signature=_float_signature),
        "floor": ConstFunctionIntr(signature=_numpy_float_unary_op_signature),
        "floor_divide": UFunc(BINARY_UFUNC),
        "fmax": UFunc(REDUCED_BINARY_UFUNC),
        "fmin": UFunc(REDUCED_BINARY_UFUNC),
        "fmod": UFunc(BINARY_UFUNC),
        "frexp": ConstFunctionIntr(),
        "frombuffer": ConstFunctionIntr(),
        "fromfunction": ConstFunctionIntr(),
        "fromiter": ConstFunctionIntr(args=("iterable", "dtype", "count"),
                                      defaults=(-1,)),
        "fromstring": ConstFunctionIntr(args=('string', 'dtype', 'count',),
                                        kwonlyargs=('sep', 'like'),
                                        defaults=(float, -1, '', None)),
        "fromfile":  FunctionIntr(args=('file', 'dtype', 'count', "sep", "offset"),
                                  defaults=(None, None, -1, None, 0),
                                  global_effects=True),
        "full": ConstFunctionIntr(signature=_numpy_ones_signature),
        "full_like": ConstFunctionIntr(signature=_numpy_ones_like_signature),
        "greater": UFunc(
            BINARY_UFUNC,
            signature=_numpy_binary_op_bool_signature,
        ),
        "greater_equal": UFunc(
            BINARY_UFUNC,
            signature=_numpy_binary_op_bool_signature,
        ),
        "heaviside": UFunc(BINARY_UFUNC),
        "hstack": ConstFunctionIntr(),
        "hypot": UFunc(BINARY_UFUNC),
        "identity": ConstFunctionIntr(),
        "imag": FunctionIntr(),
        "indices": ConstFunctionIntr(),
        "inf": ConstantIntr(),
        "Inf": ConstantIntr(),
        "inner": ConstFunctionIntr(requires_blas=True),
        "insert": ConstFunctionIntr(),
        "interp": ConstFunctionIntr(),
        "intersect1d": ConstFunctionIntr(),
        "int16": ConstFunctionIntr(signature=_int_signature),
        "int32": ConstFunctionIntr(signature=_int_signature),
        "int64": ConstFunctionIntr(signature=_int_signature),
        "int8": ConstFunctionIntr(signature=_int_signature),
        "intc": ConstFunctionIntr(signature=_int_signature),
        "intp": ConstFunctionIntr(signature=_int_signature),
        "invert": ConstFunctionIntr(),
        "isclose": ConstFunctionIntr(),
        "iscomplex": ConstFunctionIntr(),
        "isfinite": ConstFunctionIntr(),
        "isinf": ConstFunctionIntr(),
        "isnan": ConstFunctionIntr(),
        "isneginf": ConstFunctionIntr(),
        "isposinf": ConstFunctionIntr(),
        "isreal": ConstFunctionIntr(),
        "isrealobj": ConstFunctionIntr(),
        "isscalar": ConstFunctionIntr(),
        "issctype": ConstFunctionIntr(),
        "ldexp": UFunc(BINARY_UFUNC),
        "left_shift": UFunc(
            BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature,
        ),
        "less": UFunc(
            BINARY_UFUNC,
            signature=_numpy_binary_op_bool_signature,
        ),
        "less_equal": UFunc(
            BINARY_UFUNC,
            signature=_numpy_binary_op_bool_signature,
        ),
        "lexsort": ConstFunctionIntr(),
        "linalg": {
            "norm": FunctionIntr(),
            "matrix_power": ConstFunctionIntr(requires_blas=True),
        },
        "linspace": ConstFunctionIntr(),
        "log": ConstFunctionIntr(),
        "log10": ConstFunctionIntr(),
        "log1p": ConstFunctionIntr(),
        "log2": ConstFunctionIntr(),
        "logaddexp": UFunc(BINARY_UFUNC),
        "logaddexp2": UFunc(BINARY_UFUNC),
        "logspace": ConstFunctionIntr(),
        "logical_and": UFunc(
            BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature
        ),
        "logical_not": ConstFunctionIntr(),
        "logical_or": UFunc(
            BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature
        ),
        "logical_xor": UFunc(
            BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature
        ),
        "longlong": ConstFunctionIntr(signature=_int_signature),
        "max": ConstMethodIntr(signature=_numpy_unary_op_axis_signature),
        "maximum": UFunc(
            REDUCED_BINARY_UFUNC,
            signature=_numpy_binary_op_signature
        ),
        "mean": ConstMethodIntr(immediate_arguments=[4]),
        "median": ConstFunctionIntr(
            signature=_numpy_unary_op_sum_axis_signature
        ),
        "min": ConstMethodIntr(signature=_numpy_unary_op_axis_signature),
        "minimum": UFunc(
            REDUCED_BINARY_UFUNC,
            signature=_numpy_binary_op_signature
        ),
        "mod": UFunc(BINARY_UFUNC),
        "multiply": UFunc(
            REDUCED_BINARY_UFUNC,
            signature=_numpy_binary_op_signature,
        ),
        "nan": ConstantIntr(),
        "nan_to_num": ConstFunctionIntr(),
        "nanargmax": ConstFunctionIntr(),
        "nanargmin": ConstFunctionIntr(),
        "nanmax": ConstFunctionIntr(),
        "nanmin": ConstFunctionIntr(),
        "nansum": ConstFunctionIntr(),
        "ndenumerate": ConstFunctionIntr(),
        "ndarray": ClassWithConstConstructor(CLASSES["ndarray"]),
        "ndindex": ConstFunctionIntr(),
        "ndim": ConstFunctionIntr(return_range=interval.positive_values),
        "negative": ConstFunctionIntr(signature=_numpy_unary_op_signature),
        "newaxis": ConstantIntr(),
        "nextafter": UFunc(BINARY_UFUNC),
        "NINF": ConstantIntr(),
        "nonzero": ConstMethodIntr(),
        "not_equal": UFunc(BINARY_UFUNC),
        "ones": ConstFunctionIntr(signature=_numpy_ones_signature),
        "ones_like": ConstFunctionIntr(signature=_numpy_ones_like_signature),
        "outer": ConstFunctionIntr(),
        "pi": ConstantIntr(),
        "place": FunctionIntr(),
        "power": UFunc(
            BINARY_UFUNC,
            signature=_numpy_binary_op_signature
        ),
        "prod": ConstMethodIntr(),
        "product": ConstFunctionIntr(
            args=("a", "axis", "dtype", "out"),
            defaults=(None, None, None),
                ),
        "ptp": ConstMethodIntr(),
        "put": MethodIntr(),
        "putmask": FunctionIntr(),
        "rad2deg": ConstFunctionIntr(
            signature=_numpy_float_unary_op_float_signature
        ),
        "radians": ConstFunctionIntr(
            signature=_numpy_float_unary_op_float_signature
        ),
        "fft": {
            "fft": FunctionIntr(global_effects=True),
            "fftn": FunctionIntr( global_effects=True),
            "ifft": FunctionIntr(global_effects=True),
            "rfft": FunctionIntr(global_effects=True),
            "irfft": FunctionIntr(global_effects=True),
            "hfft": FunctionIntr(global_effects=True),
            "ihfft": FunctionIntr(global_effects=True),
        },
        # Note: numpy.random signatures were improved upstream, so we only need
        # to add args/defaults for older versions.
        "random": {
            "binomial": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('n', 'p', 'size'))),
            "bytes": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('length',))),
            "chisquare": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('df', 'size',))),
            "choice": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('a', 'size', 'replace', 'p'))),
            "dirichlet": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('alpha', 'size',))),
            "exponential": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('scale', 'size',),
                                      defaults=(1.0, None,))),
            "f": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('dfnum', 'dfden', 'size'))),
            "gamma": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('shape', 'scale', 'size',),
                                      defaults=(0.0, 1.0, None,))),
            "geometric": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('p', 'size',))),
            "pareto": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('a', 'size',))),
            "gumbel": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('loc', 'scale', 'size',),
                                      defaults=(0.0, 1.0, None,))),
            "poisson": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('lam', 'size',),
                                      defaults=(1.0, None,))),
            "negative_binomial": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('n', 'p', 'size',))),
            "normal": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('loc', 'scale', 'size',),
                                      defaults=(0.0, 1.0, None,))),
            "laplace": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('loc', 'scale', 'size',),
                                      defaults=(0.0, 1.0, None,))),
            "logistic": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('loc', 'scale', 'size',),
                                      defaults=(0.0, 1.0, None,))),
            "lognormal": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('mean', 'sigma', 'size',),
                                      defaults=(0.0, 1.0, None,))),
            "logseries": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('p', 'size',))),
            "power": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('a', 'size',))),
            "rand": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=())),
            "ranf": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('size',), force=True)),
            "randint": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=("low", "high", "size"),
                                      defaults=(None, None))),
            "randn": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=())),
            "random": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('size',))),
            "random_integers": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=("low", "high", "size"))),
            "random_sample": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('size',))),
            "rayleigh": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('scale', 'size',),
                                      defaults=(1.0, None,))),
            "sample": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('size',), force=True)),
            "seed": FunctionIntr(global_effects=True),
            "shuffle": FunctionIntr(global_effects=True),
            "standard_exponential": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('size',))),
            "standard_gamma": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('shape', 'size',))),
            "standard_normal": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('size',))),
            "uniform": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('low', 'high', 'size',),
                                      defaults=(0.0, 1.0, None,))),
            "weibull": FunctionIntr(global_effects=True,
                **expand_numpy_2_args(args=('a', 'size',))),
        },
        # We currently don't accurately model the fact that the return of ravel
        # shares memory with its argument, but not the type.
        # Use a global effect to prevent the optimizer from removing it.
        "ravel": ConstMethodIntr(
            #return_alias=lambda args: {args[0]},
            global_effects=True,
        ),
        "real": FunctionIntr(),
        "reciprocal": ConstFunctionIntr(),
        "remainder": UFunc(BINARY_UFUNC),
        "repeat": ConstMethodIntr(),
        "resize": ConstMethodIntr(),
        "right_shift": UFunc(
            BINARY_UFUNC,
            signature=_numpy_int_binary_op_signature,
        ),
        "rint": ConstFunctionIntr(),
        "roll": ConstFunctionIntr(),
        "rollaxis": ConstFunctionIntr(),
        "rot90": ConstFunctionIntr(),
        "round": ConstMethodIntr(),
        "round_": ConstMethodIntr(),
        "searchsorted": ConstMethodIntr(),
        "select": ConstFunctionIntr(),
        "setdiff1d": ConstFunctionIntr(),
        "shape": ConstFunctionIntr(),
        "short_": ConstFunctionIntr(signature=_int_signature),
        "sign": ConstFunctionIntr(),
        "signbit": ConstFunctionIntr(),
        "sin": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "sinh": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "size": ConstFunctionIntr(return_range=interval.positive_values),
        "sometrue": ConstFunctionIntr(
            args=("a", "axis"),
            defaults=(None,)
                ),
        "sort": ConstFunctionIntr(),
        "sort_complex": ConstFunctionIntr(),
        "spacing": ConstFunctionIntr(),
        "split": ConstFunctionIntr(),
        "sqrt": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "square": ConstFunctionIntr(),
        "stack": ConstFunctionIntr(),
        "std": ConstMethodIntr(),
        "subtract": UFunc(
            BINARY_UFUNC,
            signature=_numpy_binary_op_signature,
        ),
        "sum": ReadOnceMethodIntr(
            signature=_numpy_unary_op_sum_axis_signature),
        "swapaxes": ConstMethodIntr(),
        "short": ConstFunctionIntr(signature=_int_signature),
        "take": ConstMethodIntr(),
        "tan": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "tanh": ConstFunctionIntr(signature=_numpy_unary_op_float_signature),
        "tile": ConstFunctionIntr(),
        "trace": ConstMethodIntr(),
        "transpose": ConstMethodIntr(),
        "tri": ConstMethodIntr(),
        "tril": ConstMethodIntr(),
        "trim_zeros": ConstMethodIntr(),
        "triu": ConstMethodIntr(),
        "true_divide": UFunc(BINARY_UFUNC),
        "trunc": ConstFunctionIntr(),
        "ubyte": ConstFunctionIntr(signature=_int_signature),
        "uint16": ConstFunctionIntr(signature=_int_signature),
        "uint32": ConstFunctionIntr(signature=_int_signature),
        "uint64": ConstFunctionIntr(signature=_int_signature),
        "uintc": ConstFunctionIntr(signature=_int_signature),
        "uintp": ConstFunctionIntr(signature=_int_signature),
        "uint8": ConstFunctionIntr(signature=_int_signature),
        "ulonglong": ConstFunctionIntr(signature=_int_signature),
        "union1d": ConstFunctionIntr(),
        "unique": ConstFunctionIntr(immediate_arguments=[1, 2, 3]),
        "unwrap": ConstFunctionIntr(),
        "unravel_index": ConstFunctionIntr(),
        "ushort": ConstFunctionIntr(signature=_int_signature),
        "var": ConstMethodIntr(),
        "vectorize": ConstFunctionIntr(),
        "vdot": ConstMethodIntr(requires_blas=True),
        "vstack": ConstFunctionIntr(),
        "where": ConstFunctionIntr(),
        "zeros": ConstFunctionIntr(args=('shape', 'dtype'),
                                   defaults=("numpy.float64",),
                                   signature=_numpy_ones_signature,
                                   ),
        "zeros_like": ConstFunctionIntr(signature=_numpy_ones_like_signature),
    },
    "time": {
        "sleep": FunctionIntr(
            signature=Fun[[float], None],
            global_effects=True
        ),
        "time": FunctionIntr(
            signature=Fun[[], float],
            global_effects=True
        ),
    },
    "math": {
        "isinf": ConstFunctionIntr(signature=Fun[[float], bool]),
        "modf": ConstFunctionIntr(signature=Fun[[float], Tuple[float, float]]),
        "frexp": ConstFunctionIntr(signature=Fun[[float], Tuple[float, int]]),
        "factorial": ConstFunctionIntr(signature=Fun[[int], int]),
        "gamma": ConstFunctionIntr(signature=Fun[[float], float]),
        "lgamma": ConstFunctionIntr(signature=Fun[[float], float]),
        "trunc": ConstFunctionIntr(signature=Fun[[float], int]),
        "erf": ConstFunctionIntr(signature=Fun[[float], float]),
        "erfc": ConstFunctionIntr(signature=Fun[[float], float]),
        "asinh": ConstFunctionIntr(signature=Fun[[float], float]),
        "atanh": ConstFunctionIntr(signature=Fun[[float], float]),
        "acosh": ConstFunctionIntr(signature=Fun[[float], float]),
        "radians": ConstFunctionIntr(signature=Fun[[float], float]),
        "degrees": ConstFunctionIntr(signature=Fun[[float], float]),
        "hypot": ConstFunctionIntr(signature=Fun[[float, float], float]),
        "tanh": ConstFunctionIntr(signature=Fun[[float], float]),
        "cosh": ConstFunctionIntr(signature=Fun[[float], float]),
        "sinh": ConstFunctionIntr(signature=Fun[[float], float]),
        "atan": ConstFunctionIntr(signature=Fun[[float], float]),
        "atan2": ConstFunctionIntr(signature=Fun[[float, float], float]),
        "asin": ConstFunctionIntr(signature=Fun[[float], float]),
        "tan": ConstFunctionIntr(signature=Fun[[float], float]),
        "log": ConstFunctionIntr(signature=Fun[[float], float]),
        "log1p": ConstFunctionIntr(signature=Fun[[float], float]),
        "expm1": ConstFunctionIntr(signature=Fun[[float], float]),
        "ldexp": ConstFunctionIntr(signature=Fun[[float, int], float]),
        "fmod": ConstFunctionIntr(signature=Fun[[float, float], float]),
        "fabs": ConstFunctionIntr(signature=Fun[[float], float]),
        "copysign": UFunc(BINARY_UFUNC),
        "acos": ConstFunctionIntr(signature=Fun[[float], float]),
        "cos": ConstFunctionIntr(signature=Fun[[float], float]),
        "sin": ConstFunctionIntr(signature=Fun[[float], float]),
        "exp": ConstFunctionIntr(signature=Fun[[float], float]),
        "sqrt": ConstFunctionIntr(signature=Fun[[float], float]),
        "log10": ConstFunctionIntr(signature=Fun[[float], float]),
        "isnan": ConstFunctionIntr(signature=Fun[[float], bool]),
        "ceil": ConstFunctionIntr(signature=Fun[[float], float]),
        "floor": ConstFunctionIntr(signature=Fun[[float], float]),
        "pow": ConstFunctionIntr(signature=Fun[[float, float], float]),
        "pi": ConstantIntr(signature=float),
        "e": ConstantIntr(signature=float),
    },
    "functools": {
        "partial": FunctionIntr(
            signature=Union[
                # no arg
                Fun[[Fun[[], T0]], Fun[[], T0]],
                # 1 arg
                Fun[[Fun[[T0], T1]], Fun[[T0], T1]],
                Fun[[Fun[[T0], T1], T0], Fun[[], T1]],
                # 2 args
                Fun[[Fun[[T0, T1], T2]], Fun[[T0, T1], T2]],
                Fun[[Fun[[T0, T1], T2], T0], Fun[[T1], T2]],
                Fun[[Fun[[T0, T1], T2], T0, T1], Fun[[], T2]],
                # 3 args
                Fun[[Fun[[T0, T1, T2], T3]], Fun[[T0, T1, T2], T3]],
                Fun[[Fun[[T0, T1, T2], T3], T0], Fun[[T1, T2], T3]],
                Fun[[Fun[[T0, T1, T2], T3], T0, T1], Fun[[T2], T3]],
                Fun[[Fun[[T0, T1, T2], T3], T0, T1, T2], Fun[[], T3]],
                # 4 args
                Fun[[Fun[[T0, T1, T2, T3], T4]], Fun[[T0, T1, T2, T3], T4]],
                Fun[[Fun[[T0, T1, T2, T3], T4], T0], Fun[[T1, T2, T3], T4]],
                Fun[[Fun[[T0, T1, T2, T3], T4], T0, T1], Fun[[T2, T3], T4]],
                Fun[[Fun[[T0, T1, T2, T3], T4], T0, T1, T2], Fun[[T3], T4]],
                Fun[[Fun[[T0, T1, T2, T3], T4], T0, T1, T2, T3], Fun[[], T4]],
                # 5 args
                Fun[[Fun[[T0, T1, T2, T3, T4], T5]],
                    Fun[[T0, T1, T2, T3, T4], T5]],
                Fun[[Fun[[T0, T1, T2, T3, T4], T5], T0],
                    Fun[[T1, T2, T3, T4], T5]],
                Fun[[Fun[[T0, T1, T2, T3, T4], T5], T0, T1],
                    Fun[[T2, T3, T4], T5]],
                Fun[[Fun[[T0, T1, T2, T3, T4], T5], T0, T1, T2],
                    Fun[[T3, T4], T5]],
                Fun[[Fun[[T0, T1, T2, T3, T4], T5], T0, T1, T2, T3],
                    Fun[[T4], T5]],
                Fun[[Fun[[T0, T1, T2, T3, T4], T5], T0, T1, T2, T3, T4],
                    Fun[[], T5]],
                # 6 args
                Fun[[Fun[[T0, T1, T2, T3, T4, T5], T6]],
                    Fun[[T0, T1, T2, T3, T4, T5], T6]],
                Fun[[Fun[[T0, T1, T2, T3, T4, T5], T6],
                     T0], Fun[[T1, T2, T3, T4, T5], T6]],
                Fun[[Fun[[T0, T1, T2, T3, T4, T5], T6],
                     T0, T1], Fun[[T2, T3, T4, T5], T6]],
                Fun[[Fun[[T0, T1, T2, T3, T4, T5], T6],
                     T0, T1, T2], Fun[[T3, T4, T5], T6]],
                Fun[[Fun[[T0, T1, T2, T3, T4, T5], T6],
                     T0, T1, T2, T3], Fun[[T4, T5], T6]],
                Fun[[Fun[[T0, T1, T2, T3, T4, T5], T6],
                     T0, T1, T2, T3, T4], Fun[[T5], T6]],
                Fun[[Fun[[T0, T1, T2, T3, T4, T5], T6],
                     T0, T1, T2, T3, T4, T5], Fun[[], T6]],
            ]
        ),
        "reduce": ReadOnceFunctionIntr(signature=_functools_reduce_signature),
    },
    "bisect": {
        "bisect_left": ConstFunctionIntr(
            signature=Union[
                Fun[[List[T0], T0], int],
                Fun[[List[T0], T0, int], int],
                Fun[[List[T0], T0, int, int], int],
            ],
            return_range=interval.positive_values
        ),
        "bisect_right": ConstFunctionIntr(
            signature=Union[
                Fun[[List[T0], T0], int],
                Fun[[List[T0], T0, int], int],
                Fun[[List[T0], T0, int, int], int],
            ],
            return_range=interval.positive_values
        ),
        "bisect": ConstFunctionIntr(
            signature=Union[
                Fun[[List[T0], T0], int],
                Fun[[List[T0], T0, int], int],
                Fun[[List[T0], T0, int, int], int],
            ],
            return_range=interval.positive_values
        ),
    },
    "cmath": {
        "cos": FunctionIntr(
            signature=Union[
                Fun[[float], complex],
                Fun[[complex], complex],
            ],
        ),
        "sin": FunctionIntr(
            signature=Union[
                Fun[[float], complex],
                Fun[[complex], complex],
            ],
        ),
        "exp": FunctionIntr(
            signature=Union[
                Fun[[float], complex],
                Fun[[complex], complex],
            ],
        ),
        "sqrt": FunctionIntr(
            signature=Union[
                Fun[[float], complex],
                Fun[[complex], complex],
            ],
        ),
        "log10": FunctionIntr(
            signature=Union[
                Fun[[float], complex],
                Fun[[complex], complex],
            ],
        ),
        "isnan": FunctionIntr(
            signature=Union[
                Fun[[float], bool],
                Fun[[complex], bool],
            ],
        ),
        "pi": ConstantIntr(signature=float),
        "e": ConstantIntr(signature=float),
    },
    'io': {
        '_io': {
            "TextIOWrapper": ClassWithConstConstructor(
                CLASSES['file'], global_effects=True)
        }
    },
    "itertools": {
        "count": ReadOnceFunctionIntr(
            signature=Union[
                Fun[[], Generator[int]],
                Fun[[int], Generator[int]],
                Fun[[int, int], Generator[int]],
            ]
        ),
        "islice": ReadOnceFunctionIntr(),
        "product": ConstFunctionIntr(
            signature=Union[
                Fun[[], Generator[T0]],
                Fun[[Iterable[T0]], Generator[Tuple[T0]]],
                Fun[[Iterable[T0], Iterable[T1]], Generator[Tuple[T0, T1]]],
                Fun[[Iterable[T0], Iterable[T1], Iterable[T2]],
                    Generator[Tuple[T0, T1, T2]]],
                Fun[[Iterable[T0], Iterable[T1], Iterable[T2], Iterable[T3]],
                    Generator[Tuple[T0, T1, T2, T3]]],
            ],
        ),
        "combinations": ConstFunctionIntr(
            signature=Fun[[Iterable[T0], int], Generator[List[T0]]]),
        "permutations": ConstFunctionIntr(
            signature=Union[
                Fun[[Iterable[T0]], Generator[List[T0]]],
                Fun[[Iterable[T0], int], Generator[List[T0]]],
            ],
            immediate_arguments=[1],
        ),
        "repeat": ConstFunctionIntr(
            signature=Union[
                Fun[[T0], Iterable[T0]],
                Fun[[T0, int], Iterable[T0]],
            ],
        ),
    },
    "random": {
        "seed": FunctionIntr(
            signature=Union[
                Fun[[], None],
                Fun[[T0], None],
            ],
            global_effects=True
        ),
        "random": FunctionIntr(
            signature=Fun[[], float],
            global_effects=True
        ),
        "randint": FunctionIntr(
            signature=Fun[[int, int], int],
            global_effects=True
        ),
        "randrange": FunctionIntr(
            signature=Union[
                Fun[[int], int],
                Fun[[int, int], int],
                Fun[[int, int, int], int]
            ],
            global_effects=True
        ),
        "gauss": FunctionIntr(
            signature=Fun[[float, float], float],
            global_effects=True
        ),
        "uniform": FunctionIntr(
            signature=Fun[[float, float], float],
            global_effects=True
        ),
        "expovariate": FunctionIntr(
            signature=Fun[[float], float],
            global_effects=True
        ),
        "sample": FunctionIntr(
            signature=Fun[[Iterable[T0], int], List[T0]],
            global_effects=True
        ),
        "choice": FunctionIntr(
            signature=Fun[[Iterable[T0]], T0],
            global_effects=True
        ),
        "shuffle": FunctionIntr(
            signature=Union[
                Fun[[List[T0]], None],
                Fun[[List[T0], Fun[[], float]], None],
            ],
            global_effects=True
        ),
    },
    "omp": {
        "set_num_threads": FunctionIntr(global_effects=True),
        "get_num_threads": FunctionIntr(global_effects=True),
        "get_max_threads": FunctionIntr(global_effects=True),
        "get_thread_num": FunctionIntr(global_effects=True),
        "get_num_procs": FunctionIntr(global_effects=True),
        "in_parallel": FunctionIntr(global_effects=True),
        "set_dynamic": FunctionIntr(global_effects=True),
        "get_dynamic": FunctionIntr(global_effects=True),
        "set_nested": FunctionIntr(global_effects=True),
        "get_nested": FunctionIntr(global_effects=True),
        "init_lock": FunctionIntr(global_effects=True),
        "destroy_lock": FunctionIntr(global_effects=True),
        "set_lock": FunctionIntr(global_effects=True),
        "unset_lock": FunctionIntr(global_effects=True),
        "test_lock": FunctionIntr(global_effects=True),
        "init_nest_lock": FunctionIntr(global_effects=True),
        "destroy_nest_lock": FunctionIntr(global_effects=True),
        "set_nest_lock": FunctionIntr(global_effects=True),
        "unset_nest_lock": FunctionIntr(global_effects=True),
        "test_nest_lock": FunctionIntr(global_effects=True),
        "get_wtime": FunctionIntr(global_effects=True),
        "get_wtick": FunctionIntr(global_effects=True),
    },
    "operator": {
        "lt": ConstFunctionIntr(signature=_operator_eq_signature),
        "le": ConstFunctionIntr(signature=_operator_eq_signature),
        "eq": ConstFunctionIntr(signature=_operator_eq_signature),
        "ne": ConstFunctionIntr(signature=_operator_eq_signature),
        "ge": ConstFunctionIntr(signature=_operator_eq_signature),
        "gt": ConstFunctionIntr(signature=_operator_eq_signature),
        "__lt__": ConstFunctionIntr(signature=_operator_eq_signature),
        "__le__": ConstFunctionIntr(signature=_operator_eq_signature),
        "__eq__": ConstFunctionIntr(signature=_operator_eq_signature),
        "__ne__": ConstFunctionIntr(signature=_operator_eq_signature),
        "__ge__": ConstFunctionIntr(signature=_operator_eq_signature),
        "__gt__": ConstFunctionIntr(signature=_operator_eq_signature),
        "not_": ConstFunctionIntr(),
        "__not__": ConstFunctionIntr(),
        "truth": ConstFunctionIntr(),
        "is_": ConstFunctionIntr(),
        "is_not": ConstFunctionIntr(),
        "abs": ConstFunctionIntr(),
        "__abs__": ConstFunctionIntr(),
        "add": ConstFunctionIntr(signature=_numpy_binary_op_signature),
        "__add__": ConstFunctionIntr(signature=_numpy_binary_op_signature),
        "and_": ConstFunctionIntr(),
        "__and__": ConstFunctionIntr(),
        "floordiv": ConstFunctionIntr(signature=_numpy_binary_op_signature),
        "__floordiv__": ConstFunctionIntr(
            signature=_numpy_binary_op_signature
        ),
        "inv": ConstFunctionIntr(),
        "invert": ConstFunctionIntr(),
        "__inv__": ConstFunctionIntr(),
        "__invert__": ConstFunctionIntr(),
        "lshift": ConstFunctionIntr(signature=_numpy_int_binary_op_signature),
        "__lshift__": ConstFunctionIntr(
            signature=_numpy_int_binary_op_signature
        ),
        "matmul": ConstFunctionIntr(signature=_operator_mul_signature,
                                    requires_blas=True),
        "__matmul__": ConstFunctionIntr(signature=_operator_mul_signature,
                                        requires_blas=True),
        "imatmul": MethodIntr(update_effects, requires_blas=True),
        "__imatmul__": MethodIntr(update_effects, requires_blas=True),
        "mod": ConstFunctionIntr(signature=_operator_mod_signature),
        "__mod__": ConstFunctionIntr(signature=_operator_mod_signature),
        "mul": ConstFunctionIntr(signature=_operator_mul_signature),
        "__mul__": ConstFunctionIntr(signature=_operator_mul_signature),
        "neg": ConstFunctionIntr(),
        "__neg__": ConstFunctionIntr(),
        "or_": ConstFunctionIntr(),
        "__or__": ConstFunctionIntr(),
        "pos": ConstFunctionIntr(signature=_numpy_unary_op_signature),
        "__pos__": ConstFunctionIntr(signature=_numpy_unary_op_signature),
        "rshift": ConstFunctionIntr(signature=_numpy_int_binary_op_signature),
        "__rshift__": ConstFunctionIntr(
            signature=_numpy_int_binary_op_signature
        ),
        "sub": ConstFunctionIntr(signature=_operator_sub_signature),
        "__sub__": ConstFunctionIntr(signature=_operator_sub_signature),
        "truediv": ConstFunctionIntr(),
        "__truediv__": ConstFunctionIntr(),
        "xor": ConstFunctionIntr(),
        "__xor__": ConstFunctionIntr(),
        "concat": ConstFunctionIntr(),
        "__concat__": ConstFunctionIntr(),
        "iadd": MethodIntr(update_effects),
        "__iadd__": MethodIntr(update_effects),
        "iand": MethodIntr(update_effects),
        "__iand__": MethodIntr(update_effects),
        "iconcat": MethodIntr(update_effects),
        "__iconcat__": MethodIntr(update_effects),
        "ifloordiv": MethodIntr(update_effects),
        "__ifloordiv__": MethodIntr(update_effects),
        "ilshift": MethodIntr(update_effects),
        "__ilshift__": MethodIntr(update_effects),
        "imod": MethodIntr(update_effects),
        "__imod__": MethodIntr(update_effects),
        "imul": MethodIntr(update_effects),
        "__imul__": MethodIntr(update_effects),
        "ior": MethodIntr(update_effects),
        "__ior__": MethodIntr(update_effects),
        "ipow": MethodIntr(update_effects),
        "__ipow__": MethodIntr(update_effects),
        "irshift": MethodIntr(update_effects),
        "__irshift__": MethodIntr(update_effects),
        "isub": MethodIntr(update_effects),
        "__isub__": MethodIntr(update_effects),
        "itruediv": MethodIntr(update_effects),
        "__itruediv__": MethodIntr(update_effects),
        "ixor": MethodIntr(update_effects),
        "__ixor__": MethodIntr(update_effects),
        "contains": MethodIntr(
            update_effects,
            signature=_operator_contains_signature
        ),
        "__contains__": MethodIntr(
            update_effects,
            signature=_operator_contains_signature
        ),
        "countOf": ConstFunctionIntr(),
        "delitem": FunctionIntr(
            argument_effects=[UpdateEffect(), ReadEffect()]),
        "__delitem__": FunctionIntr(
            argument_effects=[UpdateEffect(), ReadEffect()]),
        "getitem": ConstFunctionIntr(signature=_operator_getitem_signature),
        "__getitem__": ConstFunctionIntr(
            signature=_operator_getitem_signature
        ),
        "indexOf": ConstFunctionIntr(),
        "__theitemgetter__": ConstFunctionIntr(),
        "itemgetter": MethodIntr(
            return_alias=lambda _: {
                MODULES['operator']['__theitemgetter__']}
        ),

    },
    "string": {
        "ascii_lowercase": ConstantIntr(signature=str),
        "ascii_uppercase": ConstantIntr(signature=str),
        "ascii_letters": ConstantIntr(signature=str),
        "digits": ConstantIntr(signature=str),
        "hexdigits": ConstantIntr(signature=str),
        "octdigits": ConstantIntr(signature=str),
    },
    "os": {
        "path": {
            "join": ConstFunctionIntr(
                signature=Union[
                    Fun[[str], str],
                    Fun[[str, str], str],
                    Fun[[str, str, str], str],
                    Fun[[str, str, str, str], str],
                    Fun[[str, str, str, str, str], str],
                ]
            ),
        }
    },
    # conflicting method names must be listed here
    "__dispatch__": {
        "append": MethodIntr(signature=Fun[[List[T0], T0], None]),
        "clear": MethodIntr(signature=Fun[[T0], None]),
        "conjugate": ConstMethodIntr(),
        "copy": ConstMethodIntr(signature=Fun[[T0], T0]),
        "count": ConstMethodIntr(
            signature=Union[
                Fun[[Iterable[T0], T0], int],
                Fun[[Iterable[T0], T0, int], int],
                Fun[[Iterable[T0], T0, int, int], int],
            ],
            return_range=interval.positive_values
        ),
        "extend": MethodIntr(signature=Fun[[List[T0], Iterable[T0]], None]),
        "index": ConstMethodIntr(
            signature=Union[
                Fun[[Iterable[T0], T0], int],
                Fun[[Iterable[T0], T0, int], int],
                Fun[[Iterable[T0], T0, int, int], int],
            ],
            return_range=interval.positive_values
        ),
        "insert": MethodIntr(signature=Fun[[List[T0], int, T0], None]),
        "pop": MethodIntr(),
        "remove": MethodIntr(),
        "reverse": MethodIntr(),
        "sort": MethodIntr(),
        "tolist": ConstMethodIntr(),
        "update": MethodIntr(update_effects),
    },
}

# PyPy doesn't seem to provide this.
if sys.implementation.name == 'pypy':
    del MODULES['array']['typecodes']

if sys.version_info < (3, 5):
    del MODULES['operator']['matmul']
    del MODULES['operator']['__matmul__']

# VMSError is only available on VMS
if 'VMSError' in sys.modules['builtins'].__dict__:
    MODULES['builtins']['VMSError'] = ConstExceptionIntr()

# WindowsError is only available on Windows
if 'WindowsError' in sys.modules['builtins'].__dict__:
    MODULES['builtins']['WindowsError'] = ConstExceptionIntr()

# detect and prune unsupported modules
for module_name in ["omp", "scipy", "scipy.special"]:
    try:
        import_module(module_name)
    except:
        logger.info(
            "Pythran support for package '{}' will be reduced: "
            "this module is not available at runtime.".format(module_name)
        )

# check and delete unimplemented numpy methods
for method in list(MODULES['numpy'].keys()):
    if not hasattr(numpy, method):
        del MODULES['numpy'][method]

# if openmp is available, check its version and populate the API accordingly
try:
    omp_version = getattr(__import__('omp'), 'VERSION', 45)
    if omp_version >= 30:
        MODULES['omp'].update({
            "set_schedule": FunctionIntr(global_effects=True),
            "get_schedule": FunctionIntr(global_effects=True),
            "get_thread_limit": FunctionIntr(global_effects=True),
            "set_max_active_levels": FunctionIntr(global_effects=True),
            "get_max_active_levels": FunctionIntr(global_effects=True),
            "get_level": FunctionIntr(global_effects=True),
            "get_ancestor_thread_num": FunctionIntr(global_effects=True),
            "get_team_size": FunctionIntr(global_effects=True),
            "get_active_level": FunctionIntr(global_effects=True),
            "in_final": FunctionIntr(global_effects=True),
        })
except ImportError:
    pass

def save_path(module_name, elements):
    """ Recursively sets path. """
    for elem, signature in elements.items():
        if isinstance(signature, dict):  # Submodule case
            save_path(module_name + (elem,), signature)
        else:
            signature.path = module_name + (elem,)

save_path((), CLASSES)
save_path((), MODULES)

def looks_like_a_forward_function(spec):
    return not spec.args and spec.varargs == 'args' and spec.varkw == 'kwargs'

# populate argument description through introspection
def save_arguments(module_name, elements):
    """ Recursively save arguments name and default value. """
    for elem, signature in elements.items():
        if isinstance(signature, dict):  # Submodule case
            save_arguments(module_name + (elem,), signature)
        else:
            # use introspection to get the Python obj
            try:
                themodule = import_module(".".join(module_name))
                obj = getattr(themodule, elem)
                while hasattr(obj, '__wrapped__'):
                    obj = obj.__wrapped__
            except (AttributeError, ImportError, TypeError, ValueError):
                continue

            # first try to gather info through getfullargspec
            try:
                spec = inspect.getfullargspec(obj)
            except:
                continue

            # some function are actually forward function, detect those
            # and accept to use our description instead.
            if looks_like_a_forward_function(spec):
                assert signature.args.args, "{} require an explicit description".format(elem)
                continue

            args = [ast.Name(arg, ast.Param(), None, None)
                    for arg in spec.args]

            # pop 'self' if we have a bound method
            if inspect.ismethod(obj):
                args = args[1:]

            defaults = list(spec.defaults or [])
            args += [ast.Name(arg, ast.Param(), None, None)
                     for arg in spec.kwonlyargs]
            if spec.kwonlydefaults:
                defaults += [spec.kwonlydefaults[kw] for kw in
                             spec.kwonlyargs[-len(spec.kwonlydefaults):]]

            # Check if we already have a pythran description for that object
            if signature.args.args:
                if module_name != ('numpy', 'random'):
                    # Skip this warning for `numpy.random`, because the
                    # signatures changed in 1.25.x and that may yet be reverted
                    # (see pythran#2139) in 1.25.3 and/or 1.26.x
                    logger.warning(
                        "Overriding pythran description with argspec "
                        "information for: {}".format(

                        ".".join(module_name + (elem,))))
                else:
                    continue

            # Avoid use of comprehension to fill "as much args/defaults" as
            # possible
            signature_args = args[:-len(defaults) or None]
            signature_defaults = []
            try:
                for arg, value in zip(args[-len(defaults):], defaults):
                    signature_args.append(arg)
                    signature_defaults.append(to_ast(value))
            except ToNotEval:
                continue

            # Only validate once we have a valid signature
            signature.args.args = signature_args
            signature.args.defaults = signature_defaults


save_arguments((), MODULES)


# Fill return_type field for constants
def fill_constants_types(module_name, elements):
    """ Recursively save arguments name and default value. """
    for elem, intrinsic in elements.items():
        if isinstance(intrinsic, dict):  # Submodule case
            fill_constants_types(module_name + (elem,), intrinsic)
        elif isinstance(intrinsic, ConstantIntr):
            # use introspection to get the Python constants types
            cst = getattr(import_module(".".join(module_name)), elem)
            intrinsic.signature = type(cst)


fill_constants_types((), MODULES)


# a method name to module binding
# {method_name : ((full module path), signature)}
methods = {}
duplicated_methods = {}


def save_method(elements, module_path):
    """ Recursively save methods with module name and signature. """
    for elem, signature in elements.items():
        if isinstance(signature, dict):  # Submodule case
            save_method(signature, module_path + (elem,))
        elif isinstance(signature, Class):
            save_method(signature.fields, module_path + (elem,))
        elif signature.ismethod():
            # in case of duplicates, there must be a __dispatch__ record
            # and it is the only recorded one
            if elem in MODULES['__dispatch__'] and module_path[0] != '__dispatch__':
                duplicated_methods.setdefault(elem, []).append((module_path,
                                                                signature))

            if elem in methods and module_path[0] != '__dispatch__':
                assert elem in MODULES['__dispatch__']
                path = ('__dispatch__',)
                methods[elem] = (path, MODULES['__dispatch__'][elem])
            else:
                methods[elem] = (module_path, signature)


save_method(MODULES, ())

# a function name to module binding
# {function_name : [((full module path), signature)]}
functions = {}


def save_function(elements, module_path):
    """ Recursively save functions with module name and signature. """
    for elem, signature in elements.items():
        if isinstance(signature, dict):  # Submodule case
            save_function(signature, module_path + (elem,))
        elif signature.isstaticfunction():
            functions.setdefault(elem, []).append((module_path, signature,))
        elif isinstance(signature, Class):
            save_function(signature.fields, module_path + (elem,))


save_function(MODULES, ())


# a attribute name to module binding
# {attribute_name : ((full module path), signature)}
attributes = {}


def save_attribute(elements, module_path):
    """ Recursively save attributes with module name and signature. """
    for elem, signature in elements.items():
        if isinstance(signature, dict):  # Submodule case
            save_attribute(signature, module_path + (elem,))
        elif signature.isattribute():
            assert elem not in attributes  # we need unicity
            attributes[elem] = (module_path, signature,)
        elif isinstance(signature, Class):
            save_attribute(signature.fields, module_path + (elem,))


save_attribute(MODULES, ())

blas_requires = set()

def save_blas_requires(elements, module_path):
    """ Recursively save attributes with module name and signature. """
    for elem, signature in elements.items():
        if isinstance(signature, dict):  # Submodule case
            save_blas_requires(signature, module_path + (elem,))
        elif signature.requires_blas:
            blas_requires.add(module_path + (elem,))
        elif isinstance(signature, Class):
            save_blas_requires(signature.fields, module_path + (elem,))

save_blas_requires(MODULES, ())


# patch beniget with pythran-specific builtins
import beniget
beniget.beniget.Builtins['builtins'] = __import__('builtins')
beniget.beniget.Builtins['__dispatch__'] = object()
for k, v in MODULES['builtins'].items():
    if k not in beniget.beniget.Builtins:
        beniget.beniget.Builtins[k] = v
