from mypy.util import unnamed_function
from typing import Dict, Any, Optional, Tuple
import sys

from typing_extensions import Final

PREFIX: Final = "CPyPy_"  # Python wrappers
NATIVE_PREFIX: Final = "CPyDef_"  # Native functions etc.
DUNDER_PREFIX: Final = "CPyDunder_"  # Wrappers for exposing dunder methods to the API
REG_PREFIX: Final = "cpy_r_"  # Registers
STATIC_PREFIX: Final = "CPyStatic_"  # Static variables (for literals etc.)
TYPE_PREFIX: Final = "CPyType_"  # Type object struct
MODULE_PREFIX: Final = "CPyModule_"  # Cached modules
ATTR_PREFIX: Final = "_"  # Attributes

ENV_ATTR_NAME: Final = "__mypyc_env__"
NEXT_LABEL_ATTR_NAME: Final = "__mypyc_next_label__"
TEMP_ATTR_NAME: Final = "__mypyc_temp__"
LAMBDA_NAME: Final = "__mypyc_lambda__"
PROPSET_PREFIX: Final = "__mypyc_setter__"
SELF_NAME: Final = "__mypyc_self__"

# Max short int we accept as a literal is based on 32-bit platforms,
# so that we can just always emit the same code.

TOP_LEVEL_NAME: Final = "__top_level__"  # Special function representing module top level

# Maximal number of subclasses for a class to trigger fast path in isinstance() checks.
FAST_ISINSTANCE_MAX_SUBCLASSES: Final = 2

IS_32_BIT_PLATFORM: Final = sys.maxsize < (1 << 31)

PLATFORM_SIZE = 4 if IS_32_BIT_PLATFORM else 8

# Python 3.5 on macOS uses a hybrid 32/64-bit build that requires some workarounds.
# The same generated C will be compiled in both 32 and 64 bit modes when building mypy
# wheels (for an unknown reason).
#
# Note that we use "in ['darwin']" because of https://github.com/mypyc/mypyc/issues/761.
IS_MIXED_32_64_BIT_BUILD: Final = sys.platform in ["darwin"] and sys.version_info < (3, 6)

# Maximum value for a short tagged integer.
MAX_SHORT_INT: Final = sys.maxsize >> 1

# Maximum value for a short tagged integer represented as a C integer literal.
#
# Note: Assume that the compiled code uses the same bit width as mypyc, except for
#       Python 3.5 on macOS.
MAX_LITERAL_SHORT_INT: Final = sys.maxsize >> 1 if not IS_MIXED_32_64_BIT_BUILD else 2 ** 30 - 1
MIN_LITERAL_SHORT_INT: Final = -MAX_LITERAL_SHORT_INT - 1

# Runtime C library files
RUNTIME_C_FILES: Final = [
    'init.c',
    'getargs.c',
    'getargsfast.c',
    'int_ops.c',
    'str_ops.c',
    'bytes_ops.c',
    'list_ops.c',
    'dict_ops.c',
    'set_ops.c',
    'tuple_ops.c',
    'exc_ops.c',
    'misc_ops.c',
    'generic_ops.c',
]


JsonDict = Dict[str, Any]


def shared_lib_name(group_name: str) -> str:
    """Given a group name, return the actual name of its extension module.

    (This just adds a suffix to the final component.)
    """
    return '{}__mypyc'.format(group_name)


def short_name(name: str) -> str:
    if name.startswith('builtins.'):
        return name[9:]
    return name


def use_fastcall(capi_version: Tuple[int, int]) -> bool:
    # We can use METH_FASTCALL for faster wrapper functions on Python 3.7+.
    return capi_version >= (3, 7)


def use_vectorcall(capi_version: Tuple[int, int]) -> bool:
    # We can use vectorcalls to make calls on Python 3.8+ (PEP 590).
    return capi_version >= (3, 8)


def use_method_vectorcall(capi_version: Tuple[int, int]) -> bool:
    # We can use a dedicated vectorcall API to call methods on Python 3.9+.
    return capi_version >= (3, 9)


def get_id_from_name(name: str, fullname: str, line: int) -> str:
    """Create a unique id for a function.

    This creates an id that is unique for any given function definition, so that it can be used as
    a dictionary key. This is usually the fullname of the function, but this is different in that
    it handles the case where the function is named '_', in which case multiple different functions
    could have the same name."""
    if unnamed_function(name):
        return "{}.{}".format(fullname, line)
    else:
        return fullname


def short_id_from_name(func_name: str, shortname: str, line: Optional[int]) -> str:
    if unnamed_function(func_name):
        assert line is not None
        partial_name = "{}.{}".format(shortname, line)
    else:
        partial_name = shortname
    return partial_name
