from __future__ import annotations

import re

# https://en.cppreference.com/w/c/keyword
_keywords = [
    'auto', 'break', 'case', 'char', 'const', 'continue', 'default', 'do', 'double',
    'else', 'enum', 'extern', 'float', 'for', 'goto', 'if', 'inline', 'int', 'long',
    'register', 'restrict', 'return', 'short', 'signed', 'sizeof', 'static', 'struct',
    'switch', 'typedef', 'union', 'unsigned', 'void', 'volatile', 'while',
    '_Alignas', '_Alignof', '_Atomic', '_Bool', '_Complex',
    '_Decimal32', '_Decimal64', '_Decimal128',
    '_Generic', '_Imaginary', '_Noreturn', '_Static_assert', '_Thread_local',
]
# These are only keyword'y when the corresponding headers are included.
# They are used as default value for c_extra_keywords.
_macroKeywords = [
    'alignas', 'alignof', 'bool', 'complex', 'imaginary', 'noreturn', 'static_assert',
    'thread_local',
]

# these are ordered by precedence
_expression_bin_ops = [
    ['||', 'or'],
    ['&&', 'and'],
    ['|', 'bitor'],
    ['^', 'xor'],
    ['&', 'bitand'],
    ['==', '!=', 'not_eq'],
    ['<=', '>=', '<', '>'],
    ['<<', '>>'],
    ['+', '-'],
    ['*', '/', '%'],
    ['.*', '->*'],
]
_expression_unary_ops = ["++", "--", "*", "&", "+", "-", "!", "not", "~", "compl"]
_expression_assignment_ops = ["=", "*=", "/=", "%=", "+=", "-=",
                              ">>=", "<<=", "&=", "and_eq", "^=", "xor_eq", "|=", "or_eq"]

_max_id = 1
_id_prefix = [None, 'c.', 'Cv2.']
# Ids are used in lookup keys which are used across pickled files,
# so when _max_id changes, make sure to update the ENV_VERSION.

_string_re = re.compile(r"[LuU8]?('([^'\\]*(?:\\.[^'\\]*)*)'"
                        r'|"([^"\\]*(?:\\.[^"\\]*)*)")', re.DOTALL)

# bool, complex, and imaginary are macro "keywords", so they are handled separately
_simple_type_specifiers_re = re.compile(r"""
    \b(
    void|_Bool
    |signed|unsigned
    |short|long
    |char
    |int
    |__uint128|__int128
    |__int(8|16|32|64|128)  # extension
    |float|double
    |_Decimal(32|64|128)
    |_Complex|_Imaginary
    |__float80|_Float64x|__float128|_Float128|__ibm128  # extension
    |__fp16  # extension
    |_Sat|_Fract|fract|_Accum|accum  # extension
    )\b
""", re.VERBOSE)
