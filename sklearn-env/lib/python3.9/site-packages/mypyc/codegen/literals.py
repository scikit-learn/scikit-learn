from typing import Dict, List, Union, Tuple, Any, cast

from typing_extensions import Final


# Supported Python literal types. All tuple items must have supported
# literal types as well, but we can't represent the type precisely.
LiteralValue = Union[str, bytes, int, bool, float, complex, Tuple[object, ...], None]


# Some literals are singletons and handled specially (None, False and True)
NUM_SINGLETONS: Final = 3


class Literals:
    """Collection of literal values used in a compilation group and related helpers."""

    def __init__(self) -> None:
        # Each dict maps value to literal index (0, 1, ...)
        self.str_literals: Dict[str, int] = {}
        self.bytes_literals: Dict[bytes, int] = {}
        self.int_literals: Dict[int, int] = {}
        self.float_literals: Dict[float, int] = {}
        self.complex_literals: Dict[complex, int] = {}
        self.tuple_literals: Dict[Tuple[object, ...], int] = {}

    def record_literal(self, value: LiteralValue) -> None:
        """Ensure that the literal value is available in generated code."""
        if value is None or value is True or value is False:
            # These are special cased and always present
            return
        if isinstance(value, str):
            str_literals = self.str_literals
            if value not in str_literals:
                str_literals[value] = len(str_literals)
        elif isinstance(value, bytes):
            bytes_literals = self.bytes_literals
            if value not in bytes_literals:
                bytes_literals[value] = len(bytes_literals)
        elif isinstance(value, int):
            int_literals = self.int_literals
            if value not in int_literals:
                int_literals[value] = len(int_literals)
        elif isinstance(value, float):
            float_literals = self.float_literals
            if value not in float_literals:
                float_literals[value] = len(float_literals)
        elif isinstance(value, complex):
            complex_literals = self.complex_literals
            if value not in complex_literals:
                complex_literals[value] = len(complex_literals)
        elif isinstance(value, tuple):
            tuple_literals = self.tuple_literals
            if value not in tuple_literals:
                for item in value:
                    self.record_literal(cast(Any, item))
                tuple_literals[value] = len(tuple_literals)
        else:
            assert False, 'invalid literal: %r' % value

    def literal_index(self, value: LiteralValue) -> int:
        """Return the index to the literals array for given value."""
        # The array contains first None and booleans, followed by all str values,
        # followed by bytes values, etc.
        if value is None:
            return 0
        elif value is False:
            return 1
        elif value is True:
            return 2
        n = NUM_SINGLETONS
        if isinstance(value, str):
            return n + self.str_literals[value]
        n += len(self.str_literals)
        if isinstance(value, bytes):
            return n + self.bytes_literals[value]
        n += len(self.bytes_literals)
        if isinstance(value, int):
            return n + self.int_literals[value]
        n += len(self.int_literals)
        if isinstance(value, float):
            return n + self.float_literals[value]
        n += len(self.float_literals)
        if isinstance(value, complex):
            return n + self.complex_literals[value]
        n += len(self.complex_literals)
        if isinstance(value, tuple):
            return n + self.tuple_literals[value]
        assert False, 'invalid literal: %r' % value

    def num_literals(self) -> int:
        # The first three are for None, True and False
        return (NUM_SINGLETONS + len(self.str_literals) + len(self.bytes_literals) +
                len(self.int_literals) + len(self.float_literals) + len(self.complex_literals) +
                len(self.tuple_literals))

    # The following methods return the C encodings of literal values
    # of different types

    def encoded_str_values(self) -> List[bytes]:
        return _encode_str_values(self.str_literals)

    def encoded_int_values(self) -> List[bytes]:
        return _encode_int_values(self.int_literals)

    def encoded_bytes_values(self) -> List[bytes]:
        return _encode_bytes_values(self.bytes_literals)

    def encoded_float_values(self) -> List[str]:
        return _encode_float_values(self.float_literals)

    def encoded_complex_values(self) -> List[str]:
        return _encode_complex_values(self.complex_literals)

    def encoded_tuple_values(self) -> List[str]:
        """Encode tuple values into a C array.

        The format of the result is like this:

           <number of tuples>
           <length of the first tuple>
           <literal index of first item>
           ...
           <literal index of last item>
           <length of the second tuple>
           ...
        """
        values = self.tuple_literals
        value_by_index = {}
        for value, index in values.items():
            value_by_index[index] = value
        result = []
        num = len(values)
        result.append(str(num))
        for i in range(num):
            value = value_by_index[i]
            result.append(str(len(value)))
            for item in value:
                index = self.literal_index(cast(Any, item))
                result.append(str(index))
        return result


def _encode_str_values(values: Dict[str, int]) -> List[bytes]:
    value_by_index = {}
    for value, index in values.items():
        value_by_index[index] = value
    result = []
    line: List[bytes] = []
    line_len = 0
    for i in range(len(values)):
        value = value_by_index[i]
        c_literal = format_str_literal(value)
        c_len = len(c_literal)
        if line_len > 0 and line_len + c_len > 70:
            result.append(format_int(len(line)) + b''.join(line))
            line = []
            line_len = 0
        line.append(c_literal)
        line_len += c_len
    if line:
        result.append(format_int(len(line)) + b''.join(line))
    result.append(b'')
    return result


def _encode_bytes_values(values: Dict[bytes, int]) -> List[bytes]:
    value_by_index = {}
    for value, index in values.items():
        value_by_index[index] = value
    result = []
    line: List[bytes] = []
    line_len = 0
    for i in range(len(values)):
        value = value_by_index[i]
        c_init = format_int(len(value))
        c_len = len(c_init) + len(value)
        if line_len > 0 and line_len + c_len > 70:
            result.append(format_int(len(line)) + b''.join(line))
            line = []
            line_len = 0
        line.append(c_init + value)
        line_len += c_len
    if line:
        result.append(format_int(len(line)) + b''.join(line))
    result.append(b'')
    return result


def format_int(n: int) -> bytes:
    """Format an integer using a variable-length binary encoding."""
    if n < 128:
        a = [n]
    else:
        a = []
        while n > 0:
            a.insert(0, n & 0x7f)
            n >>= 7
        for i in range(len(a) - 1):
            # If the highest bit is set, more 7-bit digits follow
            a[i] |= 0x80
    return bytes(a)


def format_str_literal(s: str) -> bytes:
    utf8 = s.encode('utf-8')
    return format_int(len(utf8)) + utf8


def _encode_int_values(values: Dict[int, int]) -> List[bytes]:
    """Encode int values into C strings.

    Values are stored in base 10 and separated by 0 bytes.
    """
    value_by_index = {}
    for value, index in values.items():
        value_by_index[index] = value
    result = []
    line: List[bytes] = []
    line_len = 0
    for i in range(len(values)):
        value = value_by_index[i]
        encoded = b'%d' % value
        if line_len > 0 and line_len + len(encoded) > 70:
            result.append(format_int(len(line)) + b'\0'.join(line))
            line = []
            line_len = 0
        line.append(encoded)
        line_len += len(encoded)
    if line:
        result.append(format_int(len(line)) + b'\0'.join(line))
    result.append(b'')
    return result


def float_to_c(x: float) -> str:
    """Return C literal representation of a float value."""
    s = str(x)
    if s == 'inf':
        return 'INFINITY'
    elif s == '-inf':
        return '-INFINITY'
    return s


def _encode_float_values(values: Dict[float, int]) -> List[str]:
    """Encode float values into a C array values.

    The result contains the number of values followed by individual values.
    """
    value_by_index = {}
    for value, index in values.items():
        value_by_index[index] = value
    result = []
    num = len(values)
    result.append(str(num))
    for i in range(num):
        value = value_by_index[i]
        result.append(float_to_c(value))
    return result


def _encode_complex_values(values: Dict[complex, int]) -> List[str]:
    """Encode float values into a C array values.

    The result contains the number of values followed by pairs of doubles
    representing complex numbers.
    """
    value_by_index = {}
    for value, index in values.items():
        value_by_index[index] = value
    result = []
    num = len(values)
    result.append(str(num))
    for i in range(num):
        value = value_by_index[i]
        result.append(float_to_c(value.real))
        result.append(float_to_c(value.imag))
    return result
