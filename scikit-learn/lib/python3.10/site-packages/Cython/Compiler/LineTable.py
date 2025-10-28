"""
Build a line table for CodeObjects, according to PEP-626 / Python 3.11.

See  https://github.com/python/cpython/blob/1054a755a3016f95fcd24b3ad20e8ed9048b7939/InternalDocs/locations.md
See  https://github.com/python/cpython/blob/1054a755a3016f95fcd24b3ad20e8ed9048b7939/Python/assemble.c#L192
"""

import cython


def build_line_table(positions, firstlineno):
    # positions is a list of four-tuples (start_lineno, end_lineno, start_col_offset, end_col_offset)
    table_bytes = []
    last_lineno = firstlineno
    for position_info in positions:
        last_lineno = encode_single_position(table_bytes, position_info, last_lineno)
    linetable = ''.join(table_bytes)

    """
    # Hacky debug helper code for the line table generation.
    code_obj = build_line_table.__code__.replace(co_linetable=linetable.encode('latin1'), co_firstlineno=firstlineno)
    print()
    print(repr(linetable))
    print(positions)
    print(list(code_obj.co_positions()))
    """

    return linetable


@cython.cfunc
def encode_single_position(table_bytes: list, position_info: tuple, last_lineno: cython.int) -> cython.int:
    start_lineno: cython.int
    end_lineno: cython.int
    start_column: cython.int
    end_column: cython.int

    start_lineno, end_lineno, start_column, end_column = position_info
    assert start_lineno >= last_lineno, f"{start_lineno} >= {last_lineno}"  # positions should be sorted

    last_lineno_delta: cython.int = start_lineno - last_lineno

    if end_lineno == start_lineno:
        # All in one line, can try short forms.
        if last_lineno_delta == 0 and start_column < 80 and 0 <= (end_column - start_column) < 16:
            # Short format (code 0-9): still on same line, small column offset
            encode_location_short(table_bytes, start_column, end_column)
            return end_lineno
        elif 0 <= last_lineno_delta < 3 and start_column < 128 and end_column < 128:
            # One line format (code 10-12): small line offsets / larger column offsets
            encode_location_oneline(table_bytes, last_lineno_delta, start_column, end_column)
            return end_lineno

    # Store in long format (code 14)
    encode_location_start(table_bytes, 14)
    # Since we sort positions, negative line deltas should never occur ==> inline encode_varint_signed()
    encode_varint(table_bytes, last_lineno_delta << 1)
    encode_varint(table_bytes, end_lineno - start_lineno)
    encode_varint(table_bytes, start_column + 1)
    encode_varint(table_bytes, end_column + 1)
    return end_lineno


@cython.exceptval(-1, check=False)
@cython.cfunc
def encode_location_start(table_bytes: list, code: cython.int) -> cython.int:
    # "Instruction" size is always 1
    # 128 | (code << 3) | (length - 1)
    table_bytes.append(chr(128 | (code << 3)))
    return 0


@cython.exceptval(-1, check=False)
@cython.cfunc
def encode_location_short(table_bytes: list, start_column: cython.int, end_column: cython.int) -> cython.int:
    low_bits: cython.int = start_column & 7
    code: cython.int = start_column >> 3
    # inlined encode_location_start()
    table_bytes.append(f"{128 | (code << 3):c}{(low_bits << 4) | (end_column - start_column):c}")
    return 0


@cython.exceptval(-1, check=False)
@cython.cfunc
def encode_location_oneline(table_bytes: list, line_delta: cython.int, start_column: cython.int, end_column: cython.int) -> cython.int:
    code: cython.int = 10 + line_delta
    # inlined encode_location_start()
    table_bytes.append(f"{128 | (code << 3):c}{start_column:c}{end_column:c}")
    return 0


"""
# Since we sort positions, negative line deltas should not occur.
@cython.cfunc
def encode_varint_signed(table_bytes: list, value: cython.int) -> cython.int:
    # (unsigned int)(-val) has undefined behavior for INT_MIN
    uval: cython.uint = cython.cast(cython.uint, value) if cython.compiled else value
    if value < 0:
        uval = ((0 - uval) << 1) | 1
    else:
        uval = uval << 1
    encode_varint(table_bytes, uval)
"""


@cython.exceptval(-1, check=False)
@cython.cfunc
def encode_varint(table_bytes: list, value: cython.uint) -> cython.int:
    assert value > 0 or value == 0
    while value >= 64:
        table_bytes.append(chr(64 | (value & 63)))
        value >>= 6
    table_bytes.append(chr(value))
    return 0
