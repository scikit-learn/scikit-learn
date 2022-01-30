"""
Formatting numeric literals.
"""
from blib2to3.pytree import Leaf


def format_hex(text: str) -> str:
    """
    Formats a hexadecimal string like "0x12B3"
    """
    before, after = text[:2], text[2:]
    return f"{before}{after.upper()}"


def format_scientific_notation(text: str) -> str:
    """Formats a numeric string utilizing scentific notation"""
    before, after = text.split("e")
    sign = ""
    if after.startswith("-"):
        after = after[1:]
        sign = "-"
    elif after.startswith("+"):
        after = after[1:]
    before = format_float_or_int_string(before)
    return f"{before}e{sign}{after}"


def format_long_or_complex_number(text: str) -> str:
    """Formats a long or complex string like `10L` or `10j`"""
    number = text[:-1]
    suffix = text[-1]
    # Capitalize in "2L" because "l" looks too similar to "1".
    if suffix == "l":
        suffix = "L"
    return f"{format_float_or_int_string(number)}{suffix}"


def format_float_or_int_string(text: str) -> str:
    """Formats a float string like "1.0"."""
    if "." not in text:
        return text

    before, after = text.split(".")
    return f"{before or 0}.{after or 0}"


def normalize_numeric_literal(leaf: Leaf) -> None:
    """Normalizes numeric (float, int, and complex) literals.

    All letters used in the representation are normalized to lowercase (except
    in Python 2 long literals).
    """
    text = leaf.value.lower()
    if text.startswith(("0o", "0b")):
        # Leave octal and binary literals alone.
        pass
    elif text.startswith("0x"):
        text = format_hex(text)
    elif "e" in text:
        text = format_scientific_notation(text)
    elif text.endswith(("j", "l")):
        text = format_long_or_complex_number(text)
    else:
        text = format_float_or_int_string(text)
    leaf.value = text
