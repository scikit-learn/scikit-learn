"""Test cases for is_subtype and is_runtime_subtype."""

import unittest

from mypyc.ir.rtypes import bit_rprimitive, bool_rprimitive, int_rprimitive
from mypyc.subtype import is_subtype
from mypyc.rt_subtype import is_runtime_subtype


class TestSubtype(unittest.TestCase):
    def test_bit(self) -> None:
        assert is_subtype(bit_rprimitive, bool_rprimitive)
        assert is_subtype(bit_rprimitive, int_rprimitive)

    def test_bool(self) -> None:
        assert not is_subtype(bool_rprimitive, bit_rprimitive)
        assert is_subtype(bool_rprimitive, int_rprimitive)


class TestRuntimeSubtype(unittest.TestCase):
    def test_bit(self) -> None:
        assert is_runtime_subtype(bit_rprimitive, bool_rprimitive)
        assert not is_runtime_subtype(bit_rprimitive, int_rprimitive)

    def test_bool(self) -> None:
        assert not is_runtime_subtype(bool_rprimitive, bit_rprimitive)
        assert not is_runtime_subtype(bool_rprimitive, int_rprimitive)
