# -*- coding: UTF-8 -*-
# pylint: disable=missing-docstring, too-few-public-methods
"""
Test the trait-type ``UseEnum``.
"""

import unittest
import enum
from ipython_genutils.py3compat import string_types
from traitlets import HasTraits, TraitError, UseEnum


# -----------------------------------------------------------------------------
# TEST SUPPORT:
# -----------------------------------------------------------------------------
class Color(enum.Enum):
    red = 1
    green = 2
    blue = 3
    yellow = 4

class OtherColor(enum.Enum):
    red = 0
    green = 1


# -----------------------------------------------------------------------------
# TESTSUITE:
# -----------------------------------------------------------------------------
class TestUseEnum(unittest.TestCase):
    # pylint: disable=invalid-name

    class Example(HasTraits):
        color = UseEnum(Color, help="Color enum")

    def test_assign_enum_value(self):
        example = self.Example()
        example.color = Color.green
        self.assertEqual(example.color, Color.green)

    def test_assign_all_enum_values(self):
        # pylint: disable=no-member
        enum_values = [value  for value in Color.__members__.values()]
        for value in enum_values:
            self.assertIsInstance(value, Color)
            example = self.Example()
            example.color = value
            self.assertEqual(example.color, value)
            self.assertIsInstance(value, Color)

    def test_assign_enum_value__with_other_enum_raises_error(self):
        example = self.Example()
        with self.assertRaises(TraitError):
            example.color = OtherColor.green

    def test_assign_enum_name_1(self):
        # -- CONVERT: string => Enum value (item)
        example = self.Example()
        example.color = "red"
        self.assertEqual(example.color, Color.red)

    def test_assign_enum_value_name(self):
        # -- CONVERT: string => Enum value (item)
        # pylint: disable=no-member
        enum_names = [enum_val.name for enum_val in Color.__members__.values()]
        for value in enum_names:
            self.assertIsInstance(value, string_types)
            example = self.Example()
            enum_value = Color.__members__.get(value)
            example.color = value
            self.assertIs(example.color, enum_value)
            self.assertEqual(example.color.name, value)

    def test_assign_scoped_enum_value_name(self):
        # -- CONVERT: string => Enum value (item)
        scoped_names = ["Color.red", "Color.green", "Color.blue", "Color.yellow"]
        for value in scoped_names:
            example = self.Example()
            example.color = value
            self.assertIsInstance(example.color, Color)
            self.assertEqual(str(example.color), value)

    def test_assign_bad_enum_value_name__raises_error(self):
        # -- CONVERT: string => Enum value (item)
        bad_enum_names = ["UNKNOWN_COLOR", "RED", "Green", "blue2"]
        for value in bad_enum_names:
            example = self.Example()
            with self.assertRaises(TraitError):
                example.color = value

    def test_assign_enum_value_number_1(self):
        # -- CONVERT: number => Enum value (item)
        example = self.Example()
        example.color = 1  # == Color.red.value
        example.color = Color.red.value
        self.assertEqual(example.color, Color.red)

    def test_assign_enum_value_number(self):
        # -- CONVERT: number => Enum value (item)
        # pylint: disable=no-member
        enum_numbers = [enum_val.value
                        for enum_val in Color.__members__.values()]
        for value in enum_numbers:
            self.assertIsInstance(value, int)
            example = self.Example()
            example.color = value
            self.assertIsInstance(example.color, Color)
            self.assertEqual(example.color.value, value)

    def test_assign_bad_enum_value_number__raises_error(self):
        # -- CONVERT: number => Enum value (item)
        bad_numbers = [-1, 0, 5]
        for value in bad_numbers:
            self.assertIsInstance(value, int)
            assert UseEnum(Color).select_by_number(value, None) is None
            example = self.Example()
            with self.assertRaises(TraitError):
                example.color = value

    def test_ctor_without_default_value(self):
        # -- IMPLICIT: default_value = Color.red (first enum-value)
        class Example2(HasTraits):
            color = UseEnum(Color)

        example = Example2()
        self.assertEqual(example.color, Color.red)

    def test_ctor_with_default_value_as_enum_value(self):
        # -- CONVERT: number => Enum value (item)
        class Example2(HasTraits):
            color = UseEnum(Color, default_value=Color.green)

        example = Example2()
        self.assertEqual(example.color, Color.green)


    def test_ctor_with_default_value_none_and_not_allow_none(self):
        # -- IMPLICIT: default_value = Color.red (first enum-value)
        class Example2(HasTraits):
            color1 = UseEnum(Color, default_value=None, allow_none=False)
            color2 = UseEnum(Color, default_value=None)
        example = Example2()
        self.assertEqual(example.color1, Color.red)
        self.assertEqual(example.color2, Color.red)

    def test_ctor_with_default_value_none_and_allow_none(self):
        class Example2(HasTraits):
            color1 = UseEnum(Color, default_value=None, allow_none=True)
            color2 = UseEnum(Color, allow_none=True)

        example = Example2()
        self.assertIs(example.color1, None)
        self.assertIs(example.color2, None)

    def test_assign_none_without_allow_none_resets_to_default_value(self):
        class Example2(HasTraits):
            color1 = UseEnum(Color, allow_none=False)
            color2 = UseEnum(Color)

        example = Example2()
        example.color1 = None
        example.color2 = None
        self.assertIs(example.color1, Color.red)
        self.assertIs(example.color2, Color.red)

    def test_assign_none_to_enum_or_none(self):
        class Example2(HasTraits):
            color = UseEnum(Color, allow_none=True)

        example = Example2()
        example.color = None
        self.assertIs(example.color, None)

    def test_assign_bad_value_with_to_enum_or_none(self):
        class Example2(HasTraits):
            color = UseEnum(Color, allow_none=True)

        example = Example2()
        with self.assertRaises(TraitError):
            example.color = "BAD_VALUE"

