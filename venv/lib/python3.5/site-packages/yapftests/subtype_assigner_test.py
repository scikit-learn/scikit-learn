# Copyright 2015 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Tests for yapf.subtype_assigner."""

import textwrap
import unittest

from yapf.yapflib import format_token
from yapf.yapflib import pytree_utils

from yapftests import yapf_test_helper


class SubtypeAssignerTest(yapf_test_helper.YAPFTest):

  def _CheckFormatTokenSubtypes(self, uwlines, list_of_expected):
    """Check that the tokens in the UnwrappedLines have the expected subtypes.

    Args:
      uwlines: list of UnwrappedLine.
      list_of_expected: list of (name, subtype) pairs. Non-semantic tokens are
        filtered out from the expected values.
    """
    actual = []
    for uwl in uwlines:
      filtered_values = [(ft.value, ft.subtypes)
                         for ft in uwl.tokens
                         if ft.name not in pytree_utils.NONSEMANTIC_TOKENS]
      if filtered_values:
        actual.append(filtered_values)

    self.assertEqual(list_of_expected, actual)

  def testFuncDefDefaultAssign(self):
    code = textwrap.dedent(r"""
        def foo(a=37, *b, **c):
          return -x[:42]
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self._CheckFormatTokenSubtypes(uwlines, [
        [('def', [format_token.Subtype.NONE]),
         ('foo', {format_token.Subtype.FUNC_DEF}),
         ('(', [format_token.Subtype.NONE]),
         ('a', {format_token.Subtype.NONE,
                format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         ('=', {format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN,
                format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         ('37', {format_token.Subtype.NONE,
                 format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         (',', {format_token.Subtype.NONE}),
         ('*', {format_token.Subtype.VARARGS_STAR,
                format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         ('b', {format_token.Subtype.NONE,
                format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         (',', {format_token.Subtype.NONE}),
         ('**', {format_token.Subtype.KWARGS_STAR_STAR,
                 format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         ('c', {format_token.Subtype.NONE,
                format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         (')', [format_token.Subtype.NONE]),
         (':', [format_token.Subtype.NONE])],
        [('return', [format_token.Subtype.NONE]),
         ('-', {format_token.Subtype.UNARY_OPERATOR}),
         ('x', [format_token.Subtype.NONE]),
         ('[', {format_token.Subtype.SUBSCRIPT_BRACKET}),
         (':', {format_token.Subtype.SUBSCRIPT_COLON}),
         ('42', [format_token.Subtype.NONE]),
         (']', {format_token.Subtype.SUBSCRIPT_BRACKET})],
    ])  # yapf: disable

  def testFuncCallWithDefaultAssign(self):
    code = textwrap.dedent(r"""
        foo(x, a='hello world')
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self._CheckFormatTokenSubtypes(uwlines, [
        [('foo', [format_token.Subtype.NONE]),
         ('(', [format_token.Subtype.NONE]),
         ('x', {format_token.Subtype.NONE,
                format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         (',', {format_token.Subtype.NONE}),
         ('a', {format_token.Subtype.NONE,
                format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST}),
         ('=', {format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN}),
         ("'hello world'", {format_token.Subtype.NONE}),
         (')', [format_token.Subtype.NONE])],
    ])  # yapf: disable

  def testSetComprehension(self):
    code = textwrap.dedent("""\
        def foo(strs):
          return {s.lower() for s in strs}
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self._CheckFormatTokenSubtypes(uwlines, [
        [('def', [format_token.Subtype.NONE]),
         ('foo', {format_token.Subtype.FUNC_DEF}),
         ('(', [format_token.Subtype.NONE]),
         ('strs', [format_token.Subtype.NONE]),
         (')', [format_token.Subtype.NONE]),
         (':', [format_token.Subtype.NONE])],
        [('return', [format_token.Subtype.NONE]),
         ('{', [format_token.Subtype.NONE]),
         ('s', {format_token.Subtype.COMP_EXPR}),
         ('.', {format_token.Subtype.COMP_EXPR}),
         ('lower', {format_token.Subtype.COMP_EXPR}),
         ('(', {format_token.Subtype.COMP_EXPR}),
         (')', {format_token.Subtype.COMP_EXPR}),
         ('for', {format_token.Subtype.DICT_SET_GENERATOR,
                  format_token.Subtype.COMP_FOR}),
         ('s', {format_token.Subtype.COMP_FOR}),
         ('in', {format_token.Subtype.COMP_FOR}),
         ('strs', {format_token.Subtype.COMP_FOR}),
         ('}', [format_token.Subtype.NONE])]
    ])  # yapf: disable

  def testUnaryNotOperator(self):
    code = textwrap.dedent("""\
        not a
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self._CheckFormatTokenSubtypes(uwlines, [
        [('not', {format_token.Subtype.UNARY_OPERATOR}),
         ('a', [format_token.Subtype.NONE])]
    ])  # yapf: disable

  def testBitwiseOperators(self):
    code = textwrap.dedent("""\
        x = ((a | (b ^ 3) & c) << 3) >> 1
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self._CheckFormatTokenSubtypes(uwlines, [
        [('x', [format_token.Subtype.NONE]),
         ('=', {format_token.Subtype.ASSIGN_OPERATOR}),
         ('(', [format_token.Subtype.NONE]),
         ('(', [format_token.Subtype.NONE]),
         ('a', [format_token.Subtype.NONE]),
         ('|', {format_token.Subtype.BINARY_OPERATOR}),
         ('(', [format_token.Subtype.NONE]),
         ('b', [format_token.Subtype.NONE]),
         ('^', {format_token.Subtype.BINARY_OPERATOR}),
         ('3', [format_token.Subtype.NONE]),
         (')', [format_token.Subtype.NONE]),
         ('&', {format_token.Subtype.BINARY_OPERATOR}),
         ('c', [format_token.Subtype.NONE]),
         (')', [format_token.Subtype.NONE]),
         ('<<', {format_token.Subtype.BINARY_OPERATOR}),
         ('3', [format_token.Subtype.NONE]),
         (')', [format_token.Subtype.NONE]),
         ('>>', {format_token.Subtype.BINARY_OPERATOR}),
         ('1', [format_token.Subtype.NONE]),],
    ])  # yapf: disable

  def testSubscriptColon(self):
    code = textwrap.dedent("""\
        x[0:42:1]
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self._CheckFormatTokenSubtypes(uwlines, [
        [
            ('x', [format_token.Subtype.NONE]),
            ('[', {format_token.Subtype.SUBSCRIPT_BRACKET}),
            ('0', [format_token.Subtype.NONE]),
            (':', {format_token.Subtype.SUBSCRIPT_COLON}),
            ('42', [format_token.Subtype.NONE]),
            (':', {format_token.Subtype.SUBSCRIPT_COLON}),
            ('1', [format_token.Subtype.NONE]),
            (']', {format_token.Subtype.SUBSCRIPT_BRACKET}),
        ],
    ])

  def testFunctionCallWithStarExpression(self):
    code = textwrap.dedent("""\
        [a, *b]
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self._CheckFormatTokenSubtypes(uwlines, [
        [('[', [format_token.Subtype.NONE]),
         ('a', [format_token.Subtype.NONE]),
         (',', [format_token.Subtype.NONE]),
         ('*', {format_token.Subtype.UNARY_OPERATOR,
                format_token.Subtype.VARARGS_STAR}),
         ('b', [format_token.Subtype.NONE]),
         (']', [format_token.Subtype.NONE]),],
    ])  # yapf: disable


if __name__ == '__main__':
  unittest.main()
