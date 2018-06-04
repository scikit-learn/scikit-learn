# Copyright 2016 Google Inc. All Rights Reserved.
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
"""Style config tests for yapf.reformatter."""

import textwrap
import unittest

from yapf.yapflib import reformatter
from yapf.yapflib import style

from yapftests import yapf_test_helper


class TestsForStyleConfig(yapf_test_helper.YAPFTest):

  def setUp(self):
    self.current_style = style.DEFAULT_STYLE

  def testSetGlobalStyle(self):
    try:
      style.SetGlobalStyle(style.CreateChromiumStyle())
      unformatted_code = textwrap.dedent(u"""\
          for i in range(5):
           print('bar')
          """)
      expected_formatted_code = textwrap.dedent(u"""\
          for i in range(5):
            print('bar')
          """)
      uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
      self.assertCodeEqual(expected_formatted_code,
                           reformatter.Reformat(uwlines))
    finally:
      style.SetGlobalStyle(style.CreatePEP8Style())
      style.DEFAULT_STYLE = self.current_style

    unformatted_code = textwrap.dedent(u"""\
        for i in range(5):
         print('bar')
        """)
    expected_formatted_code = textwrap.dedent(u"""\
        for i in range(5):
            print('bar')
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testOperatorStyle(self):
    try:
      sympy_style = style.CreatePEP8Style()
      sympy_style['NO_SPACES_AROUND_SELECTED_BINARY_OPERATORS'] = \
        style._StringSetConverter('*,/')
      style.SetGlobalStyle(sympy_style)
      unformatted_code = textwrap.dedent("""\
          a = 1+2 * 3 - 4 / 5
          """)
      expected_formatted_code = textwrap.dedent("""\
          a = 1 + 2*3 - 4/5
          """)

      uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
      self.assertCodeEqual(expected_formatted_code,
                           reformatter.Reformat(uwlines))
    finally:
      style.SetGlobalStyle(style.CreatePEP8Style())
      style.DEFAULT_STYLE = self.current_style


if __name__ == '__main__':
  unittest.main()
