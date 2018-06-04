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
"""Tests for yapf.reformatter."""

import textwrap
import unittest

from yapf.yapflib import py3compat
from yapf.yapflib import reformatter
from yapf.yapflib import style
from yapf.yapflib import verifier

from yapftests import yapf_test_helper


@unittest.skipIf(py3compat.PY3, 'Requires Python 2')
class TestVerifyNoVerify(yapf_test_helper.YAPFTest):

  @classmethod
  def setUpClass(cls):
    style.SetGlobalStyle(style.CreatePEP8Style())

  def testVerifyException(self):
    unformatted_code = textwrap.dedent("""\
        class ABC(metaclass=type):
          pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    with self.assertRaises(verifier.InternalError):
      reformatter.Reformat(uwlines, verify=True)
    reformatter.Reformat(uwlines)  # verify should be False by default.

  def testNoVerify(self):
    unformatted_code = textwrap.dedent("""\
        class ABC(metaclass=type):
          pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        class ABC(metaclass=type):
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code,
                         reformatter.Reformat(uwlines, verify=False))

  def testVerifyFutureImport(self):
    unformatted_code = textwrap.dedent("""\
        from __future__ import print_function

        def call_my_function(the_function):
          the_function("hi")

        if __name__ == "__main__":
          call_my_function(print)
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    with self.assertRaises(verifier.InternalError):
      reformatter.Reformat(uwlines, verify=True)

    expected_formatted_code = textwrap.dedent("""\
        from __future__ import print_function


        def call_my_function(the_function):
            the_function("hi")


        if __name__ == "__main__":
            call_my_function(print)
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code,
                         reformatter.Reformat(uwlines, verify=False))

  def testContinuationLineShouldBeDistinguished(self):
    unformatted_code = textwrap.dedent("""\
        class Foo(object):

            def bar(self):
                if self.solo_generator_that_is_long is None and len(
                        self.generators + self.next_batch) == 1:
                    pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        class Foo(object):
            def bar(self):
                if self.solo_generator_that_is_long is None and len(
                        self.generators + self.next_batch) == 1:
                    pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code,
                         reformatter.Reformat(uwlines, verify=False))


if __name__ == '__main__':
  unittest.main()
