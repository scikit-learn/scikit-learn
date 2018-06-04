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
"""Support module for tests for yapf."""

import difflib
import sys
import unittest

from yapf.yapflib import blank_line_calculator
from yapf.yapflib import comment_splicer
from yapf.yapflib import continuation_splicer
from yapf.yapflib import identify_container
from yapf.yapflib import pytree_unwrapper
from yapf.yapflib import pytree_utils
from yapf.yapflib import pytree_visitor
from yapf.yapflib import split_penalty
from yapf.yapflib import style
from yapf.yapflib import subtype_assigner


class YAPFTest(unittest.TestCase):

  def assertCodeEqual(self, expected_code, code):
    if code != expected_code:
      msg = ['Code format mismatch:', 'Expected:']
      linelen = style.Get('COLUMN_LIMIT')
      for l in expected_code.splitlines():
        if len(l) > linelen:
          msg.append('!> %s' % l)
        else:
          msg.append(' > %s' % l)
      msg.append('Actual:')
      for l in code.splitlines():
        if len(l) > linelen:
          msg.append('!> %s' % l)
        else:
          msg.append(' > %s' % l)
      msg.append('Diff:')
      msg.extend(
          difflib.unified_diff(
              code.splitlines(),
              expected_code.splitlines(),
              fromfile='actual',
              tofile='expected',
              lineterm=''))
      self.fail('\n'.join(msg))


def ParseAndUnwrap(code, dumptree=False):
  """Produces unwrapped lines from the given code.

  Parses the code into a tree, performs comment splicing and runs the
  unwrapper.

  Arguments:
    code: code to parse as a string
    dumptree: if True, the parsed pytree (after comment splicing) is dumped
              to stderr. Useful for debugging.

  Returns:
    List of unwrapped lines.
  """
  tree = pytree_utils.ParseCodeToTree(code)
  comment_splicer.SpliceComments(tree)
  continuation_splicer.SpliceContinuations(tree)
  subtype_assigner.AssignSubtypes(tree)
  identify_container.IdentifyContainers(tree)
  split_penalty.ComputeSplitPenalties(tree)
  blank_line_calculator.CalculateBlankLines(tree)

  if dumptree:
    pytree_visitor.DumpPyTree(tree, target_stream=sys.stderr)

  uwlines = pytree_unwrapper.UnwrapPyTree(tree)
  for uwl in uwlines:
    uwl.CalculateFormattingInformation()

  return uwlines
