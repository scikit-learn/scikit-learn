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
"""Tests for yapf.split_penalty."""

import sys
import textwrap
import unittest

from lib2to3 import pytree

from yapf.yapflib import pytree_utils
from yapf.yapflib import pytree_visitor
from yapf.yapflib import split_penalty

UNBREAKABLE = split_penalty.UNBREAKABLE
VERY_STRONGLY_CONNECTED = split_penalty.VERY_STRONGLY_CONNECTED
DOTTED_NAME = split_penalty.DOTTED_NAME
STRONGLY_CONNECTED = split_penalty.STRONGLY_CONNECTED


class SplitPenaltyTest(unittest.TestCase):

  def _ParseAndComputePenalties(self, code, dumptree=False):
    """Parses the code and computes split penalties.

    Arguments:
      code: code to parse as a string
      dumptree: if True, the parsed pytree (after penalty assignment) is dumped
        to stderr. Useful for debugging.

    Returns:
      Parse tree.
    """
    tree = pytree_utils.ParseCodeToTree(code)
    split_penalty.ComputeSplitPenalties(tree)
    if dumptree:
      pytree_visitor.DumpPyTree(tree, target_stream=sys.stderr)
    return tree

  def _CheckPenalties(self, tree, list_of_expected):
    """Check that the tokens in the tree have the correct penalties.

    Args:
      tree: the pytree.
      list_of_expected: list of (name, penalty) pairs. Non-semantic tokens are
        filtered out from the expected values.
    """

    def FlattenRec(tree):
      if pytree_utils.NodeName(tree) in pytree_utils.NONSEMANTIC_TOKENS:
        return []
      if isinstance(tree, pytree.Leaf):
        return [(tree.value,
                 pytree_utils.GetNodeAnnotation(
                     tree, pytree_utils.Annotation.SPLIT_PENALTY))]
      nodes = []
      for node in tree.children:
        nodes += FlattenRec(node)
      return nodes

    self.assertEqual(list_of_expected, FlattenRec(tree))

  def testUnbreakable(self):
    # Test function definitions.
    code = textwrap.dedent(r"""
      def foo(x):
        pass
      """)
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('def', None),
        ('foo', UNBREAKABLE),
        ('(', UNBREAKABLE),
        ('x', None),
        (')', STRONGLY_CONNECTED),
        (':', UNBREAKABLE),
        ('pass', None),
    ])

    # Test function definition with trailing comment.
    code = textwrap.dedent(r"""
      def foo(x):  # trailing comment
        pass
      """)
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('def', None),
        ('foo', UNBREAKABLE),
        ('(', UNBREAKABLE),
        ('x', None),
        (')', STRONGLY_CONNECTED),
        (':', UNBREAKABLE),
        ('pass', None),
    ])

    # Test class definitions.
    code = textwrap.dedent(r"""
      class A:
        pass
      class B(A):
        pass
      """)
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('class', None),
        ('A', UNBREAKABLE),
        (':', UNBREAKABLE),
        ('pass', None),
        ('class', None),
        ('B', UNBREAKABLE),
        ('(', UNBREAKABLE),
        ('A', None),
        (')', None),
        (':', UNBREAKABLE),
        ('pass', None),
    ])

    # Test lambda definitions.
    code = textwrap.dedent(r"""
      lambda a, b: None
      """)
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('lambda', None),
        ('a', UNBREAKABLE),
        (',', UNBREAKABLE),
        ('b', UNBREAKABLE),
        (':', UNBREAKABLE),
        ('None', UNBREAKABLE),
    ])

    # Test dotted names.
    code = textwrap.dedent(r"""
      import a.b.c
      """)
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('import', None),
        ('a', None),
        ('.', UNBREAKABLE),
        ('b', UNBREAKABLE),
        ('.', UNBREAKABLE),
        ('c', UNBREAKABLE),
    ])

  def testStronglyConnected(self):
    # Test dictionary keys.
    code = textwrap.dedent(r"""
      a = {
          'x': 42,
          y(lambda a: 23): 37,
      }
      """)
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('a', None),
        ('=', None),
        ('{', None),
        ("'x'", None),
        (':', STRONGLY_CONNECTED),
        ('42', None),
        (',', None),
        ('y', None),
        ('(', UNBREAKABLE),
        ('lambda', STRONGLY_CONNECTED),
        ('a', UNBREAKABLE),
        (':', UNBREAKABLE),
        ('23', UNBREAKABLE),
        (')', VERY_STRONGLY_CONNECTED),
        (':', STRONGLY_CONNECTED),
        ('37', None),
        (',', None),
        ('}', None),
    ])

    # Test list comprehension.
    code = textwrap.dedent(r"""
      [a for a in foo if a.x == 37]
      """)
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('[', None),
        ('a', None),
        ('for', 0),
        ('a', STRONGLY_CONNECTED),
        ('in', STRONGLY_CONNECTED),
        ('foo', STRONGLY_CONNECTED),
        ('if', 0),
        ('a', STRONGLY_CONNECTED),
        ('.', UNBREAKABLE),
        ('x', DOTTED_NAME),
        ('==', STRONGLY_CONNECTED),
        ('37', STRONGLY_CONNECTED),
        (']', None),
    ])

  def testFuncCalls(self):
    code = 'foo(1, 2, 3)\n'
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('foo', None),
        ('(', UNBREAKABLE),
        ('1', None),
        (',', UNBREAKABLE),
        ('2', None),
        (',', UNBREAKABLE),
        ('3', None),
        (')', VERY_STRONGLY_CONNECTED),
    ])

    # Now a method call, which has more than one trailer
    code = 'foo.bar.baz(1, 2, 3)\n'
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('foo', None),
        ('.', UNBREAKABLE),
        ('bar', DOTTED_NAME),
        ('.', STRONGLY_CONNECTED),
        ('baz', DOTTED_NAME),
        ('(', STRONGLY_CONNECTED),
        ('1', None),
        (',', UNBREAKABLE),
        ('2', None),
        (',', UNBREAKABLE),
        ('3', None),
        (')', VERY_STRONGLY_CONNECTED),
    ])

    # Test single generator argument.
    code = 'max(i for i in xrange(10))\n'
    tree = self._ParseAndComputePenalties(code)
    self._CheckPenalties(tree, [
        ('max', None),
        ('(', UNBREAKABLE),
        ('i', 0),
        ('for', 0),
        ('i', STRONGLY_CONNECTED),
        ('in', STRONGLY_CONNECTED),
        ('xrange', STRONGLY_CONNECTED),
        ('(', UNBREAKABLE),
        ('10', STRONGLY_CONNECTED),
        (')', VERY_STRONGLY_CONNECTED),
        (')', VERY_STRONGLY_CONNECTED),
    ])


if __name__ == '__main__':
  unittest.main()
