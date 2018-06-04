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
"""Tests for yapf.pytree_utils."""

import unittest

from lib2to3 import pygram
from lib2to3 import pytree
from lib2to3.pgen2 import token

from yapf.yapflib import pytree_utils

# More direct access to the symbol->number mapping living within the grammar
# module.
_GRAMMAR_SYMBOL2NUMBER = pygram.python_grammar.symbol2number

_FOO = 'foo'
_FOO1 = 'foo1'
_FOO2 = 'foo2'
_FOO3 = 'foo3'
_FOO4 = 'foo4'
_FOO5 = 'foo5'


class NodeNameTest(unittest.TestCase):

  def testNodeNameForLeaf(self):
    leaf = pytree.Leaf(token.LPAR, '(')
    self.assertEqual('LPAR', pytree_utils.NodeName(leaf))

  def testNodeNameForNode(self):
    leaf = pytree.Leaf(token.LPAR, '(')
    node = pytree.Node(pygram.python_grammar.symbol2number['suite'], [leaf])
    self.assertEqual('suite', pytree_utils.NodeName(node))


class ParseCodeToTreeTest(unittest.TestCase):

  def testParseCodeToTree(self):
    # Since ParseCodeToTree is a thin wrapper around underlying lib2to3
    # functionality, only a sanity test here...
    tree = pytree_utils.ParseCodeToTree('foo = 2\n')
    self.assertEqual('file_input', pytree_utils.NodeName(tree))
    self.assertEqual(2, len(tree.children))
    self.assertEqual('simple_stmt', pytree_utils.NodeName(tree.children[0]))

  def testPrintFunctionToTree(self):
    tree = pytree_utils.ParseCodeToTree(
        'print("hello world", file=sys.stderr)\n')
    self.assertEqual('file_input', pytree_utils.NodeName(tree))
    self.assertEqual(2, len(tree.children))
    self.assertEqual('simple_stmt', pytree_utils.NodeName(tree.children[0]))

  def testPrintStatementToTree(self):
    tree = pytree_utils.ParseCodeToTree('print "hello world"\n')
    self.assertEqual('file_input', pytree_utils.NodeName(tree))
    self.assertEqual(2, len(tree.children))
    self.assertEqual('simple_stmt', pytree_utils.NodeName(tree.children[0]))

  def testClassNotLocal(self):
    tree = pytree_utils.ParseCodeToTree('class nonlocal: pass\n')
    self.assertEqual('file_input', pytree_utils.NodeName(tree))
    self.assertEqual(2, len(tree.children))
    self.assertEqual('classdef', pytree_utils.NodeName(tree.children[0]))


class InsertNodesBeforeAfterTest(unittest.TestCase):

  def _BuildSimpleTree(self):
    # Builds a simple tree we can play with in the tests.
    # The tree looks like this:
    #
    #   suite:
    #     LPAR
    #     LPAR
    #     simple_stmt:
    #       NAME('foo')
    #
    lpar1 = pytree.Leaf(token.LPAR, '(')
    lpar2 = pytree.Leaf(token.LPAR, '(')
    simple_stmt = pytree.Node(_GRAMMAR_SYMBOL2NUMBER['simple_stmt'],
                              [pytree.Leaf(token.NAME, 'foo')])
    return pytree.Node(_GRAMMAR_SYMBOL2NUMBER['suite'],
                       [lpar1, lpar2, simple_stmt])

  def _MakeNewNodeRPAR(self):
    return pytree.Leaf(token.RPAR, ')')

  def setUp(self):
    self._simple_tree = self._BuildSimpleTree()

  def testInsertNodesBefore(self):
    # Insert before simple_stmt and make sure it went to the right place
    pytree_utils.InsertNodesBefore([self._MakeNewNodeRPAR()],
                                   self._simple_tree.children[2])
    self.assertEqual(4, len(self._simple_tree.children))
    self.assertEqual('RPAR',
                     pytree_utils.NodeName(self._simple_tree.children[2]))
    self.assertEqual('simple_stmt',
                     pytree_utils.NodeName(self._simple_tree.children[3]))

  def testInsertNodesBeforeFirstChild(self):
    # Insert before the first child of its parent
    simple_stmt = self._simple_tree.children[2]
    foo_child = simple_stmt.children[0]
    pytree_utils.InsertNodesBefore([self._MakeNewNodeRPAR()], foo_child)
    self.assertEqual(3, len(self._simple_tree.children))
    self.assertEqual(2, len(simple_stmt.children))
    self.assertEqual('RPAR', pytree_utils.NodeName(simple_stmt.children[0]))
    self.assertEqual('NAME', pytree_utils.NodeName(simple_stmt.children[1]))

  def testInsertNodesAfter(self):
    # Insert after and make sure it went to the right place
    pytree_utils.InsertNodesAfter([self._MakeNewNodeRPAR()],
                                  self._simple_tree.children[2])
    self.assertEqual(4, len(self._simple_tree.children))
    self.assertEqual('simple_stmt',
                     pytree_utils.NodeName(self._simple_tree.children[2]))
    self.assertEqual('RPAR',
                     pytree_utils.NodeName(self._simple_tree.children[3]))

  def testInsertNodesAfterLastChild(self):
    # Insert after the last child of its parent
    simple_stmt = self._simple_tree.children[2]
    foo_child = simple_stmt.children[0]
    pytree_utils.InsertNodesAfter([self._MakeNewNodeRPAR()], foo_child)
    self.assertEqual(3, len(self._simple_tree.children))
    self.assertEqual(2, len(simple_stmt.children))
    self.assertEqual('NAME', pytree_utils.NodeName(simple_stmt.children[0]))
    self.assertEqual('RPAR', pytree_utils.NodeName(simple_stmt.children[1]))

  def testInsertNodesWhichHasParent(self):
    # Try to insert an existing tree node into another place and fail.
    with self.assertRaises(RuntimeError):
      pytree_utils.InsertNodesAfter([self._simple_tree.children[1]],
                                    self._simple_tree.children[0])


class AnnotationsTest(unittest.TestCase):

  def setUp(self):
    self._leaf = pytree.Leaf(token.LPAR, '(')
    self._node = pytree.Node(_GRAMMAR_SYMBOL2NUMBER['simple_stmt'],
                             [pytree.Leaf(token.NAME, 'foo')])

  def testGetWhenNone(self):
    self.assertIsNone(pytree_utils.GetNodeAnnotation(self._leaf, _FOO))

  def testSetWhenNone(self):
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO, 20)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO), 20)

  def testSetAgain(self):
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO, 20)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO), 20)
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO, 30)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO), 30)

  def testMultiple(self):
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO, 20)
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO1, 1)
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO2, 2)
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO3, 3)
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO4, 4)
    pytree_utils.SetNodeAnnotation(self._leaf, _FOO5, 5)

    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO), 20)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO1), 1)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO2), 2)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO3), 3)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO4), 4)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._leaf, _FOO5), 5)

  def testSubtype(self):
    pytree_utils.AppendNodeAnnotation(self._leaf,
                                      pytree_utils.Annotation.SUBTYPE, _FOO)

    self.assertSetEqual(
        pytree_utils.GetNodeAnnotation(self._leaf,
                                       pytree_utils.Annotation.SUBTYPE), {_FOO})

    pytree_utils.RemoveSubtypeAnnotation(self._leaf, _FOO)

    self.assertSetEqual(
        pytree_utils.GetNodeAnnotation(self._leaf,
                                       pytree_utils.Annotation.SUBTYPE), set())

  def testSetOnNode(self):
    pytree_utils.SetNodeAnnotation(self._node, _FOO, 20)
    self.assertEqual(pytree_utils.GetNodeAnnotation(self._node, _FOO), 20)


if __name__ == '__main__':
  unittest.main()
