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
"""Tests for yapf.comment_splicer."""

import textwrap
import unittest

from yapf.yapflib import comment_splicer
from yapf.yapflib import py3compat
from yapf.yapflib import pytree_utils


class CommentSplicerTest(unittest.TestCase):

  def _AssertNodeType(self, expected_type, node):
    self.assertEqual(expected_type, pytree_utils.NodeName(node))

  def _AssertNodeIsComment(self, node, text_in_comment=None):
    if pytree_utils.NodeName(node) == 'simple_stmt':
      self._AssertNodeType('COMMENT', node.children[0])
      node_value = node.children[0].value
    else:
      self._AssertNodeType('COMMENT', node)
      node_value = node.value
    if text_in_comment is not None:
      self.assertIn(text_in_comment, node_value)

  def _FindNthChildNamed(self, node, name, n=1):
    for i, child in enumerate(
        py3compat.ifilter(lambda c: pytree_utils.NodeName(c) == name,
                          node.pre_order())):
      if i == n - 1:
        return child
    raise RuntimeError('No Nth child for n={0}'.format(n))

  def testSimpleInline(self):
    code = 'foo = 1 # and a comment\n'
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    expr = tree.children[0].children[0]
    # Check that the expected node is still expr_stmt, but now it has 4 children
    # (before comment splicing it had 3), the last child being the comment.
    self._AssertNodeType('expr_stmt', expr)
    self.assertEqual(4, len(expr.children))
    comment_node = expr.children[3]
    self._AssertNodeIsComment(comment_node, '# and a comment')

  def testSimpleSeparateLine(self):
    code = textwrap.dedent(r'''
      foo = 1
      # first comment
      bar = 2
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # The comment should've been added to the root's children (now 4, including
    # the ENDMARKER in the end.
    self.assertEqual(4, len(tree.children))
    comment_node = tree.children[1]
    self._AssertNodeIsComment(comment_node)

  def testTwoLineComment(self):
    code = textwrap.dedent(r'''
      foo = 1
      # first comment
      # second comment
      bar = 2
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # This is similar to the single-line standalone comment.
    self.assertEqual(4, len(tree.children))
    self._AssertNodeIsComment(tree.children[1])

  def testCommentIsFirstChildInCompound(self):
    code = textwrap.dedent(r'''
      if x:
        # a comment
        foo = 1
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # Look into the suite node under the 'if'. We don't care about the NEWLINE
    # leaf but the new COMMENT must be a child of the suite and before the
    # actual code leaf.
    if_suite = tree.children[0].children[3]
    self._AssertNodeType('NEWLINE', if_suite.children[0])
    self._AssertNodeIsComment(if_suite.children[1])

  def testCommentIsLastChildInCompound(self):
    code = textwrap.dedent(r'''
      if x:
        foo = 1
        # a comment
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # Look into the suite node under the 'if'. We don't care about the DEDENT
    # leaf but the new COMMENT must be a child of the suite and after the
    # actual code leaf.
    if_suite = tree.children[0].children[3]
    self._AssertNodeType('DEDENT', if_suite.children[-1])
    self._AssertNodeIsComment(if_suite.children[-2])

  def testInlineAfterSeparateLine(self):
    code = textwrap.dedent(r'''
      bar = 1
      # line comment
      foo = 1 # inline comment
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # The separate line comment should become a child of the root, while
    # the inline comment remains within its simple_node.
    sep_comment_node = tree.children[1]
    self._AssertNodeIsComment(sep_comment_node, '# line comment')

    expr = tree.children[2].children[0]
    inline_comment_node = expr.children[-1]
    self._AssertNodeIsComment(inline_comment_node, '# inline comment')

  def testSeparateLineAfterInline(self):
    code = textwrap.dedent(r'''
      bar = 1
      foo = 1 # inline comment
      # line comment
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # The separate line comment should become a child of the root, while
    # the inline comment remains within its simple_node.
    sep_comment_node = tree.children[-2]
    self._AssertNodeIsComment(sep_comment_node, '# line comment')

    expr = tree.children[1].children[0]
    inline_comment_node = expr.children[-1]
    self._AssertNodeIsComment(inline_comment_node, '# inline comment')

  def testCommentBeforeDedent(self):
    code = textwrap.dedent(r'''
      if bar:
        z = 1
      # a comment
      j = 2
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # The comment should go under the tree root, not under the 'if'.
    self._AssertNodeIsComment(tree.children[1])
    if_suite = tree.children[0].children[3]
    self._AssertNodeType('DEDENT', if_suite.children[-1])

  def testCommentBeforeDedentTwoLevel(self):
    code = textwrap.dedent(r'''
      if foo:
        if bar:
          z = 1
        # a comment
      y = 1
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    if_suite = tree.children[0].children[3]
    # The comment is in the first if_suite, not the nested if under it. It's
    # right before the DEDENT
    self._AssertNodeIsComment(if_suite.children[-2])
    self._AssertNodeType('DEDENT', if_suite.children[-1])

  def testCommentBeforeDedentTwoLevelImproperlyIndented(self):
    code = textwrap.dedent(r'''
      if foo:
        if bar:
          z = 1
         # comment 2
      y = 1
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # The comment here is indented by 3 spaces, which is unlike any of the
    # surrounding statement indentation levels. The splicer attaches it to the
    # "closest" parent with smaller indentation.
    if_suite = tree.children[0].children[3]
    # The comment is in the first if_suite, not the nested if under it. It's
    # right before the DEDENT
    self._AssertNodeIsComment(if_suite.children[-2])
    self._AssertNodeType('DEDENT', if_suite.children[-1])

  def testCommentBeforeDedentThreeLevel(self):
    code = textwrap.dedent(r'''
      if foo:
        if bar:
          z = 1
          # comment 2
        # comment 1
      # comment 0
      j = 2
      ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    # comment 0 should go under the tree root
    self._AssertNodeIsComment(tree.children[1], '# comment 0')

    # comment 1 is in the first if_suite, right before the DEDENT
    if_suite_1 = self._FindNthChildNamed(tree, 'suite', n=1)
    self._AssertNodeIsComment(if_suite_1.children[-2], '# comment 1')
    self._AssertNodeType('DEDENT', if_suite_1.children[-1])

    # comment 2 is in if_suite nested under the first if suite,
    # right before the DEDENT
    if_suite_2 = self._FindNthChildNamed(tree, 'suite', n=2)
    self._AssertNodeIsComment(if_suite_2.children[-2], '# comment 2')
    self._AssertNodeType('DEDENT', if_suite_2.children[-1])

  def testCommentsInClass(self):
    code = textwrap.dedent(r'''
      class Foo:
        """docstring abc..."""
        # top-level comment
        def foo(): pass
        # another comment
      ''')

    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    class_suite = tree.children[0].children[3]
    another_comment = class_suite.children[-2]
    self._AssertNodeIsComment(another_comment, '# another')

    # It's OK for the comment to be a child of funcdef, as long as it's
    # the first child and thus comes before the 'def'.
    funcdef = class_suite.children[3]
    toplevel_comment = funcdef.children[0]
    self._AssertNodeIsComment(toplevel_comment, '# top-level')

  def testMultipleBlockComments(self):
    code = textwrap.dedent(r'''
        # Block comment number 1

        # Block comment number 2
        def f():
          pass
        ''')

    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    funcdef = tree.children[0]
    block_comment_1 = funcdef.children[0]
    self._AssertNodeIsComment(block_comment_1, '# Block comment number 1')

    block_comment_2 = funcdef.children[1]
    self._AssertNodeIsComment(block_comment_2, '# Block comment number 2')

  def testCommentsOnDedents(self):
    code = textwrap.dedent(r'''
        class Foo(object):
          # A comment for qux.
          def qux(self):
            pass

          # Interim comment.

          def mux(self):
            pass
        ''')

    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    classdef = tree.children[0]
    class_suite = classdef.children[6]
    qux_comment = class_suite.children[1]
    self._AssertNodeIsComment(qux_comment, '# A comment for qux.')

    interim_comment = class_suite.children[4]
    self._AssertNodeIsComment(interim_comment, '# Interim comment.')

  def testExprComments(self):
    code = textwrap.dedent(r'''
      foo( # Request fractions of an hour.
        948.0/3600, 20)
    ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    trailer = self._FindNthChildNamed(tree, 'trailer', 1)
    comment = trailer.children[1]
    self._AssertNodeIsComment(comment, '# Request fractions of an hour.')

  def testMultipleCommentsInOneExpr(self):
    code = textwrap.dedent(r'''
      foo( # com 1
        948.0/3600, # com 2
        20 + 12 # com 3
        )
    ''')
    tree = pytree_utils.ParseCodeToTree(code)
    comment_splicer.SpliceComments(tree)

    trailer = self._FindNthChildNamed(tree, 'trailer', 1)
    self._AssertNodeIsComment(trailer.children[1], '# com 1')

    arglist = self._FindNthChildNamed(tree, 'arglist', 1)
    self._AssertNodeIsComment(arglist.children[2], '# com 2')

    arith_expr = self._FindNthChildNamed(tree, 'arith_expr', 1)
    self._AssertNodeIsComment(arith_expr.children[-1], '# com 3')


if __name__ == '__main__':
  unittest.main()
