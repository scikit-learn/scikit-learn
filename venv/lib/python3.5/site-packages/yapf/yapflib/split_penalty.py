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
"""Computation of split penalties before/between tokens."""

import re

from lib2to3 import pytree

from yapf.yapflib import format_token
from yapf.yapflib import py3compat
from yapf.yapflib import pytree_utils
from yapf.yapflib import pytree_visitor
from yapf.yapflib import style

# TODO(morbo): Document the annotations in a centralized place. E.g., the
# README file.
UNBREAKABLE = 1000 * 1000
NAMED_ASSIGN = 11000
DOTTED_NAME = 4000
VERY_STRONGLY_CONNECTED = 3500
STRONGLY_CONNECTED = 3000
CONNECTED = 500

OR_TEST = 1000
AND_TEST = 1100
NOT_TEST = 1200
COMPARISON = 1300
STAR_EXPR = 1300
EXPR = 1400
XOR_EXPR = 1500
AND_EXPR = 1700
SHIFT_EXPR = 1800
ARITH_EXPR = 1900
TERM = 2000
FACTOR = 2100
POWER = 2200
ATOM = 2300
ONE_ELEMENT_ARGUMENT = 2500


def ComputeSplitPenalties(tree):
  """Compute split penalties on tokens in the given parse tree.

  Arguments:
    tree: the top-level pytree node to annotate with penalties.
  """
  _SplitPenaltyAssigner().Visit(tree)


class _SplitPenaltyAssigner(pytree_visitor.PyTreeVisitor):
  """Assigns split penalties to tokens, based on parse tree structure.

  Split penalties are attached as annotations to tokens.
  """

  def Visit_import_as_names(self, node):  # pyline: disable=invalid-name
    # import_as_names ::= import_as_name (',' import_as_name)* [',']
    self.DefaultNodeVisit(node)
    prev_child = None
    for child in node.children:
      if (prev_child and isinstance(prev_child, pytree.Leaf) and
          prev_child.value == ','):
        _SetSplitPenalty(child, style.Get('SPLIT_PENALTY_IMPORT_NAMES'))
      prev_child = child

  def Visit_classdef(self, node):  # pylint: disable=invalid-name
    # classdef ::= 'class' NAME ['(' [arglist] ')'] ':' suite
    #
    # NAME
    _SetUnbreakable(node.children[1])
    if len(node.children) > 4:
      # opening '('
      _SetUnbreakable(node.children[2])
    # ':'
    _SetUnbreakable(node.children[-2])
    self.DefaultNodeVisit(node)

  def Visit_funcdef(self, node):  # pylint: disable=invalid-name
    # funcdef ::= 'def' NAME parameters ['->' test] ':' suite
    #
    # Can't break before the function name and before the colon. The parameters
    # are handled by child iteration.
    colon_idx = 1
    while pytree_utils.NodeName(node.children[colon_idx]) == 'simple_stmt':
      colon_idx += 1
    _SetUnbreakable(node.children[colon_idx])
    arrow_idx = -1
    while colon_idx < len(node.children):
      if isinstance(node.children[colon_idx], pytree.Leaf):
        if node.children[colon_idx].value == ':':
          break
        if node.children[colon_idx].value == '->':
          arrow_idx = colon_idx
      colon_idx += 1
    _SetUnbreakable(node.children[colon_idx])
    self.DefaultNodeVisit(node)
    if arrow_idx > 0:
      _SetSplitPenalty(
          pytree_utils.LastLeafNode(node.children[arrow_idx - 1]), 0)
      _SetUnbreakable(node.children[arrow_idx])
      _SetStronglyConnected(node.children[arrow_idx + 1])

  def Visit_lambdef(self, node):  # pylint: disable=invalid-name
    # lambdef ::= 'lambda' [varargslist] ':' test
    # Loop over the lambda up to and including the colon.
    allow_multiline_lambdas = style.Get('ALLOW_MULTILINE_LAMBDAS')
    if not allow_multiline_lambdas:
      for child in node.children:
        if pytree_utils.NodeName(child) == 'COMMENT':
          if re.search(r'pylint:.*disable=.*\bg-long-lambda', child.value):
            allow_multiline_lambdas = True
            break

    if allow_multiline_lambdas:
      _SetStronglyConnected(node)
    else:
      self._SetUnbreakableOnChildren(node)

  def Visit_parameters(self, node):  # pylint: disable=invalid-name
    # parameters ::= '(' [typedargslist] ')'
    self.DefaultNodeVisit(node)

    # Can't break before the opening paren of a parameter list.
    _SetUnbreakable(node.children[0])
    if not style.Get('DEDENT_CLOSING_BRACKETS'):
      _SetStronglyConnected(node.children[-1])

  def Visit_arglist(self, node):  # pylint: disable=invalid-name
    # arglist ::= argument (',' argument)* [',']
    self.DefaultNodeVisit(node)
    index = 1
    while index < len(node.children):
      child = node.children[index]
      if isinstance(child, pytree.Leaf) and child.value == ',':
        _SetUnbreakable(child)
      index += 1
    for child in node.children:
      if pytree_utils.NodeName(child) == 'atom':
        _IncreasePenalty(child, CONNECTED)

  def Visit_argument(self, node):  # pylint: disable=invalid-name
    # argument ::= test [comp_for] | test '=' test  # Really [keyword '='] test
    self.DefaultNodeVisit(node)
    index = 1
    while index < len(node.children) - 1:
      child = node.children[index]
      if isinstance(child, pytree.Leaf) and child.value == '=':
        _SetSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index]), NAMED_ASSIGN)
        _SetSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index + 1]), NAMED_ASSIGN)
      index += 1

  def Visit_tname(self, node):  # pylint: disable=invalid-name
    # tname ::= NAME [':' test]
    self.DefaultNodeVisit(node)
    index = 1
    while index < len(node.children) - 1:
      child = node.children[index]
      if isinstance(child, pytree.Leaf) and child.value == ':':
        _SetSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index]), NAMED_ASSIGN)
        _SetSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index + 1]), NAMED_ASSIGN)
      index += 1

  def Visit_dotted_name(self, node):  # pylint: disable=invalid-name
    # dotted_name ::= NAME ('.' NAME)*
    self._SetUnbreakableOnChildren(node)

  def Visit_dictsetmaker(self, node):  # pylint: disable=invalid-name
    # dictsetmaker ::= ( (test ':' test
    #                      (comp_for | (',' test ':' test)* [','])) |
    #                    (test (comp_for | (',' test)* [','])) )
    for child in node.children:
      self.Visit(child)
      if pytree_utils.NodeName(child) == 'COLON':
        # This is a key to a dictionary. We don't want to split the key if at
        # all possible.
        _SetStronglyConnected(child)

  def Visit_trailer(self, node):  # pylint: disable=invalid-name
    # trailer ::= '(' [arglist] ')' | '[' subscriptlist ']' | '.' NAME
    if node.children[0].value == '.':
      self._SetUnbreakableOnChildren(node)
      _SetSplitPenalty(node.children[1], DOTTED_NAME)
    elif len(node.children) == 2:
      # Don't split an empty argument list if at all possible.
      _SetSplitPenalty(node.children[1], VERY_STRONGLY_CONNECTED)
    elif len(node.children) == 3:
      name = pytree_utils.NodeName(node.children[1])
      if name in {'argument', 'comparison'}:
        # Don't split an argument list with one element if at all possible.
        _SetStronglyConnected(node.children[1])
        if (len(node.children[1].children) > 1 and
            pytree_utils.NodeName(node.children[1].children[1]) == 'comp_for'):
          # Don't penalize splitting before a comp_for expression.
          _SetSplitPenalty(pytree_utils.FirstLeafNode(node.children[1]), 0)
        else:
          _SetSplitPenalty(
              pytree_utils.FirstLeafNode(node.children[1]),
              ONE_ELEMENT_ARGUMENT)
      elif (pytree_utils.NodeName(node.children[0]) == 'LSQB' and
            len(node.children[1].children) > 2 and
            (name.endswith('_test') or name.endswith('_expr'))):
        _SetStronglyConnected(node.children[1].children[0])
        _SetStronglyConnected(node.children[1].children[2])

        # Still allow splitting around the operator.
        split_before = ((name.endswith('_test') and
                         style.Get('SPLIT_BEFORE_LOGICAL_OPERATOR')) or
                        (name.endswith('_expr') and
                         style.Get('SPLIT_BEFORE_BITWISE_OPERATOR')))
        if split_before:
          _SetSplitPenalty(
              pytree_utils.LastLeafNode(node.children[1].children[1]), 0)
        else:
          _SetSplitPenalty(
              pytree_utils.FirstLeafNode(node.children[1].children[2]), 0)

        # Don't split the ending bracket of a subscript list.
        _SetVeryStronglyConnected(node.children[-1])
      elif name not in {
          'arglist', 'argument', 'term', 'or_test', 'and_test', 'comparison',
          'atom', 'power'
      }:
        # Don't split an argument list with one element if at all possible.
        _SetStronglyConnected(node.children[1], node.children[2])

      if name == 'arglist':
        _SetStronglyConnected(node.children[-1])

    self.DefaultNodeVisit(node)

  def Visit_power(self, node):  # pylint: disable=invalid-name,missing-docstring
    # power ::= atom trailer* ['**' factor]
    self.DefaultNodeVisit(node)

    # When atom is followed by a trailer, we can not break between them.
    # E.g. arr[idx] - no break allowed between 'arr' and '['.
    if (len(node.children) > 1 and
        pytree_utils.NodeName(node.children[1]) == 'trailer'):
      # children[1] itself is a whole trailer: we don't want to
      # mark all of it as unbreakable, only its first token: (, [ or .
      _SetUnbreakable(node.children[1].children[0])

      # A special case when there are more trailers in the sequence. Given:
      #   atom tr1 tr2
      # The last token of tr1 and the first token of tr2 comprise an unbreakable
      # region. For example: foo.bar.baz(1)
      # We can't put breaks between either of the '.', '(', or '[' and the names
      # *preceding* them.
      prev_trailer_idx = 1
      while prev_trailer_idx < len(node.children) - 1:
        cur_trailer_idx = prev_trailer_idx + 1
        cur_trailer = node.children[cur_trailer_idx]
        if pytree_utils.NodeName(cur_trailer) == 'trailer':
          # Now we know we have two trailers one after the other
          prev_trailer = node.children[prev_trailer_idx]
          if prev_trailer.children[-1].value != ')':
            # Set the previous node unbreakable if it's not a function call:
            #   atom tr1() tr2
            # It may be necessary (though undesirable) to split up a previous
            # function call's parentheses to the next line.
            _SetStronglyConnected(prev_trailer.children[-1])
          _SetStronglyConnected(cur_trailer.children[0])
          prev_trailer_idx = cur_trailer_idx
        else:
          break

    # We don't want to split before the last ')' of a function call. This also
    # takes care of the special case of:
    #   atom tr1 tr2 ... trn
    # where the 'tr#' are trailers that may end in a ')'.
    for trailer in node.children[1:]:
      if pytree_utils.NodeName(trailer) != 'trailer':
        break
      if trailer.children[0].value in '([':
        if len(trailer.children) > 2:
          subtypes = pytree_utils.GetNodeAnnotation(
              trailer.children[0], pytree_utils.Annotation.SUBTYPE)
          if subtypes and format_token.Subtype.SUBSCRIPT_BRACKET in subtypes:
            _SetStronglyConnected(
                pytree_utils.FirstLeafNode(trailer.children[1]))

          last_child_node = pytree_utils.LastLeafNode(trailer)
          if last_child_node.value.strip().startswith('#'):
            last_child_node = last_child_node.prev_sibling
          if not style.Get('DEDENT_CLOSING_BRACKETS'):
            last = pytree_utils.LastLeafNode(last_child_node.prev_sibling)
            if last.value != ',':
              if last_child_node.value == ']':
                _SetUnbreakable(last_child_node)
              else:
                _SetSplitPenalty(last_child_node, VERY_STRONGLY_CONNECTED)
        else:
          # If the trailer's children are '()', then make it a strongly
          # connected region.  It's sometimes necessary, though undesirable, to
          # split the two.
          _SetStronglyConnected(trailer.children[-1])

    # If the original source has a "builder" style calls, then we should allow
    # the reformatter to retain that.
    _AllowBuilderStyleCalls(node)

  def Visit_subscript(self, node):  # pylint: disable=invalid-name
    # subscript ::= test | [test] ':' [test] [sliceop]
    _SetStronglyConnected(*node.children)
    self.DefaultNodeVisit(node)

  def Visit_comp_for(self, node):  # pylint: disable=invalid-name
    # comp_for ::= 'for' exprlist 'in' testlist_safe [comp_iter]
    _SetSplitPenalty(pytree_utils.FirstLeafNode(node), 0)
    _SetStronglyConnected(*node.children[1:])
    self.DefaultNodeVisit(node)

  def Visit_comp_if(self, node):  # pylint: disable=invalid-name
    # comp_if ::= 'if' old_test [comp_iter]
    _SetSplitPenalty(node.children[0],
                     style.Get('SPLIT_PENALTY_BEFORE_IF_EXPR'))
    _SetStronglyConnected(*node.children[1:])
    self.DefaultNodeVisit(node)

  def Visit_or_test(self, node):  # pylint: disable=invalid-name
    # or_test ::= and_test ('or' and_test)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, OR_TEST)
    index = 1
    while index + 1 < len(node.children):
      if style.Get('SPLIT_BEFORE_LOGICAL_OPERATOR'):
        _DecrementSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index]), OR_TEST)
      else:
        _DecrementSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index + 1]), OR_TEST)
      index += 2

  def Visit_and_test(self, node):  # pylint: disable=invalid-name
    # and_test ::= not_test ('and' not_test)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, AND_TEST)
    index = 1
    while index + 1 < len(node.children):
      if style.Get('SPLIT_BEFORE_LOGICAL_OPERATOR'):
        _DecrementSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index]), AND_TEST)
      else:
        _DecrementSplitPenalty(
            pytree_utils.FirstLeafNode(node.children[index + 1]), AND_TEST)
      index += 2

  def Visit_not_test(self, node):  # pylint: disable=invalid-name
    # not_test ::= 'not' not_test | comparison
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, NOT_TEST)

  def Visit_comparison(self, node):  # pylint: disable=invalid-name
    # comparison ::= expr (comp_op expr)*
    self.DefaultNodeVisit(node)
    if len(node.children) == 3 and _StronglyConnectedCompOp(node):
      _SetSplitPenalty(
          pytree_utils.FirstLeafNode(node.children[1]), STRONGLY_CONNECTED)
      _SetSplitPenalty(
          pytree_utils.FirstLeafNode(node.children[2]), STRONGLY_CONNECTED)
    else:
      _IncreasePenalty(node, COMPARISON)

  def Visit_star_expr(self, node):  # pylint: disable=invalid-name
    # star_expr ::= '*' expr
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, STAR_EXPR)

  def Visit_expr(self, node):  # pylint: disable=invalid-name
    # expr ::= xor_expr ('|' xor_expr)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, EXPR)
    index = 1
    while index < len(node.children) - 1:
      child = node.children[index]
      if isinstance(child, pytree.Leaf) and child.value == '|':
        if style.Get('SPLIT_BEFORE_BITWISE_OPERATOR'):
          _SetSplitPenalty(child, style.Get('SPLIT_PENALTY_BITWISE_OPERATOR'))
        else:
          _SetSplitPenalty(
              pytree_utils.FirstLeafNode(node.children[index + 1]),
              style.Get('SPLIT_PENALTY_BITWISE_OPERATOR'))
      index += 1

  def Visit_xor_expr(self, node):  # pylint: disable=invalid-name
    # xor_expr ::= and_expr ('^' and_expr)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, XOR_EXPR)

  def Visit_and_expr(self, node):  # pylint: disable=invalid-name
    # and_expr ::= shift_expr ('&' shift_expr)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, AND_EXPR)

  def Visit_shift_expr(self, node):  # pylint: disable=invalid-name
    # shift_expr ::= arith_expr (('<<'|'>>') arith_expr)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, SHIFT_EXPR)

  _ARITH_OPS = frozenset({'PLUS', 'MINUS'})

  def Visit_arith_expr(self, node):  # pylint: disable=invalid-name
    # arith_expr ::= term (('+'|'-') term)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, ARITH_EXPR)

    index = 1
    while index < len(node.children) - 1:
      child = node.children[index]
      if pytree_utils.NodeName(child) in self._ARITH_OPS:
        next_node = pytree_utils.FirstLeafNode(node.children[index + 1])
        _SetSplitPenalty(next_node, ARITH_EXPR)
      index += 1

  _TERM_OPS = frozenset({'STAR', 'AT', 'SLASH', 'PERCENT', 'DOUBLESLASH'})

  def Visit_term(self, node):  # pylint: disable=invalid-name
    # term ::= factor (('*'|'@'|'/'|'%'|'//') factor)*
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, TERM)

    index = 1
    while index < len(node.children) - 1:
      child = node.children[index]
      if pytree_utils.NodeName(child) in self._TERM_OPS:
        next_node = pytree_utils.FirstLeafNode(node.children[index + 1])
        _SetSplitPenalty(next_node, TERM)
      index += 1

  def Visit_factor(self, node):  # pyline: disable=invalid-name
    # factor ::= ('+'|'-'|'~') factor | power
    self.DefaultNodeVisit(node)
    _IncreasePenalty(node, FACTOR)

  def Visit_atom(self, node):  # pylint: disable=invalid-name
    # atom ::= ('(' [yield_expr|testlist_gexp] ')'
    #           '[' [listmaker] ']' |
    #           '{' [dictsetmaker] '}')
    self.DefaultNodeVisit(node)
    if node.children[0].value == '(':
      if node.children[-1].value == ')':
        if pytree_utils.NodeName(node.parent) == 'if_stmt':
          _SetSplitPenalty(node.children[-1], STRONGLY_CONNECTED)
        else:
          if len(node.children) > 2:
            _SetSplitPenalty(pytree_utils.FirstLeafNode(node.children[1]), EXPR)
          _SetSplitPenalty(node.children[-1], ATOM)
    elif node.children[0].value in '[{' and len(node.children) == 2:
      # Keep empty containers together if we can.
      _SetUnbreakable(node.children[-1])

  def Visit_testlist_gexp(self, node):  # pylint: disable=invalid-name
    self.DefaultNodeVisit(node)
    prev_was_comma = False
    for child in node.children:
      if isinstance(child, pytree.Leaf) and child.value == ',':
        _SetUnbreakable(child)
        prev_was_comma = True
      else:
        if prev_was_comma:
          _SetSplitPenalty(pytree_utils.FirstLeafNode(child), 0)
        prev_was_comma = False

  ############################################################################
  # Helper methods that set the annotations.

  def _SetUnbreakableOnChildren(self, node):
    """Set an UNBREAKABLE penalty annotation on children of node."""
    for child in node.children:
      self.Visit(child)
    start = 2 if hasattr(node.children[0], 'is_pseudo') else 1
    for i in py3compat.range(start, len(node.children)):
      _SetUnbreakable(node.children[i])


def _SetUnbreakable(node):
  """Set an UNBREAKABLE penalty annotation for the given node."""
  _RecAnnotate(node, pytree_utils.Annotation.SPLIT_PENALTY, UNBREAKABLE)


def _SetStronglyConnected(*nodes):
  """Set a STRONGLY_CONNECTED penalty annotation for the given nodes."""
  for node in nodes:
    _RecAnnotate(node, pytree_utils.Annotation.SPLIT_PENALTY,
                 STRONGLY_CONNECTED)


def _SetVeryStronglyConnected(*nodes):
  """Set a VERY_STRONGLY_CONNECTED penalty annotation for the given nodes."""
  for node in nodes:
    _RecAnnotate(node, pytree_utils.Annotation.SPLIT_PENALTY,
                 VERY_STRONGLY_CONNECTED)


def _SetExpressionPenalty(node, penalty):
  """Set a penalty annotation on children nodes."""

  def RecExpression(node, first_child_leaf):
    if node is first_child_leaf:
      return

    if isinstance(node, pytree.Leaf):
      if node.value in {'(', 'for', 'if'}:
        return
      penalty_annotation = pytree_utils.GetNodeAnnotation(
          node, pytree_utils.Annotation.SPLIT_PENALTY, default=0)
      if penalty_annotation < penalty:
        _SetSplitPenalty(node, penalty)
    else:
      for child in node.children:
        RecExpression(child, first_child_leaf)

  RecExpression(node, pytree_utils.FirstLeafNode(node))


def _IncreasePenalty(node, amt):
  """Increase a penalty annotation on children nodes."""

  def RecExpression(node, first_child_leaf):
    if node is first_child_leaf:
      return

    if isinstance(node, pytree.Leaf):
      if node.value in {'(', 'for', 'if'}:
        return
      penalty = pytree_utils.GetNodeAnnotation(
          node, pytree_utils.Annotation.SPLIT_PENALTY, default=0)
      _SetSplitPenalty(node, penalty + amt)
    else:
      for child in node.children:
        RecExpression(child, first_child_leaf)

  RecExpression(node, pytree_utils.FirstLeafNode(node))


def _RecAnnotate(tree, annotate_name, annotate_value):
  """Recursively set the given annotation on all leafs of the subtree.

  Takes care to only increase the penalty. If the node already has a higher
  or equal penalty associated with it, this is a no-op.

  Args:
    tree: subtree to annotate
    annotate_name: name of the annotation to set
    annotate_value: value of the annotation to set
  """
  for child in tree.children:
    _RecAnnotate(child, annotate_name, annotate_value)
  if isinstance(tree, pytree.Leaf):
    cur_annotate = pytree_utils.GetNodeAnnotation(
        tree, annotate_name, default=0)
    if cur_annotate < annotate_value:
      pytree_utils.SetNodeAnnotation(tree, annotate_name, annotate_value)


def _StronglyConnectedCompOp(op):
  if (len(op.children[1].children) == 2 and
      pytree_utils.NodeName(op.children[1]) == 'comp_op' and
      pytree_utils.FirstLeafNode(op.children[1]).value == 'not' and
      pytree_utils.LastLeafNode(op.children[1]).value == 'in'):
    return True
  if (isinstance(op.children[1], pytree.Leaf) and
      op.children[1].value in {'==', 'in'}):
    return True
  return False


def _DecrementSplitPenalty(node, amt):
  penalty = pytree_utils.GetNodeAnnotation(
      node, pytree_utils.Annotation.SPLIT_PENALTY, default=amt)
  penalty = penalty - amt if amt < penalty else 0
  _SetSplitPenalty(node, penalty)


def _SetSplitPenalty(node, penalty):
  pytree_utils.SetNodeAnnotation(node, pytree_utils.Annotation.SPLIT_PENALTY,
                                 penalty)


def _AllowBuilderStyleCalls(node):
  """Allow splitting before '.' if it's a builder style function call."""

  def RecGetLeaves(node):
    if isinstance(node, pytree.Leaf):
      return [node]
    children = []
    for child in node.children:
      children += RecGetLeaves(child)
    return children

  list_of_children = RecGetLeaves(node)
  prev_child = None
  for child in list_of_children:
    if child.value == '.':
      if prev_child.lineno != child.lineno:
        _SetSplitPenalty(child, 0)
    prev_child = child
