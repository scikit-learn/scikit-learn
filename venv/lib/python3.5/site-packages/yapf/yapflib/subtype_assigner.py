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
"""Subtype assigner for lib2to3 trees.

This module assigns extra type information to the lib2to3 trees. This
information is more specific than whether something is an operator or an
identifier. For instance, it can specify if a node in the tree is part of a
subscript.

  AssignSubtypes(): the main function exported by this module.

Annotations:
  subtype: The subtype of a pytree token. See 'format_token' module for a list
      of subtypes.
"""

from lib2to3 import pytree
from lib2to3.pgen2 import token
from lib2to3.pygram import python_symbols as syms

from yapf.yapflib import format_token
from yapf.yapflib import pytree_utils
from yapf.yapflib import pytree_visitor
from yapf.yapflib import style


def AssignSubtypes(tree):
  """Run the subtype assigner visitor over the tree, modifying it in place.

  Arguments:
    tree: the top-level pytree node to annotate with subtypes.
  """
  subtype_assigner = _SubtypeAssigner()
  subtype_assigner.Visit(tree)


# Map tokens in argument lists to their respective subtype.
_ARGLIST_TOKEN_TO_SUBTYPE = {
    '=': format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN,
    ':': format_token.Subtype.TYPED_NAME,
    '*': format_token.Subtype.VARARGS_STAR,
    '**': format_token.Subtype.KWARGS_STAR_STAR,
}


class _SubtypeAssigner(pytree_visitor.PyTreeVisitor):
  """_SubtypeAssigner - see file-level docstring for detailed description.

  The subtype is added as an annotation to the pytree token.
  """

  def Visit_dictsetmaker(self, node):  # pylint: disable=invalid-name
    # dictsetmaker ::= (test ':' test (comp_for |
    #                                   (',' test ':' test)* [','])) |
    #                  (test (comp_for | (',' test)* [',']))
    for child in node.children:
      self.Visit(child)

    comp_for = False
    dict_maker = False

    for child in node.children:
      if pytree_utils.NodeName(child) == 'comp_for':
        comp_for = True
        _AppendFirstLeafTokenSubtype(child,
                                     format_token.Subtype.DICT_SET_GENERATOR)
      elif pytree_utils.NodeName(child) in ('COLON', 'DOUBLESTAR'):
        dict_maker = True

    if not comp_for and dict_maker:
      last_was_colon = False
      unpacking = False
      for child in node.children:
        if pytree_utils.NodeName(child) == 'DOUBLESTAR':
          _AppendFirstLeafTokenSubtype(child,
                                       format_token.Subtype.KWARGS_STAR_STAR)
        if last_was_colon:
          if style.Get('INDENT_DICTIONARY_VALUE'):
            _InsertPseudoParentheses(child)
          else:
            _AppendFirstLeafTokenSubtype(child,
                                         format_token.Subtype.DICTIONARY_VALUE)
        elif (isinstance(child, pytree.Node) or
              (not child.value.startswith('#') and child.value not in '{:,')):
          # Mark the first leaf of a key entry as a DICTIONARY_KEY. We
          # normally want to split before them if the dictionary cannot exist
          # on a single line.
          if not unpacking or pytree_utils.FirstLeafNode(child).value == '**':
            _AppendFirstLeafTokenSubtype(child,
                                         format_token.Subtype.DICTIONARY_KEY)
          _AppendSubtypeRec(child, format_token.Subtype.DICTIONARY_KEY_PART)
        last_was_colon = pytree_utils.NodeName(child) == 'COLON'
        if pytree_utils.NodeName(child) == 'DOUBLESTAR':
          unpacking = True
        elif last_was_colon:
          unpacking = False

  def Visit_expr_stmt(self, node):  # pylint: disable=invalid-name
    # expr_stmt ::= testlist_star_expr (augassign (yield_expr|testlist)
    #               | ('=' (yield_expr|testlist_star_expr))*)
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == '=':
        _AppendTokenSubtype(child, format_token.Subtype.ASSIGN_OPERATOR)

  def Visit_or_test(self, node):  # pylint: disable=invalid-name
    # or_test ::= and_test ('or' and_test)*
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == 'or':
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_and_test(self, node):  # pylint: disable=invalid-name
    # and_test ::= not_test ('and' not_test)*
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == 'and':
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_not_test(self, node):  # pylint: disable=invalid-name
    # not_test ::= 'not' not_test | comparison
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == 'not':
        _AppendTokenSubtype(child, format_token.Subtype.UNARY_OPERATOR)

  def Visit_comparison(self, node):  # pylint: disable=invalid-name
    # comparison ::= expr (comp_op expr)*
    # comp_op ::= '<'|'>'|'=='|'>='|'<='|'<>'|'!='|'in'|'not in'|'is'|'is not'
    for child in node.children:
      self.Visit(child)
      if (isinstance(child, pytree.Leaf) and
          child.value in {'<', '>', '==', '>=', '<=', '<>', '!=', 'in', 'is'}):
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)
      elif pytree_utils.NodeName(child) == 'comp_op':
        for grandchild in child.children:
          _AppendTokenSubtype(grandchild, format_token.Subtype.BINARY_OPERATOR)

  def Visit_star_expr(self, node):  # pylint: disable=invalid-name
    # star_expr ::= '*' expr
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == '*':
        _AppendTokenSubtype(child, format_token.Subtype.UNARY_OPERATOR)
        _AppendTokenSubtype(child, format_token.Subtype.VARARGS_STAR)

  def Visit_expr(self, node):  # pylint: disable=invalid-name
    # expr ::= xor_expr ('|' xor_expr)*
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == '|':
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_xor_expr(self, node):  # pylint: disable=invalid-name
    # xor_expr ::= and_expr ('^' and_expr)*
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == '^':
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_and_expr(self, node):  # pylint: disable=invalid-name
    # and_expr ::= shift_expr ('&' shift_expr)*
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == '&':
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_shift_expr(self, node):  # pylint: disable=invalid-name
    # shift_expr ::= arith_expr (('<<'|'>>') arith_expr)*
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value in {'<<', '>>'}:
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_arith_expr(self, node):  # pylint: disable=invalid-name
    # arith_expr ::= term (('+'|'-') term)*
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value in '+-':
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_term(self, node):  # pylint: disable=invalid-name
    # term ::= factor (('*'|'/'|'%'|'//') factor)*
    for child in node.children:
      self.Visit(child)
      if (isinstance(child, pytree.Leaf) and
          child.value in {'*', '/', '%', '//'}):
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_factor(self, node):  # pylint: disable=invalid-name
    # factor ::= ('+'|'-'|'~') factor | power
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value in '+-~':
        _AppendTokenSubtype(child, format_token.Subtype.UNARY_OPERATOR)

  def Visit_power(self, node):  # pylint: disable=invalid-name
    # power ::= atom trailer* ['**' factor]
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == '**':
        _AppendTokenSubtype(child, format_token.Subtype.BINARY_OPERATOR)

  def Visit_trailer(self, node):  # pylint: disable=invalid-name
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value in '[]':
        _AppendTokenSubtype(child, format_token.Subtype.SUBSCRIPT_BRACKET)

  def Visit_subscript(self, node):  # pylint: disable=invalid-name
    # subscript ::= test | [test] ':' [test] [sliceop]
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == ':':
        _AppendTokenSubtype(child, format_token.Subtype.SUBSCRIPT_COLON)

  def Visit_sliceop(self, node):  # pylint: disable=invalid-name
    # sliceop ::= ':' [test]
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == ':':
        _AppendTokenSubtype(child, format_token.Subtype.SUBSCRIPT_COLON)

  def Visit_argument(self, node):  # pylint: disable=invalid-name
    # argument ::=
    #     test [comp_for] | test '=' test
    self._ProcessArgLists(node)

  def Visit_arglist(self, node):  # pylint: disable=invalid-name
    # arglist ::=
    #     (argument ',')* (argument [',']
    #                     | '*' test (',' argument)* [',' '**' test]
    #                     | '**' test)
    self._ProcessArgLists(node)
    _SetArgListSubtype(node, format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN,
                       format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST)

  def Visit_tname(self, node):  # pylint: disable=invalid-name
    self._ProcessArgLists(node)
    _SetArgListSubtype(node, format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN,
                       format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST)

  def Visit_decorator(self, node):  # pylint: disable=invalid-name
    # decorator ::=
    #     '@' dotted_name [ '(' [arglist] ')' ] NEWLINE
    for child in node.children:
      if isinstance(child, pytree.Leaf) and child.value == '@':
        _AppendTokenSubtype(child, subtype=format_token.Subtype.DECORATOR)
      self.Visit(child)

  def Visit_funcdef(self, node):  # pylint: disable=invalid-name
    # funcdef ::=
    #     'def' NAME parameters ['->' test] ':' suite
    for child in node.children:
      if pytree_utils.NodeName(child) == 'NAME' and child.value != 'def':
        _AppendTokenSubtype(child, format_token.Subtype.FUNC_DEF)
        break
    for child in node.children:
      self.Visit(child)

  def Visit_typedargslist(self, node):  # pylint: disable=invalid-name
    # typedargslist ::=
    #     ((tfpdef ['=' test] ',')*
    #          ('*' [tname] (',' tname ['=' test])* [',' '**' tname]
    #           | '**' tname)
    #     | tfpdef ['=' test] (',' tfpdef ['=' test])* [','])
    self._ProcessArgLists(node)
    _SetArgListSubtype(node, format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN,
                       format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST)
    tname = False
    for child in node.children:
      if pytree_utils.NodeName(child) == 'tname':
        tname = True
        _SetArgListSubtype(child, format_token.Subtype.TYPED_NAME,
                           format_token.Subtype.TYPED_NAME_ARG_LIST)
      if not isinstance(child, pytree.Leaf):
        continue
      if child.value == ',':
        tname = False
      elif child.value == '=' and tname:
        _AppendTokenSubtype(child, subtype=format_token.Subtype.TYPED_NAME)
        tname = False

  def Visit_varargslist(self, node):  # pylint: disable=invalid-name
    # varargslist ::=
    #     ((vfpdef ['=' test] ',')*
    #          ('*' [vname] (',' vname ['=' test])*  [',' '**' vname]
    #           | '**' vname)
    #      | vfpdef ['=' test] (',' vfpdef ['=' test])* [','])
    self._ProcessArgLists(node)
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf) and child.value == '=':
        _AppendTokenSubtype(child, format_token.Subtype.VARARGS_LIST)

  def Visit_comp_for(self, node):  # pylint: disable=invalid-name
    # comp_for ::= 'for' exprlist 'in' testlist_safe [comp_iter]
    _AppendSubtypeRec(node, format_token.Subtype.COMP_FOR)
    # Mark the previous node as COMP_EXPR unless this is a nested comprehension
    # as these will have the outer comprehension as their previous node.
    attr = pytree_utils.GetNodeAnnotation(node.parent,
                                          pytree_utils.Annotation.SUBTYPE)
    if not attr or format_token.Subtype.COMP_FOR not in attr:
      _AppendSubtypeRec(node.parent.children[0], format_token.Subtype.COMP_EXPR)
    self.DefaultNodeVisit(node)

  def Visit_comp_if(self, node):  # pylint: disable=invalid-name
    # comp_if ::= 'if' old_test [comp_iter]
    _AppendSubtypeRec(node, format_token.Subtype.COMP_IF)
    self.DefaultNodeVisit(node)

  def _ProcessArgLists(self, node):
    """Common method for processing argument lists."""
    for child in node.children:
      self.Visit(child)
      if isinstance(child, pytree.Leaf):
        _AppendTokenSubtype(
            child,
            subtype=_ARGLIST_TOKEN_TO_SUBTYPE.get(child.value,
                                                  format_token.Subtype.NONE))


def _SetArgListSubtype(node, node_subtype, list_subtype):
  """Set named assign subtype on elements in a arg list."""

  def HasSubtype(node):
    """Return True if the arg list has a named assign subtype."""
    if isinstance(node, pytree.Leaf):
      if node_subtype in pytree_utils.GetNodeAnnotation(
          node, pytree_utils.Annotation.SUBTYPE, set()):
        return True
      return False
    has_subtype = False
    for child in node.children:
      if pytree_utils.NodeName(child) != 'arglist':
        has_subtype |= HasSubtype(child)
    return has_subtype

  if HasSubtype(node):
    for child in node.children:
      if pytree_utils.NodeName(child) != 'COMMA':
        _AppendFirstLeafTokenSubtype(child, list_subtype)


def _AppendTokenSubtype(node, subtype):
  """Append the token's subtype only if it's not already set."""
  pytree_utils.AppendNodeAnnotation(node, pytree_utils.Annotation.SUBTYPE,
                                    subtype)


def _AppendFirstLeafTokenSubtype(node, subtype):
  """Append the first leaf token's subtypes."""
  if isinstance(node, pytree.Leaf):
    _AppendTokenSubtype(node, subtype)
    return
  _AppendFirstLeafTokenSubtype(node.children[0], subtype)


def _AppendSubtypeRec(node, subtype, force=True):
  """Append the leafs in the node to the given subtype."""
  if isinstance(node, pytree.Leaf):
    _AppendTokenSubtype(node, subtype)
    return
  for child in node.children:
    _AppendSubtypeRec(child, subtype, force=force)


def _InsertPseudoParentheses(node):
  """Insert pseudo parentheses so that dicts can be formatted correctly."""
  comment_node = None
  if isinstance(node, pytree.Node):
    if node.children[-1].type == token.COMMENT:
      comment_node = node.children[-1].clone()
      node.children[-1].remove()

  first = pytree_utils.FirstLeafNode(node)
  last = pytree_utils.LastLeafNode(node)

  if first == last and first.type == token.COMMENT:
    # A comment was inserted before the value, which is a pytree.Leaf.
    # Encompass the dictionary's value into an ATOM node.
    last = first.next_sibling
    new_node = pytree.Node(syms.atom, [first.clone(), last.clone()])
    node.replace(new_node)
    node = new_node
    last.remove()

    first = pytree_utils.FirstLeafNode(node)
    last = pytree_utils.LastLeafNode(node)

  lparen = pytree.Leaf(
      token.LPAR, u'(', context=('', (first.get_lineno(), first.column - 1)))
  last_lineno = last.get_lineno()
  if last.type == token.STRING and '\n' in last.value:
    last_lineno += last.value.count('\n')

  if last.type == token.STRING and '\n' in last.value:
    last_column = len(last.value.split('\n')[-1]) + 1
  else:
    last_column = last.column + len(last.value) + 1
  rparen = pytree.Leaf(
      token.RPAR, u')', context=('', (last_lineno, last_column)))

  lparen.is_pseudo = True
  rparen.is_pseudo = True

  if isinstance(node, pytree.Node):
    node.insert_child(0, lparen)
    node.append_child(rparen)
    if comment_node:
      node.append_child(comment_node)
    _AppendFirstLeafTokenSubtype(node, format_token.Subtype.DICTIONARY_VALUE)
  else:
    clone = node.clone()
    new_node = pytree.Node(syms.atom, [lparen, clone, rparen])
    node.replace(new_node)
    _AppendFirstLeafTokenSubtype(clone, format_token.Subtype.DICTIONARY_VALUE)
