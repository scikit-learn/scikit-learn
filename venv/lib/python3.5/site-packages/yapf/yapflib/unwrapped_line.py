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
"""UnwrappedLine primitive for formatting.

An unwrapped line is the containing data structure produced by the parser. It
collects all nodes (stored in FormatToken objects) that could appear on a
single line if there were no line length restrictions. It's then used by the
parser to perform the wrapping required to comply with the style guide.
"""

from yapf.yapflib import format_token
from yapf.yapflib import py3compat
from yapf.yapflib import pytree_utils
from yapf.yapflib import split_penalty
from yapf.yapflib import style


class UnwrappedLine(object):
  """Represents a single unwrapped line in the output.

  Attributes:
    depth: indentation depth of this line. This is just a numeric value used to
      distinguish lines that are more deeply nested than others. It is not the
      actual amount of spaces, which is style-dependent.
  """

  def __init__(self, depth, tokens=None):
    """Constructor.

    Creates a new unwrapped line with the given depth an initial list of tokens.
    Constructs the doubly-linked lists for format tokens using their built-in
    next_token and previous_token attributes.

    Arguments:
      depth: indentation depth of this line
      tokens: initial list of tokens
    """
    self.depth = depth
    self._tokens = tokens or []
    self.disable = False

    if self._tokens:
      # Set up a doubly linked list.
      for index, tok in enumerate(self._tokens[1:]):
        # Note, 'index' is the index to the previous token.
        tok.previous_token = self._tokens[index]
        self._tokens[index].next_token = tok

  def CalculateFormattingInformation(self):
    """Calculate the split penalty and total length for the tokens."""
    # Say that the first token in the line should have a space before it. This
    # means only that if this unwrapped line is joined with a predecessor line,
    # then there will be a space between them.
    self.first.spaces_required_before = 1
    self.first.total_length = len(self.first.value)

    prev_token = self.first
    prev_length = self.first.total_length
    for token in self._tokens[1:]:
      if (token.spaces_required_before == 0 and
          _SpaceRequiredBetween(prev_token, token)):
        token.spaces_required_before = 1

      tok_len = len(token.value) if not token.is_pseudo_paren else 0
      token.total_length = prev_length + tok_len + token.spaces_required_before

      # The split penalty has to be computed before {must|can}_break_before,
      # because these may use it for their decision.
      token.split_penalty += _SplitPenalty(prev_token, token)
      token.must_break_before = _MustBreakBefore(prev_token, token)
      token.can_break_before = (
          token.must_break_before or _CanBreakBefore(prev_token, token))

      prev_length = token.total_length
      prev_token = token

  def Split(self):
    """Split the line at semicolons."""
    if not self.has_semicolon or self.disable:
      return [self]

    uwlines = []
    uwline = UnwrappedLine(self.depth)
    for tok in self._tokens:
      if tok.value == ';':
        uwlines.append(uwline)
        uwline = UnwrappedLine(self.depth)
      else:
        uwline.AppendToken(tok)

    if uwline.tokens:
      uwlines.append(uwline)

    for uwline in uwlines:
      pytree_utils.SetNodeAnnotation(uwline.first.node,
                                     pytree_utils.Annotation.MUST_SPLIT, True)
      uwline.first.previous_token = None
      uwline.last.next_token = None

    return uwlines

  ############################################################################
  # Token Access and Manipulation Methods                                    #
  ############################################################################

  def AppendToken(self, token):
    """Append a new FormatToken to the tokens contained in this line."""
    if self._tokens:
      token.previous_token = self.last
      self.last.next_token = token
    self._tokens.append(token)

  def AppendNode(self, node):
    """Convenience method to append a pytree node directly.

    Wraps the node with a FormatToken.

    Arguments:
      node: the node to append
    """
    self.AppendToken(format_token.FormatToken(node))

  @property
  def first(self):
    """Returns the first non-whitespace token."""
    return self._tokens[0]

  @property
  def last(self):
    """Returns the last non-whitespace token."""
    return self._tokens[-1]

  ############################################################################
  # Token -> String Methods                                                  #
  ############################################################################

  def AsCode(self, indent_per_depth=2):
    """Return a "code" representation of this line.

    The code representation shows how the line would be printed out as code.

    TODO(eliben): for now this is rudimentary for debugging - once we add
    formatting capabilities, this method will have other uses (not all tokens
    have spaces around them, for example).

    Arguments:
      indent_per_depth: how much spaces to indend per depth level.

    Returns:
      A string representing the line as code.
    """
    indent = ' ' * indent_per_depth * self.depth
    tokens_str = ' '.join(tok.value for tok in self._tokens)
    return indent + tokens_str

  def __str__(self):  # pragma: no cover
    return self.AsCode()

  def __repr__(self):  # pragma: no cover
    tokens_repr = ','.join(
        ['{0}({1!r})'.format(tok.name, tok.value) for tok in self._tokens])
    return 'UnwrappedLine(depth={0}, tokens=[{1}])'.format(
        self.depth, tokens_repr)

  ############################################################################
  # Properties                                                               #
  ############################################################################

  @property
  def tokens(self):
    """Access the tokens contained within this line.

    The caller must not modify the tokens list returned by this method.

    Returns:
      List of tokens in this line.
    """
    return self._tokens

  @property
  def lineno(self):
    """Return the line number of this unwrapped line.

    Returns:
      The line number of the first token in this unwrapped line.
    """
    return self.first.lineno

  @property
  def is_comment(self):
    return self.first.is_comment

  @property
  def has_semicolon(self):
    return any(tok.value == ';' for tok in self._tokens)


def _IsIdNumberStringToken(tok):
  return tok.is_keyword or tok.is_name or tok.is_number or tok.is_string


def _IsUnaryOperator(tok):
  return format_token.Subtype.UNARY_OPERATOR in tok.subtypes


def _SpaceRequiredBetween(left, right):
  """Return True if a space is required between the left and right token."""
  lval = left.value
  rval = right.value
  if (left.is_pseudo_paren and _IsIdNumberStringToken(right) and
      left.previous_token and _IsIdNumberStringToken(left.previous_token)):
    # Space between keyword... tokens and pseudo parens.
    return True
  if left.is_pseudo_paren or right.is_pseudo_paren:
    # There should be a space after the ':' in a dictionary.
    if left.OpensScope():
      return True
    # The closing pseudo-paren shouldn't affect spacing.
    return False
  if left.is_continuation or right.is_continuation:
    # The continuation node's value has all of the spaces it needs.
    return False
  if right.name in pytree_utils.NONSEMANTIC_TOKENS:
    # No space before a non-semantic token.
    return False
  if _IsIdNumberStringToken(left) and _IsIdNumberStringToken(right):
    # Spaces between keyword, string, number, and identifier tokens.
    return True
  if lval == ',' and rval == ':':
    # We do want a space between a comma and colon.
    return True
  if rval in ':,':
    # Otherwise, we never want a space before a colon or comma.
    return False
  if lval == ',' and rval in ']})':
    # Add a space between ending ',' and closing bracket if requested.
    return style.Get('SPACE_BETWEEN_ENDING_COMMA_AND_CLOSING_BRACKET')
  if lval == ',':
    # We want a space after a comma.
    return True
  if lval == 'from' and rval == '.':
    # Space before the '.' in an import statement.
    return True
  if lval == '.' and rval == 'import':
    # Space after the '.' in an import statement.
    return True
  if (lval == '=' and rval == '.' and
      format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN not in left.subtypes):
    # Space between equal and '.' as in "X = ...".
    return True
  if ((right.is_keyword or right.is_name) and
      (left.is_keyword or left.is_name)):
    # Don't merge two keywords/identifiers.
    return True
  if (format_token.Subtype.SUBSCRIPT_COLON in left.subtypes or
      format_token.Subtype.SUBSCRIPT_COLON in right.subtypes):
    # A subscript shouldn't have spaces separating its colons.
    return False
  if (format_token.Subtype.TYPED_NAME in left.subtypes or
      format_token.Subtype.TYPED_NAME in right.subtypes):
    # A typed argument should have a space after the colon.
    return True
  if left.is_string:
    if (rval == '=' and
        format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST in right.subtypes
       ):
      # If there is a type hint, then we don't want to add a space between the
      # equal sign and the hint.
      return False
    if rval not in '[)]}.':
      # A string followed by something other than a subscript, closing bracket,
      # or dot should have a space after it.
      return True
  if left.is_binary_op and lval != '**' and _IsUnaryOperator(right):
    # Space between the binary operator and the unary operator.
    return True
  if left.is_keyword and _IsUnaryOperator(right):
    # Handle things like "not -3 < x".
    return True
  if _IsUnaryOperator(left) and _IsUnaryOperator(right):
    # No space between two unary operators.
    return False
  if left.is_binary_op or right.is_binary_op:
    if lval == '**' or rval == '**':
      # Space around the "power" operator.
      return style.Get('SPACES_AROUND_POWER_OPERATOR')
    # Enforce spaces around binary operators except the blacklisted ones.
    blacklist = style.Get('NO_SPACES_AROUND_SELECTED_BINARY_OPERATORS')
    return lval not in blacklist and rval not in blacklist
  if (_IsUnaryOperator(left) and lval != 'not' and
      (right.is_name or right.is_number or rval == '(')):
    # The previous token was a unary op. No space is desired between it and
    # the current token.
    return False
  if (format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN in left.subtypes or
      format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN in right.subtypes):
    # A named argument or default parameter shouldn't have spaces around it.
    return style.Get('SPACES_AROUND_DEFAULT_OR_NAMED_ASSIGN')
  if (format_token.Subtype.VARARGS_LIST in left.subtypes or
      format_token.Subtype.VARARGS_LIST in right.subtypes):
    return False
  if (format_token.Subtype.VARARGS_STAR in left.subtypes or
      format_token.Subtype.KWARGS_STAR_STAR in left.subtypes):
    # Don't add a space after a vararg's star or a keyword's star-star.
    return False
  if lval == '@' and format_token.Subtype.DECORATOR in left.subtypes:
    # Decorators shouldn't be separated from the 'at' sign.
    return False
  if left.is_keyword and rval == '.' or lval == '.' and right.is_keyword:
    # Add space between keywords and dots.
    return lval != 'None'
  if lval == '.' or rval == '.':
    # Don't place spaces between dots.
    return False
  if ((lval == '(' and rval == ')') or (lval == '[' and rval == ']') or
      (lval == '{' and rval == '}')):
    # Empty objects shouldn't be separated by spaces.
    return False
  if (lval in pytree_utils.OPENING_BRACKETS and
      rval in pytree_utils.OPENING_BRACKETS):
    # Nested objects' opening brackets shouldn't be separated.
    return False
  if (lval in pytree_utils.CLOSING_BRACKETS and
      rval in pytree_utils.CLOSING_BRACKETS):
    # Nested objects' closing brackets shouldn't be separated.
    return False
  if lval in pytree_utils.CLOSING_BRACKETS and rval in '([':
    # A call, set, dictionary, or subscript that has a call or subscript after
    # it shouldn't have a space between them.
    return False
  if lval in pytree_utils.OPENING_BRACKETS and _IsIdNumberStringToken(right):
    # Don't separate the opening bracket from the first item.
    return False
  if left.is_name and rval in '([':
    # Don't separate a call or array access from the name.
    return False
  if rval in pytree_utils.CLOSING_BRACKETS:
    # Don't separate the closing bracket from the last item.
    # FIXME(morbo): This might be too permissive.
    return False
  if lval == 'print' and rval == '(':
    # Special support for the 'print' function.
    return False
  if lval in pytree_utils.OPENING_BRACKETS and _IsUnaryOperator(right):
    # Don't separate a unary operator from the opening bracket.
    return False
  if (lval in pytree_utils.OPENING_BRACKETS and
      (format_token.Subtype.VARARGS_STAR in right.subtypes or
       format_token.Subtype.KWARGS_STAR_STAR in right.subtypes)):
    # Don't separate a '*' or '**' from the opening bracket.
    return False
  if rval == ';':
    # Avoid spaces before a semicolon. (Why is there a semicolon?!)
    return False
  if lval == '(' and rval == 'await':
    # Special support for the 'await' keyword. Don't separate the 'await'
    # keyword from an opening paren.
    return False
  return True


def _MustBreakBefore(prev_token, cur_token):
  """Return True if a line break is required before the current token."""
  if prev_token.is_comment or (prev_token.previous_token and
                               prev_token.is_pseudo_paren and
                               prev_token.previous_token.is_comment):
    # Must break if the previous token was a comment.
    return True
  if (cur_token.is_string and prev_token.is_string and
      IsSurroundedByBrackets(cur_token)):
    # We want consecutive strings to be on separate lines. This is a
    # reasonable assumption, because otherwise they should have written them
    # all on the same line, or with a '+'.
    return True
  return pytree_utils.GetNodeAnnotation(
      cur_token.node, pytree_utils.Annotation.MUST_SPLIT, default=False)


def _CanBreakBefore(prev_token, cur_token):
  """Return True if a line break may occur before the current token."""
  pval = prev_token.value
  cval = cur_token.value
  if py3compat.PY3:
    if pval == 'yield' and cval == 'from':
      # Don't break before a yield argument.
      return False
    if pval in {'async', 'await'} and cval in {'def', 'with', 'for'}:
      # Don't break after sync keywords.
      return False
  if cur_token.split_penalty >= split_penalty.UNBREAKABLE:
    return False
  if pval == '@':
    # Don't break right after the beginning of a decorator.
    return False
  if cval == ':':
    # Don't break before the start of a block of code.
    return False
  if cval == ',':
    # Don't break before a comma.
    return False
  if prev_token.is_name and cval == '(':
    # Don't break in the middle of a function definition or call.
    return False
  if prev_token.is_name and cval == '[':
    # Don't break in the middle of an array dereference.
    return False
  if prev_token.is_name and cval == '.':
    # Don't break before the '.' in a dotted name.
    return False
  if cur_token.is_comment and prev_token.lineno == cur_token.lineno:
    # Don't break a comment at the end of the line.
    return False
  if format_token.Subtype.UNARY_OPERATOR in prev_token.subtypes:
    # Don't break after a unary token.
    return False
  return True


def IsSurroundedByBrackets(tok):
  """Return True if the token is surrounded by brackets."""
  paren_count = 0
  brace_count = 0
  sq_bracket_count = 0
  previous_token = tok.previous_token
  while previous_token:
    if previous_token.value == ')':
      paren_count -= 1
    elif previous_token.value == '}':
      brace_count -= 1
    elif previous_token.value == ']':
      sq_bracket_count -= 1

    if previous_token.value == '(':
      if paren_count == 0:
        return previous_token
      paren_count += 1
    elif previous_token.value == '{':
      if brace_count == 0:
        return previous_token
      brace_count += 1
    elif previous_token.value == '[':
      if sq_bracket_count == 0:
        return previous_token
      sq_bracket_count += 1

    previous_token = previous_token.previous_token
  return None


_LOGICAL_OPERATORS = frozenset({'and', 'or'})
_BITWISE_OPERATORS = frozenset({'&', '|', '^'})
_TERM_OPERATORS = frozenset({'*', '/', '%', '//'})


def _SplitPenalty(prev_token, cur_token):
  """Return the penalty for breaking the line before the current token."""
  pval = prev_token.value
  cval = cur_token.value
  if pval == 'not':
    return split_penalty.UNBREAKABLE

  if cur_token.node_split_penalty > 0:
    return cur_token.node_split_penalty

  if style.Get('SPLIT_BEFORE_LOGICAL_OPERATOR'):
    # Prefer to split before 'and' and 'or'.
    if pval in _LOGICAL_OPERATORS:
      return style.Get('SPLIT_PENALTY_LOGICAL_OPERATOR')
    if cval in _LOGICAL_OPERATORS:
      return 0
  else:
    # Prefer to split after 'and' and 'or'.
    if pval in _LOGICAL_OPERATORS:
      return 0
    if cval in _LOGICAL_OPERATORS:
      return style.Get('SPLIT_PENALTY_LOGICAL_OPERATOR')

  if style.Get('SPLIT_BEFORE_BITWISE_OPERATOR'):
    # Prefer to split before '&', '|', and '^'.
    if pval in _BITWISE_OPERATORS:
      return style.Get('SPLIT_PENALTY_BITWISE_OPERATOR')
    if cval in _BITWISE_OPERATORS:
      return 0
  else:
    # Prefer to split after '&', '|', and '^'.
    if pval in _BITWISE_OPERATORS:
      return 0
    if cval in _BITWISE_OPERATORS:
      return style.Get('SPLIT_PENALTY_BITWISE_OPERATOR')

  if (format_token.Subtype.COMP_FOR in cur_token.subtypes or
      format_token.Subtype.COMP_IF in cur_token.subtypes):
    # We don't mind breaking before the 'for' or 'if' of a list comprehension.
    return 0
  if format_token.Subtype.UNARY_OPERATOR in prev_token.subtypes:
    # Try not to break after a unary operator.
    return style.Get('SPLIT_PENALTY_AFTER_UNARY_OPERATOR')
  if pval == ',':
    # Breaking after a comma is fine, if need be.
    return 0
  if prev_token.is_binary_op:
    # We would rather not split after an equality operator.
    return 20
  if (format_token.Subtype.VARARGS_STAR in prev_token.subtypes or
      format_token.Subtype.KWARGS_STAR_STAR in prev_token.subtypes):
    # Don't split after a varargs * or kwargs **.
    return split_penalty.UNBREAKABLE
  if prev_token.OpensScope() and cval != '(':
    # Slightly prefer
    return style.Get('SPLIT_PENALTY_AFTER_OPENING_BRACKET')
  if cval == ':':
    # Don't split before a colon.
    return split_penalty.UNBREAKABLE
  if cval == '=':
    # Don't split before an assignment.
    return split_penalty.UNBREAKABLE
  if (format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN in prev_token.subtypes or
      format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN in cur_token.subtypes):
    # Don't break before or after an default or named assignment.
    return split_penalty.UNBREAKABLE
  if cval == '==':
    # We would rather not split before an equality operator.
    return split_penalty.STRONGLY_CONNECTED
  if cur_token.ClosesScope():
    # Give a slight penalty for splitting before the closing scope.
    return 100
  if pval in _TERM_OPERATORS or cval in _TERM_OPERATORS:
    return 50
  return 0
