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
"""Implements a format decision state object that manages whitespace decisions.

Each token is processed one at a time, at which point its whitespace formatting
decisions are made. A graph of potential whitespace formattings is created,
where each node in the graph is a format decision state object. The heuristic
tries formatting the token with and without a newline before it to determine
which one has the least penalty. Therefore, the format decision state object for
each decision needs to be its own unique copy.

Once the heuristic determines the best formatting, it makes a non-dry run pass
through the code to commit the whitespace formatting.

  FormatDecisionState: main class exported by this module.
"""

from yapf.yapflib import format_token
from yapf.yapflib import object_state
from yapf.yapflib import split_penalty
from yapf.yapflib import style
from yapf.yapflib import unwrapped_line


class FormatDecisionState(object):
  """The current state when indenting an unwrapped line.

  The FormatDecisionState object is meant to be copied instead of referenced.

  Attributes:
    first_indent: The indent of the first token.
    column: The number of used columns in the current line.
    next_token: The next token to be formatted.
    paren_level: The level of nesting inside (), [], and {}.
    lowest_level_on_line: The lowest paren_level on the current line.
    newline: Indicates if a newline is added along the edge to this format
      decision state node.
    previous: The previous format decision state in the decision tree.
    stack: A stack (of _ParenState) keeping track of properties applying to
      parenthesis levels.
    comp_stack: A stack (of ComprehensionState) keeping track of properties
      applying to comprehensions.
    ignore_stack_for_comparison: Ignore the stack of _ParenState for state
      comparison.
  """

  def __init__(self, line, first_indent):
    """Initializer.

    Initializes to the state after placing the first token from 'line' at
    'first_indent'.

    Arguments:
      line: (UnwrappedLine) The unwrapped line we're currently processing.
      first_indent: (int) The indent of the first token.
    """
    self.next_token = line.first
    self.column = first_indent
    self.line = line
    self.paren_level = 0
    self.lowest_level_on_line = 0
    self.ignore_stack_for_comparison = False
    self.stack = [_ParenState(first_indent, first_indent)]
    self.comp_stack = []
    self.first_indent = first_indent
    self.newline = False
    self.previous = None
    self.column_limit = style.Get('COLUMN_LIMIT')

  def Clone(self):
    """Clones a FormatDecisionState object."""
    new = FormatDecisionState(self.line, self.first_indent)
    new.next_token = self.next_token
    new.column = self.column
    new.line = self.line
    new.paren_level = self.paren_level
    new.line.depth = self.line.depth
    new.lowest_level_on_line = self.lowest_level_on_line
    new.ignore_stack_for_comparison = self.ignore_stack_for_comparison
    new.first_indent = self.first_indent
    new.newline = self.newline
    new.previous = self.previous
    new.stack = [state.Clone() for state in self.stack]
    new.comp_stack = [state.Clone() for state in self.comp_stack]
    return new

  def __eq__(self, other):
    # Note: 'first_indent' is implicit in the stack. Also, we ignore 'previous',
    # because it shouldn't have a bearing on this comparison. (I.e., it will
    # report equal if 'next_token' does.)
    return (self.next_token == other.next_token and
            self.column == other.column and
            self.paren_level == other.paren_level and
            self.line.depth == other.line.depth and
            self.lowest_level_on_line == other.lowest_level_on_line and
            (self.ignore_stack_for_comparison or
             other.ignore_stack_for_comparison or
             self.stack == other.stack and self.comp_stack == other.comp_stack))

  def __ne__(self, other):
    return not self == other

  def __hash__(self):
    return hash((self.next_token, self.column, self.paren_level,
                 self.line.depth, self.lowest_level_on_line))

  def __repr__(self):
    return ('column::%d, next_token::%s, paren_level::%d, stack::[\n\t%s' %
            (self.column, repr(self.next_token), self.paren_level,
             '\n\t'.join(repr(s) for s in self.stack) + ']'))

  def CanSplit(self, must_split):
    """Determine if we can split before the next token.

    Arguments:
      must_split: (bool) A newline was required before this token.

    Returns:
      True if the line can be split before the next token.
    """
    current = self.next_token
    previous = current.previous_token

    if current.is_pseudo_paren:
      return False

    if (not must_split and
        format_token.Subtype.DICTIONARY_KEY_PART in current.subtypes and
        format_token.Subtype.DICTIONARY_KEY not in current.subtypes and
        not style.Get('ALLOW_MULTILINE_DICTIONARY_KEYS')):
      # In some situations, a dictionary may be multiline, but pylint doesn't
      # like it. So don't allow it unless forced to.
      return False

    if (not must_split and
        format_token.Subtype.DICTIONARY_VALUE in current.subtypes and
        not style.Get('ALLOW_SPLIT_BEFORE_DICT_VALUE')):
      return False

    if previous and previous.value == '(' and current.value == ')':
      # Don't split an empty function call list if we aren't splitting before
      # dict values.
      token = previous.previous_token
      while token:
        prev = token.previous_token
        if not prev or prev.name not in {'NAME', 'DOT'}:
          break
        token = token.previous_token
      if token and format_token.Subtype.DICTIONARY_VALUE in token.subtypes:
        if not style.Get('ALLOW_SPLIT_BEFORE_DICT_VALUE'):
          return False

    if previous and previous.value == '.' and current.value == '.':
      return False

    return current.can_break_before

  def MustSplit(self):
    """Returns True if the line must split before the next token."""
    current = self.next_token
    previous = current.previous_token

    if current.is_pseudo_paren:
      return False

    if current.must_break_before:
      return True

    if not previous:
      return False

    if style.Get('SPLIT_ALL_COMMA_SEPARATED_VALUES') and previous.value == ',':
      return True

    if (self.stack[-1].split_before_closing_bracket and
        current.value in '}]' and style.Get('SPLIT_BEFORE_CLOSING_BRACKET')):
      # Split before the closing bracket if we can.
      return current.node_split_penalty != split_penalty.UNBREAKABLE

    if (current.value == ')' and previous.value == ',' and
        not _IsSingleElementTuple(current.matching_bracket)):
      return True

    # Prevent splitting before the first argument in compound statements
    # with the exception of function declarations.
    if (style.Get('SPLIT_BEFORE_FIRST_ARGUMENT') and
        _IsCompoundStatement(self.line.first) and
        not _IsFunctionDef(self.line.first)):
      return False

    ###########################################################################
    # List Splitting
    if (style.Get('DEDENT_CLOSING_BRACKETS') or
        style.Get('SPLIT_BEFORE_FIRST_ARGUMENT')):
      bracket = current if current.ClosesScope() else previous
      if format_token.Subtype.SUBSCRIPT_BRACKET not in bracket.subtypes:
        if bracket.OpensScope():
          if style.Get('COALESCE_BRACKETS'):
            if current.OpensScope():
              # Prefer to keep all opening brackets together.
              return False

          if (not _IsLastScopeInLine(bracket) or
              unwrapped_line.IsSurroundedByBrackets(bracket)):
            last_token = bracket.matching_bracket
          else:
            last_token = _LastTokenInLine(bracket.matching_bracket)

          if not self._FitsOnLine(bracket, last_token):
            # Split before the first element if the whole list can't fit on a
            # single line.
            self.stack[-1].split_before_closing_bracket = True
            return True

        elif style.Get('DEDENT_CLOSING_BRACKETS') and current.ClosesScope():
          # Split before and dedent the closing bracket.
          return self.stack[-1].split_before_closing_bracket

    if (style.Get('SPLIT_BEFORE_EXPRESSION_AFTER_OPENING_PAREN') and
        current.is_name):
      # An expression that's surrounded by parens gets split after the opening
      # parenthesis.
      def SurroundedByParens(token):
        """Check if it's an expression surrounded by parentheses."""
        while token:
          if token.value == ',':
            return False
          if token.value == ')':
            return not token.next_token
          if token.OpensScope():
            token = token.matching_bracket.next_token
          else:
            token = token.next_token
        return False

      if (previous.value == '(' and not previous.is_pseudo_paren and
          not unwrapped_line.IsSurroundedByBrackets(previous)):
        pptoken = previous.previous_token
        if (pptoken and not pptoken.is_name and not pptoken.is_keyword and
            SurroundedByParens(current)):
          return True

    if (current.is_name or current.is_string) and previous.value == ',':
      # If the list has function calls in it and the full list itself cannot
      # fit on the line, then we want to split. Otherwise, we'll get something
      # like this:
      #
      #     X = [
      #         Bar(xxx='some string',
      #             yyy='another long string',
      #             zzz='a third long string'), Bar(
      #                 xxx='some string',
      #                 yyy='another long string',
      #                 zzz='a third long string')
      #     ]
      #
      # or when a string formatting syntax.
      func_call_or_string_format = False
      tok = current.next_token
      if current.is_name:
        while tok and (tok.is_name or tok.value == '.'):
          tok = tok.next_token
        func_call_or_string_format = tok and tok.value == '('
      elif current.is_string:
        while tok and tok.is_string:
          tok = tok.next_token
        func_call_or_string_format = tok and tok.value == '%'
      if func_call_or_string_format:
        open_bracket = unwrapped_line.IsSurroundedByBrackets(current)
        if open_bracket:
          if open_bracket.value in '[{':
            if not self._FitsOnLine(open_bracket,
                                    open_bracket.matching_bracket):
              return True
          elif tok.value == '(':
            if not self._FitsOnLine(current, tok.matching_bracket):
              return True

    ###########################################################################
    # Dict/Set Splitting
    if (style.Get('EACH_DICT_ENTRY_ON_SEPARATE_LINE') and
        format_token.Subtype.DICTIONARY_KEY in current.subtypes and
        not current.is_comment):
      # Place each dictionary entry onto its own line.
      if previous.value == '{' and previous.previous_token:
        opening = _GetOpeningBracket(previous.previous_token)
        if (opening and opening.value == '(' and opening.previous_token and
            opening.previous_token.is_name):
          # This is a dictionary that's an argument to a function.
          if (self._FitsOnLine(previous, previous.matching_bracket) and
              previous.matching_bracket.next_token and
              (not opening.matching_bracket.next_token or
               opening.matching_bracket.next_token.value != '.') and
              _ScopeHasNoCommas(previous)):
            # Don't split before the key if:
            #   - The dictionary fits on a line, and
            #   - The function call isn't part of a builder-style call and
            #   - The dictionary has one entry and no trailing comma
            return False
      return True

    if (style.Get('SPLIT_BEFORE_DICT_SET_GENERATOR') and
        format_token.Subtype.DICT_SET_GENERATOR in current.subtypes):
      # Split before a dict/set generator.
      return True

    if (format_token.Subtype.DICTIONARY_VALUE in current.subtypes or
        (previous.is_pseudo_paren and previous.value == '(' and
         not current.is_comment)):
      # Split before the dictionary value if we can't fit every dictionary
      # entry on its own line.
      if not current.OpensScope():
        opening = _GetOpeningBracket(current)
        if not self._EachDictEntryFitsOnOneLine(opening):
          return style.Get('ALLOW_SPLIT_BEFORE_DICT_VALUE')

    if previous.value == '{':
      # Split if the dict/set cannot fit on one line and ends in a comma.
      closing = previous.matching_bracket
      if (not self._FitsOnLine(previous, closing) and
          closing.previous_token.value == ','):
        self.stack[-1].split_before_closing_bracket = True
        return True

    ###########################################################################
    # Argument List Splitting
    if (style.Get('SPLIT_BEFORE_NAMED_ASSIGNS') and not current.is_comment and
        format_token.Subtype.DEFAULT_OR_NAMED_ASSIGN_ARG_LIST in
        current.subtypes):
      if (previous.value not in {'=', ':', '*', '**'} and
          current.value not in ':=,)' and not _IsFunctionDefinition(previous)):
        # If we're going to split the lines because of named arguments, then we
        # want to split after the opening bracket as well. But not when this is
        # part of a function definition.
        if previous.value == '(':
          # Make sure we don't split after the opening bracket if the
          # continuation indent is greater than the opening bracket:
          #
          #  a(
          #      b=1,
          #      c=2)
          if (self._FitsOnLine(previous, previous.matching_bracket) and
              unwrapped_line.IsSurroundedByBrackets(previous)):
            # An argument to a function is a function call with named
            # assigns.
            return False

          column = self.column - self.stack[-1].last_space
          return column > style.Get('CONTINUATION_INDENT_WIDTH')

        opening = _GetOpeningBracket(current)
        if opening:
          arglist_length = (
              opening.matching_bracket.total_length - opening.total_length +
              self.stack[-1].indent)
          return arglist_length > self.column_limit

    if (current.value not in '{)' and previous.value == '(' and
        self._ArgumentListHasDictionaryEntry(current)):
      return True

    if style.Get('SPLIT_ARGUMENTS_WHEN_COMMA_TERMINATED'):
      # Split before arguments in a function call or definition if the
      # arguments are terminated by a comma.
      opening = _GetOpeningBracket(current)
      if opening and opening.previous_token and opening.previous_token.is_name:
        if previous.value in '(,':
          if opening.matching_bracket.previous_token.value == ',':
            return True

    if ((current.is_name or current.value in {'*', '**'}) and
        previous.value == ','):
      # If we have a function call within an argument list and it won't fit on
      # the remaining line, but it will fit on a line by itself, then go ahead
      # and split before the call.
      opening = _GetOpeningBracket(current)
      if (opening and opening.value == '(' and opening.previous_token and
          (opening.previous_token.is_name or
           opening.previous_token.value in {'*', '**'})):
        is_func_call = False
        opening = current
        while opening:
          if opening.value == '(':
            is_func_call = True
            break
          if (not (opening.is_name or opening.value in {'*', '**'}) and
              opening.value != '.'):
            break
          opening = opening.next_token

        if is_func_call:
          if (not self._FitsOnLine(current, opening.matching_bracket) or
              (opening.matching_bracket.next_token and
               opening.matching_bracket.next_token.value != ',' and
               not opening.matching_bracket.next_token.ClosesScope())):
            return True

    pprevious = previous.previous_token
    if (current.is_name and pprevious and pprevious.is_name and
        previous.value == '('):
      if (not self._FitsOnLine(previous, previous.matching_bracket) and
          _IsFunctionCallWithArguments(current)):
        # There is a function call, with more than 1 argument, where the first
        # argument is itself a function call with arguments.  In this specific
        # case, if we split after the first argument's opening '(', then the
        # formatting will look bad for the rest of the arguments. E.g.:
        #
        #     outer_function_call(inner_function_call(
        #         inner_arg1, inner_arg2),
        #                         outer_arg1, outer_arg2)
        #
        # Instead, enforce a split before that argument to keep things looking
        # good.
        return True

    if (previous.OpensScope() and not current.OpensScope() and
        not current.is_comment and
        format_token.Subtype.SUBSCRIPT_BRACKET not in previous.subtypes):
      if pprevious and not pprevious.is_keyword and not pprevious.is_name:
        # We want to split if there's a comment in the container.
        token = current
        while token != previous.matching_bracket:
          if token.is_comment:
            return True
          token = token.next_token
      if previous.value == '(':
        pptoken = previous.previous_token
        if not pptoken or not pptoken.is_name:
          # Split after the opening of a tuple if it doesn't fit on the current
          # line and it's not a function call.
          if self._FitsOnLine(previous, previous.matching_bracket):
            return False
        elif not self._FitsOnLine(previous, previous.matching_bracket):
          if len(previous.container_elements) == 1:
            return False

          elements = previous.container_elements + [previous.matching_bracket]
          i = 1
          while i < len(elements):
            if (not elements[i - 1].OpensScope() and
                not self._FitsOnLine(elements[i - 1], elements[i])):
              return True
            i += 1

          if (self.column_limit - self.column) / float(self.column_limit) < 0.3:
            # Try not to squish all of the arguments off to the right.
            return True
      else:
        # Split after the opening of a container if it doesn't fit on the
        # current line.
        if not self._FitsOnLine(previous, previous.matching_bracket):
          return True

    ###########################################################################
    # Original Formatting Splitting
    # These checks rely upon the original formatting. This is in order to
    # attempt to keep hand-written code in the same condition as it was before.
    # However, this may cause the formatter to fail to be idempotent.
    if (style.Get('SPLIT_BEFORE_BITWISE_OPERATOR') and current.value in '&|' and
        previous.lineno < current.lineno):
      # Retain the split before a bitwise operator.
      return True

    if (current.is_comment and
        previous.lineno < current.lineno - current.value.count('\n')):
      # If a comment comes in the middle of an unwrapped line (like an if
      # conditional with comments interspersed), then we want to split if the
      # original comments were on a separate line.
      return True

    return False

  def AddTokenToState(self, newline, dry_run, must_split=False):
    """Add a token to the format decision state.

    Allow the heuristic to try out adding the token with and without a newline.
    Later on, the algorithm will determine which one has the lowest penalty.

    Arguments:
      newline: (bool) Add the token on a new line if True.
      dry_run: (bool) Don't commit whitespace changes to the FormatToken if
        True.
      must_split: (bool) A newline was required before this token.

    Returns:
      The penalty of splitting after the current token.
    """
    penalty = 0
    if newline:
      penalty = self._AddTokenOnNewline(dry_run, must_split)
    else:
      self._AddTokenOnCurrentLine(dry_run)

    penalty += self._CalculateComprehensionState(newline)

    return self.MoveStateToNextToken() + penalty

  def _AddTokenOnCurrentLine(self, dry_run):
    """Puts the token on the current line.

    Appends the next token to the state and updates information necessary for
    indentation.

    Arguments:
      dry_run: (bool) Commit whitespace changes to the FormatToken if True.
    """
    current = self.next_token
    previous = current.previous_token

    spaces = current.spaces_required_before
    if not dry_run:
      current.AddWhitespacePrefix(newlines_before=0, spaces=spaces)

    if previous.OpensScope():
      if not current.is_comment:
        # Align closing scopes that are on a newline with the opening scope:
        #
        #     foo = [a,
        #            b,
        #           ]
        self.stack[-1].closing_scope_indent = self.column - 1
        if style.Get('ALIGN_CLOSING_BRACKET_WITH_VISUAL_INDENT'):
          self.stack[-1].closing_scope_indent += 1
        self.stack[-1].indent = self.column + spaces
      else:
        self.stack[-1].closing_scope_indent = (
            self.stack[-1].indent - style.Get('CONTINUATION_INDENT_WIDTH'))

    self.column += spaces

  def _AddTokenOnNewline(self, dry_run, must_split):
    """Adds a line break and necessary indentation.

    Appends the next token to the state and updates information necessary for
    indentation.

    Arguments:
      dry_run: (bool) Don't commit whitespace changes to the FormatToken if
        True.
      must_split: (bool) A newline was required before this token.

    Returns:
      The split penalty for splitting after the current state.
    """
    current = self.next_token
    previous = current.previous_token

    self.column = self._GetNewlineColumn()

    if not dry_run:
      indent_level = self.line.depth
      spaces = self.column
      if spaces:
        spaces -= indent_level * style.Get('INDENT_WIDTH')
      current.AddWhitespacePrefix(
          newlines_before=1, spaces=spaces, indent_level=indent_level)

    if not current.is_comment:
      self.stack[-1].last_space = self.column
    self.lowest_level_on_line = self.paren_level

    if (previous.OpensScope() or
        (previous.is_comment and previous.previous_token is not None and
         previous.previous_token.OpensScope())):
      self.stack[-1].closing_scope_indent = max(
          0, self.stack[-1].indent - style.Get('CONTINUATION_INDENT_WIDTH'))

      self.stack[-1].split_before_closing_bracket = True

    # Calculate the split penalty.
    penalty = current.split_penalty

    if must_split:
      # Don't penalize for a must split.
      return penalty

    if previous.is_pseudo_paren and previous.value == '(':
      # Small penalty for splitting after a pseudo paren.
      penalty += 50

    # Add a penalty for each increasing newline we add, but don't penalize for
    # splitting before an if-expression or list comprehension.
    if current.value not in {'if', 'for'}:
      last = self.stack[-1]
      last.num_line_splits += 1
      penalty += (
          style.Get('SPLIT_PENALTY_FOR_ADDED_LINE_SPLIT') *
          last.num_line_splits)

    if current.OpensScope() and previous.OpensScope():
      # Prefer to keep opening brackets coalesced (unless it's at the beginning
      # of a function call).
      pprev = previous.previous_token
      if not pprev or not pprev.is_name:
        penalty += 10

    return penalty + 10

  def MoveStateToNextToken(self):
    """Calculate format decision state information and move onto the next token.

    Before moving onto the next token, we first calculate the format decision
    state given the current token and its formatting decisions. Then the format
    decision state is set up so that the next token can be added.

    Returns:
      The penalty for the number of characters over the column limit.
    """
    current = self.next_token
    if not current.OpensScope() and not current.ClosesScope():
      self.lowest_level_on_line = min(self.lowest_level_on_line,
                                      self.paren_level)

    # If we encounter an opening bracket, we add a level to our stack to prepare
    # for the subsequent tokens.
    if current.OpensScope():
      last = self.stack[-1]
      new_indent = style.Get('CONTINUATION_INDENT_WIDTH') + last.last_space

      self.stack.append(_ParenState(new_indent, self.stack[-1].last_space))
      self.paren_level += 1

    # If we encounter a closing bracket, we can remove a level from our
    # parenthesis stack.
    if len(self.stack) > 1 and current.ClosesScope():
      if format_token.Subtype.DICTIONARY_KEY_PART in current.subtypes:
        self.stack[-2].last_space = self.stack[-2].indent
      else:
        self.stack[-2].last_space = self.stack[-1].last_space
      self.stack.pop()
      self.paren_level -= 1

    is_multiline_string = current.is_string and '\n' in current.value
    if is_multiline_string:
      # This is a multiline string. Only look at the first line.
      self.column += len(current.value.split('\n')[0])
    elif not current.is_pseudo_paren:
      self.column += len(current.value)

    self.next_token = self.next_token.next_token

    # Calculate the penalty for overflowing the column limit.
    penalty = 0
    if not current.is_pylint_comment and self.column > self.column_limit:
      excess_characters = self.column - self.column_limit
      penalty += style.Get('SPLIT_PENALTY_EXCESS_CHARACTER') * excess_characters

    if is_multiline_string:
      # If this is a multiline string, the column is actually the
      # end of the last line in the string.
      self.column = len(current.value.split('\n')[-1])

    return penalty

  def _CalculateComprehensionState(self, newline):
    """Makes required changes to comprehension state.

    Args:
      newline: Whether the current token is to be added on a newline.

    Returns:
      The penalty for the token-newline combination given the current
      comprehension state.
    """
    current = self.next_token
    previous = current.previous_token
    top_of_stack = self.comp_stack[-1] if self.comp_stack else None
    penalty = 0

    if top_of_stack is not None:
      # Check if the token terminates the current comprehension.
      if current == top_of_stack.closing_bracket:
        last = self.comp_stack.pop()
        # Lightly penalize comprehensions that are split across multiple lines.
        if last.has_interior_split:
          penalty += style.Get('SPLIT_PENALTY_COMPREHENSION')

        return penalty

      if newline:
        top_of_stack.has_interior_split = True

    if (format_token.Subtype.COMP_EXPR in current.subtypes and
        format_token.Subtype.COMP_EXPR not in previous.subtypes):
      self.comp_stack.append(object_state.ComprehensionState(current))
      return penalty

    if (current.value == 'for' and
        format_token.Subtype.COMP_FOR in current.subtypes):
      if top_of_stack.for_token is not None:
        # Treat nested comprehensions like normal comp_if expressions.
        # Example:
        #     my_comp = [
        #         a.qux + b.qux
        #         for a in foo
        #   -->   for b in bar   <--
        #         if a.zut + b.zut
        #     ]
        if (style.Get('SPLIT_COMPLEX_COMPREHENSION') and
            top_of_stack.has_split_at_for != newline and
            (top_of_stack.has_split_at_for or
             not top_of_stack.HasTrivialExpr())):
          penalty += split_penalty.UNBREAKABLE
      else:
        top_of_stack.for_token = current
        top_of_stack.has_split_at_for = newline

        # Try to keep trivial expressions on the same line as the comp_for.
        if (style.Get('SPLIT_COMPLEX_COMPREHENSION') and newline and
            top_of_stack.HasTrivialExpr()):
          penalty += split_penalty.CONNECTED

    if (format_token.Subtype.COMP_IF in current.subtypes and
        format_token.Subtype.COMP_IF not in previous.subtypes):
      # Penalize breaking at comp_if when it doesn't match the newline structure
      # in the rest of the comprehension.
      if (style.Get('SPLIT_COMPLEX_COMPREHENSION') and
          top_of_stack.has_split_at_for != newline and
          (top_of_stack.has_split_at_for or not top_of_stack.HasTrivialExpr())):
        penalty += split_penalty.UNBREAKABLE

    return penalty

  def _GetNewlineColumn(self):
    """Return the new column on the newline."""
    current = self.next_token
    previous = current.previous_token
    top_of_stack = self.stack[-1]

    if current.spaces_required_before > 2 or self.line.disable:
      return current.spaces_required_before

    if current.OpensScope():
      return top_of_stack.indent if self.paren_level else self.first_indent

    if current.ClosesScope():
      if (previous.OpensScope() or
          (previous.is_comment and previous.previous_token is not None and
           previous.previous_token.OpensScope())):
        return max(0,
                   top_of_stack.indent - style.Get('CONTINUATION_INDENT_WIDTH'))
      return top_of_stack.closing_scope_indent

    if (previous and previous.is_string and current.is_string and
        format_token.Subtype.DICTIONARY_VALUE in current.subtypes):
      return previous.column

    if style.Get('INDENT_DICTIONARY_VALUE'):
      if previous and (previous.value == ':' or previous.is_pseudo_paren):
        if format_token.Subtype.DICTIONARY_VALUE in current.subtypes:
          return top_of_stack.indent

    if (_IsCompoundStatement(self.line.first) and
        (not style.Get('DEDENT_CLOSING_BRACKETS') or
         style.Get('SPLIT_BEFORE_FIRST_ARGUMENT'))):
      token_indent = (
          len(self.line.first.whitespace_prefix.split('\n')[-1]) +
          style.Get('INDENT_WIDTH'))
      if token_indent == top_of_stack.indent:
        return top_of_stack.indent + style.Get('CONTINUATION_INDENT_WIDTH')

    return top_of_stack.indent

  def _FitsOnLine(self, start, end):
    """Determines if line between start and end can fit on the current line."""
    length = end.total_length - start.total_length
    if not start.is_pseudo_paren:
      length += len(start.value)
    return length + self.column <= self.column_limit

  def _EachDictEntryFitsOnOneLine(self, opening):
    """Determine if each dict elems can fit on one line."""

    def PreviousNonCommentToken(tok):
      tok = tok.previous_token
      while tok.is_comment:
        tok = tok.previous_token
      return tok

    def ImplicitStringConcatenation(tok):
      num_strings = 0
      if tok.is_pseudo_paren:
        tok = tok.next_token
      while tok.is_string:
        num_strings += 1
        tok = tok.next_token
      return num_strings > 1

    closing = opening.matching_bracket
    entry_start = opening.next_token
    current = opening.next_token.next_token

    while current and current != closing:
      if format_token.Subtype.DICTIONARY_KEY in current.subtypes:
        prev = PreviousNonCommentToken(current)
        length = prev.total_length - entry_start.total_length
        length += len(entry_start.value)
        if length + self.stack[-2].indent >= self.column_limit:
          return False
        entry_start = current
      if current.OpensScope():
        if ((current.value == '{' or
             (current.is_pseudo_paren and current.next_token.value == '{') and
             format_token.Subtype.DICTIONARY_VALUE in current.subtypes) or
            ImplicitStringConcatenation(current)):
          # A dictionary entry that cannot fit on a single line shouldn't matter
          # to this calculation. If it can't fit on a single line, then the
          # opening should be on the same line as the key and the rest on
          # newlines after it. But the other entries should be on single lines
          # if possible.
          if current.matching_bracket:
            current = current.matching_bracket
          while current:
            if current == closing:
              return True
            if format_token.Subtype.DICTIONARY_KEY in current.subtypes:
              entry_start = current
              break
            current = current.next_token
        else:
          current = current.matching_bracket
      else:
        current = current.next_token

    # At this point, current is the closing bracket. Go back one to get the the
    # end of the dictionary entry.
    current = PreviousNonCommentToken(current)
    length = current.total_length - entry_start.total_length
    length += len(entry_start.value)
    return length + self.stack[-2].indent <= self.column_limit

  def _ArgumentListHasDictionaryEntry(self, token):
    """Check if the function argument list has a dictionary as an arg."""
    if _IsArgumentToFunction(token):
      while token:
        if token.value == '{':
          length = token.matching_bracket.total_length - token.total_length
          return length + self.stack[-2].indent > self.column_limit
        if token.ClosesScope():
          break
        if token.OpensScope():
          token = token.matching_bracket
        token = token.next_token
    return False


_COMPOUND_STMTS = frozenset(
    {'for', 'while', 'if', 'elif', 'with', 'except', 'def', 'class'})


def _IsCompoundStatement(token):
  if token.value == 'async':
    token = token.next_token
  return token.value in _COMPOUND_STMTS


def _IsFunctionDef(token):
  if token.value == 'async':
    token = token.next_token
  return token.value == 'def'


def _IsFunctionCallWithArguments(token):
  while token:
    if token.value == '(':
      token = token.next_token
      return token and token.value != ')'
    elif token.name not in {'NAME', 'DOT', 'EQUAL'}:
      break
    token = token.next_token
  return False


def _IsArgumentToFunction(token):
  bracket = unwrapped_line.IsSurroundedByBrackets(token)
  if not bracket or bracket.value != '(':
    return False
  previous = bracket.previous_token
  return previous and previous.is_name


def _GetLengthOfSubtype(token, subtype, exclude=None):
  current = token
  while (current.next_token and subtype in current.subtypes and
         (exclude is None or exclude not in current.subtypes)):
    current = current.next_token
  return current.total_length - token.total_length + 1


def _GetOpeningBracket(current):
  """Get the opening bracket containing the current token."""
  if current.matching_bracket and not current.is_pseudo_paren:
    return current.matching_bracket
  while current:
    if current.ClosesScope():
      current = current.matching_bracket
    elif current.is_pseudo_paren:
      current = current.previous_token
    elif current.OpensScope():
      return current
    current = current.previous_token
  return None


def _LastTokenInLine(current):
  while not current.is_comment and current.next_token:
    current = current.next_token
  return current


def _IsFunctionDefinition(current):
  prev = current.previous_token
  return (current.value == '(' and prev and
          format_token.Subtype.FUNC_DEF in prev.subtypes)


def _IsLastScopeInLine(current):
  while current:
    current = current.next_token
    if current and current.OpensScope():
      return False
  return True


def _IsSingleElementTuple(token):
  """Check if it's a single-element tuple."""
  close = token.matching_bracket
  token = token.next_token
  num_commas = 0
  while token != close:
    if token.value == ',':
      num_commas += 1
    if token.OpensScope():
      token = token.matching_bracket
    else:
      token = token.next_token
  return num_commas == 1


def _ScopeHasNoCommas(token):
  """Check if the scope has no commas."""
  close = token.matching_bracket
  token = token.next_token
  while token != close:
    if token.value == ',':
      return False
    if token.OpensScope():
      token = token.matching_bracket
    else:
      token = token.next_token
  return True


class _ParenState(object):
  """Maintains the state of the bracket enclosures.

  A stack of _ParenState objects are kept so that we know how to indent relative
  to the brackets.

  Attributes:
    indent: The column position to which a specified parenthesis level needs to
      be indented.
    last_space: The column position of the last space on each level.
    split_before_closing_bracket: Whether a newline needs to be inserted before
      the closing bracket. We only want to insert a newline before the closing
      bracket if there also was a newline after the beginning left bracket.
    num_line_splits: Number of line splits this _ParenState contains already.
      Each subsequent line split gets an increasing penalty.
  """

  # TODO(morbo): This doesn't track "bin packing."

  def __init__(self, indent, last_space):
    self.indent = indent
    self.last_space = last_space
    self.closing_scope_indent = 0
    self.split_before_closing_bracket = False
    self.num_line_splits = 0

  def Clone(self):
    state = _ParenState(self.indent, self.last_space)
    state.closing_scope_indent = self.closing_scope_indent
    state.split_before_closing_bracket = self.split_before_closing_bracket
    state.num_line_splits = self.num_line_splits
    return state

  def __repr__(self):
    return '[indent::%d, last_space::%d, closing_scope_indent::%d]' % (
        self.indent, self.last_space, self.closing_scope_indent)

  def __eq__(self, other):
    return hash(self) == hash(other)

  def __ne__(self, other):
    return not self == other

  def __hash__(self, *args, **kwargs):
    return hash((self.indent, self.last_space, self.closing_scope_indent,
                 self.split_before_closing_bracket, self.num_line_splits))
