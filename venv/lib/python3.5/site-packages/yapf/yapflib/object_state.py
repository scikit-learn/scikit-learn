# Copyright 2017 Google Inc. All Rights Reserved.
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
"""Represents the state of Python objects being formatted.

Objects (e.g., list comprehensions, dictionaries, etc.) have specific
requirements on how they're formatted. These state objects keep track of these
requirements.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


class ComprehensionState(object):
  """Maintains the state of list comprehension formatting decisions.

  A stack of ComprehensionState objects are kept to ensure that list
  comprehensions are wrapped with well-defined rules.

  Attributes:
    expr_token: The first token in the comprehension.
    for_token: The first 'for' token of the comprehension.
    has_split_at_for: Whether there is a newline immediately before the
        for_token.
    has_interior_split: Whether there is a newline within the comprehension.
        That is, a split somewhere after expr_token or before closing_bracket.
  """

  def __init__(self, expr_token):
    self.expr_token = expr_token
    self.for_token = None
    self.has_split_at_for = False
    self.has_interior_split = False

  def HasTrivialExpr(self):
    """Returns whether the comp_expr is "trivial" i.e. is a single token."""
    return self.expr_token.next_token.value == 'for'

  @property
  def opening_bracket(self):
    return self.expr_token.previous_token

  @property
  def closing_bracket(self):
    return self.opening_bracket.matching_bracket

  def Clone(self):
    clone = ComprehensionState(self.expr_token)
    clone.for_token = self.for_token
    clone.has_split_at_for = self.has_split_at_for
    clone.has_interior_split = self.has_interior_split
    return clone

  def __repr__(self):
    return ('[opening_bracket::%s, for_token::%s, has_split_at_for::%s,'
            ' has_interior_split::%s, has_trivial_expr::%s]' %
            (self.opening_bracket, self.for_token, self.has_split_at_for,
             self.has_interior_split, self.HasTrivialExpr()))

  def __eq__(self, other):
    return hash(self) == hash(other)

  def __ne__(self, other):
    return not self == other

  def __hash__(self, *args, **kwargs):
    return hash((self.expr_token, self.for_token, self.has_split_at_for,
                 self.has_interior_split))
