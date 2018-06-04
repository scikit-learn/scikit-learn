# Copyright 2018 Google Inc. All Rights Reserved.
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
"""Identify containers for lib2to3 trees.

This module identifies containers and the elements in them. Each element points
to the opening bracket and vice-versa.

  IdentifyContainers(): the main function exported by this module.
"""

from yapf.yapflib import pytree_utils
from yapf.yapflib import pytree_visitor


def IdentifyContainers(tree):
  """Run the identify containers visitor over the tree, modifying it in place.

  Arguments:
    tree: the top-level pytree node to annotate with subtypes.
  """
  identify_containers = _IdentifyContainers()
  identify_containers.Visit(tree)


class _IdentifyContainers(pytree_visitor.PyTreeVisitor):
  """_IdentifyContainers - see file-level docstring for detailed description."""

  def Visit_trailer(self, node):  # pylint: disable=invalid-name
    for child in node.children:
      self.Visit(child)

    if len(node.children) != 3:
      return
    if pytree_utils.NodeName(node.children[0]) != 'LPAR':
      return

    if pytree_utils.NodeName(node.children[1]) == 'arglist':
      for child in node.children[1].children:
        pytree_utils.SetOpeningBracket(
            pytree_utils.FirstLeafNode(child), node.children[0])
    else:
      pytree_utils.SetOpeningBracket(
          pytree_utils.FirstLeafNode(node.children[1]), node.children[0])

  def Visit_atom(self, node):  # pylint: disable=invalid-name
    for child in node.children:
      self.Visit(child)

    if len(node.children) != 3:
      return
    if pytree_utils.NodeName(node.children[0]) != 'LPAR':
      return

    for child in node.children[1].children:
      pytree_utils.SetOpeningBracket(
          pytree_utils.FirstLeafNode(child), node.children[0])
