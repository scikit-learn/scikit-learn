"""
Unit tests for the stack formatting utilities
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2010 Gael Varoquaux
# License: BSD Style, 3 clauses.

import nose

from ..format_stack import safe_repr


###############################################################################

class Vicious(object):
    def __repr__(self):
        raise ValueError


def test_safe_repr():
    safe_repr(Vicious())
