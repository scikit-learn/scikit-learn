"""
Test the format_stack module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2010 Gael Varoquaux
# License: BSD Style, 3 clauses.

from ..format_stack import safe_repr

import nose

################################################################################
class Nasty(object):
    def __repr__(self):
        raise ValueError

    __str__ = __repr__

################################################################################
# Test safe_repr 

def test_safe_repr():
    """ Smoke test safe_repr on a nasty class.
    """
    nasty = Nasty()
    safe_repr(nasty)

